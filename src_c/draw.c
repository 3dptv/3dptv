/****************************************************************************

Author/Copyright:		Jochen Willneff

Address:				Institute of Geodesy and Photogrammetry
						ETH - Hoenggerberg
						CH - 8093 Zurich

Creation Date:			end of 97
	
Description:			different drawing function to interact
						with Tcl/Tk
	
Routines contained: 	drawcross, drawvector, draw_pnr, mark_detections
						mark_correspondences, mark_corr, mark_track_c

****************************************************************************/


/* --- Improved zooming functionality,  ad holten 2013 --- */ 
#include "ptv.h"

Marker *markers[4];		/* lists with objects, like crosses, that are to be redrawn after zooming */
int    nmarkers[4], maxnmarkers[4];
BOOL   skip_adding;		// if set, draw objects without saving

void expand_markersbuffer(int imgnr)
{
	int i;
	Marker *oldmarkers = markers[imgnr];
	
	maxnmarkers[imgnr] = nmarkers[imgnr] + 1000;
	markers[imgnr] = (Marker*) malloc(maxnmarkers[imgnr]*sizeof(Marker));
	if (markers[imgnr]>0) { 
		for (i=0; i<nmarkers[imgnr]; i++)
			markers[imgnr][i] = oldmarkers[i];
		free(oldmarkers);
	}
}


void clear_drawnobjectslist()
{
	int i_img;
	for (i_img=0; i_img<n_img; i_img++) {
		// if (nmarkers[i_img] > 0 && markers[i_img] != NULL)
		if (markers[i_img] != NULL)
			free(markers[i_img]);
		nmarkers[i_img] = 0;
		markers[i_img]  = NULL;
		maxnmarkers[i_img] = 0;
	}
}

void add_mark(int x0, int y0, int imgnr)
{
	Marker* pm;
	int i;

    i = nmarkers[imgnr];
    if (i >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

	pm = &markers[imgnr][i];
    pm->x = x0;
    pm->y = y0;
    pm->type = OBJ_MARK;
    nmarkers[imgnr]++;
}

void add_cross(int x0, int y0, int size, int imgnr, char* color)
{
	Marker* pm;
    if (nmarkers[imgnr] >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

    pm = &markers[imgnr][nmarkers[imgnr]];
    pm->x = x0;
    pm->y = y0;
    pm->size = size;
    pm->type = OBJ_CROSS;
    strcpy(pm->color, color);
    nmarkers[imgnr]++;
}

void add_oval(int x0, int y0, int size, int imgnr, char* color)
{
	Marker* pm;
    if (nmarkers[imgnr] >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

    pm = &markers[imgnr][nmarkers[imgnr]];
    pm->x = x0;
    pm->y = y0;
    pm->size = size;
    pm->type = OBJ_OVAL;
    strcpy(pm->color, color);
    nmarkers[imgnr]++;
}

void add_pnr(int x, int y, int pnr, int imgnr, char* color, int minzoomlevel)
{
	Marker* pm;
    if (nmarkers[imgnr] >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

    pm = &markers[imgnr][nmarkers[imgnr]];
    pm->x = x;
    pm->y = y;
    pm->size = pnr;
	pm->minzoom = minzoomlevel;
    pm->type = OBJ_PNR;
    strcpy(pm->color, color);
    nmarkers[imgnr]++;
}
 
void add_value(int x, int y, double value, int imgnr, char* color, int minzoomlevel)
{
	Marker* pm;
    if (nmarkers[imgnr] >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

    pm = &markers[imgnr][nmarkers[imgnr]];
    pm->x = x;
    pm->y = y;
    pm->value = value;
	pm->minzoom = minzoomlevel;
    pm->type = OBJ_VALUE;
    strcpy(pm->color, color);
    nmarkers[imgnr]++;
}

void add_vector(int imgnr, int x0, int y0, int x1, int y1, int width, char *color)
{
	Marker* pm;
    if (nmarkers[imgnr] >= maxnmarkers[imgnr])
		expand_markersbuffer(imgnr);

    pm = &markers[imgnr][nmarkers[imgnr]];
    pm->x  = x0;
    pm->y  = y0;
    pm->x1 = x1;
    pm->y1 = y1;
    pm->size = width;
    pm->type = OBJ_VECTOR;
    strcpy(pm->color, color);
    nmarkers[imgnr]++;
}

int clearmarkers_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	clear_drawnobjectslist();
	return TCL_OK;
}

void draw_objects(Tcl_Interp* interp, int imgnr)
{
	// (re)drawing all the object in the buffer
    int i;
	Zoompar zoompar;

	get_tclzoomparms(interp, &zoompar, imgnr);
	skip_adding = TRUE;
    for (i=0; i<nmarkers[imgnr]; i++) {
        Marker* pm = &markers[imgnr][i];
        switch (pm->type) {
            case OBJ_MARK:
				drawmarker(interp, *pm, imgnr);
                break;
            case OBJ_CROSS:
				drawcross(interp, pm->x, pm->y, pm->size, imgnr, pm->color);
				break;
            case OBJ_OVAL:
				drawoval(interp, pm->x, pm->y, pm->size, imgnr, pm->color);
				break;
            case OBJ_PNR:
				if (zoompar.fac >= pm->minzoom)
					draw_pnr(interp, pm->x, pm->y, pm->size, imgnr, pm->color);
                break;
            case OBJ_VALUE:
				if (zoompar.fac >= pm->minzoom)
					draw_value(interp, pm->x, pm->y, pm->value, imgnr, pm->color);
                break;
            case OBJ_VECTOR:
				drawvector(interp, pm->x, pm->y, pm->x1, pm->y1, pm->size, imgnr, pm->color);
                break;
        }
    }
	skip_adding = FALSE;        // restore the adding of objects
}



int drawmarker(Tcl_Interp* interp, Marker mark, int imgnr)
{
    char val[256];
	if (!skip_adding) 		// add mark to the objects to be redrawn 
		add_mark(mark.x, mark.y, imgnr);

	img_to_view_coordinates(&mark.x, &mark.y, (double)mark.x, (double)mark.y, imgnr);
    sprintf(val, "mark .cam%d.pic %d %d %d", imgnr+1, mark.x, mark.y, imgnr+1);
    Tcl_Eval(interp, val);
    return TCL_OK;
}

int drawcross (Tcl_Interp* interp, int x0, int y0, int size, int imgnr, char color[])
{
	char val[256];
	if (!skip_adding) 		// add cross to the objects to be redrawn 
		add_cross(x0,y0,size,imgnr,color);
	
	img_to_view_coordinates(&x0, &y0, (double)x0, (double)y0, imgnr);
	sprintf(val, "markparticle %d %d %d %d %s", x0, y0, size, imgnr+1, color);
	Tcl_Eval(interp, val);
	return TCL_OK;
}

int drawoval (Tcl_Interp* interp, int x0, int y0, int size, int imgnr, char color[])
{
	char val[256];
	if (!skip_adding) 		// add cross to the objects to be redrawn 
		add_oval(x0,y0,size,imgnr,color);

	img_to_view_coordinates(&x0, &y0, (double)x0, (double)y0, imgnr);
	sprintf(val, "draw_oval %d %d %d %d %s", x0, y0, size, imgnr+1, color);
	Tcl_Eval(interp, val);
	return TCL_OK;
}

int drawvector (Tcl_Interp* interp, int x0, int y0, int x1, int y1, int width, int imgnr, char color[])
{
	char val[256];
	if (!skip_adding) 		// add vector to the objects to be redrawn 
		add_vector(imgnr,x0,y0,x1,y1,width,color);

	img_to_view_coordinates(&x0, &y0, (double)x0, (double)y0, imgnr);
	img_to_view_coordinates(&x1, &y1, (double)x1, (double)y1, imgnr);
	sprintf(val, "drawline %d %d %d %d %d %d %s", x0, y0, x1, y1, width, imgnr+1, color);
	Tcl_Eval(interp, val );
	return TCL_OK;
}

int draw_pnr (Tcl_Interp* interp, int x, int y, int pnr, int imgnr, char color[])
{
	// By calling the draw_pnr_autohide() with -999 as min. zoom level,
	// the point number will always be redrawn.
	return draw_pnr_autohide (interp, x, y, pnr, imgnr, color, -999, 0);
}

int draw_pnr_autohide (Tcl_Interp* interp, int x, int y, int pnr, 
	int imgnr, char color[], int minzoomlevel, int zoomlevel)
{
	char buf[256], val[256];
	if (!skip_adding) 		// add pnr to the objects to be redrawn 
		add_pnr(x,y,pnr,imgnr,color,minzoomlevel);

	if (zoomlevel >= minzoomlevel) {
		img_to_view_coordinates(&x, &y, (double)x, (double)y, imgnr);
		sprintf (val, "%d", pnr);
		sprintf(buf, "drawtext %d %d %s %d %s", x, y, val, imgnr+1, color); 
		Tcl_Eval(interp, buf);
	}
	return TCL_OK; 
}

int draw_value (Tcl_Interp* interp, int x, int y, double pnr, int imgnr, char color[])
{
	// By calling the draw_value_autohide() with -999 as min. zoom level,
	// the value will always be redrawn.
	return draw_value_autohide (interp, x, y, pnr, imgnr, color, -999, 0);
}

int draw_value_autohide (Tcl_Interp* interp, int x, int y, double pnr,
	int imgnr, char color[], int minzoomlevel, int zoomlevel)
{
	char buf[256], val[256];
	if (!skip_adding) 		// add value to the objects to be redrawn 
		add_value(x,y,pnr,imgnr,color,minzoomlevel);

	if (zoomlevel >= minzoomlevel) {
		img_to_view_coordinates(&x, &y, (double)x, (double)y, imgnr);
		sprintf (val, "%5.3f", pnr);
		sprintf(buf, "drawtext %d %d %s %d %s", x, y, val, imgnr+1, color); 
		Tcl_Eval(interp, buf);
	}
	return TCL_OK; 
}

#ifdef EVER_CALLED		// Unused function, ad holten 04-2013
void mark_detections (Tcl_Interp* interp, int nr)
/* draws crosses for detected points in a displayed image */
{
	int i,limx, limy, intx, inty;

	if (num[nr] == 0) {
		printf ("No points detected");
		return;
	}
	limy = imy/(2*zoom_f[nr]);
	limx = imx/(2*zoom_f[nr]);
	for (i=0; i<num[nr]; i++) {
		if (   (fabs(pix[nr][i].x-zoom_x[nr]) < limx)
			&& (fabs(pix[nr][i].y-zoom_y[nr]) < limy))
		{
			intx = (int)(imx/2+zoom_f[nr]*(pix[nr][i].x-zoom_x[nr]));
			inty = (int)(imy/2+zoom_f[nr]*(pix[nr][i].y-zoom_y[nr]));
			drawcross (interp, intx, inty, cr_sz , nr, "blue");
		}
	}
}
#endif 

#ifdef EVER_CALLED		// Unused function, ad holten 04-2013
void mark_correspondences (Tcl_Interp* interp, int nr)
/* draws crosses and numbers for corresponding points in a displayed window */
{
	int    i,j, pnr, lim, sum, intx, inty;
	double x, y;

	if (match == 0) return;

	lim = imx/(2*zoom_f[nr]);

	for (i=0; i<match; i++) {
		pnr = geo[nr][con[i].p[nr]].pnr;
		if (pnr < 0 || con[i].p[nr] < 0)  continue;

		x = pix[nr][pnr].x;  y = pix[nr][pnr].y;
		if ((fabs (x-zoom_x[nr]) < lim) && (fabs (y-zoom_y[nr]) < lim)) {
			intx = (int) ( imx/2 + zoom_f[nr] * (x-zoom_x[nr]));
			inty = (int) ( imy/2 + zoom_f[nr] * (y-zoom_y[nr]));

			/* check whether quadruplet, triplet or pair -> select color */
			for (j=0, sum=0; j<4; j++)
				if (con[i].p[j] > 0) sum++; 
			if (sum == 2) sprintf(buf ,"yellow");
			if (sum == 3) sprintf(buf ,"green");
			if (sum == 4) sprintf(buf ,"red");

			drawcross (interp, intx, inty, cr_sz, nr, buf);

			/* draw point number */
			// if (examine && zoom_f[nr] > 2) {
			//	  /* number of established correspondence */
			//	  draw_pnr (interp, intx+3 , inty+3, i, nr, "white");
			// } 
		}
	}
}
#endif

#ifdef EVER_CALLED		// Unused function, ad holten 04-2013
void mark_corr (Tcl_Interp* interp, int nr)
/* draws crosses and numbers for corresponding points in a displayed window */
{
	int    i,j, pnr, sum, intx, inty;
	double x, y;
  
	if (match == 0) return;
  
	for (i=0; i<match; i++) {
		pnr = geo[nr][con[i].p[nr]].pnr;
		if (pnr < 0|| con[i].p[nr] < 0) continue;

		x = pix[nr][pnr].x;  y = pix[nr][pnr].y;

		intx = (int) (imx/2 + zoom_f[nr] * (x-zoom_x[nr]));
		inty = (int) (imy/2 + zoom_f[nr] * (y-zoom_y[nr]));

		/* check whether quadruplet, triplet or pair -> select color */
		for (j=0, sum=0; j<4; j++)
			if (con[i].p[j] > 0) sum++;
		if (sum == 2) sprintf(buf ,"yellow");
		if (sum == 3) sprintf(buf ,"green");
		if (sum == 4) sprintf(buf ,"red");

		/* i is the number of the established correspondence */
		draw_pnr (interp, intx+5 , inty, i, nr, "white"); 
		// hilfi = geo[nr][con[i].p[nr]].pnr;
		// drawcross (interp, intx, inty, cr_sz, nr, buf);
		// /* hilfi is the number of the detected point, see in targetlist */
		// draw_pnr (interp, intx+25 , inty-15, hilfi, nr, "red"); 
	}
}
#endif

BOOL get_tclzoomparms(Tcl_Interp* interp, Zoompar *zoompar, int i_img)
{
	char val[256];
	int  ic = i_img+1;

	sprintf(val, "get_viewsize %d", ic);
	Tcl_Eval(interp, val);
	sprintf(val, "xc%d", ic);
	zoompar->xc = atof( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );

	sprintf(val, "yc%d", ic);
	zoompar->yc = atof( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );

	sprintf(val, "vwx%d", ic);
	zoompar->vwx = atoi( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );

	sprintf(val, "vwy%d", ic);
	zoompar->vwy = atoi( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );

	sprintf(val, "fac%d", ic);
	zoompar->fac = atoi( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );

	sprintf(val, "fixed%d", ic);
	zoompar->fixed = atoi( Tcl_GetVar2(interp, "zoom", val,  TCL_GLOBAL_ONLY) );
	return TRUE;
}

BOOL isinview(double x, double y, int i_img)
{
	// With scrollbars, all the parts of the image can be scrolled in view.
	// so, checking if the point falls within the image area, is enough.
	return x >= 0 && x < imx && y >= 0 && y < imy;
}

int mark_track_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
/* draws crosses for detected points in a displayed image */
{
	char seq_name[4][128];
	int  i_img, i_seq, h, intx, inty;
	Zoompar zoompar[4];

	cr_sz = atoi(Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY));

	fpp = fopen_rp ("parameters/sequence.par");			// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;
	for (i_img=0; i_img<4; i_img++)
		fscanf (fpp, "%s\n", seq_name[i_img]); 

	/* name of sequence */
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);

	sprintf (buf, "Show detected particles "); puts (buf);
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 2");
	Tcl_Eval(interp, ".text insert 2 $tbuf");


	for (i_img=0; i_img<n_img; i_img++)					// added, ad holten, 04-2013
		get_tclzoomparms(interp, &zoompar[i_img], i_img);
	clear_drawnobjectslist();

	/* track sequence */
	for (i_seq=seq_first; i_seq<=seq_last; i_seq++) {
		read_ascii_data(i_seq);
		/* treat the cameras one after the other */
		for (i_img=0; i_img<n_img; i_img++) {
			for (h=0; h<nt4[3][i_img]; h++) {
				// changed the next test with isinview(), ad holten, 04-2013
				// if ( ( fabs(t4[3][i_img][h].x-zoom_x[i_img]) < imx/(2*zoom_f[i_img]))
				//   && ( fabs(t4[3][i_img][h].y-zoom_y[i_img]) < imy/(2*zoom_f[i_img])) )
				if (isinview(t4[3][i_img][h].x, t4[3][i_img][h].y, i_img))
				{
					// replaced the next conversion by img_to_view_coordinates(), ad holten, 04-2013
					// intx = (int)(imx/2+zoom_f[i_img]*(t4[3][i_img][h].x-zoom_x[i_img]));
					// inty = (int)(imy/2+zoom_f[i_img]*(t4[3][i_img][h].y-zoom_y[i_img]));		  
					img_to_view_coordinates(&intx, &inty, 
						(int)t4[3][i_img][h].x, (int)t4[3][i_img][h].y, i_img);
					if (t4[3][i_img][h].tnr > -1) { 
						drawcross (interp, intx, inty, cr_sz+1, i_img, "green");
						if (zoompar[i_img].fac >= 6) {
							draw_pnr (interp, intx, inty+10, i_seq, i_img, "orange");
							draw_pnr (interp, intx, inty, t4[3][i_img][h].tnr, i_img, "green");
						}
					}
					else
						drawcross (interp, intx, inty, cr_sz, i_img, "blue");
				}
			}
			Tcl_Eval(interp, "update idletasks");		   
		}
	}
	sprintf(val, "...done");
	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text insert 3 $tbuf");

	return TCL_OK;	
}

int trajectories_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
/* draws crosses for detected points in a displayed image */
{
	int k, intx1, inty1, intx2, inty2;
	int i, anz1, anz2, m, j;
	FILE *fp1;
	char val[256];
	vector *line1, *line2;
	double color;
	coord_2d p1[4], p2[4];

	cr_sz = atoi(Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY));

	fpp = fopen_rp ("parameters/sequence.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	/* name of sequence */
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);

	sprintf (buf, "Show trajectories "); puts (buf);
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 2");
	Tcl_Eval(interp, ".text insert 2 $tbuf");

	line1 = line2 = NULL;		// added, ad holten, 12-2012
	for (i=seq_first; i<seq_last;i++) {
		// ad holten, 12-2012, replaced next lines
		//		 if      (i < 10)  sprintf (val, "res/ptv_is.%1d", i);
		//		 else if (i < 100) sprintf (val, "res/ptv_is.%2d", i);
		//		 else			   sprintf (val, "res/ptv_is.%3d", i);
		sprintf(val, "res/ptv_is.%d", i);

		fp1 = fopen_rp(val);		// replaced fopen, ad holten 12-2012			
		if (!fp1) break;
	  
		color = ((double)(i-seq_first))/((double)(seq_last-2-seq_first));
		fscanf (fp1,"%d\n", &anz1);
	  
		line1 = (vector*) malloc(anz1 * sizeof (vector));
		for (j=0; j<anz1; j++) {
			fscanf (fp1, "%d\n", &line1[j].p);
			fscanf (fp1, "%d\n", &line1[j].n);
			fscanf (fp1, "%lf\n", &line1[j].x1);
			fscanf (fp1, "%lf\n", &line1[j].y1);
			fscanf (fp1, "%lf\n", &line1[j].z1);
		}
		strcpy(val, "");	 
		fclose (fp1);

		/* read next time step */	  
		// ad holten, 12-2012, replaced next lines,
		//		if (i+1 < 10)       sprintf (val, "res/ptv_is.%1d", i+1);
		//		else if (i+1 < 100) sprintf (val, "res/ptv_is.%2d", i+1);
		//		else				sprintf (val, "res/ptv_is.%3d", i+1);
		sprintf (val, "res/ptv_is.%d", i+1);

		fp1 = fopen_rp(val);			// replaced fopen, ad holten 12-2012
		if (!fp1) break;

		fscanf (fp1,"%d\n", &anz2);
		line2 = (vector*) calloc(anz2, sizeof (vector));

		for (j=0; j<anz2; j++) {
			fscanf (fp1, "%d\n", &line2[j].p);
			fscanf (fp1, "%d\n", &line2[j].n);
			fscanf (fp1, "%lf\n", &line2[j].x1);
			fscanf (fp1, "%lf\n", &line2[j].y1);
			fscanf (fp1, "%lf\n", &line2[j].z1);
		}
		fclose (fp1);
	  
		for(j=0; j<anz1; j++) { 	
			m = line1[j].n;

			if (m >= 0) {	   
				for (k=0; k<n_img; k++) {
					// code replaced, ad holten, 04-2013
					//   img_coord (line2[m].x1, line2[m].y1, line2[m].z1, Ex[k],I[k], G[k], ap[k], mmp, &p2[k].x, &p2[k].y);
					//   metric_to_pixel (p2[k].x, p2[k].y, imx,imy, pix_x,pix_y, &p2[k].x, &p2[k].y, chfield); 
					pixelcoord_from_3Dpnt(&p2[k].x, &p2[k].y, k, line2[m].x1, line2[m].y1, line2[m].z1);
		  
					// replaced next by isinview(), ad holten 04-2013
					//if ( fabs( p2[k].x-zoom_x[k]) < imx/(2*zoom_f[k])
					//   && ( fabs(p2[k].y-zoom_y[k]) < imy/(2*zoom_f[k])) )
					if (isinview(p2[k].x, p2[k].y, k))
					{	 
						// replaced next by img_to_view_coordinates(), ad holten 04-2013
						//intx1 = (int)(imx/2+zoom_f[k]*(p1[k].x-zoom_x[k]));
						//inty1 = (int)(imy/2+zoom_f[k]*(p1[k].y-zoom_y[k]));
						//intx2 = (int)(imx/2+zoom_f[k]*(p2[k].x-zoom_x[k]));
						//inty2 = (int)(imy/2+zoom_f[k]*(p2[k].y-zoom_y[k]));
						img_to_view_coordinates(&intx1, &inty1, (int)p1[k].x, (int)p1[k].y, k);
						img_to_view_coordinates(&intx2, &inty2, (int)p2[k].x, (int)p2[k].y, k);

						drawcross ( interp, intx1, inty1, cr_sz+1, k, "blue");	   
						drawcross ( interp, intx2, inty2, cr_sz+1, k, "red");
						drawvector (interp, intx1, inty1, intx2, inty2, 2, k, "green");
					}	 
				}
			}
		}
		Tcl_Eval(interp, "update idletasks");		   

		strcpy(val, "");
		free(line1); free(line2);
		line1 = line2 = NULL;			// added, ad holten, 12-2012
	}  /* end of sequence loop */
 
	if (line1 != NULL) free(line1);		// added, ad holten, 12-2012
	if (line2 != NULL) free(line2);		// ie. call free() if not deallocated yet
  
	sprintf(val, "...done");
	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text insert 3 $tbuf"); 
  
	return TCL_OK;
}


void clicked_to_imgcoordinates(int *px, int *py, double x, double y, int n)
{
	// Tcl will do the conversion
	*px = (int)x;
	*py = (int)y;
}

void img_to_view_coordinates(int *px, int *py, double x, double y, int n)
{
	// Tcl will do the conversion
	*px = (int)x;
	*py = (int)y;
}
