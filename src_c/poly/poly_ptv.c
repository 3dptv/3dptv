
#include "../ptv.h"

int map_method;
int examine3;						/* flag for selecting debugging output */

unsigned char *img[4];				/* image data */
unsigned char *img_mask[4]; 		/* mask data */
unsigned char *img_new[4];			/* image data for reducing mask */
unsigned char *img0[4]; 			/* image data for filtering etc */
unsigned char *zoomimg; 			/* zoom image data */

BOOL copyparmsfiles(char* savefolder, int loop);

int  cmp_dist(const void *arg1, const void *arg2);
void ReadDataFromDistribFile(char* pathname, int ic, int* pnelem, target* pixels, coord_3d* fixed);

int sortgrid_poly_mult(Tcl_Interp* interp, int iplane);
int sortgrid_poly_body(Tcl_Interp* interp);
int sortgrid_file_poly_body(Tcl_Interp* interp);
int calibration_poly_body(Tcl_Interp* interp);
int calibration_poly_mult(Tcl_Interp* interp);
int save_sorting_results(char *img_name, coord_3d *fix, target *pix, int nfix);


void build_indexedname(char* pathname, char* basename, int index)
{
    char drive[_MAX_DRIVE], dir[_MAX_DIR], fname[_MAX_FNAME], ext[_MAX_EXT];
    char buf[64];
    _splitpath(basename, drive, dir, fname, ext);
    sprintf(buf, "%d%s", index, ext);
    _makepath(pathname, drive,dir,fname,buf);
}

void reset_target(target* ptarget)
{
    ptarget->pnr = (int)(ptarget->x = ptarget->y = -999);
    ptarget->n = ptarget->nx = ptarget->ny = ptarget->sumg = 0;
}

void savedataset(char* filename, coord_3d* fixpnts, coord_2d* pixpnts, int n)
{
    int i;
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Can't open debug output file %s\n", filename);
        return;
    }
    for (i=0; i<n; i++) {
        fprintf(fp, "%10.3f %10.3f  %10.3f %10.3f %10.3f\n", 
            pixpnts[i].x, pixpnts[i].y, fixpnts[i].x, fixpnts[i].y, fixpnts[i].z); 
    }
    fclose(fp);
}

void tcl_text(Tcl_Interp* interp, char* title, char* val)
{
    Tcl_SetVar(interp, "tbuf", title, TCL_GLOBAL_ONLY);
    Tcl_Eval(interp, ".text delete 2");
    Tcl_Eval(interp, ".text insert 2 $tbuf");
    Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
    Tcl_Eval(interp, ".text delete 3");
    Tcl_Eval(interp, ".text insert 3 $tbuf");
}

const char* getfilename(char* pathname)
{
    char   drive[_MAX_DRIVE], dir[_MAX_DIR], /*fname[_MAX_FNAME],*/ ext[_MAX_EXT];
    static char fname[_MAX_FNAME];
    _splitpath(pathname, drive, dir, fname, ext);
    // sprintf(fname, "%s%s", index, ext);
    strcat(fname, ext);
    return fname;
}


// void draw_fitted_fixpoints(Tcl_Interp* interp, int ic, int loop, FitInfo fit, coord_3d* fix, int nfix);
void draw_fiterror(Tcl_Interp* interp, double x1, double y1, int ic, FitInfo fit, coord_3d fix);
void find_grid_points(coord_3d *fixori, coord_2d *pixori, int *pcnt, double eps, 
                      FitInfo fit, coord_3d* fix, int nfix, target* pix, int npix);



/* --- polynomial version of the calib_cmd -------------------------------------------------- */

int calibration_poly_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
    int   sel, retval, sup, ncam, nori, nplanes, deg3d2d[2], deg2depi[2];
    int   rv, ic, i, j, iplane, ne, ori[4][MaxOri];
    char  filename[256], val[256], buf[256];
    FILE  *fp;

	/* read support of unsharp mask */
	fp = fopen ("parameters/unsharp_mask.par", "r");
	if (! fp)
		sup = 12;
	else {
		fscanf (fp, "%d\n", &sup); 
		fclose (fp);
	}

    /* Get Selection value from TclTk */
    sel = atoi( Tcl_GetVar(interp, "sel",  TCL_GLOBAL_ONLY) );
    *val = *buf = 0;

	rv = parse_polyptv_par("parameters/poly_ptv.par");

	if (usingZplanes) {
		/* Get the current plane index */
		iplane = atoi( Tcl_GetVar(interp, "cp(plane_id)", TCL_GLOBAL_ONLY) );
		iplane--;	// index start from 1 in dialog

		rv = parse_mult_calori_par("parameters/mult_calori.par", 
				iplane, &ncam, &nplanes, img_name, fixp_name, deg3d2d, deg2depi);
	}
	else
		rv = parse_body_calori_par("parameters/body_calori.par", 
				img_name, fixp_name, deg3d2d, deg2depi);
	if (!rv) return TCL_OK;
	
	// clearing all objects in the objects list
	if (sel == 1 || sel == 2 || sel == 6)
		clear_drawnobjectslist();

	retval = TCL_OK;
    switch (sel)
	{
    case 1: // --- read and show the calibration images ---
        // Tcl_Eval(interp,  "clearcam");
        for (i=0; i<n_img; i++) {
            read_image(interp, img_name[i], img[i]);

            sprintf(val, "camcanvas %d", i+1);
            Tcl_Eval(interp, val);
			sprintf(val, "newimage %d %f %f %d %d", i+1, 0.5, 0.5, 1, 0);
            Tcl_Eval(interp, val);
        }
        break;

    case 2:  puts("Detection procedure"); strcpy(val,"");		// same as in jw_ptv.c
        
        // Highpass Filtering
        pre_processing_c(clientData, interp, argc, argv);

		// copy the images because the target recognition will set greyvalues to zero
        for (i=0; i<n_img; i++)
            copy_images (img[i], img0[i]);

		/* target recognition */
		for (i=0; i<n_img; i++) {
			targ_rec (interp, img[i], img0[i], "parameters/detect_plate.par",
					  0, imx, 1, imy, pix[i], i, &num[i]);

			sprintf (buf,"image %d: %d,  ", i+1, num[i]);
			strcat(val, buf);

			if (num[i] > nmax) {
				printf("Aborted, too many targets detected!\n");	// ad holten, 12-2012
				exit(1);
			}
		}

		/* save pixel coord as approx. for template matching */
		if (examine) {
			for (i=0; i<n_img; i++) {
				sprintf (filename, "%s_pix", img_name[i]);
				fp1 = fopen (filename, "w");
				for (j=0; j<num[i]; j++)
					fprintf (fp1, "%4d	%8.3f  %8.3f\n",
						pix[i][j].pnr, pix[i][j].x, pix[i][j].y);

				fclose (fp1);
			}
		}

        sprintf(buf,"Number of detected targets, interaction enabled");
        Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
        Tcl_Eval(interp, ".text delete 2");
        Tcl_Eval(interp, ".text insert 2 $tbuf");
        Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
        Tcl_Eval(interp, ".text delete 3");
        Tcl_Eval(interp, ".text insert 3 $tbuf");
        break;

    case 3: // --- Manual Orientation ---
        pp1=0;  pp2=0;  pp3=0;  pp4=0;
        for (i=0; i<n_img; i++) {
            sprintf (buf, "%d targets remain", num[i]);
            puts (buf);
        }
		if (usingZplanes)
			rv = parse_multi_manori_par("parameters/mult_manori.par", iplane, ori, &ncam, &nori);
		else
			rv = parse_body_manori_par("parameters/body_manori.par", ori, &ncam, &nori);
		if (!rv) return TCL_OK;

        for (i=0; i<n_img; i++) {
			if (nori == 4)
				sprintf(val, "measure %d %d %d %d %d", ori[i][0], ori[i][1], ori[i][2], ori[i][3], i+1);
			else
				sprintf(val, "measure2 %d %d %d %d %d %d %d",
					ori[i][0], ori[i][1], ori[i][2], ori[i][3], ori[i][4], ori[i][5], i+1);
            Tcl_Eval(interp, val);
            for (j=0; j<nori; j++) {
                sprintf(val, "px%d", j);
                pix0[i][j].x = atoi( Tcl_GetVar(interp, val, TCL_GLOBAL_ONLY) );
                sprintf(val, "py%d", j);
                pix0[i][j].y = atoi( Tcl_GetVar(interp, val, TCL_GLOBAL_ONLY) );
            }
        }

		/* write measured coordinates to file for next trial */
		if (usingZplanes)
			sprintf(filename,"mult_manori%d.dat", iplane+1);
		else
			strcpy(filename, "body_manori.dat");

        fp1 = fopen (filename, "w");
        for (i=0; i<n_img; i++)
            for (j=0; j<nori; j++)
                fprintf (fp1, "%f %f\n", pix0[i][j].x, pix0[i][j].y);
        fclose (fp1);
        break;

    case 4:  /* read pixel coordinates of older pre-orientation */

		/* read point numbers of pre-clicked points */
		if (usingZplanes) {
			rv = parse_multi_manori_par("parameters/mult_manori.par", iplane, ori, &ncam, &nori);
			sprintf(filename,"mult_manori%d.dat",iplane+1);
		}
		else {
			rv = parse_body_manori_par("parameters/body_manori.par", ori, &ncam, &nori);
			strcpy(filename, "body_manori.dat");
		}
		// read coordinates of pre-clicked points 
        fp = fopen_rp(filename);
        if (!rv || !fp) break;
		ne = 0;
        for (ic=0; ic<n_img; ic++)
            for (i=0; i<nori; i++) {
                ne += fscanf(fp, "%lf%lf\n", &pix0[ic][i].x, &pix0[ic][i].y);
				drawcross(interp, (int) pix0[ic][i].x, (int) pix0[ic][i].y, cr_sz+2, ic, "red");
				draw_pnr (interp, (int) pix0[ic][i].x, (int) pix0[ic][i].y,	ori[ic][i], ic, "red");
            }
        fclose (fp);
        if (ne < n_img*nori*2)
            puts("Insufficient number of clicked points found !!!");
        else
			puts("Clicked points from file loaded.");
        break;


    case 5: // --- Sort grid points --- 
		if (usingZplanes) rv = sortgrid_poly_mult(interp, iplane);
		else              rv = sortgrid_poly_body(interp);
		if (!rv) return TCL_ERROR;
		break;

	case 6: // --- Calibration with grid points --- 
		if (usingZplanes) rv = calibration_poly_mult(interp); 
		else              rv = calibration_poly_body(interp);
		if (!rv) return TCL_ERROR;
        break;

	case 7: checkpoint_poly (interp);
		sprintf(val,"blue: planimetry,	 yellow: height");
		Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text delete 2");
		Tcl_Eval(interp, ".text insert 2 $tbuf");
		Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text delete 3");
		Tcl_Eval(interp, ".text insert 3 $tbuf");
		break;


	case 16: // -- Sort grid points using file --
		if (usingZplanes)
			puts ("For multi planes, sort grid points using files is not implemented.");
		else {
			puts ("Sort grid points using files");
			if (! sortgrid_file_poly_body(interp) )
				return TCL_ERROR;
		}
		break;

	case 10: // -- calibration from particles --
		puts ("Calibration from particles"); strcpy(buf, "");
		poly_shaking_c(clientData, interp, argc, argv);
		break;

	case 15: puts ("Show number on detected points");
		for (i=0; i<n_img; i++) {
			for (j=0; j<num[i]; j++) { 
				draw_pnr (interp, (int)pix[i][j].x, (int)pix[i][j].y, j, i, "blue");
			}
		}
		break; //Beat and Debashish Jan 2011 

	}
 	return TCL_OK;
}


int sortgrid_poly_body(Tcl_Interp* interp)
{
	FILE *fp;
	double Rsort;
	int i, j, k, ic, nident, degree[2], order, crossorder, loop, ix, iy, rv, ncam, ori[4][MaxOri];
    int deg3d2d[2], deg2depi[2];
    coord_3d *fixori;
    coord_2d *pixori;
    FitInfo  fit;
    target *orgpix;

	rv  = parse_body_calori_par("parameters/body_calori.par", img_name, fixp_name, deg3d2d, deg2depi);
	rv &= parse_body_manori_par("parameters/body_manori.par", ori, &ncam, &nori);
	if (!rv)
		return FALSE;

	puts("\nSort grid points");

	/* reading the 3D-coordinates of the calibration target */
	fp1 = fopen_rp(fixp_name);
	if (!fp1) return FALSE;
	k = 0;
	while (fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
			&fix[k].x, &fix[k].y, &fix[k].z) != EOF)  k++;
	fclose (fp1);
	nfix = k;

    // defining the search radius 
	Rsort = 10.0;
	fp = fopen("parameters/sortgrid.par", "r");
	if (fp) {
		if (fscanf (fp, "%lf", &Rsort) == 1)
			printf ("Sortgrid search radius: %.1f pixel (from sortgrid.par)\n", Rsort);
		fclose(fp);
	}
	else
        printf ("No parameters/sortgrid.par found, search radius set to %.1f.\n", Rsort);

    // construct and fill arrays used in the fitting process
	fixori = (coord_3d*) malloc(nfix * sizeof(coord_3d));
	pixori = (coord_2d*) malloc(nfix * sizeof(coord_2d));

	for (ic=0; ic<n_img; ic++) {
		printf("Camera: %d\n", ic+1);

		// find the 3D-coordinates of the clicked points
		for (i=0; i<nori; i++)
			for (k=0; k<nfix; k++)
				if (fix[k].pnr == ori[ic][i]) {
					fixori[i]   = fix[k];
					pixori[i].x = pix0[ic][i].x;
					pixori[i].y = pix0[ic][i].y;
				}

		// 1: fit the clicked 3D-points with 2D pixel coordinates
		// 2: using the fit found, try to find more grid points
		// 3: redo the cycle using higher order terms
		nident = nori;
		for (loop=0; loop<=1; loop++) {
			switch (loop) {
				case 0: order = 1; crossorder = 0; break;    // raw fit
				case 1: order = 1; crossorder = 2; break;    // better fit
				// Using higher orders result in better fits for the points found,
				// but leads to bad results if extrapolating to other areas. 
			}
			// Test if we have enough points to fit
			k = getNumParms3D(order, crossorder);
			if (nident < k) {
				printf("Found only %d grid points, %d needed for the next step.\n"
						"Please try more and/or different ori-points\n", nident, k);
				free(fixori);
				free(pixori);
				return FALSE;
			}
			degree[0] = order;
			degree[1] = crossorder;
			Fit3D_to_2D(&fit, degree, fixori, pixori, nident);

			if (examine3 == 50) {    // Debugging: save current results
				if (createfolder("debug\\sortgrid")) {
					sprintf(filename, "debug\\sortgrid\\coord%d_%d.txt", ic+1, loop+1);
					fp1 = fopen(filename, "w");
					for (i=0; i<nident; i++)
						fprintf(fp1, "%3d  %3d %8.3f %8.3f  %8.3f    %3d %8.3f  %8.3f\n",
							i, fixori[i].pnr, fixori[i].x, fixori[i].y, fixori[i].z, pixori[i].pnr, pixori[i].x, pixori[i].y);
					fclose(fp1);

					sprintf(filename, "debug\\sortgrid\\loop%d_%d.txt", ic+1, loop+1);
					saveparms_3D_to_pix(filename, "w", &fit);
				}
			}

			find_grid_points(fixori, pixori, &nident, Rsort, fit, fix, nfix, pix[ic], num[ic]);
                
			if (examine3 == 50) {    // Debugging: Save the points found
				sprintf(filename, "debug\\sortgrid\\ident%d_%d.txt", ic+1, loop+1);
				fp1 = fopen(filename, "w");
				for (i=0; i<nident; i++)
					fprintf(fp1, "%3d  %3d %8.3f %8.3f  %8.3f    %3d %8.3f  %8.3f\n",
						i, fixori[i].pnr, fixori[i].x, fixori[i].y, fixori[i].z, pixori[i].pnr, pixori[i].x, pixori[i].y);
				fclose(fp1);
			}
		}
		// copy and re-initialize pixel data before sorting
		orgpix = (target*) malloc(num[ic] * sizeof(target));
		for (i=0; i<num[ic]; i++) {
			target *p = &pix[ic][i];
			orgpix[i] = *p;
			p->pnr = (int)(p->x = p->y = -999);
			p->n = p->nx = p->ny = p->sumg = 0;
		}

		// adapt # of detected points
		num[ic] = nfix;
		for (i=0; i<nfix; i++)
			for (j=0; j<nident; j++)
				if (fix[i].pnr == fixori[j].pnr) {
					k = pixori[j].pnr;
					pix[ic][i] = orgpix[k];
					pix[ic][i].pnr = fix[i].pnr;
				}
		free(orgpix);

		// -- Save the sorting results --
		// These files will be read during the calibration procedure
		save_sorting_results(img_name[ic], fix, pix[ic], nfix);

		// -- draw the grid points found --
		for (i=0; i<nfix; i++) {
			if (pix[ic][i].pnr <= 0)
				continue;
			ix = (int)pix[ic][i].x;
			iy = (int)pix[ic][i].y;
			drawcross(interp, ix, iy, cr_sz, ic, "white");
			draw_pnr(interp, ix, iy, fix[i].pnr, ic, "white");
			draw_fiterror(interp, pix[ic][i].x, pix[ic][i].y, ic, fit, fix[i]);
		}
	}
	free(fixori);
	free(pixori);
	return TRUE;
}


int sortgrid_poly_mult(Tcl_Interp* interp, int iplane)
{
    coord_2d *pixori, *fixori;
	coord_3d *fix;
    target   *orgpix;
    double   r, rsearch, rstep, xp, yp, xfx, yfx, dmin, Rsort;
    int      i,j,k, nprev, icrd, ic, ip, nident, ix, iy;
    FitInfo  fit;
	sort_point *Pfix;
	int      rv, ncam, nplanes, nori, ori[4][MaxOri], deg3d2d[2], deg2depi[2];
	int      degree[2];

	printf("\nSorting grid points for plane %d\n", iplane+1);

	/* read the names of the calibration images */
	rv  = parse_mult_calori_par("parameters/mult_calori.par", iplane,
			&ncam, &nplanes, img_name, fixp_name, deg3d2d, deg2depi);
	rv &= parse_multi_manori_par("parameters/mult_manori.par", iplane,
			ori, &ncam, &nori);
	if (!rv) return FALSE;

	/* allocate fix[] and read the 3D-coordinates of the calibration target (plane) */
	if (! read_target_crds(fixp_name, &fix, &nfix))
		return FALSE;

    /* reading the search radius */
    Rsort = get_sortradius("parameters/sortgrid.par", 10.0);
    printf ("Sortgrid search radius: %.1f pixels\n", Rsort);

	for (ic=0; ic<n_img; ic++) {
        printf("Camera: %d\n", ic+1);

        // construct and fill arrays used in the fitting process
        fixori  = (coord_2d*) malloc(nfix * sizeof(coord_2d));    // 2d-coordinates in the current plane
        pixori  = (coord_2d*) malloc(nfix * sizeof(coord_2d));    // pixel coordinates

		/* Use the four clicked points for the first raw fit */  
		// fill the fixori and pixori arrays with the coordinates of the four points
        for (j=0; j<4; j++) {						// assuming 4 points are clicked!
            for (k=0; k<nfix; k++)
                if (fix[k].pnr == ori[ic][j]) {		// ie k is a clicked point
					fixori[j].x = fix[k].x;      fixori[j].y = fix[k].y; fixori[j].pnr = fix[k].pnr;
                    pixori[j].x = pix0[ic][j].x; pixori[j].y = pix0[ic][j].y;
                }
        }

        // 1: fit the fixed2d points against pixel coordinates
		//    using a0 + a1*x + a2*y + a3*x*y polynomials
		degree[0] = 1; degree[1] = 2;
        Fit2D_to_2D(&fit, degree, fixori, pixori, 4);

        // 2: start the search at the first clicked point.
		//    first sort all the 3D-crd by distance to this starting point.
		Pfix = (sort_point*) malloc(nfix *sizeof(sort_point));   // container for the sorted list

        for (j=0; j<nfix; j++) {
            Pfix[j].x = fix[j].x; Pfix[j].y = fix[j].y; Pfix[j].pnr = fix[j].pnr;
            Pfix[j].r = distance(Pfix[j].x, Pfix[j].y, fixori[0].x, fixori[0].y);
        }
		// sortdistance(Pfix, nfix, fixori[0]);
        qsort((void*)Pfix, nfix, sizeof(sort_point), cmp_dist);

        // 3: search the area using annuli, centered at the starting point, 
        //    increazing in size. After each annulus has been searched, new points found
		//	  are used to improve the fit before the next search is started.
		//    For this, a0 + a1*x + a2*x^2 + a3*y + a4*y^2 + a5*x*y polynomials are used.

		rstep   = Pfix[nfix-1].r/10;	// step size annulus radius = max.distance/10
        rsearch = rstep;				// the initial outer radius

        icrd = 0;
        nident = 0;         // the manori point will also be detected and added to the list
        nprev = nident;
        while (icrd < nfix) {
            while (icrd < nfix && Pfix[icrd].r <= rsearch) {
                xfx = Pfix[icrd].x;
                yfx = Pfix[icrd].y;
                xp  = fitdata_Map2D_x(xfx, yfx, &fit);
                yp  = fitdata_Map2D_y(xfx, yfx, &fit);
                
                dmin = 1e6;
                ip   = -999;
                for (j=0; j<num[ic]; j++) {
                    r = distance(xp,yp,  pix[ic][j].x, pix[ic][j].y);
                    if (r<Rsort && r<dmin) {
                        ip = j;
                        dmin = r;
                    }
                }
                if (ip > 0) {   // add the point to the collection to fit
                    fixori[nident].pnr  = Pfix[icrd].pnr;
                    fixori[nident].x    = Pfix[icrd].x;
                    fixori[nident].y    = Pfix[icrd].y;
                    pixori[nident].pnr  = ip;
                    pixori[nident].x    = pix[ic][ip].x;
                    pixori[nident].y    = pix[ic][ip].y;
                    nident++;
                }
                icrd = icrd+1;
            }
            if (nident>6 && nident>nprev) {		// calculate more accurate fit
				freefitparms(&fit);				// deallocate old parameters first
				degree[0] = 2; degree[1] = 2;
                Fit2D_to_2D(&fit, degree, fixori, pixori, nident);
                nprev = nident;
            }
            rsearch += rstep;
        }
        draw_fitted_fixpoints2d(interp, ic, fit, Pfix, nfix);

        // not identified targets are cleared from the list pix
        // First, the pixel data in pix is copied to orgpix, then pix is cleared
        // and filled agiain with the targets identified as grid points
        orgpix = (target*) malloc(num[ic] * sizeof(target));
        for (i=0; i<num[ic]; i++)
            orgpix[i] = pix[ic][i];

        // adapt # of detected points
        num[ic] = nfix;
        for (i=0; i<nfix; i++) {
			pix[ic][i].pnr =  (int)(pix[ic][i].x = pix[ic][i].y = -999);
            for (j=0; j<nident; j++)
                if (fixori[j].pnr == fix[i].pnr) {
                    pix[ic][i] = orgpix[ pixori[j].pnr ];
                    pix[ic][i].pnr = fix[i].pnr;
                    break;
                }
        }

		// -- Save the sorting results --
		// These files will be read during the calibration procedure
		save_sorting_results(img_name[ic], fix, pix[ic], nfix);

		// -- draw the grid points found --
        for (i=0; i<nfix; i++) {
            if (pix[ic][i].pnr < 0)
                continue;
            ix = (int)pix[ic][i].x;
            iy = (int)pix[ic][i].y;
            drawcross(interp, ix, iy, cr_sz, ic, "white");
            draw_pnr(interp, ix, iy, fix[i].pnr, ic, "white");
        }

		freefitparms(&fit);		// cleaning up
        free(orgpix); 
        free(Pfix);
        free(fixori);
        free(pixori);
    }
	free(fix);
    return TRUE;
}


int save_sorting_results(char *img_name, coord_3d *fix, target *pix, int nfix)
{
	FILE *fp;
	int i;
	char filename[256];

 	// 3D-coordinates
	sprintf(filename, "%s.fix", img_name);
    fp = fopen(filename, "w");
    for (i=0; i<nfix; i++)
        fprintf(fp, "%3d %10.4f %10.4f %10.4f\n", i+1, fix[i].x, fix[i].y, fix[i].z);
    fclose(fp);

    // image coordinates (pixels)
    sprintf(filename, "%s.crd", img_name);
    fp = fopen(filename, "w");
    for (i=0; i<nfix; i++)
		fprintf(fp, "%3d  %9.5f  %9.5f\n", i+1, pix[i].x, pix[i].y);
    fclose(fp);
	return TRUE;
}


int sortgrid_file_poly_body(Tcl_Interp* interp)
{
	int    ic, npix, i, k, dummy, *detection_pnr, ix, iy;
	target *old;
	char   file_sort[256];

	for (ic=0; ic<n_img; ic++) {
		// -- read the target indices from the sortgrid file
		sprintf (file_sort, "for_sortgrid.%1d", ic+1);
		fpp = fopen_rp (file_sort);
		if (!fpp) return FALSE;

		npix = num[ic];
		old           = (target*) malloc(npix * sizeof(target));
		detection_pnr = (int*)    malloc(nfix * sizeof(int));
		k = 0;
		while (k<nfix && fscanf (fpp, "%d %lf %lf %lf %d", &fix[k].pnr,
					&fix[k].x, &fix[k].y, &fix[k].z, &dummy) == 5) {
			detection_pnr[k++] = dummy;
		}
		fclose (fpp);

		// copy the targets list pix[] to old[]
		for (i=0; i<npix; i++)
			old[i] = pix[ic][i];

		// reset pix[]
		for (i=0; i<nfix; i++)
			reset_target(&pix[ic][i]);

		// fill pix[] with the information from the file
		for (i=0; i<nfix; i++) {
			if (detection_pnr[i] !=-999) {
				pix[ic][i] = old[detection_pnr[i]];
				pix[ic][i].pnr = fix[i].pnr;
			}
		}
		free (old);
		free (detection_pnr);

		// -- Save the sorting results --
		// These files will be read during the calibration procedure
		save_sorting_results(img_name[ic], fix, pix[ic], nfix);

		for (i=0; i<nfix; i++) {
			if (pix[ic][i].pnr <= 0)
				continue;
			ix = (int)pix[ic][i].x;
			iy = (int)pix[ic][i].y;
			drawcross(interp, ix, iy, cr_sz, ic, "white");
			draw_pnr(interp, ix, iy, fix[i].pnr, ic, "white");
		}
	}
	return TRUE;
}


int calibration_poly_body(Tcl_Interp* interp)
{
	// After the sortgrid function has been carried out successfully, 
    // this function will calculate the fit parameters 3D->2D and 2D->lines
	int      rv, i, ic, n, nzlev, ok, nfix;
	int      degree3d2d[2], degree2depi[2], order, crossorder;
    double   *Zlev;
	char     filename[256], imgname[4][256], fixp_name[256];
    coord_3d *fixpnts;
    coord_2d *pixpnts;
    FitInfo  fit, *planesfit, linefit_x, linefit_y;
    Rect rect;

    puts ("\nCalibration with grid points");
	rv = parse_body_calori_par("parameters/body_calori.par", imgname, fixp_name, degree3d2d, degree2depi);
	if (!rv) return FALSE;

	nzlev = 10;
	planesfit = (FitInfo*) malloc(nzlev * sizeof(FitInfo));
	Zlev      =  (double*) malloc(nzlev * sizeof(double));

	ok = 1; 
    for (ic=0; ic<n_img; ic++) {
		printf("\nProcessing camera %d\n", ic+1);

		/* --- collect the necessary data for this camera --- */
		printf("- Reading 3D and image coordinates.\n");
		sprintf(filename, "%s.fix", imgname[ic]);
		ok &= read_3Dcoordinates(filename, &fixpnts, &nfix);
		sprintf(filename, "%s.crd", imgname[ic]);
		ok &= read_2Dcoordinates(filename, &pixpnts, &nfix);
		if (!ok) break;

		// 1. Find the fit: 3D-points to pixel coordinates
		printf("- Fitting 3D-points and pixel positions.\n");
        // first remove points not detected during the sorting process
        n = 0;
        for (i=0; i<nfix; i++) {
            if (pix[ic][i].x >= 0) {	// = -999 if not found
                fixpnts[n] = fixpnts[i];
                pixpnts[n] = pixpnts[i]; 
                n++;
            }
        }
		nfix = n;

		// fitting the 3D-points and pixel coordinates
        order      = degree3d2d[0];
        crossorder = degree3d2d[1];
        Fit3D_to_2D_body(&fit, fixpnts, pixpnts, n);
		saveparms_3D_to_pix(img_cal[ic], "w", &fit);

		// 2.   Find the fit for pixels to pixel epilpolar lines
		printf("- Fitting pixel coordinates and epipolar lines.\n");
		// building a new data set for the whole volume
		// temp. solution
		// Get the fit of pixels to 3D-points in a series of z-planes.
		GetFitPixels_Planes(planesfit, Zlev, nzlev, &rect, &fit);

		// With this result, get the fit between pixels coordinates
		// and the line equation of the corresponding epilpolar lines.
		GetFit_Pixels_Lines(&linefit_x, &linefit_y, planesfit, Zlev, nzlev, &rect);    // fixed order/cross 3 3
		saveparms_2D_to_rays(img_cal[ic], "a", &linefit_x, &linefit_y);

		// cleaning up
		free(fixpnts);
		free(pixpnts);
		freefitparms(&fit);
		freefitparms(&linefit_x);
		freefitparms(&linefit_y);
		for (i=0; i<nzlev; i++)
			freefitparms(&planesfit[i]);
		puts("");
	}
    free(planesfit);
    free(Zlev);

    //sprintf(buf, "%f %f %f %f", fit[0].result, fit[1].result, fit[2].result, fit[3].result);
    //Tcl_Eval(interp, ".text delete 3");
    //Tcl_Eval(interp, ".text delete 1");
    //Tcl_Eval(interp, ".text insert 1 \"Calibration 3d -> 3d and 2d -> epipolar lines\"");
    //Tcl_Eval(interp, ".text delete 2");
    //Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each image -> \"");
    //Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
    //Tcl_Eval(interp, ".text insert 3 $tbuf");

	return TRUE;
}

int calibration_poly_mult(Tcl_Interp* interp)
{
	// If, for all the palnes, the sortgrid function has been carried out 
	// successfully, this function will calculate the fit parameters 3D->2D and 2D->lines

    int ic, iplane, ncam, nplanes, i, ok;
	char *imgname, filename[256];
    int rv, k, ii, np, *cnt;
    double *Zlev;
    coord_3d *fixpnts, **fixed;
    coord_2d *pixpnts, **pixcrd, *planepix, *planefix;
    FitInfo fit, *planesfit, linefit_x, linefit_y;
    Rect rect;
	int degree3d2d[2], degree2depi[2];

    puts ("Calibration with grid points");

	rv = parse_mult_calori_par("parameters/mult_calori.par", 0, 
			&ncam, &nplanes, img_name, fixp_name, degree3d2d, degree2depi);
	if (!rv) return FALSE;

    cnt    = (int*) malloc(nplanes * sizeof(int));
	fixed  = (coord_3d**) malloc(nplanes * sizeof(coord_3d*));
	pixcrd = (coord_2d**) malloc(nplanes * sizeof(coord_2d*));
	ok = 1;
	for (ic=0; ic<ncam; ic++) {
		/* --- collect the necessary data for this camera --- */
		printf("\nProcessing camera %d\n", ic+1);
		printf("- Reading 3D and image coordinates.\n");
		for (iplane=0; ok && iplane<nplanes; iplane++) {
			parse_mult_calori_par("parameters/mult_calori.par", iplane, 
					&ncam, &nplanes, img_name, fixp_name, degree3d2d, degree2depi);

			imgname = img_name[ic]; //val[3 + iplane*5];
			sprintf(filename, "%s.fix", imgname);
			ok &= read_3Dcoordinates(filename, &fixed[iplane], &cnt[iplane]);
			sprintf(filename, "%s.crd", imgname);
			ok &= read_2Dcoordinates(filename, &pixcrd[iplane], &cnt[iplane]);
		}
		if (!ok) break;

		// 1.   Find the fit for 3D-points to pixel coordinates
		printf("- Fitting grid points and pixel positions.\n");
		// Put all the valid grid points of all planes in one data set.
		nfix = 0;
		for (ii=0; ii<nplanes; ii++) nfix += cnt[ii];

		fixpnts = (coord_3d*) malloc(nfix * sizeof(coord_3d));
		pixpnts = (coord_2d*) malloc(nfix * sizeof(coord_2d));
		np = 0;
		for (iplane=0; iplane<nplanes; iplane++) {
			for (i=0; i<cnt[iplane]; i++) {
				if (pixcrd[iplane][i].x >= 0) {			// = -999 if not found
					fixpnts[np] = fixed[iplane][i];
					pixpnts[np] = pixcrd[iplane][i]; 
					np++;
				}
			}
		}
		// fitting the 3D-points to pixel coordinates
		Fit3D_to_2D(&fit, degree3d2d, fixpnts, pixpnts, np);
		saveparms_3D_to_pix(img_cal[ic], "w", &fit);

		// cleaning up
		freefitparms(&fit);
		free(fixpnts);
		free(pixpnts);

		// Fitting pixel coordinates and lines in 3D space
		// 1: Using the calibration grid data, get the coordinates of plane points 
		//    as function of pixel coordinates.
		// 2: create an artificial grid in pixel coordinates
		// 3: calculate for each grid point the 3D-point in each plane to find the
		//    intersection points of the path of light illuminating this grid point.
		// 4: by linear regression, find for each grid point a line going through these points
		//    find polynomials for the next functions to describe such a light ray
		//    Xo = fx(xp,yp): x- and y-coordinate of a point in the Z=0 plane
		//    Yo = fy(xp,yp)
		//    Xr = gx(xp,yp): direction coefficients of the line through (Xo,Yo,0)
		//    Yr = gy(xp,yp)

		printf("- Fitting pixel coordinates and epipolar lines.\n");

		// step 1:
		planesfit = (FitInfo*) malloc(nplanes * sizeof(FitInfo));
		Zlev      = (double*)  malloc(nplanes * sizeof(double));
		for (iplane=0; iplane<nplanes; iplane++) {
			// building the data set to fit
			planefix = (coord_2d*) malloc(cnt[iplane] * sizeof(coord_2d));
			planepix = (coord_2d*) malloc(cnt[iplane] * sizeof(coord_2d));
			k = 0;
			for (i=0; i<cnt[iplane]; i++) {
				if (pix[iplane][i].pnr < 0)
					continue;
				planefix[k].pnr = fixed[iplane][i].pnr;
				planefix[k].x   = fixed[iplane][i].x;
				planefix[k].y   = fixed[iplane][i].y;
				planepix[k++]   = pixcrd[iplane][i];
			}
			boundingrect(&rect, planepix, k, iplane);  // iplane=0 set rect, iplane>0 adjust rect
			Zlev[iplane] = fixed[iplane][0].z;

			Fit2D_to_2D(&planesfit[iplane], degree2depi, planepix, planefix, k);

			free (planefix);
			free (planepix);
		}
		for (iplane=0; iplane<nplanes; iplane++) {
			free(fixed[iplane]);
			free(pixcrd[iplane]);
		}

		// The steps 2, 3 and 4 are carried out by the next function
		// The function can be tested with the examine3 flag 64
		GetFit_Pixels_Lines(&linefit_x, &linefit_y, planesfit, Zlev, nplanes, &rect);

		// Append the resulting coefficients to the calibration file
		saveparms_2D_to_rays(img_cal[ic], "a", &linefit_x, &linefit_y);
		
		// cleaning up
		freefitparms(&linefit_x);
		freefitparms(&linefit_y);
		for (iplane=0; iplane<nplanes; iplane++)
			freefitparms(&planesfit[iplane]);
		free (planesfit);
		free (Zlev);

		//sprintf(buf, "%f %f %f %f", fit[0].result, fit[1].result, fit[2].result, fit[3].result);
		//Tcl_Eval(interp, ".text delete 3");
		//Tcl_Eval(interp, ".text delete 1");
		//Tcl_Eval(interp, ".text insert 1 \"Calibration 3d -> 3d and 2d -> epipolar lines\"");
		//Tcl_Eval(interp, ".text delete 2");
		//Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each image -> \"");
		//Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		//Tcl_Eval(interp, ".text insert 3 $tbuf");
 	}
    free(cnt);
	free(fixed);
	free(pixcrd);
	return ok;
}

int search_tracks_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
    int  rv, loop;
    char cmdbuf[1024], destnm[256];
    int  ac = 2;
    char *av[2];

    loop = 1;
    if (!createfolder("../data"))
        rv = TCL_ERROR;
    else {
        puts("Collecting links started");
        _flushall();
        av[0] = "search tracks";
        sprintf(destnm, "../data/3Dlinks_%d.txt", loop);
        sprintf(cmdbuf, "%s %s %d %d %d %s",
            "D:/vs2008/jemil/searchtracks/1.2/SearchTracksCon/Release/SearchTracksCon.exe",
            "res/ptv_is.%d", seq_first, seq_last, 2, destnm);
        rv = system(cmdbuf);
    }

    if (rv != TCL_OK) {
        sprintf (buf,"The %s command failed, process aborted", av[0]);
        puts(buf);
        Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
        Tcl_Eval(interp, ".text delete 2");
        Tcl_Eval(interp, ".text insert 2 $tbuf");
        return TCL_ERROR;
    }
    return TCL_OK;
}


BOOL copyparmsfiles(char* savefolder, int loop)        // copy the current calibration files to a save place
{
    char pathname[256], destfolder[256];
    int  ic;
    FILE *fp1, *fp2;
    BOOL ok;

    sprintf(destfolder, "%s/cycle%d", savefolder, loop);
    ok = createfolder(destfolder);

    for (ic=0; ok && ic<n_img; ic++) {
        // sprintf(pathname, "%s/%s", destfolder, getfilename(fnm_3dfit[ic])); 
        sprintf(pathname, "%s/%s", destfolder, getfilename(img_cal[ic])); 
        fp1 = fopen(img_cal[ic], "r");
        fp2 = fopen(pathname, "w");
        ok  = fp1!=NULL && fp2!=NULL;
        while (ok && fgets(buf, 256, fp1))
            fputs(buf, fp2);
        fclose(fp1);
        fclose(fp2);
    }
    return ok;
}
