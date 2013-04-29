// tstoutput.cpp : Defines the entry point for the console application.
//

#include "../ptv.h"

extern char safety[4][128];
extern char safety_addpar[4][128];

char* getbasename(char* folder, char* basename);
BOOL  read_ptv_is(char* pathname, ptv_is **plist, int *pcount);
BOOL  parse_ptv_is(char* pathname, ptv_is **plist, int *pcount);
BOOL  read_rt_is(char* pathname, rt_is **plist, int *pcount);
void  RevFitFunction3D_newest(double index, double *iodata, int iosize);
void  FitIncOrder_XYZ_uv(FitInfo* pfi, int order, int crossorder, FitData *fdata, int ndata);

void FreeFitInfoParms(FitInfo* pfit)
{
    free(pfit->parms1);
    free(pfit->parms2);
}

int poly_shaking_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
    char    pathname[256];
    FitInfo fit, linefit_x, linefit_y, fitZplane[8];
    double  x1,y1,z1, x2,y2,z2, z, dz, Zlev[8];
    Rect    rect = {0,0,0,0};
    int     ok, i, iz, ic, nstep, nx=11, ny=11, nz=8, cnt, trgcnt, a[4];
	int		seq_first, seq_last, max_shake_points, max_shake_frames, filenum;
	int		frameCount, currentFrame, step_shake, k, ipnt;
	int		degree32[2] = {3,3}, degree23[2] = {3,3};
    ptv_is  *ptv;
    rt_is   *rtis;
    target	*targets;
	coord_2d *targ, planepnt[11*11], pixelpnt[11*11];
	coord_3d *fix;
	critparameters critpar;
	BOOL byall, islink, bythis;

	ok  = parse_polycriteria_par2("parameters/poly_criteria.par", &critpar);
	ok &= parse_sequence_par("parameters/sequence.par", seq_name, &seq_first, &seq_last);
	ok &= parse_shaking_par("parameters/shaking.par", &seq_first, &seq_last,
				&max_shake_points, &max_shake_frames);
	if (!ok) return FALSE;

	fix  = (coord_3d*) malloc(max_shake_points * sizeof(coord_3d));
	targ = (coord_2d*) malloc(max_shake_points * sizeof(coord_2d));

	nstep = seq_last - seq_first + 1;
	step_shake = (int)((double)nstep/max_shake_frames + 0.5);
	printf("\nShaking's frame step size is: %d\n", step_shake);

	for (ic=0; ok && ic<n_img; ic++) {
		ipnt         = 0;
		frameCount   = 0;
		currentFrame = 0;
		printf("Processing camera %d\n", ic+1);
		for (filenum=seq_first+2; ok && filenum<seq_last+1-2; filenum+=step_shake)
		{
			printf("Reading 3D-positions and target data, file number %d\n", filenum);
			sprintf(pathname, "res/ptv_is.%d", filenum);
			ok = parse_ptv_is(pathname, &ptv, &cnt);	// allocates and fills ptv

			sprintf(pathname, "res/rt_is.%d", filenum);
			ok &= parse_rt_is(pathname, &rtis, &cnt);	// allocates and fills rtis
			
			/* read targets of this camera */
			sprintf(pathname, "%s%d_targets", seq_name[ic], filenum);
            ok &= parse_targetfile(pathname, &targets, &trgcnt);

			k = 0;
			for (i=0; ok && i<cnt; i++) {
				for (k=0; k<4; k++) a[k] = rtis[i].tnr[k];
				
				byall  = a[0]>-1 && a[1]>-1 && a[2]>-1 && a[3]>-1;		// seen by all
				islink = ptv[i].next >-1 &&  ptv[i].prev >-1;
				bythis = a[ic]>-1;
				if ( islink && ipnt<max_shake_points && 
					 ((bythis && frameCount<max_shake_frames+1) || (byall && filenum==seq_first+2)) ) 
				{
					targ[ipnt].x   = targets[a[ic]].x;
					targ[ipnt].y   = targets[a[ic]].y;
					targ[ipnt].pnr = i; 
					fix[ipnt].x    = rtis[i].X;
					fix[ipnt].y    = rtis[i].Y;
					fix[ipnt].z    = rtis[i].Z;
					fix[ipnt].pnr  = i;
						   
					nfix = ++ipnt;
					if (currentFrame < filenum) {
						currentFrame = filenum; 
						frameCount++;
					}
				}
			}
			// cleaning up
			free (ptv);
			free (rtis);
            free (targets);
		}
		if (!ok) break;

		// end of loop through seq, collected all necessary data for this camera.
		printf("Using %d linked points of %d frame%c for camera %d\n", 
			nfix, frameCount, (frameCount==1) ? 0 : 's', ic+1);

		// -- fitting -- 
		if (nfix < 40) {
            printf("Number of valid points found: %d. NOT ENOUGH TO PROCEED\n", nfix);
			return FALSE;
        }

		// Fitting fixed points with target coordinates.
		printf("Fitting 3D fixed points with 2D target coordinates\n");
		Fit3D_to_2D(&fit, degree32, fix, targ, nfix);

		// Fitting target coordinates with epipolar lines
		// - Get the fit of target coordinates to 3D-points in a series of z-planes.
		// - With this result, get the fit between pixels coordinates
		//     and the line equation of the corresponding epilpolar lines.

		// building a new data set over the whole volume
		x1 = critpar.xmin; x2 = critpar.xmax;
		y1 = critpar.ymin; y2 = critpar.ymax;
		z1 = critpar.zmin; z2 = critpar.zmax;

		dz = (z2-z1)/(nz-1);
		rect.left   = imx; rect.right = 0;
		rect.bottom = imy; rect.top   = 0;
		for (iz=0; iz<nz; iz++) {
			z = z1 + iz*dz;
			CreateZplaneFitdata(pixelpnt, planepnt, &fit, z, x1,x2,nx, y1,y2,ny);
			ok = Fit2D_to_2D(&fitZplane[iz], degree23, pixelpnt, planepnt, nx*ny);
            
			Zlev[iz] = z;
			// Calculation of the bounding rect of the projection of the
			// measuring volume on the CCD
			for (i=0; i<nx+ny; i++) {
				if (rect.left   > pixelpnt[i].x) rect.left   = pixelpnt[i].x;
				if (rect.right  < pixelpnt[i].x) rect.right  = pixelpnt[i].x;
				if (rect.bottom > pixelpnt[i].y) rect.bottom = pixelpnt[i].y;
				if (rect.top    < pixelpnt[i].y) rect.top    = pixelpnt[i].y;
			}
		}
		GetFit_Pixels_Lines(&linefit_x, &linefit_y, fitZplane, Zlev, nz, &rect);

		if (ok) {
			printf("Saving the new fit parameters.\n");
			saveparms_3D_to_pix(img_cal[ic],  "w", &fit);
			saveparms_2D_to_rays(img_cal[ic], "a", &linefit_x, &linefit_y);
			
			// cleaning up
			FreeFitInfoParms(&fit);
			FreeFitInfoParms(&linefit_x);
			FreeFitInfoParms(&linefit_y);
			for (iz=0; iz<nz; iz++)
				FreeFitInfoParms(&fitZplane[iz]);
		}
		else
			printf("Shaking failed for this camera!\nThe new fit parameters are NOT saved.\n\n");
	}
	free (fix);
	free (targ);
    return TRUE;
}

#ifdef UNDER_CONSTRUCTION	// ad holten, 04-2013
int orientation_from_particles_discarding_bad_3d_points(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int    ok, i, j, frameCount, currentFrame, filenumber, step_shake, max_shake_points, max_shake_frames;
	int    i_img, a[4], prev, next, count_inner, n, m;
	double d_inner=0., av_dist=0., x, y;
	double X[4],Y[4],Z[4],aa[4],b[4],c[4],dist,X_pos,Y_pos,Z_pos,XX,YY,ZZ;
	double pX[6],pY[6],pZ[6];
	double stdX,stdY,stdZ,rmsX[10000],av_rmsX;
	int    count_inner=0,count_outer=0,pair_count=0,count_dist=0;
	char   filein[256], filein_ptv[256], filein_T[256];
	FILE   *FILEIN, *FILEIN_ptv, *FILEIN_T;
	int    good[10000];
	Line3D line[4];

	puts ("Orientation from particles, discarding bad 3d-points"); strcpy(buf, "");
	
	for (i=0; i<4; i++) {
		sprintf(safety[i], "safety_%d.ori", i);
		sprintf(safety_addpar[i], "safety_%d.addpar", i);
	}

	// ok  = parse_polycriteria_par2("parameters/poly_criteria.par", &critpar);
	ok  = parse_sequence_par("parameters/sequence.par", seq_name, &seq_first, &seq_last);
	ok &= parse_shaking_par("parameters/shaking.par", &seq_first, &seq_last,
				&max_shake_points, &max_shake_frames);
	if (!ok) return FALSE;

	i = 0;
	frameCount = 0;
	currentFrame = 0;
	step_shake = (int)((double)(seq_last-seq_first+1)/(double)max_shake_frames+0.5);
	printf("\nframe step size for each camera is %d\n", step_shake);
	for (filenumber=seq_first+2; filenumber<seq_last+1-2; filenumber=filenumber+step_shake) { //changed by Beat Feb 08

		/* read targets of each camera */
		// to only use quadruplets for shaking that can be linked

		sprintf (filein, "res/rt_is.%d", filenumber);
		FILEIN = fopen_rp(filein);				// replaced fopen, ad holten, 12-2012
		if (! FILEIN) return TCL_OK;

		sprintf (filein_ptv, "res/ptv_is.%d", filenumber);
		FILEIN_ptv = fopen_rp (filein_ptv);		// replaced fopen, ad holten, 12-2012
		if (! FILEIN_ptv) return TCL_OK;

		/* read targets of each camera */
		for (i_img=0; i_img<n_img; i_img++) {
			nt4[3][i]=0;
			compose_name_plus_nr_str (seq_name[i_img], "_targets",filenumber, filein_T);
				
			// replaced next lines, ad holten
			//	FILEIN_T= fopen (filein_T, "r");
			//	if (! FILEIN_T) printf("Can't open ascii file: %s\n", filein_T);
			FILEIN_T = fopen_rp (filein_T);		// replaced fopen, ad holten, 12-2012
			if (! FILEIN_T) return TCL_OK;

			fscanf (FILEIN_T, "%d\n", &nt4[3][i_img]);
			for (j=0; j<nt4[3][i_img]; j++) {
				fscanf (FILEIN_T, "%4d %lf %lf %d %d %d %d %d\n",
					&t4[3][i_img][j].pnr, &t4[3][i_img][j].x,
					&t4[3][i_img][j].y, &t4[3][i_img][j].n ,
					&t4[3][i_img][j].nx ,&t4[3][i_img][j].ny,
					&t4[3][i_img][j].sumg, &t4[3][i_img][j].tnr);
			}
			fclose (FILEIN_T);
		}
				
		fscanf(FILEIN,     "%*d\n");	// skip # of 3D points
		fscanf(FILEIN_ptv, "%*d\n");	// skip # of 3D points
		do {
			/* read dataset row by row, x,y,z and correspondences */
			a[0] = a[1] = a[2] = a[3] = -1;
			fscanf(FILEIN, "%*d %lf %lf %lf %d %d %d %d\n", 
				&fix[i].x, &fix[i].y, &fix[i].z, &a[0], &a[1], &a[2], &a[3]);
			for (j=3; j>=n_img; j--) a[j] = -1;

			fscanf(FILEIN_ptv, "%d %d %*lf %*lf %*lf\n", &prev, &next);

			// read pixel data according a0,a1,a2,a3 
			if ( (a[0]>-1 && a[1]>-1 && a[2]>-1 && a[3]>-1  
					&&  next>-1 && prev>-1 && filenumber==seq_first+2 ) )		// OR ALLE QUADRUPLETS
			{
				for (i_img=0; i_img<n_img; i_img++) {
					pix[i_img][i].x = t4[3][i_img][a[i_img]].x;
					pix[i_img][i].y = t4[3][i_img][a[i_img]].y;
					pix[i_img][i].pnr = i; 
					fix[i].pnr = i;
				}		 
				nfix = ++i;
				if (currentFrame<filenumber) {
					currentFrame=filenumber;
					frameCount++;
				}
			}		   
		} while (!feof(FILEIN));
		fclose(FILEIN);
 
	} //end of loop through seq
  
	//for (i=0; i<nfix ; i++) {
	//	for (i_img=0; i_img<n_img; i_img++) {
	//		pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
	//						imx,imy, pix_x, pix_y,
	//						&crd[i_img][i].x, &crd[i_img][i].y,
	//						chfield);
	//		crd[i_img][i].pnr = pix[i_img][i].pnr;
	//	}
	//}
    // Keep using pixel coordinates,  ad holten, 04-2013
    for (i_img=0; i_img<n_img; i_img++) {
        for (i=0; i<num[i_img]; i++) {
            crd[i_img][i].x   = pix[i_img][i].x;
            crd[i_img][i].y   = pix[i_img][i].y,
            crd[i_img][i].pnr = pix[i_img][i].pnr;
        }
    }
	for (i=0; i<nfix ; i++) {
		// here comes the actual check!!!!!
		for (j=0; j<4; j++) {
			if (crd[j][i].x>-999) {
				//  x = crd[j][i].x - I[j].xh;
				//  y = crd[j][i].y - I[j].yh;
				//  // correct_brown_affin (x, y, ap[0], &x, &y);
				//  ray_tracing_v2 (x,y, Ex[j], I[j], G[j], mmp, &X[j], &Y[j], &Z[j], &aa[j], &b[j], &c[j]);
				ray_tracing_poly(x, y, j, &line[j]);
			}
		}

		count_inner=0;
		X_pos=0.; Y_pos=0.; Z_pos=0.;
		for (n=0; n<n_img; n++) {
			for(m=n+1;m<n_img;m++){
				if(crd[n][i].x>-999 && crd[m][i].x>-999){
					mid_point_poly(line[n], line[m], &dist, &XX, &YY, &ZZ);
					// mid_point(X[n],Y[n],Z[n],aa[n],b[n],c[n],X[m],Y[m],Z[m],aa[m],b[m],c[m],&dist,&XX,&YY,&ZZ);
					pX[count_inner] = XX;
					pY[count_inner] = YY;
					pZ[count_inner] = ZZ;
					count_inner++;
					d_inner += dist;
					X_pos += XX;  Y_pos += YY;  Z_pos += ZZ;
				}
			}
		}
		d_inner/=(double)count_inner;
		X_pos/=(double)count_inner; Y_pos/=(double)count_inner; Z_pos/=(double)count_inner;
		stdX=0; stdY=0; stdZ=0;
		for (n=0; n<count_inner; n++) {
			stdX += pow(pX[n]-X_pos,2);
			stdY += pow(pY[n]-Y_pos,2);
			stdZ += pow(pZ[n]-Z_pos,2);
		}
		rmsX[i] = pow(stdX+stdY+stdZ,0.5);
		//end of new det_lsq
		//end of actual check
	}
	av_rmsX=0;
	for (i=0; i<nfix ; i++)
		av_rmsX += rmsX[i];

	av_rmsX = av_rmsX/(double)nfix;
	for (i=0; i<nfix ; i++) {
		if (rmsX[i]<2*av_rmsX)
			good[i] = 1;
		else
			good[i] = 0;
	}
	// end of check for each fix -point.
		
	for (i_img=0; i_img<n_img; i_img++) {
	   
		orient_v6 (interp, Ex[i_img], I[i_img], G[i_img], ap[i_img], mmp,
				nfix, fix, good, crd[i_img],
				&Ex[i_img], &I[i_img], &G[i_img], &ap[i_img], i_img);

		/* save orientation and additional parameters */
		//make safety copy of ori files

		fp1 = fopen(img_ori[i_img], "r" );
		if (fp1 != NULL) {
			fclose(fp1);
			read_ori (&sEx[i_img], &sI[i_img], &sG[i_img], img_ori[i_img]);
			fp1 = fopen (img_addpar0[i_img], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					&sap[i_img].k1,&sap[i_img].k2,&sap[i_img].k3,
					&sap[i_img].p1,&sap[i_img].p2,
					&sap[i_img].scx,&sap[i_img].she);
				fclose (fp1);} 
			else {
				printf("no addpar.raw\n");
				sap[i_img].k1=sap[i_img].k2=sap[i_img].k3=sap[i_img].p1=sap[i_img].p2=sap[i_img].she=0.0;
				sap[i_img].scx=1.0;
			}

			write_ori (sEx[i_img], sI[i_img], sG[i_img], safety[i_img]);
			fp1 = fopen (safety_addpar[i_img], "w");
			fprintf (fp1, "%f %f %f %f %f %f %f",
				sap[i_img].k1, sap[i_img].k2, sap[i_img].k3,
				sap[i_img].p1, sap[i_img].p2,
				sap[i_img].scx, sap[i_img].she);
			fclose (fp1);
		}
		else {
			write_ori (Ex[i_img], I[i_img], G[i_img], safety[i_img]);
			fp1 = fopen (safety_addpar[i_img], "w");
			fprintf (fp1, "%f %f %f %f %f %f %f",
				ap[i_img].k1, ap[i_img].k2, ap[i_img].k3,
				ap[i_img].p1, ap[i_img].p2,
				ap[i_img].scx, ap[i_img].she);
			fclose (fp1);
		}
		write_ori (Ex[i_img], I[i_img], G[i_img], img_ori[i_img]);
		fp1 = fopen (img_addpar[i_img], "w");
		fprintf (fp1, "%f %f %f %f %f %f %f",
				ap[i_img].k1, ap[i_img].k2, ap[i_img].k3,
				ap[i_img].p1, ap[i_img].p2,
				ap[i_img].scx, ap[i_img].she);
		fclose (fp1);
	} //end of loop through images

	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text delete 1");
	Tcl_Eval(interp, ".text insert 1 \"Orientation from particles \"");
	Tcl_Eval(interp, ".text delete 2");
	if (examine != 4)
		Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each image -> \"");
	if (examine == 4 && multi==0)
		Tcl_Eval(interp, ".text insert 2 \"resection data written to disk \"");
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text insert 3 $tbuf");

	return TCL_OK;
}
#endif


char* getbasename(char* folder, char* basename)
{
    static char buf[256];

    sprintf(buf, "%s/%s", folder, basename);
    return buf;
}

