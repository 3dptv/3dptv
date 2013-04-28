/****************************************************************************

Routine:			checkpoints.c

Author/Copyright:	Hans-Gerd Maas

Address:			Institute of Geodesy and Photogrammetry
					ETH - Hoenggerberg
					CH - 8093 Zurich

Creation Date:		1988
	
Description:		calibration quality control from check points
					on calibration plate

Routines contained:

****************************************************************************/
/*
Copyright (c) 1990-2011 ETH Zurich

See the file license.txt for copying permission.
*/

/* --- original version modfified for the use with polynomials --- */
/*      ad holten, 04-2013                                         */

#include "../ptv.h"


void checkpoint_poly (Tcl_Interp* interp)
{
	int 	  ic, ok, i, j, iplane, count1=0, useflag;
	int 	  intx1, inty1, intx2, inty2;
	double	  X, Y, Z, xp[4], yp[4];
	double	  sigmax1=0, sigmay1=0, sigmaz1=0;
	FILE	  *fp;
	coord_3d  *fix;
	coord_2d  *pix[4];
	int ncam, nplanes, deg3d2d[2], deg2depi[2];
	char	  img_name[4][256], fixp_name[256];
	
	puts ("check points");

	if (n_img < 1) {
		printf("No image data available\n");
		return;
	}

	// read the current calibration parameters -------------
	if (usingZplanes) {
		// --- Get the current plane index ---
		iplane = atoi( Tcl_GetVar(interp, "cp(plane_id)", TCL_GLOBAL_ONLY) );
		iplane--;	// index start from 1 in dialog

		ok =  parse_mult_calori_par("parameters/mult_calori.par", 
				iplane, &ncam, &nplanes, img_name, fixp_name, deg3d2d, deg2depi);
	}
	else
		ok = parse_body_calori_par("parameters/body_calori.par", 
				img_name, fixp_name, deg3d2d, deg2depi);
	if (!ok) return;


	// --- update the values of the fit parameters ----------
    if (fit3dpix != NULL) {
		destroy_POLYFIT_list(fit3dpix);
		free (fit3dpix);
	}  
	if (fitpixepi != NULL) {
		destroy_POLYFIT_list(fitpixepi);
		free (fitpixepi);
	}
    fit3dpix  = (POLYFIT*) malloc(sizeof(POLYFIT) * n_img);
    fitpixepi = (POLYFIT*) malloc(sizeof(POLYFIT) * n_img);
    for (i=0; i<n_img; i++)
        parsepolycalibfile(img_cal[i], &fit3dpix[i], &fitpixepi[i]);

	// --- collect the necessary fit data -------------------
	printf("- Reading 3D and image coordinates.\n");
	sprintf(filename, "%s.fix", img_name[0]);
	ok = read_3Dcoordinates(filename, &fix, &nfix);
	for (ic=0; ok && ic<n_img; ic++) {
		sprintf(filename, "%s.crd", img_name[ic]);
		ok &= read_2Dcoordinates(filename, &pix[ic], &nfix);
	}
	if (!ok) return;

	/* read, which points shall be used  (those not used for orientation) */
	fp = fopen_rp ("parameters/orient.par");			// replaced fopen_r(), ad holten, 12-2012
	if (!fp) return;
	fscanf (fp,"%d", &useflag);
	fclose (fp);  

	rmsX = 0;  rmsY = 0;  rmsZ = 0;  mean_sigma0 = 0;

	for (i=0; i<nfix; i++) {	
		/* use (all) (even) (odd) point numbers as check points  */
		switch (useflag) {
			case 0: break;
			case 1: if ((fix[i].pnr % 2) != 0) continue; break;
			case 2: if ((fix[i].pnr % 2) == 0) continue; break;
			case 3: if ((fix[i].pnr % 3) != 0) continue; break;
		}
	  
		/* safety check */
		ok = 1;
		for (ic=0; ok && ic<n_img; ic++)
			ok = pix[ic][i].pnr == fix[i].pnr;
		if (!ok) continue;

		// --- find the 3D-point from the pixel coordinates ---
		for (j=0; j<n_img; j++) {
			xp[j] = pix[j][i].x;
			yp[j] = pix[j][i].y;
		}
		det_lsq_3d_poly(xp, yp, 4, &X, &Y, &Z);
		   
		/*	compute deviations	*/
		sigmax1 += (X - fix[i].x) * (X - fix[i].x);
		sigmay1 += (Y - fix[i].y) * (Y - fix[i].y);
		sigmaz1 += (Z - fix[i].z) * (Z - fix[i].z);
		printf ("%3d  -   %7.3f    %7.3f   %7.3f   |   %6.3f   %6.3f  %6.3f\n",
			fix[i].pnr, X,Y,Z, X-fix[i].x, Y-fix[i].y, Z-fix[i].z);
	  
		/* draw residual vectors into img0 */
		ic    = 0;
		intx1 = (int) pix[0][i].x;
		inty1 = (int) pix[0][i].y;
		intx2 = (int)(intx1
		         + 5*(Map3D_x(fix[i].x, fix[i].y, fix[i].z, &fit3dpix[ic]) -
			          Map3D_x(X, Y, Z, &fit3dpix[ic])) );
		inty2 = (int)(inty1
		         + 5*(Map3D_y(fix[i].x, fix[i].y, fix[i].z, &fit3dpix[ic]) -
			          Map3D_y(X, Y, Z, &fit3dpix[ic])) );
		//intz1 = inty1;
		//intz2 = (int)(intz1
		//         + 5*(Map3D_x(fix[i].x, fix[i].y, fix[i].z, &fit3dpix[ic]) -
		//	          Map3D_y(X, Y, Z, &fit3dpix[ic])) );

		drawvector (interp,intx1, inty1, intx2, inty2, 1, 0, "yellow");
		//inty2 = (int)(inty1 + (fix[i].z - Z)*5*Ex[0].z0/I[0].cc);
		//drawvector (interp, intx1, inty1, intx1, inty2, 1, 0, "blue");

		count1++;
	}

	/*	compute a priori RMS  */
	rmsX = sqrt(rmsX/count1);
	rmsY = sqrt(rmsY/count1);
	rmsZ = sqrt(rmsZ/count1);
	mean_sigma0 = sqrt (mean_sigma0/count1);	
	printf ("RMS from %d checkpoints:\n", count1);
	printf ("a priori => sigma0 = %4.2f micron, RMS = %6.3f/%6.3f/%6.3f mm",
		mean_sigma0*1000, rmsX, rmsY, rmsZ);

	/*	compute a posteriori RMS from check points	*/
	puts ("\n accuracies from check points:\n");
	sigmax1 = sqrt (sigmax1/count1);
	sigmay1 = sqrt (sigmay1/count1);
	sigmaz1 = sqrt (sigmaz1/count1);
	printf ("a posteriori	 => 			   RMS = %6.3f/%6.3f/%6.3f mm\n",
		sigmax1, sigmay1, sigmaz1);

	sprintf (buf, "%d check points => mx=%6.3f mm, my=%6.3f mm, mz=%6.3f mm",
		count1, sigmax1, sigmay1, sigmaz1);
	puts (buf);

	free(fix);
	for (ic=0; ic<n_img; ic++) {
		free(pix[ic]);
	}
}
