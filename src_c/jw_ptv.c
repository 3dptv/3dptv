/****************************************************************************

Author/Copyright:	Hans-Gerd Maas / Jochen Willneff

Address:			Institute of Geodesy and Photogrammetry
					ETH - Hoenggerberg
					CH - 8093 Zurich

Creation Date:		took a longer time ...

Description:		target detection, correspondences and
					positioning with tclTk display
					-> 4 camera version

Routines contained: many ...

****************************************************************************/

/*
Copyright (c) 1990-2011 ETH Zurich

See the file license.txt for copying permission.
*/


#include "ptv.h"

#define nmax 20240

/*	global declarations for ptv  */
/*-------------------------------------------------------------------------*/

int		n_img;						/* no of images */
int		hp_flag=0;					/* flag for highpass */
int		allCam_flag=0;				/* flag for using all cams for points */
int		tiff_flag=0; 				/* flag for tiff header */
int		pair_flag=0; 				/*flag for accept pair */
int		chfield; 					/* flag for field mode */
int		nfix;						/* no. of control points */
int		num[4];						/* no. of targets per image */
int		numc[4]; 					/* no. of targets in current image */
int		nump[4]; 					/* no. of targets in previous image */
int		numn[4]; 					/* no. of targets in next image */
int		n_trac[4];					/* no. of tracks */
int		match=0; 					/* no. of matches */
int		match2=0;					/* no. of matches in 2nd pass */
int		nr[4][4];					/* point numbers for man. ori */
int		imx, imy, imgsize;			   /* image size */
int		zoom_x[4],zoom_y[4],zoom_f[4];  /* zoom parameters */
int		pp1=0, pp2=0, pp3=0, pp4=0,pp5=0;  /* for man. orientation */
int		seq_first, seq_last; 			  /* 1. and last img of seq */
int		max_shake_points, max_shake_frames, step_shake;
int		demo_nr; 					/* for demo purposes */
int		examine = 0; 				/* for more detailed output */
int		dump_for_rdb;				/* # of dumpfiles for rdb */
int		cr_sz;						/* size of crosses */
int		display; 					/* display flag */
int		corp, corc, corn;
int		m[4];
int		trackallocflag = 0;			/* checkflag if mega, c4, t4 already allocated */
int		mask;						/*checkmark for subtract mask*/

double	pix_x, pix_y;							/* pixel size */
double	ro;										/* 200/pi */
double	cn, cnx, cny, csumg, eps0, corrmin;		/* correspondences par */
double	rmsX, rmsY, rmsZ, mean_sigma0;			/* a priori rms */
double	X_lay[2], Zmin_lay[2], Zmax_lay[2];		/* illu. layer data */
double	db_scale;								/*dumbbell length, Beat Mai 2010*/	

FILE	*fp1, *fp2, *fp3, *fp4, *fpp;

char	img_name[4][256];			/* original image names */
char	img_lp_name[4][256];		/* lowpass image names */
char	img_hp_name[4][256];		/* highpass image names */
char	img_cal[4][128];			/* calibrayion image names */
char	img_ori[4][128];			/* image orientation data */
char	img_ori0[4][128];			/* orientation approx. values */
char	img_addpar[4][128]; 		/* image additional parameters */
char	safety[4][128];
char	safety_addpar[4][128];
char	img_addpar0[4][128];		/* ap approx. values */
char	seq_name[4][128];			/* sequence names */
char	img_mask_name[4][256];		/* mask image names*/
char	img_mask_path[256];
char	track_dir[128]; 			/* directory with dap track data */
char	fixp_name[128];
char	res_name[128];				/* result destination */
char	filename[128];				/* for general use */
char	buf[256], val[256]; 		/* buffer */
char	name[128];					//Beat Dez 08
double	xp, yp; 					//Beat Dez 08

unsigned char *img[4];				/* image data */
unsigned char *img_mask[4]; 		/* mask data */
unsigned char *img_new[4];			/* image data for reducing mask */
unsigned char *img0[4]; 			/* image data for filtering etc */
unsigned char *zoomimg; 			/* zoom image data */

Exterior Ex[4],sEx[4];				/* exterior orientation */
Interior I[4],sI[4];				/* interior orientation */
Glass	 G[4],sG[4];				/* glass orientation */
ap_52	 ap[4],sap[4];				/* add. parameters k1,k2,k3,p1,p2,scx,she */
mm_np	 mmp;						/* n-media parameters */
target	 pix[4][nmax];				/* target pixel data */
target	 pix0[4][12];				/* pixel data for man_ori points */

target	 *t4[4][4];
int 	 nt4[4][4];

coord_2d crd[4][nmax];				/* (distorted) metric coordinates */
coord_2d geo[4][nmax];				/* corrected metric coordinates */
coord_3d fix[20096];				/* testfield points coordinates */ //Beat changed it on 090325
n_tupel  con[nmax]; 				/* list of correspondences */

corres	 *c4[4];
trackparameters tpar;				/* tracking parameters */

mm_LUT	 mmLUT[4];					/* LUT for multimedia radial displacement */
coord_3d *p_c3d;
P		 *mega[4];

/***************************************************************************/

int init_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int  i;
	const char *valp;

	puts ("\n Multimedia Particle Positioning and Tracking Software \n");

	valp = Tcl_GetVar(interp, "examine",  TCL_GLOBAL_ONLY);
	examine = atoi (valp);

	ro = 200/M_PI;

	fpp = fopen ("parameters/pft_version.par", "r");
	if (! fpp) {
		fpp = fopen ("parameters/pft_version.par", "w");
		fprintf(fpp,"%d\n", 0);
		fclose(fpp);
	}

	/*	read from main parameter file  */
	fpp = fopen_rp ("parameters/ptv.par");			// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%d\n", &n_img);

	for (i=0; i<4; i++) {
		fscanf (fpp, "%s\n", img_name[i]);
		fscanf (fpp, "%s\n", img_cal[i]);
	}
	fscanf (fpp, "%d\n",  &hp_flag);
	fscanf (fpp, "%d\n",  &allCam_flag);
	fscanf (fpp, "%d\n",  &tiff_flag);
	fscanf (fpp, "%d\n",  &imx);
	fscanf (fpp, "%d\n",  &imy);
	fscanf (fpp, "%lf\n", &pix_x);
	fscanf (fpp, "%lf\n", &pix_y);
	fscanf (fpp, "%d\n",  &chfield);
	fscanf (fpp, "%lf\n", &mmp.n1);
	fscanf (fpp, "%lf\n", &mmp.n2[0]);
	fscanf (fpp, "%lf\n", &mmp.n3);
	fscanf (fpp, "%lf\n", &mmp.d[0]);
	fclose (fpp);

	/* read illuminated layer data */
	fpp = fopen_rp ("parameters/criteria.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%lf\n", &X_lay[0]);
	fscanf (fpp, "%lf\n", &Zmin_lay[0]);
	fscanf (fpp, "%lf\n", &Zmax_lay[0]);
	fscanf (fpp, "%lf\n", &X_lay[1]);
	fscanf (fpp, "%lf\n", &Zmin_lay[1]);
	fscanf (fpp, "%lf\n", &Zmax_lay[1]);
	fscanf (fpp, "%lf", &cnx);
	fscanf (fpp, "%lf", &cny);
	fscanf (fpp, "%lf", &cn);
	fscanf (fpp, "%lf", &csumg);
	fscanf (fpp, "%lf", &corrmin);
	fscanf (fpp, "%lf", &eps0);
	fclose (fpp);

	mmp.nlay = 1;

	/* read sequence parameters (needed for some demos) */

	fpp = fopen_rp ("parameters/sequence.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	for (i=0; i<4; i++) 
		fscanf (fpp, "%s\n", seq_name[i]);
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);

	/* initialize zoom parameters and image positions */
	for (i=0; i<n_img; i++) {
		num[i] = 0;
		zoom_x[i] = imx/2; zoom_y[i] = imy/2; zoom_f[i] = 1;
	}
	imgsize = imx*imy;

	/* allocate memory for images */
	for (i=0; i<n_img; i++) {
		img[i] = (unsigned char*) calloc(imgsize, 1);
		if (! img[i]) {
			printf ("calloc for img%d --> error\n", i);
			exit (1);
		}
	}

	for (i=0; i<n_img; i++) {
		img_mask[i] = (unsigned char*) calloc(imgsize, 1);
		if (! img_mask[i]) {
			printf ("calloc for img_mask%d --> error\n", i);
			exit (1);
		}
	}

	for (i=0; i<n_img; i++) {
		img0[i] = (unsigned char*) calloc(imgsize, 1);
		if (! img0[i]) {
			printf ("calloc for img0%d --> error\n", i);
			exit (1);
		}
	}

	for (i=0; i<n_img; i++) {
		img_new[i] = (unsigned char*) calloc(imgsize, 1);
		if (! img_new[i]) {
			printf ("calloc for img_new%d --> error\n", i);
			exit (1);
		}
	}

	zoomimg = (unsigned char*) calloc(imgsize, 1);
	if (! zoomimg) {
		printf ("calloc for zoomimg --> error\n");
		return TCL_ERROR;
	}

	parameter_panel_init(interp);
	cr_sz = atoi(Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY));

	display = 1;
	return TCL_OK;
}


int start_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int  i, k;

	/*	read from main parameter file  */
	fpp = fopen_rp ("parameters/ptv.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%d\n", &n_img);
	for (i=0; i<4; i++) {
		fscanf (fpp, "%s\n", img_name[i]);
		fscanf (fpp, "%s\n", img_cal[i]);
	}
	fscanf (fpp, "%d\n",  &hp_flag);
	fscanf (fpp, "%d\n",  &allCam_flag);  
	fscanf (fpp, "%d\n",  &tiff_flag);
	fscanf (fpp, "%d\n",  &imx);
	fscanf (fpp, "%d\n",  &imy);
	fscanf (fpp, "%lf\n", &pix_x);
	fscanf (fpp, "%lf\n", &pix_y);
	fscanf (fpp, "%d\n",  &chfield);
	fscanf (fpp, "%lf\n", &mmp.n1);
	fscanf (fpp, "%lf\n", &mmp.n2[0]);
	fscanf (fpp, "%lf\n", &mmp.n3);
	fscanf (fpp, "%lf\n", &mmp.d[0]);
	fclose (fpp);

	if (imgsize < imx*imy) {	// added, ad holten 12-2012
		printf("The allocated image buffers are to small for the\n"
			   "image size, defined in the calibration parameters dialog.\n."
			   "Please restart the program.\n");
		return TCL_ERROR;
	}

	/* read illuminated layer data */
	fpp = fopen_rp ("parameters/criteria.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%lf\n", &X_lay[0]);
	fscanf (fpp, "%lf\n", &Zmin_lay[0]);
	fscanf (fpp, "%lf\n", &Zmax_lay[0]);
	fscanf (fpp, "%lf\n", &X_lay[1]);
	fscanf (fpp, "%lf\n", &Zmin_lay[1]);
	fscanf (fpp, "%lf\n", &Zmax_lay[1]);
	fscanf (fpp, "%lf", &cnx);
	fscanf (fpp, "%lf", &cny);
	fscanf (fpp, "%lf", &cn);
	fscanf (fpp, "%lf", &csumg);
	fscanf (fpp, "%lf", &corrmin);
	fscanf (fpp, "%lf", &eps0);
	fclose (fpp);

	mmp.nlay = 1;

	/* read sequence parameters (needed for some demos) */

	fpp = fopen_rp ("parameters/sequence.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	for (i=0; i<4; i++)
		fscanf (fpp, "%s\n", seq_name[i]);
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);

	/*	create file names  */
	for (i=0; i<n_img; i++) {
		strcpy (img_lp_name[i], img_name[i]); strcat (img_lp_name[i], "_lp");
		strcpy (img_hp_name[i], img_name[i]); strcat (img_hp_name[i], "_hp");
		strcpy (img_ori[i], 	img_cal[i]);  strcat (img_ori[i], ".ori");
		strcpy (img_addpar[i],	img_cal[i]);  strcat (img_addpar[i],".addpar");
	}

	/*	read orientation and additional parameters	*/
	for (i=0; i<n_img; i++) {
		if (!read_ori (&Ex[i], &I[i], &G[i], img_ori[i]))
			return TCL_OK;						// added, ad holten, 12-2012
		rotation_matrix (Ex[i], Ex[i].dm);

		fp1 = fopen_rp (img_addpar[i]);			// replaced fopen_r, ad holten 12-2012
		if (!fp1) return TCL_OK;
		fscanf (fp1,"%lf %lf %lf %lf %lf %lf %lf",
			&ap[i].k1, &ap[i].k2, &ap[i].k3, &ap[i].p1, &ap[i].p2,
			&ap[i].scx, &ap[i].she);
		fclose (fp1);
	}

	/* read and display original images */
	for (i=0; i<n_img; i++) {
		/* reading */
		sprintf(val, "camcanvas %d", i+1);
		Tcl_Eval(interp, val);

		if (!read_image (interp, img_name[i], img[i]))
			return TCL_OK;							// added, ad holten, 12-2012
		sprintf(val, "newimage %d", i+1);

		Tcl_Eval(interp, val);
		sprintf(val, "keepori %d", i+1);
		Tcl_Eval(interp, val);
	}

	if (!trackallocflag) {
		for (i=0; i<4; i++) {
			mega[i] = (P*) calloc(sizeof(P),M);
			c4[i]	= (corres*) calloc(sizeof(corres),M);
			for (k=0; k<4; k++)
				t4[i][k] = (target*) calloc(sizeof (target),M);
		}
		trackallocflag = 1;
	}
	return TCL_OK;
}

int pre_processing_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int i_img, sup, i;

	Tk_PhotoHandle img_handle;
	Tk_PhotoImageBlock img_block;

	sprintf(val, "Filtering with Highpass");
	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 2");
	Tcl_Eval(interp, ".text insert 2 $tbuf");

	/* read support of unsharp mask */
	fpp = fopen ("parameters/unsharp_mask.par", "r");
	if (fpp == 0)
		sup = 12;
	else {
		fscanf (fpp, "%d\n", &sup);
		fclose (fpp);
	}

	//_____________________Matthias subtract mask__________________________
	/* Matthias JULI 08 read checkmark for masks and create mask names*/

	fpp = fopen_rp ("parameters/targ_rec.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;

	for (i=0; i<14; i++)
		fscanf (fpp, "%d", &mask);		/*checkmark for subtract mask */
	fscanf (fpp, "%s\n", img_mask_path);
	fclose (fpp);

	/*read mask names*/
	strcpy (img_mask_name[0], img_mask_path); strcat (img_mask_name[0], ".0");
	strcpy (img_mask_name[1], img_mask_path); strcat (img_mask_name[1], ".1");
	strcpy (img_mask_name[2], img_mask_path); strcat (img_mask_name[2], ".2");
	strcpy (img_mask_name[3], img_mask_path); strcat (img_mask_name[3], ".3");

	/* if the checkmark is set, read mask-image and subtract it from the filtered-original image.*/
	// - read mask
	// - highpass original image
	// - subtract mask from original image
	// - copy subtracted imgage on the original image
	if (mask == 1) {	// read mask image
		for (i_img=0; i_img<n_img; i_img++) {
			read_image (interp, img_mask_name[i_img], img_mask[i_img]); 				
			highpass (img_name[i_img], img[i_img], img[i_img], sup, 0, chfield, i_img); 
			subtract_mask (img[i_img], img_mask[i_img], img_new[i_img]);				
			copy_images (img_new[i_img], img[i_img]);									

			if (display) {
				img_handle = Tk_FindPhoto( interp, "temp");
				Tk_PhotoGetImage (img_handle, &img_block);
				tclimg2cimg (interp, img[i_img], &img_block);

				sprintf(val, "newimage %d", i_img+1);
				Tcl_GlobalEval(interp, val);
			}
		}
	}

	if (mask == 2) {	//Beat April 090402 was ==0
		for (i_img=0; i_img<n_img; i_img++) {
			//highpass original image
			highpass (img_name[i_img], img[i_img], img[i_img], sup, 0, chfield, i_img);

			if (display) {
				img_handle = Tk_FindPhoto( interp, "temp");
				Tk_PhotoGetImage (img_handle, &img_block);
				tclimg2cimg (interp, img[i_img], &img_block);

				sprintf(val, "newimage %d", i_img+1);
				Tcl_GlobalEval(interp, val);
			}
		}
	}

	/*------------------------------------------------------------*/
	sprintf(val, "...done");
	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text insert 3 $tbuf");

	return TCL_OK;
}


int detection_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
{
	int  i, i_img, j;
	int  xmin, pft_version=3;
	char val[256];
	char filename[256];
	FILE *FILEIN;

	Tk_PhotoHandle img_handle;
	Tk_PhotoImageBlock img_block;

	/* process info */
	sprintf(val, "Detection of Particles");
	Tcl_Eval(interp, ".text delete 2");
	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text insert 2 $tbuf");

	if (display) {
		for (i_img=0; i_img<n_img; i_img++) {
			img_handle = Tk_FindPhoto( interp, "temp");
			Tk_PhotoGetImage (img_handle, &img_block);
			tclimg2cimg (interp, img[i_img], &img_block);
			sprintf(val, "newimage %d", i_img+1);
			Tcl_Eval(interp, val);
		}
	}
	strcpy(val, "");

	xmin = 0;

	/*	read pft version  */
	fpp = fopen ("parameters/pft_version.par", "r");
	if (fpp) {
		fscanf (fpp, "%d\n", &pft_version);
		pft_version = pft_version+3;
		fclose (fpp);
	}
	else{
		fpp = fopen ("parameters/pft_version.par", "w");
		fprintf(fpp,"%d\n", 0);
		fclose(fpp);
	}

	/* reset zoom values */
	for (i_img=0; i_img<n_img; i_img++) {
		zoom_x[i_img] = imx/2; zoom_y[i_img] = imy/2;  zoom_f[i_img] = 1;
	}

	/* copy images because the target recognition will set greyvalues to 0 */
	for (i_img=0; i_img<n_img; i_img++)
		copy_images (img[i_img], img0[i_img]);

	/* target recognition */
	for (i_img=0; i_img<n_img; i_img++) {
		switch (pft_version)
		{
		case 3: 	/* pft with profile and distance check */
			/* newest version */
			xmin=0; /* vertical line restriction */
			num[i_img] = peak_fit_new (interp, img[i_img], "parameters/targ_rec.par",
							           xmin, imx, 1, imy, pix[i_img], i_img);
			break;

		case 0:    /* without peak fitting technique */
			simple_connectivity (interp, img[i_img], img0[i_img], "parameters/targ_rec.par",
				                       xmin, imx, 1, imy, pix[i_img], i_img, &num[i_img]);
			break;

		case 1:    /* with old (but fast) peak fitting technique */
			targ_rec (interp, img[i_img], img0[i_img], "parameters/targ_rec.par",
					  xmin, imx, 1, imy, pix[i_img], i_img, &num[i_img]);
			break;
	
		case 4: 	/* new option for external image processing routines */
			/* added by Alex, 19.04.10 */
			/* this works here only for the pre-processing stage, see img_name[i_img] is not from a sequence */

			sprintf (filename, "%s%s", img_name[i_img],"_targets");
			/* read targets of each camera */
			nt4[3][i_img]=0;

			// FILEIN= fopen (filename, "r");
			// if (! FILEIN) printf("Can't open ascii file: %s\n", filename);
			FILEIN = fopen_rp (filename);		// replaced fopen, ad holten, 12-2012
			if (! FILEIN) return TCL_OK;

			fscanf (FILEIN, "%d\n", &nt4[3][i_img]);
			for (j=0; j<nt4[3][i_img]; j++) {
				fscanf (FILEIN, "%4d %lf %lf %d %d %d %d %d\n",
					&pix[i_img][j].pnr, &pix[i_img][j].x,
					&pix[i_img][j].y, &pix[i_img][j].n ,
					&pix[i_img][j].nx ,&pix[i_img][j].ny,
					&pix[i_img][j].sumg, &pix[i_img][j].tnr);
			}
			fclose (FILEIN);
			num[i_img] = nt4[3][i_img];
		  
			if (display) {
				for (j=0; j<num[i_img]; j++)
					drawcross (interp, (int)pix[i_img][j].x, (int)pix[i_img][j].y, cr_sz, i_img, "blue");
			}
			break;
		}
		sprintf (buf,"%d: %d,  ", i_img+1, num[i_img]);
		strcat(val, buf);

		/* proper sort of targets in y-direction for later binary search */
		/* and for dimitris' tracking */
		quicksort_target_y (pix[i_img], num[i_img]);

		/* reorganize target numbers */
		for (i=0; i<num[i_img]; i++)
			pix[i_img][i].pnr = i;
	}

	sprintf (buf, "Number of detected particles per image");
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 2");
	Tcl_Eval(interp, ".text insert 2 $tbuf");

	Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text insert 3 $tbuf");

	printf("%s\n", val);
	return TCL_OK;
}



int correspondences_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int    i, i_img;
	double x,y;

	puts ("\nTransformation to metric coordinates\n");

	/* rearrange point numbers after manual deletion of points */
	for (i_img=0; i_img<n_img; i_img++)
		for (i=0; i<num[i_img]; i++)
			pix[i_img][i].pnr = i;
	/* transformations pixel coordinates -> metric coordinates */
	/* transformations metric coordinates -> corrected metric coordinates */
	for (i_img=0; i_img<n_img; i_img++) {
		for (i=0; i<num[i_img]; i++) {
			pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
							 imx,imy, pix_x, pix_y,
							 &crd[i_img][i].x, &crd[i_img][i].y, chfield);
			crd[i_img][i].pnr = pix[i_img][i].pnr;

			x = crd[i_img][i].x - I[i_img].xh;
			y = crd[i_img][i].y - I[i_img].yh;
			correct_brown_affin (x, y, ap[i_img], &geo[i_img][i].x, &geo[i_img][i].y);

			geo[i_img][i].pnr = crd[i_img][i].pnr;
		}
	}

	/* sort coordinates for binary search in correspondences_proc */
	for (i_img=0; i_img<n_img; i_img++)
		quicksort_coord2d_x (geo[i_img], num[i_img]);

	/* init multimedia radial displacement LUTs */
	/* ======================================== */

	if ( !mmp.lut && (mmp.n1 != 1 || mmp.n2[0] != 1 || mmp.n3 != 1)) {
		puts ("Init multimedia displacement LUTs");
		for (i_img=0; i_img<n_img; i_img++) init_mmLUT(i_img);
		mmp.lut = 1;
	}

	correspondences_4 (interp, argv);

	/* ------ save pixel coords for tracking  ------- */
	for (i_img=0; i_img<n_img; i_img++) {
		sprintf (filename, "%s_targets", img_name[i_img]);
		fp1 = fopen (filename, "w");

		fprintf(fp1,"%d\n", num[i_img]);
		for (i=0; i<num[i_img]; i++) {
			fprintf (fp1, "%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n", 
				pix[i_img][i].pnr, pix[i_img][i].x,
				pix[i_img][i].y, pix[i_img][i].n ,
				pix[i_img][i].nx ,pix[i_img][i].ny,
				pix[i_img][i].sumg, pix[i_img][i].tnr);
		}
		fclose (fp1);
	}
	return TCL_OK;
}


int determination_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int    i, j, n,dummy;
	int    p[4];
	double x[4], y[4], X,Y,Z;
	double Zlo = 1e20, Zhi = -1e20;
	int    dumbbell=0, i1, i2;
	double x1,y1,z1,x2,y2,z2;
	int a1[4],a2[4],checksum_1,checksum_2;

	puts ("Determinate");

	sprintf (buf, "Point positioning (mid_point in 3d)");
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 2");
	Tcl_Eval(interp, ".text insert 2 $tbuf");

	/* Beat Mai 2007 to set the variable examine for mulit-plane calibration */
	fp1 = fopen_rp ("parameters/examine.par");		// replaced fopen_r, ad holten 12-2012
	if (!fp1) return TCL_OK;
	fscanf (fp1,"%d\n", &dummy);
	fclose (fp1);
	if (dummy == 1)
		examine=4;
	else
		examine=0;

	fp1 = fopen (res_name, "w");
	if (! fp1) {
		sprintf(res_name,"res/dt_lsq");
		fp1 = fopen (res_name, "w");
	}
	if (! fp1) {
		printf ("cannot find dir: res,	data written to dt_lsq in same dir\n");
		sprintf (res_name, "dt_lsq");
		fp1 = fopen (res_name, "w");
	}
	/* create dump file for rdb */
	if (examine == 4) {
		/* create filename for dumped dataset */
		sprintf (res_name, "dump_for_rdb");
		printf ("dataset dumped into %s\n", res_name);
		fp2 = fopen (res_name, "w");

		/* write # of points to file */
		fprintf (fp2, "%d\n", match);
	}
	/* first line to be updated in res_name file */
	fprintf (fp1, "%4d\n", match);
	/* least squares determination for triplets */

	rmsX = 0; rmsY = 0; rmsZ = 0; mean_sigma0 = 0;

	for (i=0; i<match; i++) {
		for (j=0; j<4; j++) {
			if (con[i].p[j] >= 0)  p[j] = geo[j][con[i].p[j]].pnr;
			else				   p[j] = -1;
		}
		for (j=0, n=0; j<4; j++) {
			if (p[j] > -1) {
				x[j] = crd[j][p[j]].x;	  y[j] = crd[j][p[j]].y;
				n++;
			}
			else {
				x[j] = -1e10;	y[j] = -1e10;
				if (p[j] == -2) n = -100;
			}
		}

		// take only points which are matched in all images
		// or triplets/quadruplets which result from object model
		// e.g.: quad -> n=4; model triplet -> n=3; model pair -> n=2;
		// unrestricted triplet -> n<0; unrestricted pair -> n<0
		/* if (n_img > 2  &&  n < 3)	continue; */

		if ((n_img > 2 && num[0]>64 && num[1]>64 && num[2]>64 && num[3]>64)
			&&	n < 3) continue;

		/* hack due to problems with approx in det_lsq: */
		X = 0.0; Y = 0.0; Z = (Zmin_lay[0]+Zmax_lay[0])/2.0;
		for (j=0; j<n_img; j++) {
			X += Ex[j].x0;
			Y += Ex[j].y0;
		}
		X /= n_img; Y /= n_img;
	  
		// det_lsq_old (Ex, I, ap, mmp,
		//	   x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], &X, &Y, &Z);

		det_lsq_3d (Ex, I, G, ap, mmp,
			x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], &X, &Y, &Z);

		/* write a sequential point number,
		   sumg, if the point was used, and the 3D coordinates */

		// next code replaced by ad holten, 12-2012 
		// 		fprintf (fp1, "%4d", i+1);

		// 		/*	if (p[0] > -1) fprintf (fp1, "	%4d", pix[0][p[0]].sumg);
		//			else		   fprintf (fp1, "	%4d", -1);
		//			if (p[1] > -1) fprintf (fp1, "	%4d", pix[1][p[1]].sumg);
		//			else		   fprintf (fp1, "	%4d", -1);
		//			if (p[2] > -1) fprintf (fp1, "	%4d", pix[2][p[2]].sumg);
		//			else		   fprintf (fp1, "	%4d", -1);
		//			if (p[3] > -1) fprintf (fp1, "	%4d", pix[3][p[3]].sumg);
		//			else		   fprintf (fp1, "	%4d", -1);
		// 		*/

		// 		fprintf (fp1, " %9.3f %9.3f %9.3f", X, Y, Z);
		// 		if (p[0] > -1) fprintf (fp1, " %4d",   pix[0][p[0]].pnr);
		// 		else		   fprintf (fp1, " %4d",   -1);
		// 		if (p[1] > -1) fprintf (fp1, " %4d",   pix[1][p[1]].pnr);
		// 		else		   fprintf (fp1, " %4d",   -1);
		// 		if (p[2] > -1) fprintf (fp1, " %4d",   pix[2][p[2]].pnr);
		// 		else		   fprintf (fp1, " %4d",   -1);
		// 		if (p[3] > -1) fprintf (fp1, " %4d\n", pix[3][p[3]].pnr);
		// 		else		   fprintf (fp1, " %4d\n", -1);

		fprintf (fp1, "%4d %9.3f %9.3f %9.3f", i+1, X, Y, Z);
		for (j=0; j<4; j++) {
			if (p[j] > -1) fprintf (fp1, " %4d", pix[j][p[j]].pnr);
			else		   fprintf (fp1, " %4d", -1);
		}
		fprintf (fp1, "\n");

		/* write data as new points to dump for rdb */
		if (examine == 4) {
			fprintf (fp2, "%d %10.3f %10.3f %10.3f	 %d    ", i, X, Y, Z, 3);
			for (j=0; j<n_img; j++) {
				if (x[j] != -1e10)
					fprintf (fp2, "%4d %8.5f %8.5f	  ", i, x[j], y[j]);
				else
					fprintf (fp2, "%4d %8.5f %8.5f	  ", -999, x[j], y[j]);
			}
			fprintf (fp2, "\n");
			fclose (fp2);
		}
		if (Z < Zlo) Zlo = Z;
		if (Z > Zhi) Zhi = Z;
	}
	fclose (fp1);

	// Beat Mai 2010: now we should open the file db_is.* again, check
	//				  if it has exactly two points, rescale them, write them again and close the file.
	if ((argv[1]) > 0) {
		if (atoi(argv[1])==3){
			dumbbell=1;
			display=0;
		}
	}
	if (dumbbell==1) {
		fpp = fopen_rp ("parameters/dumbbell.par");			// replaced fopen, ad holten, 12-2012
		if (!fpp) return TCL_OK;
		fscanf (fpp, "%lf", &eps0);
		fscanf (fpp, "%lf", &db_scale);
		fclose (fpp);

		fpp = fopen_rp (res_name);							// replaced fopen, ad holten, 12-2012
		if (!fpp) return TCL_OK;
		fscanf (fpp, "%d\n", &match);
		if (match==2) {
			fscanf(fpp, "%d %lf %lf %lf %d %d %d %d\n",
				&i1, &x1, &y1, &z1, &a1[0], &a1[1], &a1[2], &a1[3]);
			fscanf(fpp, "%d %lf %lf %lf %d %d %d %d\n",
				&i2, &x2, &y2, &z2, &a2[0], &a2[1], &a2[2], &a2[3]);
			// // now adapt x,y,z
			// /* dist=pow(pow(x2-x1,2.)+pow(y2-y1,2.)+pow(z2-z1,2.),0.5);
			//	  mx=0.5*(x1+x2);
			//	  my=0.5*(y1+y2);
			//	  mz=0.5*(z1+z2);
			//	  nx=(x2-x1)/dist;
			//	  ny=(y2-y1)/dist;
			//	  nz=(z2-z1)/dist;
			//	  x1=mx-0.5*db_scale*nx;
			//	  x2=mx+0.5*db_scale*nx;
			//	  y1=my-0.5*db_scale*ny;
			//	  y2=my+0.5*db_scale*ny;
			//	  z1=mz-0.5*db_scale*nz;
			//	  z2=mz+0.5*db_scale*nz;
			// */
			
			// // check if reasonable
			// /* dist=pow(pow(x2-x1,2.)+pow(y2-y1,2.)+pow(z2-z1,2.),0.5);
			// if (fabs(dist-38)>1){
			//	  match=0;
			// }*/
		 
			//check if all quadruplets or triplets
			checksum_1=0;
			checksum_2=0;
			for(j=0; j<4; j++) {
				if (a1[1] < 0)
					checksum_1++;
				if (a2[1] < 0)
					checksum_2++;
			}
			if (checksum_1>1 || checksum_2>1)
				match=0;
			// end of check if all quadruplets or triplets
		}
		else
			match=0;
		
		fclose (fpp);
		fpp = fopen (res_name, "w");
		if (match==2) {
			fprintf (fpp, "%4d\n", match);
			fprintf (fpp, " %4d %9.3f %9.3f %9.3f %4d %4d %4d %4d\n", i1,x1,y1,z1,a1[0],a1[1],a1[2],a1[3]);
			fprintf (fpp, " %4d %9.3f %9.3f %9.3f %4d %4d %4d %4d\n", i2,x2,y2,z2,a2[0],a2[1],a2[2],a2[3]);
		}
		else
			fprintf (fpp, "%4d\n", 0);
		fclose (fpp);
	}
	//end of dumbbell treatment

	rmsX = sqrt(rmsX/match); rmsY = sqrt(rmsY/match); rmsZ = sqrt(rmsZ/match);
	mean_sigma0 = sqrt (mean_sigma0/match);

	sprintf (buf, "Match: %d, => rms = %4.2f micron, rms_x,y,z = %5.3f/%5.3f/%5.3f mm", match, mean_sigma0*1000, rmsX, rmsY, rmsZ);
	puts (buf);
	Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	Tcl_Eval(interp, ".text delete 3");
	Tcl_Eval(interp, ".text insert 3 $tbuf");

	/* sort coordinates for binary search in epi line segment drawing */
	for (i=0; i<n_img; i++)
		quicksort_coord2d_x (geo[0], num[0]);

	puts ("Determinate done\n");
	return TCL_OK;
}


int sequence_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int    i, j, ok, k, nslices=19, slicepos=0, pft_version = 3;
	char   seq_ch[128], seq_name[4][128];
	Tk_PhotoHandle img_handle;
	Tk_PhotoImageBlock img_block;
	double slice_step;
	double slicethickness;
	double zdim;						// z_cen_slice[19], removed ad holten, 12-2012
	int dumbbell=0,step_shake=1;
	double dummy;

	fpp = fopen_rp ("parameters/sequence.par");		// replaced fopen_r, ad holten, 12-2012
	if (!fpp) return TCL_OK;

	for (i=0; i<4; i++)
		fscanf (fpp, "%s\n", seq_name[i]);	   /* name of sequence */
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);

	display = atoi(argv[1]); 
	// Beat Mai 2010 for dumbbell
	if (atoi(argv[1])==3) {
		dumbbell=1;
		display=0;
	}

	/* scanning ptv ************** */
	printf("\nObject volume is scanned in %d slices!\n", nslices);
	slicepos=0;
	/* read illuminated Volume */
	fpp = fopen_rp ("parameters/criteria.par");		// replaced fopen_r, ad holten, 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%lf\n", &X_lay[0]);
	fscanf (fpp, "%lf\n", &Zmin_lay[0]);
	fscanf (fpp, "%lf\n", &Zmax_lay[0]);
	fscanf (fpp, "%lf\n", &X_lay[1]);
	fscanf (fpp, "%lf\n", &Zmin_lay[1]);
	fscanf (fpp, "%lf\n", &Zmax_lay[1]);
	fscanf (fpp, "%lf", &cnx);
	fscanf (fpp, "%lf", &cny);
	fscanf (fpp, "%lf", &cn);
	fscanf (fpp, "%lf", &csumg);
	fscanf (fpp, "%lf", &corrmin);
	fscanf (fpp, "%lf", &eps0);
	fclose (fpp);

	/* read illuminated layer data */
	if (dumbbell==1) {
		fpp = fopen ("parameters/dumbbell.par", "r");
		if (fpp){
			fscanf (fpp, "%lf", &eps0);
			fscanf (fpp, "%lf", &dummy);
			fscanf (fpp, "%lf", &dummy);
			fscanf (fpp, "%lf", &dummy);
			fscanf (fpp, "%d", &step_shake);
			fclose (fpp);
		}
		else {
			fpp = fopen ("parameters/dumbbell.par", "w");
			fprintf(fpp,"%lf\n", 5.0);
			fprintf(fpp,"%lf\n", 46.5);
			fprintf(fpp,"%lf\n", 0.5);
			fprintf(fpp,"%lf\n", 2.);
			fprintf(fpp,"%d\n", 2);
			fprintf(fpp,"%d\n",500);
			fprintf(fpp,"\n\n");
			fprintf(fpp,"explanation for parameters:\n");
			fprintf(fpp,"\n");
			fprintf(fpp,"1: eps (mm)\n");
			fprintf(fpp,"2: dunbbell scale\n");
			fprintf(fpp,"3: gradient descent factor\n");
			fprintf(fpp,"4: weight for dumbbell penalty\n");
			fprintf(fpp,"5: step size through sequence\n");
			fprintf(fpp,"6: num iterations per click\n");
			fclose(fpp);
			eps0 = 10;
		}
	}

	mmp.nlay = 1;

	zdim = Zmax_lay[0]-Zmin_lay[0];
	slice_step = zdim/nslices;
	slicethickness = 5.0;

	printf("\nzdim: %f, max: %f, min: %f, st: %f\n", zdim,Zmax_lay[0], Zmin_lay[0], slice_step);

	// changed by Beat July 2010 to allow for faster proc. when stepping through data
	for (i=seq_first; i<seq_last+1; i+=step_shake) {

		printf("\nstep: %d, zslice[j]: %f, slicepos: %d\n", i);

		//		Zmax_lay[0]= z_cen_slice[slicepos] - slicethickness/2.0;
		//		Zmin_lay[0]= z_cen_slice[slicepos] + slicethickness/2.0;
		//		Zmax_lay[1]= z_cen_slice[slicepos] - slicethickness/2.0;
		//		Zmin_lay[1]= z_cen_slice[slicepos] + slicethickness/2.0;
		//printf("in sequence zslice[j]: %f, zmin0: %f, zmax0: %f\n",
		//z_cen_slice[slicepos], Zmax_lay[0],Zmin_lay[0] );

		//slicepos++; if (slicepos==nslices) {slicepos=0;}

		// if      (i < 10)  sprintf (seq_ch, "%1d", i);
		// else if (i < 100) sprintf (seq_ch, "%2d", i);
		// else              sprintf (seq_ch, "%3d", i);
		// changed, ad holten 12-2012
		sprintf (seq_ch, "%d", i);

		for (j=0; j<n_img; j++) {
			sprintf (img_name[j],    "%s%s",    seq_name[j], seq_ch);
			sprintf (img_lp_name[j], "%s%s_lp", seq_name[j], seq_ch);
			sprintf (img_hp_name[j], "%s%s_hp", seq_name[j], seq_ch);
		}

		// Beat Mai 2010 for dumbbell
		if (dumbbell==0) {
			if (chfield == 0) sprintf (res_name, "res/rt_is.%s", seq_ch);
			else			  sprintf (res_name, "res/rt_is.%s_%1d", seq_ch, chfield);
		}
		else {
			if (chfield == 0) sprintf (res_name, "res/db_is.%s", seq_ch);
			else			  sprintf (res_name, "res/db_is.%s_%1d", seq_ch, chfield);
		}
		sprintf (buf, "\nImages:");
		for (j=0; j<n_img; j++) 
			sprintf (buf, "%s  %s", buf, img_name[j]);
		puts (buf);

		/* calling function for each sequence-n-tupel */
		/* read and display original images */

		for (k=0; k<n_img; k++) {
			/* reading */
			read_image (interp, img_name[k], img[k]);

			if (display) {
				img_handle = Tk_FindPhoto( interp, "temp");
				Tk_PhotoGetImage (img_handle, &img_block);
				tclimg2cimg (interp, img[k], &img_block);
				sprintf(buf, "newimage %d", k+1);
				Tcl_Eval(interp, buf);
			}
		}

		/* read pft version  */
		/* added by Alex for external detection procedure, 19.4.10 */
		fpp = fopen ("parameters/pft_version.par", "r");
		if (fpp) {
			fscanf (fpp, "%d\n", &pft_version);
			pft_version=pft_version+3;
			fclose (fpp);
		}

		if (hp_flag) {
			pre_processing_c (clientData, interp, argc, argv);
			puts("\nHighpass switched on\n");
		} 
		else
			puts("\nHighpass switched off\n");

		if (display)
			Tcl_Eval(interp, "update idletasks");

		/**************************************************************************/
		/* pft_version = 4 means external detection, Alex, 19.4.10 */

		if (pft_version == 4) { 
			for (k=0; k<n_img; k++) {
				read_targets(k,i,&num[k]);
				// sprintf (buf,"%d: %d,  ", k+1, num[k]);
				// strcat(val, buf);
				/* proper sort of targets in y-direction for later binary search */
				/* and for dimitris' tracking */
				quicksort_target_y (pix[k], num[k]);
				/* reorganize target numbers */
				for (j=0; j<num[k]; j++)  pix[k][j].pnr = j;
			}
			/*
			sprintf (buf, "Number of detected particles per image");
			Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
			Tcl_Eval(interp, ".text delete 2");
			Tcl_Eval(interp, ".text insert 2 $tbuf");

			Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
			Tcl_Eval(interp, ".text delete 3");
			Tcl_Eval(interp, ".text insert 3 $tbuf");

			printf("%s\n", val);
			return TCL_OK;
			*/
		} 
		/**************************************************************************/
		else {
			// added i to the detection_proc_c to get 'filenumber' for external API, Alex, 19.04.10
			detection_proc_c (clientData, interp, argc, argv);
		}

		if (display)  Tcl_Eval(interp, "update idletasks");
		correspondences_proc_c (clientData, interp, argc, argv);
		if (display)  Tcl_Eval(interp, "update idletasks");
		if(n_img>1)
			determination_proc_c (clientData, interp, argc, argv);

		/* delete unneeded files */
		for (j=0; j<n_img; j++) {
			ok = remove (img_lp_name[j]);
			ok = remove (img_hp_name[j]);
		}
	}
	/* reset of display flag */
	display = 1;

	return TCL_OK;
}

int restore_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int i_img,i;

	fpp = fopen_rp("parameters/ptv.par");		// replaced fopen_r, ad holten, 12-2012
	if (!fpp) return TCL_OK;

	fscanf (fpp, "%d\n", &n_img);
	for (i=0; i<n_img; i++) {
		fscanf (fpp, "%s\n", img_name[i]);
		fscanf (fpp, "%s\n", img_cal[i]);
	}
	fclose (fpp);

	/*	create file names  */
	for (i=0; i<n_img; i++) {
		strcpy (img_ori[i], img_cal[i]);	strcat (img_ori[i], ".ori");
		strcpy (img_addpar[i], img_cal[i]); strcat (img_addpar[i],".addpar");
	}

	// -- changed syntax, ad holten 12-2012
	//		strcpy (safety[0], "safety_0");
	//		strcat (safety[0], ".ori");
	//		strcpy (safety[1], "safety_1");
	//		strcat (safety[1], ".ori");
	//		strcpy (safety[2], "safety_2");
	//		strcat (safety[2], ".ori");
	//		strcpy (safety[3], "safety_3");
	//		strcat (safety[3], ".ori");
	//		strcpy (safety_addpar[0], "safety_0");
	//		strcat (safety_addpar[0], ".addpar");
	//		strcpy (safety_addpar[1], "safety_1");
	//		strcat (safety_addpar[1], ".addpar");
	//		strcpy (safety_addpar[2], "safety_2");
	//		strcat (safety_addpar[2], ".addpar");
	//		strcpy (safety_addpar[3], "safety_3");
	//		strcat (safety_addpar[3], ".addpar");
	for (i=0; i<4; i++) {
		sprintf(safety[i],        "safety_%d.ori", i);
		sprintf(safety_addpar[i], "safety_%d.addpar", i);
	}

	for (i_img=0; i_img<n_img; i_img++) {
		read_ori (&Ex[i_img], &I[i_img], &G[i_img], safety[i_img]);
		fp1 = fopen (safety_addpar[i_img], "r");
		if (! fp1)	fp1 = fopen ("addpar.raw", "r");

		if (fp1) {
			fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
				&ap[i_img].k1, &ap[i_img].k2, &ap[i_img].k3,
				&ap[i_img].p1, &ap[i_img].p2, &ap[i_img].scx, &ap[i_img].she);
			fclose (fp1);}
		else {
			printf("no addpar.raw\n");
			ap[i_img].k1=ap[i_img].k2=ap[i_img].k3=ap[i_img].p1=ap[i_img].p2=ap[i_img].she=0.0;
			ap[i_img].scx=1.0;
		}
	
		write_ori (Ex[i_img], I[i_img], G[i_img], img_ori[i_img]);
		fp1 = fopen (img_addpar[i_img], "w");
		fprintf (fp1, "%f %f %f %f %f %f %f",
			ap[i_img].k1, ap[i_img].k2, ap[i_img].k3,
			ap[i_img].p1, ap[i_img].p2, ap[i_img].scx, ap[i_img].she);
		fclose (fp1);
	}
	puts ("\nOrientation files restored");

	return TCL_OK;
}

int calibration_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int 	 i, j, sel, i_img, k, n, sup, dummy, multi,planes, prev, next; 
	double	 e1, e2, e3, varx, vary, cx, cy, nxp, nyp;
	int 	 intx1, inty1, intx2, inty2;
	coord_2d apfig1[11][11];	/* regular grid for ap figures */
	coord_2d apfig2[11][11];	/* ap figures */
	coord_3d fix4[4];			/* object points for preorientation */
	coord_2d crd0[4][12];		/* image points for preorientation */
	char	 multi_filename[10][256],filename[256], val[256],filename2[256];
	const char *valp;

	FILE	*FILEIN;
	char	filein[256];
	FILE	*FILEIN_ptv;
	char	filein_ptv[256];
	FILE	*FILEIN_T;
	char	filein_T[256];
	int 	filenumber, dumy, frameCount, currentFrame;
	int 	a[4], a1[4], a2[4], success=1;
	int 	good[10000];

	double	d_inner=0., av_dist=0., x, y;
	double	X[4],Y[4],Z[4],aa[4],b[4],c[4],dist,X_pos,Y_pos,Z_pos,XX,YY,ZZ;
	double	pX[6],pY[6],pZ[6];
	double	stdX,stdY,stdZ,rmsX[10000],av_rmsX;
	int 	count_inner=0,count_outer=0,pair_count=0,count_dist=0,m;

	// removed unreferenced variables, ad holten, 12-2012
	// Y2, Y1, residual, Z1, temp, xb, xa, X2, dist_error, ya, X1, Z2, dummy_float, yb

	Tk_PhotoHandle img_handle;
	Tk_PhotoImageBlock img_block;

	/* read support of unsharp mask */
	fp1 = fopen ("parameters/unsharp_mask.par", "r");
	if (! fp1)
		sup = 12;
	else {
		fscanf (fp1, "%d\n", &sup); 
		fclose (fp1);
	}

	/* Get Selection value from TclTk */
	valp = Tcl_GetVar(interp, "sel",  TCL_GLOBAL_ONLY);
	sel = atoi (valp);

	/* Beat Mai 2007 to set the variable examine for mulit-plane calibration*/
	fp1 = fopen_rp ("parameters/examine.par");		// replaced fopen_r, ad holten, 12-2012
	if (!fp1) return TCL_OK;
	fscanf (fp1,"%d\n", &dummy);
	fscanf (fp1,"%d\n", &multi);
	fclose (fp1);
	if (dummy==1) 
		examine = 4;
	else
		examine=0;

	/* Oswald Juni 2008 accept pairs------------------------------- */

	fp1 = fopen_rp ("parameters/cal_ori.par");		// replaced fopen_r, ad holten, 12-2012
	if (!fp1) return TCL_OK;

	fscanf (fp1, "%s\n", fixp_name);
	for (i=0; i<4; i++) {
		fscanf (fp1, "%s\n", img_name[i]);
		fscanf (fp1, "%s\n", img_ori0[i]);
	}
	fscanf (fp1, "%d\n", &tiff_flag);		// ad holten 12-2012, 2x bug repaired, was fpp 
	fscanf (fp1, "%d\n", &pair_flag);
	fscanf (fp1, "%d\n", &chfield);			// added, ad holten 12-2012
	fclose (fp1);
	
	//if (pair_flag==1) {
	//	  int OSWALDDUMMY=1;
	//}
	//else {
	//	  int OSWALDDUMMY=0;
	//}

	///////////////////////////////////////////////////////////////////////////////

	switch (sel)
	{
	case 1: 	/*	read calibration parameter file  */
		// 10 lines commented out, file already parsed, ad holten 12-2012
		// 		fp1 = fopen_r ("parameters/cal_ori.par");
		// 		fscanf (fp1,"%s\n", fixp_name);
		// 		for (i=0; i<4; i++) {
		//			fscanf (fp1, "%s\n", img_name[i]);
		// 			fscanf (fp1, "%s\n", img_ori0[i]);
		// 		}
		// 		fscanf (fpp, "%d\n", &tiff_flag);
		// 		fscanf (fpp, "%d\n", &pair_flag);
		// 		fscanf (fp1, "%d\n", &chfield);
		// 		fclose (fp1);

		/*	create file names  */
		// replaced by sprintf, ad holten 12-2012 
		//		for (i=0; i<n_img; i++) {
		//			strcpy (img_ori[i], img_name[i]);
		//			strcat (img_ori[i], ".ori");
		//			strcpy (img_addpar0[i], img_name[i]);
		//			strcat (img_addpar0[i], ".addpar0");
		//			strcpy (img_addpar[i], img_name[i]);
		//			strcat (img_addpar[i], ".addpar");
		//			strcpy (img_hp_name[i], img_name[i]);
		//			strcat (img_hp_name[i], "_hp");
		//		}
		//		strcpy (safety[0], "safety_0");
		//		strcat (safety[0], ".ori");
		//		strcpy (safety[1], "safety_1");
		//		strcat (safety[1], ".ori");
		//		strcpy (safety[2], "safety_2");
		//		strcat (safety[2], ".ori");
		//		strcpy (safety[3], "safety_3");
		//		strcat (safety[3], ".ori");
		//		strcpy (safety_addpar[0], "safety_0");
		//		strcat (safety_addpar[0], ".addpar");
		//		strcpy (safety_addpar[1], "safety_1");
		//		strcat (safety_addpar[1], ".addpar");
		//		strcpy (safety_addpar[2], "safety_2");
		//		strcat (safety_addpar[2], ".addpar");
		//		strcpy (safety_addpar[3], "safety_3");
		//		strcat (safety_addpar[3], ".addpar");

		for (i=0; i<n_img; i++) {
			sprintf(img_ori[i],       "%s.ori",     img_name[i]);
			sprintf(img_addpar0[i],   "%s.addpar0", img_name[i]);
			sprintf(img_addpar[i],    "%s.addpar",  img_name[i]);
			sprintf(img_hp_name[i],   "%s_hp",      img_name[i]);
			sprintf(safety[i],        "safety_%d.ori",   i);
			sprintf(safety_addpar[i], "safety_%d.addpa", i);
		}

		for (i=0; i<n_img; i++) {
			zoom_x[i] = imx/2, zoom_y[i] = imy/2, zoom_f[i] = 1;

			read_image (interp, img_name[i], img[i]);

			sprintf(val, "camcanvas %d", i+1);
			Tcl_Eval(interp, val);

			img_handle = Tk_FindPhoto( interp, "temp");
			Tk_PhotoGetImage (img_handle, &img_block);
			tclimg2cimg (interp, img[i], &img_block);

			sprintf(val, "newimage %d", i+1);
			Tcl_Eval(interp, val);
		}
		break;


	case 2: puts ("Detection procedure"); strcpy(val,"");

		/* Highpass Filtering */
		pre_processing_c (clientData, interp, argc, argv);

		/* reset zoom values */
		for (i=0; i<n_img; i++) {
			zoom_x[i] = imx/2; zoom_y[i] = imy/2; zoom_f[i] = 1;
		}

		/* copy images because the target recognition will set greyvalues to zero */
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


	case 3:    pp1=0;	 pp2=0;    pp3=0;	 pp4=0;

		for (i=0; i<n_img; i++) {
			sprintf (buf, "%d targets remain", num[i]);
			puts (buf);
		}
		fp1 = fopen_rp ("parameters/man_ori.par");		// replaced fopen_r, ad holten, 12-2012
		if (!fp1) break;

		for (i=0; i<n_img; i++)
			fscanf (fp1, "%d %d %d %d \n", &nr[i][0], &nr[i][1], &nr[i][2], &nr[i][3]);
		fclose (fp1);

		for (i=0; i<n_img; i++) {
			sprintf(val, "measure %d %d %d %d %d", nr[i][0], nr[i][1], nr[i][2], nr[i][3], i+1);
			Tcl_Eval(interp, val);
			// replaced by a for loop, ad holten 12-2012
			//valp = Tcl_GetVar(interp, "px0",  TCL_GLOBAL_ONLY);
			//pix0[i][0].x = atoi (valp);
			//valp = Tcl_GetVar(interp, "py0",  TCL_GLOBAL_ONLY);
			//pix0[i][0].y = atoi (valp);
			//valp = Tcl_GetVar(interp, "px1",  TCL_GLOBAL_ONLY);
			//pix0[i][1].x = atoi (valp);
			//valp = Tcl_GetVar(interp, "py1",  TCL_GLOBAL_ONLY);
			//pix0[i][1].y = atoi (valp);
			//valp = Tcl_GetVar(interp, "px2",  TCL_GLOBAL_ONLY);
			//pix0[i][2].x = atoi (valp);
			//valp = Tcl_GetVar(interp, "py2",  TCL_GLOBAL_ONLY);
			//pix0[i][2].y = atoi (valp);
			//valp = Tcl_GetVar(interp, "px3",  TCL_GLOBAL_ONLY);
			//pix0[i][3].x = atoi (valp);
			//valp = Tcl_GetVar(interp, "py3",  TCL_GLOBAL_ONLY);
			//pix0[i][3].y = atoi (valp);
			for (j=0; j<4; j++) {
				sprintf(val, "px%d", j);
				pix0[i][j].x = atoi( Tcl_GetVar(interp, val,  TCL_GLOBAL_ONLY) );
				sprintf(val, "py%d", j);
				pix0[i][j].y = atoi( Tcl_GetVar(interp, val,  TCL_GLOBAL_ONLY) );
			}
		}

		/* write measured coordinates to file for next trial */
		fp1 = fopen ("man_ori.dat", "w");
		for (i=0; i<n_img; i++)
			for (j=0; j<4; j++)
				fprintf (fp1, "%f %f\n", pix0[i][j].x, pix0[i][j].y);
		fclose (fp1);

		break;


	case 4: /* read pixel coordinates of older pre-orientation */

		/* read point numbers of pre-clicked points */
		fp1 = fopen_rp ("parameters/man_ori.par");		// replaced fopen_r, ad holten, 12-2012
		if (!fp1) break;
		for (i=0; i<n_img; i++)
			fscanf (fp1, "%d %d %d %d \n", &nr[i][0], &nr[i][1], &nr[i][2], &nr[i][3]);
		fclose (fp1);

		/* read coordinates of pre-clicked points */
		fp1 = fopen_rp("man_ori.dat");					// replaced fopen_r, ad holten, 12-2012
		if (!fp1) break;
		for (i_img=0; i_img<n_img; i_img++)
			for (i=0; i<4; i++)
		{
			fscanf (fp1, "%lf %lf\n", &pix0[i_img][i].x, &pix0[i_img][i].y);
			drawcross (interp,	(int) pix0[i_img][i].x,
				(int) pix0[i_img][i].y, cr_sz+2, i_img, "red");
			draw_pnr (interp, (int) pix0[i_img][i].x, (int) pix0[i_img][i].y,
				nr[i_img][i], i_img, "red");
		}
		fclose (fp1);

		break;


	case 5: puts ("Sort grid points");
		for (i=0; i<n_img; i++) {
			/* read control point coordinates for man_ori points */
			fp1 = fopen_rp(fixp_name);					// replaced fopen_r, ad holten, 12-2012
			if (!fp1) return TCL_OK;
			k = 0;
			while (fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
					&fix[k].x, &fix[k].y, &fix[k].z) != EOF)  k++;
			fclose (fp1);
			nfix = k;

			/* take clicked points from control point data set */
			for (j=0; j<4; j++)
				for (k=0; k<nfix; k++)
					if (fix[k].pnr == nr[i][j])  fix4[j] = fix[k];

			/* get approx for orientation and ap */
			read_ori (&Ex[i], &I[i], &G[i], img_ori0[i]);
			fp1 = fopen (img_addpar0[i], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					&ap[i].k1, &ap[i].k2, &ap[i].k3,
					&ap[i].p1, &ap[i].p2, &ap[i].scx,&ap[i].she);
				fclose (fp1);
			}
			else {
				printf("no addpar.raw\n");
				ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
				ap[i].scx=1.0;
			}


			/* transform clicked points */
			for (j=0; j<4; j++) {
				pixel_to_metric (pix0[i][j].x, pix0[i][j].y, imx,imy, pix_x, pix_y,
					&crd0[i][j].x, &crd0[i][j].y, chfield);
				correct_brown_affin (crd0[i][j].x, crd0[i][j].y, ap[i],
					&crd0[i][j].x, &crd0[i][j].y);
			}

			/* raw orientation with 4 points */
			raw_orient_v3 (Ex[i], I[i], G[i], ap[i], mmp, 4, fix4, crd0[i], &Ex[i],&G[i],0); //Beat Nov 2008
			sprintf (filename, "raw%d.ori", i);
			write_ori (Ex[i], I[i], G[i], filename);
	 
			/* sorting of detected points by back-projection */
			sortgrid_man (interp, Ex[i], I[i], G[i], ap[i], mmp, imx, imy,
				pix_x,pix_y, nfix, fix, num[i], pix[i], chfield, i);

			/* adapt # of detected points */
			num[i] = nfix;

			for (j=0; j<nfix; j++) {
				if (pix[i][j].pnr < 0)	  continue;
				intx1 = (int) pix[i][j].x ;
				inty1 = (int) pix[i][j].y ;

				drawcross (interp, intx1, inty1, cr_sz, i, "white");
				draw_pnr (interp, intx1, inty1, fix[j].pnr, i, "white");
			}
		}

		/* dump dataset for rdb */
		if (examine == 4) {
			/* create filename for dumped dataset */
			sprintf (filename, "dump_for_rdb");
			fp1 = fopen (filename, "w");

			/* write # of points to file */
			fprintf (fp1, "%d\n", nfix);

			/* write point and image coord to file */
			for (i=0; i<nfix; i++) {
				fprintf (fp1, "%4d %10.3f %10.3f %10.3f   %d	",
					fix[i].pnr, fix[i].x, fix[i].y, fix[i].z, 0);
				for (i_img=0; i_img<n_img; i_img++) {
					if (pix[i_img][i].pnr >= 0) {
						/* transform pixel coord to metric */
						pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y, imx,imy, 
							pix_x, pix_y, &crd[i_img][i].x, &crd[i_img][i].y, chfield);
						fprintf (fp1, "%4d %8.5f %8.5f	  ", pix[i_img][i].pnr,
							crd[i_img][i].x, crd[i_img][i].y);
					}
					else {
						fprintf (fp1, "%4d %8.5f %8.5f	  ", pix[i_img][i].pnr, 0.0, 0.0);
					}
				}
				fprintf (fp1, "\n");
			}
			fclose (fp1);
			printf ("dataset dumped into %s\n", filename);
		}
		break;


	case 16: puts ("Sort grid points using files"); //Beat Jan 2011
		for (i=0; i<n_img; i++) {
			/* read control point coordinates for man_ori points */
			fp1 = fopen_rp(fixp_name);					// replaced fopen_r, ad holten, 12-2012
			if (!fp1) return TCL_OK;
			k = 0;
			while ( fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
					&fix[k].x, &fix[k].y, &fix[k].z) != EOF)  k++;
			fclose (fp1);
			nfix = k;

			/* take clicked points from control point data set */
			for (j=0; j<4; j++)
				for (k=0; k<nfix; k++)
					if (fix[k].pnr == nr[i][j])    fix4[j] = fix[k];

			/* get approx for orientation and ap */
			read_ori (&Ex[i], &I[i], &G[i], img_ori0[i]);
			fp1 = fopen (img_addpar0[i], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf", &ap[i].k1,&ap[i].k2,&ap[i].k3,
					&ap[i].p1,&ap[i].p2, &ap[i].scx,&ap[i].she);
				fclose (fp1);
			}
			else {
				printf("no addpar.raw\n");
				ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
				ap[i].scx=1.0;
			}

			/* transform clicked points */
			for (j=0; j<4; j++) {
				pixel_to_metric (pix0[i][j].x, pix0[i][j].y, imx,imy, pix_x, pix_y,
				   &crd0[i][j].x, &crd0[i][j].y, chfield);
				correct_brown_affin (crd0[i][j].x, crd0[i][j].y, ap[i],
				   &crd0[i][j].x, &crd0[i][j].y);
			}

			/* raw orientation with 4 points */
			raw_orient_v3 (Ex[i], I[i], G[i], ap[i], mmp, 4, fix4, crd0[i], &Ex[i],&G[i],0); //Beat Nov 2008
			sprintf (filename, "raw%d.ori", i);
			write_ori (Ex[i], I[i], G[i], filename);
	 
			/* sorting of detected points by back-projection */
			sortgrid_file (interp, Ex[i], I[i], G[i], ap[i], mmp, imx,imy, pix_x,pix_y,
				nfix, fix, num[i], pix[i], chfield, i);

			/* adapt # of detected points */
			num[i] = nfix;

			for (j=0; j<nfix; j++) {
				if (pix[i][j].pnr < 0)	  continue;
				intx1 = (int) pix[i][j].x ;
				inty1 = (int) pix[i][j].y ;

				drawcross (interp, intx1, inty1, cr_sz, i, "white");
				draw_pnr (interp, intx1, inty1, fix[j].pnr, i, "white");
			}
		}

		/* dump dataset for rdb */
		if (examine == 4) {
			/* create filename for dumped dataset */
			sprintf (filename, "dump_for_rdb");
			fp1 = fopen (filename, "w");

			/* write # of points to file */
			fprintf (fp1, "%d\n", nfix);

			/* write point and image coord to file */
			for (i=0; i<nfix; i++) {
				fprintf (fp1, "%4d %10.3f %10.3f %10.3f   %d	",
					fix[i].pnr, fix[i].x, fix[i].y, fix[i].z, 0);
				for (i_img=0; i_img<n_img; i_img++) {
					if (pix[i_img][i].pnr >= 0) {
						/* transform pixel coord to metric */
						pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y, imx, imy,
							pix_x, pix_y, &crd[i_img][i].x, &crd[i_img][i].y, chfield);
						fprintf (fp1, "%4d %8.5f %8.5f	  ", pix[i_img][i].pnr,
							crd[i_img][i].x, crd[i_img][i].y);
					}
					else {
						fprintf (fp1, "%4d %8.5f %8.5f	  ", pix[i_img][i].pnr, 0.0, 0.0);
					}
				}
				fprintf (fp1, "\n");
			}
			fclose (fp1);
			printf ("dataset dumped into %s\n", filename);
		}
		break;


	case 6: puts ("Orientation"); strcpy(buf, "");
		// replaced next lines by for loop, ad holten 12-2012
		// strcpy (safety[0], "safety_0");
		// strcat (safety[0], ".ori");
		// strcpy (safety[1], "safety_1");
		// strcat (safety[1], ".ori");
		// strcpy (safety[2], "safety_2");
		// strcat (safety[2], ".ori");
		// strcpy (safety[3], "safety_3");
		// strcat (safety[3], ".ori");
		// strcpy (safety_addpar[0], "safety_0");
		// strcat (safety_addpar[0], ".addpar");
		// strcpy (safety_addpar[1], "safety_1");
		// strcat (safety_addpar[1], ".addpar");
		// strcpy (safety_addpar[2], "safety_2");
		// strcat (safety_addpar[2], ".addpar");
		// strcpy (safety_addpar[3], "safety_3");
		// strcat (safety_addpar[3], ".addpar");
		for (i=0; i<4; i++) {
			sprintf(safety[i], "safety_%d.ori", i);
			sprintf(safety_addpar[i], "safety_%d.addpar", i);
		}

		for (i_img=0; i_img<n_img; i_img++) {
			for (i=0; i<nfix ; i++) {
				pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y, imx,imy, 
					pix_x, pix_y, &crd[i_img][i].x, &crd[i_img][i].y, chfield);
				crd[i_img][i].pnr = pix[i_img][i].pnr;
			}

			/* save data for special use of resection routine */
			if (examine == 4 && multi == 0) {
				printf ("try write resection data to disk\n");
				/* point coordinates */
				//sprintf (filename, "resect_%s.fix", img_name[i_img]);
				sprintf (filename, "%s.fix", img_name[i_img]);
				write_ori (Ex[i_img], I[i_img], G[i_img], img_ori[i_img]);
				fp1 = fopen (filename, "w");
				for (i=0; i<nfix; i++)
					fprintf (fp1, "%3d	%10.5f	%10.5f	%10.5f\n",
						fix[i].pnr, fix[i].x, fix[i].y, fix[i].z);
				fclose (fp1);

				/* metric image coordinates */
				//sprintf (filename, "resect_%s.crd", img_name[i_img]);
				sprintf (filename, "%s.crd", img_name[i_img]);
				fp1 = fopen (filename, "w");
				for (i=0; i<nfix; i++)
					fprintf (fp1, "%3d	%9.5f  %9.5f\n", crd[i_img][i].pnr,
						crd[i_img][i].x, crd[i_img][i].y);
				fclose (fp1);

				/* orientation and calibration approx data */
				write_ori (Ex[i_img], I[i_img], G[i_img], "resect.ori0");
				fp1 = fopen ("resect.ap0", "w");
				fprintf (fp1, "%f %f %f %f %f %f %f",
					ap[i_img].k1, ap[i_img].k2, ap[i_img].k3, 
					ap[i_img].p1, ap[i_img].p2, ap[i_img].scx, ap[i_img].she);
				fclose (fp1);
				printf ("resection data written to disk\n");
			}

			/* resection routine */
			if (examine != 4)
				orient_v3 (interp, Ex[i_img], I[i_img], G[i_img], ap[i_img], mmp,
					nfix, fix, crd[i_img], &Ex[i_img], &I[i_img], &G[i_img], &ap[i_img], i_img);

			/* resection with dumped datasets */
			if (examine == 4) {

				//printf("Resection with dumped datasets? (y/n)");
				//scanf("%s",buf);
				//if (buf[0] != 'y')	continue;
				//strcpy (buf, "");
				if (multi==0)  continue;

				/* read calibration frame datasets */
				//sprintf (multi_filename[0],"img/calib_a_cam");
				//sprintf (multi_filename[1],"img/calib_b_cam");

				fp1 = fopen_rp ("parameters/multi_planes.par");			// replaced fopen_r, ad holten, 12-2012
				if (!fp1) return TCL_OK;
				fscanf (fp1,"%d\n", &planes);
				for(i=0; i<planes; i++) 
					fscanf (fp1,"%s\n", &multi_filename[i]);
				fclose(fp1);

				for (n=0, nfix=0, dump_for_rdb=0; n<planes; n++) {		// replaced 10 by 'planes', ad holten, 01-2013
					//sprintf (filename, "resect.fix%d", n);

					sprintf (filename, "%s%d.tif.fix", multi_filename[n],i_img+1);
					fp1 = fopen_rp (filename);							// replaced fopen by fopen_rp, and
					if (! fp1) break;									// continue by break, ad holten, 01-2013

					printf("reading file: %s\n", filename);
					printf("reading dumped resect data #%d\n", n);
					k = 0;
					while (fscanf (fp1, "%d %lf %lf %lf", &fix[nfix+k].pnr, 
								   &fix[nfix+k].x, &fix[nfix+k].y, &fix[nfix+k].z)
						   != EOF) k++;
					fclose (fp1);
					/* read metric image coordinates */
					//sprintf (filename, "resect_%d.crd%d", i_img, n);
					sprintf (filename, "%s%d.tif.crd", multi_filename[n], i_img+1);
					printf("reading file: %s\n", filename);
					fp1 = fopen_rp (filename);							// replaced fopen by fopen_rp and
					if (! fp1) break;									// added this line, ad holten, 01-2013
					for (i=nfix; i<nfix+k; i++)
						fscanf (fp1, "%d %lf %lf", &crd[i_img][i].pnr,
								&crd[i_img][i].x, &crd[i_img][i].y);
					nfix += k;
					fclose (fp1);
				}
				if (n<planes) {											// added, ad holten, 01-2013
					printf("Action aborted. (Error in multi_planes.par ?)\n");
					return TCL_OK;
				}

				/* resection */
				/*Beat Mai 2007*/
				sprintf (filename, "raw%d.ori", i_img);
				read_ori (&Ex[i_img], &I[i_img], &G[i_img], filename);
				fp1 = fopen ("addpar.raw", "r");

				if (fp1) {
					fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
						&ap[i_img].k1,&ap[i_img].k2,&ap[i_img].k3,
						&ap[i_img].p1,&ap[i_img].p2,
						&ap[i_img].scx,&ap[i_img].she);
					fclose (fp1);
				}
				else {
					printf("no addpar.raw\n");
					ap[i_img].k1=ap[i_img].k2=ap[i_img].k3=ap[i_img].p1=ap[i_img].p2=ap[i_img].she=0.0;
					ap[i_img].scx=1.0;
				}
				////////////////////////////////////////

				/* markus 14.05.2007 show coordinates combined */
				for (i=0; i<nfix; i++) {
					/* first crd->pix */
					metric_to_pixel (crd[i_img][i].x, crd[i_img][i].y, imx,imy, pix_x,pix_y, &pix[i_img][i].x, &pix[i_img][i].y, chfield);
					/*then draw crosses*/
					//intx1 = (int) pix[i_img][i].x;
					//inty1 = (int) pix[i_img][i].y;
					//drawcross (interp, intx1, inty1, 3, i_img, "orange");
				}
				orient_v3 (interp, Ex[i_img], I[i_img], G[i_img], ap[i_img], mmp,
					nfix, fix, crd[i_img], &Ex[i_img], &I[i_img], &G[i_img], &ap[i_img], i_img);

				///////////////////////////////////////////
			}

			/* save orientation and additional parameters */
			write_ori (Ex[i_img], I[i_img], G[i_img], img_ori[i_img]);
			fp1 = fopen( img_ori[i_img], "r" );
			if(fp1 != NULL) {
				fclose(fp1);
				read_ori (&sEx[i_img], &sI[i_img], &sG[i_img], img_ori[i_img]);
				fp1 = fopen (img_addpar0[i_img], "r");
				if (! fp1)	fp1 = fopen ("addpar.raw", "r");

				if (fp1) {
					fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					    &sap[i_img].k1, &sap[i_img].k2,&sap[i_img].k3,
					    &sap[i_img].p1, &sap[i_img].p2,
					    &sap[i_img].scx, &sap[i_img].she);
					fclose (fp1);
				} 
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
			else{
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
		}
		Tcl_Eval(interp, ".text delete 3");
		Tcl_Eval(interp, ".text delete 1");
		Tcl_Eval(interp, ".text insert 1 \"Orientation and self calibration \"");
		Tcl_Eval(interp, ".text delete 2");
		if (examine != 4)
			Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each image -> \"");
		if (examine == 4 && multi==0)
			Tcl_Eval(interp, ".text insert 2 \"resection data written to disk \"");
		Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text insert 3 $tbuf");

		break;

	case 7: checkpoint_proc (interp);
		sprintf(val,"blue: planimetry,	 yellow: height");
		Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text delete 2");
		Tcl_Eval(interp, ".text insert 2 $tbuf");
		Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text delete 3");
		Tcl_Eval(interp, ".text insert 3 $tbuf");
		break;


	case 8: /* draw additional parameter figures */

		Tcl_Eval(interp, "clearcam");

		/*	read orientation and additional parameters	*/
		for (i=0; i<n_img; i++) {
			fp1 = fopen_rp (img_addpar[i]);				// replaced fopen_r, ad holten, 12-2012
			if (!fp1 || !read_ori (&Ex[i], &I[i], &G[i], img_ori[i]))
				return TCL_OK;

			fscanf (fp1,"%lf %lf %lf %lf %lf %lf %lf",
				&ap[i].k1, &ap[i].k2, &ap[i].k3,
				&ap[i].p1, &ap[i].p2, &ap[i].scx, &ap[i].she);
			fclose (fp1);
		}
		for (i_img=0; i_img<n_img; i_img++) {
			/* create undistorted grid */
			for (i=0; i<11; i++)
				for (j=0; j<11; j++)
			{
				apfig1[i][j].x = i * imx/10;
				apfig1[i][j].y = j * imy/10;
			}
			/* draw undistorted grid */
			for (i=0; i<10; i++)
				for (j=0; j<10; j++)
			{
				intx1 = (int) apfig1[i][j].x;
				inty1 = (int) apfig1[i][j].y;
				intx2 = (int) apfig1[i+1][j].x;
				inty2 = (int) apfig1[i][j+1].y;
				drawvector (interp, intx1, inty1, intx2, inty1, 1, i_img, "black");
				drawvector (interp, intx1, inty1, intx1, inty2, 1, i_img, "black");
			}
			for (j=0; j<10; j++) {
				intx1 = (int) apfig1[10][j].x;
				inty1 = (int) apfig1[10][j].y;
				inty2 = (int) apfig1[10][j+1].y;
				drawvector (interp, intx1, inty1, intx1, inty2, 1, i_img, "black");
			}
			for (i=0; i<10; i++) {
				intx1 = (int) apfig1[i][10].x;
				inty1 = (int) apfig1[i][10].y;
				intx2 = (int) apfig1[i+1][10].x;
				drawvector (interp, intx1, inty1, intx2, inty1, 1, i_img, "black");
			}
			/* distort grid */
			for (i=0; i<11; i++)
				for (j=0; j<11; j++)
			{
				/* transform to metric, distort and re-transform */
				pixel_to_metric (apfig1[i][j].x, apfig1[i][j].y, imx,imy, pix_x,pix_y,
					&apfig2[i][j].x, &apfig2[i][j].y, chfield);
				distort_brown_affin (apfig2[i][j].x, apfig2[i][j].y,
					ap[i_img], &apfig2[i][j].x, &apfig2[i][j].y);
				metric_to_pixel (apfig2[i][j].x, apfig2[i][j].y,
				   imx,imy, pix_x,pix_y, &apfig2[i][j].x, &apfig2[i][j].y, chfield);
				/* exaggerate distortion by factor 5 */
				apfig2[i][j].x = 5*apfig2[i][j].x - 4*apfig1[i][j].x;
				apfig2[i][j].y = 5*apfig2[i][j].y - 4*apfig1[i][j].y;
			}
			/* draw distorted grid */
			for (i=0; i<10; i++)
				for (j=0; j<10; j++)
			{
				intx1 = (int) apfig2[i][j].x;
				inty1 = (int) apfig2[i][j].y;
				intx2 = (int) apfig2[i+1][j].x;
				inty2 = (int) apfig2[i+1][j].y;
				drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
				intx2 = (int) apfig2[i][j+1].x ;
				inty2 = (int) apfig2[i][j+1].y ;
				drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
			}
			for (j=0; j<10; j++) {
				intx1 = (int) apfig2[10][j].x;
				inty1 = (int) apfig2[10][j].y;
				intx2 = (int) apfig2[10][j+1].x;
				inty2 = (int) apfig2[10][j+1].y;
				drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
			}
			for (i=0; i<10; i++) {
				intx1 = (int) apfig2[i][10].x;
				inty1 = (int) apfig2[i][10].y;
				intx2 = (int) apfig2[i+1][10].x;
				inty2 = (int) apfig2[i+1][10].y ;
				drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
			}
		}
		// case 8: no break needed ???? adh, 10-12-2012
	case 9: puts ("Plot initial guess");
		
		for (i=0; i<n_img; i++) {
			/* read calblock points  */
			fp1 = fopen_rp(fixp_name);				// replaced fopen_r, ad holten, 12-2012
			if (!fp1) return TCL_OK;
			k = 0;
			while ( fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
				&fix[k].x, &fix[k].y, &fix[k].z) != EOF) k++;
			fclose (fp1);
			nfix = k;

			// /* take clicked points from control point data set */
			// /*for (j=0; j<4; j++)	for (k=0; k<nfix; k++)
			// {
			//	 if (fix[k].pnr == nr[i][j])	fix4[j] = fix[k];
			// }*/

			/* read orientation orientation and ap */
			read_ori (&Ex[i], &I[i], &G[i], img_ori0[i]);
			fp1 = fopen (img_addpar0[i], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					&ap[i].k1,&ap[i].k2,&ap[i].k3,
					&ap[i].p1,&ap[i].p2, &ap[i].scx,&ap[i].she);
				fclose (fp1);
			}
			else {
				printf("no addpar.raw\n");
					ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
					ap[i].scx=1.0;
			}
			// /* transform clicked points */
			// /*for (j=0; j<4; j++)
			//	{
			//	  pixel_to_metric (pix0[i][j].x, pix0[i][j].y,
			//			   imx,imy, pix_x, pix_y,
			//			   &crd0[i][j].x, &crd0[i][j].y,
			//			   chfield);
			//	  correct_brown_affin (crd0[i][j].x, crd0[i][j].y, ap[i],
			//			   &crd0[i][j].x, &crd0[i][j].y);
			// }*/
			// /* raw orientation with 4 points */
			// /*raw_orient_v3 (Ex[i], I[i], G[i], ap[i], mmp, 4, fix4, crd0[i], &Ex[i],&G[i],1);*/
	  
			/* just plot the stuff */
			just_plot (interp, Ex[i], I[i], G[i], ap[i], mmp,
				imx,imy, pix_x,pix_y, nfix, fix,  chfield, i);

			/*write artifical images*/
		}
		break;

	case 30: puts ("map mm to pixel");
		for (i=0; i<n_img; i++) {

			/* read calblock points  */
			// fp1 = fopen ("points4mapping.txt", "r");	// commented out, ad holten
			fp1 = fopen_rp(fixp_name);					// replaced fopen, ad holten, 12-2012
			if (!fp1) return TCL_OK;

			k = 0;
			while ( fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
					&fix[k].x, &fix[k].y, &fix[k].z) != EOF) k++;
			fclose (fp1);
			nfix = k;	   

			/* read orientation orientation and ap */
			read_ori (&Ex[i], &I[i], &G[i], img_ori0[i]);
			fp1 = fopen (img_addpar0[i], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					&ap[i].k1,&ap[i].k2,&ap[i].k3,
					&ap[i].p1,&ap[i].p2, &ap[i].scx,&ap[i].she);
				fclose (fp1);} else {
				printf("no addpar.raw\n");
					ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
					ap[i].scx=1.0;
			}

			/* reproject all calibration plate points into pixel space, varying all co*/
			sprintf (filename2, "mm2pixel_mapping_cam%d.txt", i);
			fp2 = fopen (filename2, "w");		  
			for (n=0; n<nfix; n++) {
				//first, compute central projection
				img_coord (fix[n].x, fix[n].y, fix[n].z,  Ex[i], I[i], G[i], ap[i], mmp, &nxp,&nyp);
				metric_to_pixel (nxp, nyp, imx,imy, pix_x,pix_y, &nxp, &nyp, chfield);
				cx=nxp;
				cy=nyp;
				//now project variation +/- 0.5mm
				varx=0;
				vary=0;
				for (e1=-0.5;e1<=0.5;e1=e1+1) {
					for (e2=-0.5;e2<=0.5;e2=e2+1) {
						for (e3=-0.5;e3<=0.5;e3=e3+1) {
							img_coord (fix[n].x+e1, fix[n].y+e2, fix[n].z+e3,  Ex[i], I[i], G[i], ap[i], mmp, &nxp,&nyp);
							metric_to_pixel (nxp, nyp, imx,imy, pix_x,pix_y, &nxp, &nyp, chfield);
							varx+=fabs(nxp-cx)/8;
							vary+=fabs(nyp-cy)/8;
						}
					}
				}
				fprintf (fp2, "%f %f\n",2*varx,2*vary); // so it write out 1mm per soandsomany pixel
			}
			fclose (fp2);
		}
		break;

case 14: puts ("Sortgrid = initial guess");

		for (i=0; i<n_img; i++)
		{
			/* read control point coordinates for man_ori points */
			fp1 = fopen_rp(fixp_name);				// replaced fopen_r, ad holten, 12-2012
			if (!fp1) return TCL_OK;
			k = 0;
			while ( fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
				&fix[k].x, &fix[k].y, &fix[k].z) != EOF) k++;
			fclose (fp1);
			nfix = k;

			/* take clicked points from control point data set */
			for (j=0; j<4; j++) 
				for (k=0; k<nfix; k++)
			{
				if (fix[k].pnr == nr[i][j])  fix4[j] = fix[k];
			}

			/* get approx for orientation and ap */
			read_ori (&Ex[i], &I[i], &G[i], img_ori0[i]);
			fp1 = fopen (img_addpar0[i], "r");
			if (! fp1)	fp1 = fopen ("addpar.raw", "r");

			if (fp1) {
				fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
					&ap[i].k1,&ap[i].k2,&ap[i].k3,
					&ap[i].p1,&ap[i].p2,
					&ap[i].scx,&ap[i].she);
				fclose (fp1);
			}
			else {
				printf("no addpar.raw\n");
				ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
				ap[i].scx=1.0;
			}

			/* transform clicked points */
			for (j=0; j<4; j++) {
				pixel_to_metric (pix0[i][j].x, pix0[i][j].y,
					imx,imy, pix_x, pix_y,
					&crd0[i][j].x, &crd0[i][j].y,
					chfield);
				correct_brown_affin (crd0[i][j].x, crd0[i][j].y, ap[i],
					&crd0[i][j].x, &crd0[i][j].y);
			}

			/* raw orientation with 4 points */
			raw_orient_v3 (Ex[i], I[i], G[i], ap[i], mmp, 4, fix4, crd0[i], &Ex[i],&G[i],1); /* correction 0 to 1 , al*/
			sprintf (filename, "raw%d.ori", i);
			write_ori (Ex[i], I[i], G[i], filename);

			/* sorting of detected points by back-projection */
			sortgrid_man (interp, Ex[i], I[i], G[i], ap[i], mmp,
				imx,imy, pix_x,pix_y,
				nfix, fix, num[i], pix[i], chfield, i);

			/* adapt # of detected points */
			num[i] = nfix;

			for (j=0; j<nfix; j++) {
				if (pix[i][j].pnr < 0)	  continue;
				intx1 = (int) pix[i][j].x;
				inty1 = (int) pix[i][j].y;

				drawcross (interp, intx1, inty1, cr_sz, i, "yellow");
				draw_pnr (interp, intx1, inty1, fix[j].pnr, i, "yellow");
			}
		}

		break;


	case 15: puts ("Show numbe on detected points");
		for (i=0; i<n_img; i++) {
			for (j=0; j<num[i]; j++) { 
				draw_pnr (interp, (int)pix[i][j].x, (int)pix[i][j].y, j, i, "blue");
			}
		}
		break; //Beat and Debashish Jan 2011 

 
	case 10: puts ("Orientation from particles"); strcpy(buf, "");
		
		// replaced next lines by for loop, ad holten 12-2012
		//		 strcpy (safety[0], "safety_0");
		//		 strcat (safety[0], ".ori");
		//		 strcpy (safety[1], "safety_1");
		//		 strcat (safety[1], ".ori");
		//		 strcpy (safety[2], "safety_2");
		//		 strcat (safety[2], ".ori");
		//		 strcpy (safety[3], "safety_3");
		//		 strcat (safety[3], ".ori");
		//		 strcpy (safety_addpar[0], "safety_0");
		//		 strcat (safety_addpar[0], ".addpar");
		//		 strcpy (safety_addpar[1], "safety_1");
		//		 strcat (safety_addpar[1], ".addpar");
		//		 strcpy (safety_addpar[2], "safety_2");
		//		 strcat (safety_addpar[2], ".addpar");
		//		 strcpy (safety_addpar[3], "safety_3");
		//		 strcat (safety_addpar[3], ".addpar");
		for (i=0; i<4; i++) {
			sprintf(safety[i], "safety_%d.ori", i);
			sprintf(safety_addpar[i], "safety_%d.addpar", i);
		}

		for (i_img=0; i_img<n_img; i_img++) {
			/* read control point coordinates for man_ori points */
   
			fpp = fopen_rp("parameters/sequence.par");				// replaced fopen_r, ad holten, 12-2012
			if (!fpp) return TCL_OK;
			for (i=0; i<4; i++)
				fscanf (fpp, "%s\n", seq_name[i]);	   /* name of sequence */
			fclose (fpp);

			fpp = fopen_rp("parameters/shaking.par");				// replaced fopen_r, ad holten, 12-2012
			if (!fpp) return TCL_OK;
			fscanf (fpp,"%d\n", &seq_first); 
			fscanf (fpp,"%d\n", &seq_last);
			fscanf (fpp,"%d\n", &max_shake_points);
			fscanf (fpp,"%d\n", &max_shake_frames);
			fclose (fpp);

			// next lines commented out, ad holten, 12-2012
			// /*	read from main parameter file  */
			// fpp = fopen_r("parameters/ptv.par");
			// fscanf (fpp, "%d\n", &n_img);
			// fclose (fpp);

			i=0;
			frameCount=0;
			currentFrame=0;
			step_shake=(int)((double)(seq_last-seq_first+1)/(double)max_shake_frames+0.5);
			printf("\nframe step size for camera %d is %d\n", i_img+1, step_shake);
			for (filenumber=seq_first+2; filenumber<seq_last+1-2; filenumber=filenumber+step_shake) { //changed by Beat Feb 08
				// replaced the next lines, ad holten 12-2012 
				//	if (filenumber < 10)		sprintf (filein, "res/rt_is.%1d", filenumber);
				//	else if (filenumber < 100)	sprintf (filein, "res/rt_is.%2d",  filenumber);
				//	else						sprintf (filein, "res/rt_is.%3d", filenumber);
				sprintf (filein, "res/rt_is.%d", filenumber);
				FILEIN = fopen_rp(filein);				// replaced fopen, ad holten, 12-2012
				if (! FILEIN) return TCL_OK;

				//	if (filenumber < 10)		sprintf (filein_ptv, "res/ptv_is.%1d", filenumber);
				//	else if (filenumber < 100)	sprintf (filein_ptv, "res/ptv_is.%2d",	filenumber);
				//	else						sprintf (filein_ptv, "res/ptv_is.%3d", filenumber);

				sprintf (filein_ptv, "res/ptv_is.%d", filenumber);
				FILEIN_ptv = fopen_rp(filein_ptv);		// replaced fopen, ad holten, 12-2012
				if (! FILEIN_ptv) return TCL_OK;

				/* read targets of each camera */
				// to only use quadruplets for shaking that can be linked

				nt4[3][i]=0;
				compose_name_plus_nr_str (seq_name[i_img], "_targets",filenumber, filein_T);
				
				FILEIN_T = fopen_rp(filein_T);		// replaced fopen, ad holten, 12-2012
				if (! FILEIN_T)	return TCL_OK;

				fscanf (FILEIN_T, "%d\n", &nt4[3][i_img]);
				for (j=0; j<nt4[3][i_img]; j++) {
					fscanf (FILEIN_T, "%4d %lf %lf %d %d %d %d %d\n",
						&t4[3][i_img][j].pnr, &t4[3][i_img][j].x,
						&t4[3][i_img][j].y, &t4[3][i_img][j].n ,
						&t4[3][i_img][j].nx ,&t4[3][i_img][j].ny,
						&t4[3][i_img][j].sumg, &t4[3][i_img][j].tnr);
				}
				fclose (FILEIN_T);

				// modified the next code, removed memory bugs, ad holten, 12-2012
				//  fscanf(FILEIN,     "%d\n", &dumy); /* read # of 3D points on dumy */
				//  fscanf(FILEIN_ptv, "%d\n", &dumy); /* read # of 3D points on dumy */
				//	do {
				//		/* read dataset row by row, x,y,z and correspondences */
				//		a[0]=-1;a[1]=-1;a[2]=-1;a[3]=-1;
				//		if (n_img==4) {
				//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
				//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
				//				&a[0], &a[1], &a[2], &a[3]);
				//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
				//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
				//		}
				//		if (n_img==3){
				//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
				//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
				//				&a[0], &a[1], &a[2]);
				//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
				//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
				//		}
				//		if (n_img==2){ // Alex's patch. 24.09.09. Working on Wesleyan data of 2 cameras only
				//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
				//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
				//				&a[0], &a[1]);
				//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
				//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
				//		}

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
					if (( a[i_img]>-1 && next>-1 && prev>-1 && i<max_shake_points && frameCount<max_shake_frames+1  ) ||
						( a[0]>-1 && a[1]>-1 && a[2]>-1 && a[3]>-1 && next>-1 && prev>-1 && filenumber==seq_first+2 ) )
					{
						pix[i_img][i].x   = t4[3][i_img][a[i_img]].x;
						pix[i_img][i].y   = t4[3][i_img][a[i_img]].y;
						pix[i_img][i].pnr = i; 
						fix[i].pnr=i;
						   
						nfix = ++i;
						if (currentFrame < filenumber) {
							currentFrame = filenumber; 
							frameCount++;
						}
					}
				} while (!feof(FILEIN));
				fclose(FILEIN);
	  
			} // end of loop through seq, but loop i_img still open
			if (frameCount==1)
				printf("Using %d linked points of %d frame for camera %d\n", nfix,frameCount, i_img+1);
			else
				printf("Using %d linked points of %d frames for camera %d\n", nfix,frameCount, i_img+1);	

			for (i=0; i<nfix; i++) {
				pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
					imx,imy, pix_x, pix_y,
					&crd[i_img][i].x, &crd[i_img][i].y,
					chfield);
				crd[i_img][i].pnr = pix[i_img][i].pnr;
			}
	  
			/* ================= */

			orient_v3 (interp, Ex[i_img], I[i_img], G[i_img], ap[i_img], mmp,
				nfix, fix, crd[i_img],
				&Ex[i_img], &I[i_img], &G[i_img], &ap[i_img], i_img);

			/* ================= */


			/* save orientation and additional parameters */
			// make safety copy of ori files

			fp1 = fopen( img_ori[i_img], "r" );
			if(fp1 != NULL) {
				fclose(fp1);
				read_ori (&sEx[i_img], &sI[i_img], &sG[i_img], img_ori[i_img]);
				fp1 = fopen (img_addpar0[i_img], "r");
				if (! fp1)	fp1 = fopen ("addpar.raw", "r");

				if (fp1) {
					fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
						&sap[i_img].k1,&sap[i_img].k2,&sap[i_img].k3,
						&sap[i_img].p1,&sap[i_img].p2,
						&sap[i_img].scx,&sap[i_img].she);
					fclose (fp1);
				} 
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
			else{
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

		break;

	case 20: puts ("Orientation from particles, discarding bad 3d-points"); strcpy(buf, "");
		// replaced next lines by for loop, ad holten 12-2012
		//		 strcpy (safety[0], "safety_0");
		//		 strcat (safety[0], ".ori");
		//		 strcpy (safety[1], "safety_1");
		//		 strcat (safety[1], ".ori");
		//		 strcpy (safety[2], "safety_2");
		//		 strcat (safety[2], ".ori");
		//		 strcpy (safety[3], "safety_3");
		//		 strcat (safety[3], ".ori");
		//		 strcpy (safety_addpar[0], "safety_0");
		//		 strcat (safety_addpar[0], ".addpar");
		//		 strcpy (safety_addpar[1], "safety_1");
		//		 strcat (safety_addpar[1], ".addpar");
		//		 strcpy (safety_addpar[2], "safety_2");
		//		 strcat (safety_addpar[2], ".addpar");
		//		 strcpy (safety_addpar[3], "safety_3");
		//		 strcat (safety_addpar[3], ".addpar");
		for (i=0; i<4; i++) {
			sprintf(safety[i], "safety_%d.ori", i);
			sprintf(safety_addpar[i], "safety_%d.addpar", i);
		}

		fpp = fopen_rp("parameters/sequence.par");		// replaced fopen_r, ad holten, 12-2012
		if (! fpp) return TCL_OK;
		for (i=0; i<4; i++)
			fscanf (fpp, "%s\n", seq_name[i]);	   /* name of sequence */
		fclose (fpp);

		fpp = fopen_rp("parameters/shaking.par");		// replaced fopen_r, ad holten, 12-2012
		if (! fpp) return TCL_OK;
		fscanf (fpp,"%d\n", &seq_first); 
		fscanf (fpp,"%d\n", &seq_last);
		fscanf (fpp,"%d\n", &max_shake_points);
		fscanf (fpp,"%d\n", &max_shake_frames);
		fclose (fpp);

		// commented out, ad holten, 12-2012
		//		/*	read from main parameter file  */
		//		fpp = fopen_r ("parameters/ptv.par");
		//		fscanf (fpp, "%d\n", &n_img);
		//		fclose (fpp);
		//			// here I will insert a check for each fix -point. Beat Nov 2011
		//			// The idea is to mark those n of nfix, which are bad, so not to use them in thenext loop
		//			// here I will insert a check for each fix -point. Beat Nov 2011
		//			// The idea is to mark those n of nfix, which are bad, so not to use them in thenext loop
		//			// here I will insert a check for each fix -point. Beat Nov 2011
		//			// The idea is to mark those n of nfix, which are bad, so not to use them in thenext loop

		i = 0;
		frameCount = 0;
		currentFrame = 0;
		step_shake = (int)((double)(seq_last-seq_first+1)/(double)max_shake_frames+0.5);
		printf("\nframe step size for each camera is %d\n", step_shake);
		for (filenumber=seq_first+2; filenumber<seq_last+1-2; filenumber=filenumber+step_shake) { //changed by Beat Feb 08
			// modified the next code, ad holten 12-2012
			//		if (filenumber < 10)		sprintf (filein, "res/rt_is.%1d", filenumber);
			//		else if (filenumber < 100)	sprintf (filein, "res/rt_is.%2d",  filenumber);
			//		else						sprintf (filein, "res/rt_is.%3d", filenumber);

			//		FILEIN = fopen (filein, "r");
			//		if (! FILEIN) printf("Can't open ascii file: %s\n", filein);
			//		/////////open target file(s)!
			//		/* read targets of each camera */

			//		if (filenumber < 10)		sprintf (filein_ptv, "res/ptv_is.%1d", filenumber);
			//		else if (filenumber < 100)	sprintf (filein_ptv, "res/ptv_is.%2d",	filenumber);
			//		else						sprintf (filein_ptv, "res/ptv_is.%3d", filenumber);
			//		
			//		// to only use quadruplets for shaking that can be linked
			//		FILEIN_ptv = fopen (filein_ptv, "r");
			//		if (! FILEIN_ptv) printf("Can't open ascii file: %s\n", filein_ptv);
			// to

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
			////////done reading target files

			// modified the next code, removed memory bugs, ad holten, 12-2012
			//   fscanf(FILEIN, "%*d\n", &dumy); /* read # of 3D points on dumy */
			//   fscanf(FILEIN_ptv, "%d\n", &dumy); /* read # of 3D points on dumy */
			//   do {
			//		/*read dataset row by row, x,y,z and correspondences */
			//		a[0]=-1;a[1]=-1;a[2]=-1;a[3]=-1;
			//		if (n_img==4){
			//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
			//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
			//				&a[0], &a[1], &a[2], &a[3]);
			//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
			//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
			//		}
			//		if (n_img==3){
			//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
			//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
			//				&a[0], &a[1], &a[2]);
			//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
			//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
			//		}
			//		if (n_img==2){ // Alex's patch. 24.09.09. Working on Wesleyan data of 2 cameras only
			//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
			//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,
			//				&a[0], &a[1]);
			//			fscanf(FILEIN_ptv, "%d %d %lf %lf %lf\n",
			//				&prev, &next, &dummy_float, &dummy_float, &dummy_float);
			//		}
				
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
  
		for (i=0; i<nfix ; i++) {
			for (i_img=0; i_img<n_img; i_img++) {
				pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
								imx,imy, pix_x, pix_y,
								&crd[i_img][i].x, &crd[i_img][i].y,
								chfield);
				crd[i_img][i].pnr = pix[i_img][i].pnr;
			}
		}
		for (i=0; i<nfix ; i++) {
			// replaced the next code, ad holten, 12-2012
			//		 ////here comes the actual check!!!!!
			//		 ////here comes the actual check!!!!!
			//		 ////here comes the actual check!!!!!
			//		 ////here comes the actual check!!!!!
			//		 ////here comes the actual check!!!!!
			//		 // new det_lsq function, bloody fast!
			//		 if(crd[0][i].x>-999){
			//			x = crd[0][i].x - I[0].xh;
			//			y = crd[0][i].y - I[0].yh;
			//			//correct_brown_affin (x, y, ap[0], &x, &y);
			//			ray_tracing_v2 (x,y, Ex[0], I[0], G[0], mmp, &X[0], &Y[0], &Z[0], &aa[0], &b[0], &c[0]);
			//		}		 
			//		if(crd[1][i].x>-999) {
			//			x = crd[1][i].x - I[1].xh;
			//			y = crd[1][i].y - I[1].yh;
			//			//correct_brown_affin (x, y, ap[1], &x, &y);
			//			ray_tracing_v2 (x,y, Ex[1], I[1], G[1], mmp, &X[1], &Y[1], &Z[1], &aa[1], &b[1], &c[1]);
			//		}		 
			//		if(crd[2][i].x>-999) {
			//			x = crd[2][i].x - I[2].xh;
			//			y = crd[2][i].y - I[2].yh;
			//			//correct_brown_affin (x, y, ap[2], &x, &y);
			//			ray_tracing_v2 (x,y, Ex[2], I[2], G[2], mmp, &X[2], &Y[2], &Z[2], &aa[2], &b[2], &c[2]);
			//		}		 
			//		if(crd[3][i].x>-999){
			//			x = crd[3][i].x - I[3].xh;
			//			y = crd[3][i].y - I[3].yh;
			//			//correct_brown_affin (x, y, ap[3], &x, &y);
			//			ray_tracing_v2 (x,y, Ex[3], I[3], G[3], mmp, &X[3], &Y[3], &Z[3], &aa[3], &b[3], &c[3]);
			//		}
			// by

			// here comes the actual check!!!!!
			for (j=0; j<4; j++) {
				if (crd[j][i].x>-999) {
					x = crd[j][i].x - I[j].xh;
					y = crd[j][i].y - I[j].yh;
					// correct_brown_affin (x, y, ap[0], &x, &y);
					ray_tracing_v2 (x,y, Ex[j], I[j], G[j], mmp, &X[j], &Y[j], &Z[j], &aa[j], &b[j], &c[j]);
				}
			}

			count_inner=0;
			X_pos=0.; Y_pos=0.; Z_pos=0.;
			for (n=0; n<n_img; n++) {
				for(m=n+1;m<n_img;m++){
					if(crd[n][i].x>-999 && crd[m][i].x>-999){
						mid_point(X[n],Y[n],Z[n],aa[n],b[n],c[n],X[m],Y[m],Z[m],aa[m],b[m],c[m],&dist,&XX,&YY,&ZZ);
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

		break;


	case 11: puts ("Orientation from dumbbells"); strcpy(buf, ""); 

		// replaced next lines by for loop, ad holten 12-2012
		// 		strcpy (safety[0], "safety_0");
		// 		strcat (safety[0], ".ori");
		// 		strcpy (safety[1], "safety_1");
		// 		strcat (safety[1], ".ori");
		// 		strcpy (safety[2], "safety_2");
		// 		strcat (safety[2], ".ori");
		// 		strcpy (safety[3], "safety_3");
		// 		strcat (safety[3], ".ori");
		// 		strcpy (safety_addpar[0], "safety_0");
		// 		strcat (safety_addpar[0], ".addpar");
		// 		strcpy (safety_addpar[1], "safety_1");
		// 		strcat (safety_addpar[1], ".addpar");
		// 		strcpy (safety_addpar[2], "safety_2");
		// 		strcat (safety_addpar[2], ".addpar");
		// 		strcpy (safety_addpar[3], "safety_3");
		// 		strcat (safety_addpar[3], ".addpar");
		for (i=0; i<4; i++) {
			sprintf(safety[i], "safety_%d.ori", i);
			sprintf(safety_addpar[i], "safety_%d.addpar", i);
		}

		for (i_img=0; i_img<n_img; i_img++) {
		
			/* read control point coordinates for man_ori points */
			fpp = fopen_rp("parameters/sequence.par");		// replaced fopen, ad holten, 12-2012
			if (!fpp) return TCL_OK;
	 
			//fscanf (fpp, "%s\n", seq_name[0]);	// changed, ad holten, 12-2012
			//fscanf (fpp, "%s\n", seq_name[1]);
			//fscanf (fpp, "%s\n", seq_name[2]);
			//fscanf (fpp, "%s\n", seq_name[3]);
			for (i=0; i<4; i++)
				fscanf (fpp, "%s\n", seq_name[i]); /* name of sequence */
			fscanf (fpp,"%d\n", &seq_first);
			fscanf (fpp,"%d\n", &seq_last);
			fclose (fpp);
  
			i=0;
			frameCount=0;
			currentFrame=0;
			step_shake=1;
			//printf("\nframe step size for camera %d is %d\n", i_img+1, step_shake);
			for (filenumber=seq_first; filenumber<seq_last+1; filenumber=filenumber+step_shake) { //changed by Beat Feb 08
				sprintf (filein, "res/db_is.%d", filenumber);

				// changed nexte lines, ad holten, 12-2012
				FILEIN = fopen_rp (filein);			// replaced fopen, ad holten, 12-2012
				if (! FILEIN) return TCL_OK;

				/////////open target file(s)!
				/* read targets of each camera */
				nt4[3][i]=0;
				compose_name_plus_nr_str (seq_name[i_img], "_targets",filenumber, filein_T);
				// changed nexte lines, ad holten, 12-2012
				FILEIN_T = fopen_rp (filein_T);		// replaced fopen, ad holten, 12-2012
				if (! FILEIN_T) return TCL_OK;

				fscanf (FILEIN_T, "%d\n", &nt4[3][i_img]);
				for (j=0; j<nt4[3][i_img]; j++){
					fscanf (FILEIN_T, "%4d %lf %lf %d %d %d %d %d\n",
						&t4[3][i_img][j].pnr, &t4[3][i_img][j].x,
						&t4[3][i_img][j].y, &t4[3][i_img][j].n ,
						&t4[3][i_img][j].nx ,&t4[3][i_img][j].ny,
						&t4[3][i_img][j].sumg, &t4[3][i_img][j].tnr);
				}
				fclose (FILEIN_T);
				////////done reading target files

				fscanf(FILEIN, "%*d\n", &dumy); /* read # of 3D points on dumy */
				if (dumy==2) {
					// modified the next code and removed memory bugs, ad holten, 12-2012
					//		/*read the two rows x,y,z and correspondences */
					//		 a1[0]=-1;a1[1]=-1;a1[2]=-1;a1[3]=-1;
					//		 a2[0]=-1;a2[1]=-1;a2[2]=-1;a2[3]=-1;
					//		 if (n_img==4){
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,&a1[0], &a1[1], &a1[2], &a1[3]);
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i+1].x, &fix[i+1].y, &fix[i+1].z,&a2[0], &a2[1], &a2[2], &a2[3]);
					//		 }
					//		 if (n_img==3){
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,&a1[0], &a1[1], &a1[2]);
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i+1].x, &fix[i+1].y, &fix[i+1].z,&a2[0], &a2[1], &a2[2]);
					//		 }
					//		 if (n_img==2){ // Alex's patch. 24.09.09. Working on Wesleyan data of 2 cameras only
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i].x, &fix[i].y, &fix[i].z,&a1[0], &a1[1]);
					//			fscanf(FILEIN, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[i+1].x, &fix[i+1].y, &fix[i+1].z,&a2[0], &a2[1]);
					//		 }

					a1[0] = a1[1] = a1[2] = a1[3] = -1;
					a2[0] = a2[1] = a2[2] = a2[3] = -1;
					fscanf(FILEIN, "%*d %lf %lf %lf %d %d %d %d\n",
						&fix[i].x, &fix[i].y, &fix[i].z,&a1[0], &a1[1], &a1[2], &a1[3]);
					fscanf(FILEIN, "%*d %lf %lf %lf %d %d %d %d\n",
						&fix[i+1].x, &fix[i+1].y, &fix[i+1].z,&a2[0], &a2[1], &a2[2], &a2[3]);

					for (j=3; j>=n_img; j--)
						a1[j] = a2[j] = -1;
					
					// read pixel data according a0,a1,a2,a3 
					pix[i_img][i].x = t4[3][i_img][a1[i_img]].x;
					pix[i_img][i].y = t4[3][i_img][a1[i_img]].y;
					pix[i_img][i].pnr = i; 
					fix[i].pnr = i;
					pix[i_img][i+1].x = t4[3][i_img][a2[i_img]].x;
					pix[i_img][i+1].y = t4[3][i_img][a2[i_img]].y;
					pix[i_img][i+1].pnr = i+1; 
					fix[i+1].pnr = i+1;

					i += 2;
					nfix = i;
				}
				fclose(FILEIN);  
			} //end of loop through seq, but loop i_img still open
			printf("Using %d points for camera %d\n", nfix, i_img+1);
			if (nfix>2) {
				for (i=0; i<nfix ; i++){
					pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
						imx,imy, pix_x, pix_y,
						&crd[i_img][i].x, &crd[i_img][i].y,
						chfield);
					crd[i_img][i].pnr = pix[i_img][i].pnr;
				}

				orient_v3 (interp, Ex[i_img], I[i_img], G[i_img], ap[i_img], mmp,
					nfix, fix, crd[i_img],
					&Ex[i_img], &I[i_img], &G[i_img], &ap[i_img], i_img);

				/* save orientation and additional parameters */
				// make safety copy of ori files

				fp1 = fopen( img_ori[i_img], "r" );
				if(fp1 != NULL) {
					fclose(fp1);
					read_ori (&sEx[i_img], &sI[i_img], &sG[i_img], img_ori[i_img]);
					fp1 = fopen (img_addpar0[i_img], "r");
					if (! fp1)	fp1 = fopen ("addpar.raw", "r");

					if (fp1) {
						fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
							&sap[i_img].k1,&sap[i_img].k2,&sap[i_img].k3,
							&sap[i_img].p1,&sap[i_img].p2,
							&sap[i_img].scx,&sap[i_img].she);
						fclose (fp1);
					} 
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
			} //end if nfix>2
			else {
				success=0;
			}
		} //loop through images

		if (success==1) {
			Tcl_Eval(interp, ".text delete 3");
			Tcl_Eval(interp, ".text delete 1");
			Tcl_Eval(interp, ".text insert 1 \"Orientation from dumbbells \"");
			Tcl_Eval(interp, ".text delete 2");
			Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each view -> \"");
			Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
			Tcl_Eval(interp, ".text insert 3 $tbuf");
		}
		else {
			Tcl_Eval(interp, ".text delete 3");
			Tcl_Eval(interp, ".text delete 1");
			Tcl_Eval(interp, ".text insert 1 \"Orientation from dummbbells \"");
			Tcl_Eval(interp, ".text delete 2");
			Tcl_Eval(interp, ".text insert 2 \"...failed \"");
			Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
			Tcl_Eval(interp, ".text insert 3 $tbuf");
		}
		break;

	case 12: puts ("Orientation from dumbbells"); strcpy(buf, "");
		
		prepare_eval(n_img,&nfix); //goes and looks up what sequence is defined and takes all cord. from rt_is
		// orient_v5 (n_img, nfix, &Ex, &I, &G, &ap);	corrected pointer errors, ad holten
		orient_v5 (n_img, nfix, Ex, I, G, ap);
		  
		for(i_img=0;i_img<n_img;i_img++){
			write_ori (Ex[i_img], I[i_img], G[i_img], img_ori[i_img]);
			fp1 = fopen (img_addpar[i_img], "w");
			fprintf (fp1, "%f %f %f %f %f %f %f",
					 ap[i_img].k1, ap[i_img].k2, ap[i_img].k3,
					 ap[i_img].p1, ap[i_img].p2,
					 ap[i_img].scx, ap[i_img].she);
			fclose (fp1);
		}
		break;
	
	}
	return TCL_OK;
}


int quit_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int i, k;

	for (i=0; i<n_img; i++) {
		free (img[i]);
		free (img0[i]);
	}
	free (zoomimg);

	/* delete unneeded files */
	for (i=0; i<n_img; i++)
		k = remove (img_lp_name[i]);
	return TCL_OK;
}
