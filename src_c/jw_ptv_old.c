/****************************************************************************
*****************************************************************************

Author/Copyright:      	Hans-Gerd Maas / Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	took a longer time ...

Description:	       	target detection, correspondences and
		       	positioning with tclTk display
		       	-> 4 camera version

Routines contained:    	many ...

****************************************************************************/
#include "ptv.h"

#define nmax 5120
#define MAXSLICES 40
#define MAXRESOLUTION 400

/*  global declarations for ptv  */
/*-------------------------------------------------------------------------*/

int	n_img;	       		      	/* no of images */
int	hp_flag=0;           	      	/* flag for highpass */
int	tiff_flag=0;           	      	/* flag for tiff header */
int	chfield;       		       	/* flag for field mode */
int	nfix;	       	       	       	/* no. of control points */
int	num[4];	       		       	/* no. of targets per image */
int numc[4];                        /* no. of targets in current image */
int nump[4];                        /* no. of targets in previous image */
int numn[4];                        /* no. of targets in next image */
int n_trac[4];	           	/* no. of tracks */
int	match=0;		      	/* no. of matches */
int	match2=0;	      	       	/* no. of matches in 2nd pass */
int	nr[4][4];		     	/* point numbers for man. ori */
int	imx, imy, imgsize;	      	/* image size */
int	zoom_x[4],zoom_y[4],zoom_f[4];  /* zoom parameters */
int	pp1=0, pp2=0, pp3=0, pp4=0;   	/* for man. orientation */
int	seq_first, seq_last;	       	/* 1. and last img of seq */
int	demo_nr;		      	/* for demo purposes */
int	examine = 0;		       	/* for more detailed output */
int	dump_for_rdb;		       	/* # of dumpfiles for rdb */
int cr_sz;                          /* size of crosses */
int display;                        /* display flag */
int corp, corc, corn;
int m[4];
int trackallocflag2 = 0;      /* checkflag if mega, c4, t4 already allocated */
int trackallocflag = 0;  
int openfile=0, closefile=0;		// flag to assemble the slices in one file


double	pix_x, pix_y;			      	/* pixel size */
double	ro;			      	        /* 200/pi */
double	cn, cnx, cny, csumg, eps0, corrmin;	/* correspondences par */
double 	rmsX, rmsY, rmsZ, mean_sigma0;		/* a priori rms */
double  X_lay[2], Zmin_lay[2], Zmax_lay[2];	/* illu. layer data */
double  Zmin0_lay[2], Zmax0_lay[2];			/* illu. layer data, t=0*/

FILE	*fp1, *fp2, *fp3, *fp4, *fpp, *fp5;

char	img_name[4][256];      	/* original image names */
char   	img_lp_name[4][256]; 	/* lowpass image names */
char   	img_hp_name[4][256];   	/* highpass image names */
char   	img_cal[4][128];       	/* calibrayion image names */
char   	img_ori[4][128];       	/* image orientation data */
char   	img_ori0[4][128];      	/* orientation approx. values */
char   	img_addpar[4][128];    	/* image additional parameters */
char   	img_addpar0[4][128];   	/* ap approx. values */
char   	seq_name[4][128];      	/* sequence names */
char   	track_dir[128];	       	/* directory with dap track data */
char    fixp_name[128];
char   	res_name[128];	      	/* result destination */
char   	filename[128];	      	/* for general use */
char   	buf[256], val[256];	       	/* buffer */

unsigned char	*img[4];      	/* image data */
unsigned char	*img0[4];      	/* image data for filtering etc */
unsigned char	*zoomimg;     	/* zoom image data */

Exterior       	Ex[4];	      	/* exterior orientation */
Interior       	I[4];	       	/* interior orientation */
ap_52	       	ap[4];	       	/* add. parameters k1,k2,k3,p1,p2,scx,she */
mm_np	       	mmp;	       	/* n-media parameters */
target	       	pix[4][nmax]; 	/* target pixel data */
target	       	pix0[4][4];    	/* pixel data for man_ori points */

target          *t4[4][4];
int             nt4[4][4];

coord_2d       	crd[4][nmax];  	/* (distorted) metric coordinates */
coord_2d       	geo[4][nmax];  	/* corrected metric coordinates */
coord_3d       	fix[4096];     	/* testfield points coordinates */
n_tupel	       	con[nmax];     	/* list of correspondences */

corres	       	*c4[4];
trackparameters tpar;           /* tracking parameters */
scanparameters  scan_par;		// scanning Parameters
mm_LUT	       	mmLUT[4];     	/* LUT for multimedia radial displacement */
coord_3d        *p_c3d;
P *mega[4];

corres          *c4_slice[4][40];
target          *t4_slice[4][4][40];
P				*mega_slice[4][40];
int				m_slice[4][40];
int				nt4_slice[4][4][40];

/***************************************************************************/

int init_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int  i;
  const char *valp;

  puts ("\nMultimedia Particle Positioning and Tracking (Scanning Mode)\n");

  valp = Tcl_GetVar(interp, "examine",  TCL_GLOBAL_ONLY);
  examine = atoi (valp);

  ro = 200/M_PI;

  /*  read from main parameter file  */
  fpp = fopen_r ("parameters/ptv.par");

  fscanf (fpp, "%d\n", &n_img);

  for (i=0; i<4; i++)
    {
      fscanf (fpp, "%s\n", img_name[i]);
      fscanf (fpp, "%s\n", img_cal[i]);
    }
  fscanf (fpp, "%d\n", &hp_flag);
  fscanf (fpp, "%d\n", &tiff_flag);
  fscanf (fpp, "%d\n", &imx);
  fscanf (fpp, "%d\n", &imy);
  fscanf (fpp, "%lf\n", &pix_x);
  fscanf (fpp, "%lf\n", &pix_y);
  fscanf (fpp, "%d\n", &chfield);
  fscanf (fpp, "%lf\n", &mmp.n1);
  fscanf (fpp, "%lf\n", &mmp.n2[0]);
  fscanf (fpp, "%lf\n", &mmp.n3);
  fscanf (fpp, "%lf\n", &mmp.d[0]);
  fclose (fpp);

  /* read illuminated layer data */
  fpp = fopen_r ("parameters/criteria.par");
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

  fpp = fopen_r ("parameters/sequence.par");

  for (i=0; i<4; i++)		fscanf (fpp, "%s\n", seq_name[i]);
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);

  /* initialize zoom parameters and image positions */
  for (i=0; i<n_img; i++)
    {
      num[i] = 0;
      zoom_x[i] = imx/2; zoom_y[i] = imy/2; zoom_f[i] = 1;
    }
  imgsize = imx*imy;

  /* allocate memory for images */
  for (i=0; i<n_img; i++)
    {
      img[i] = (unsigned char *) calloc (imgsize, 1);
      if ( ! img[i])
	{
	  printf ("calloc for img%d --> error\n", i);
	  exit (1);
	}
    }

  for (i=0; i<n_img; i++)
    {
      img0[i] = (unsigned char *) calloc (imgsize, 1);
      if ( ! img0[i])
	{
	  printf ("calloc for img0%d --> error\n", i);
	  exit (1);
	}
    }

  zoomimg = (unsigned char *) calloc (imgsize, 1);
  if ( ! zoomimg)
    {
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
  
  /*  read from main parameter file  */
  fpp = fopen_r ("parameters/ptv.par");

  fscanf (fpp, "%d\n", &n_img);

  for (i=0; i<4; i++)
    {
      fscanf (fpp, "%s\n", img_name[i]);
      fscanf (fpp, "%s\n", img_cal[i]);
    }
  fscanf (fpp, "%d\n", &hp_flag);
  fscanf (fpp, "%d\n", &tiff_flag);
  fscanf (fpp, "%d\n", &imx);
  fscanf (fpp, "%d\n", &imy);
  fscanf (fpp, "%lf\n", &pix_x);
  fscanf (fpp, "%lf\n", &pix_y);
  fscanf (fpp, "%d\n", &chfield);
  fscanf (fpp, "%lf\n", &mmp.n1);
  fscanf (fpp, "%lf\n", &mmp.n2[0]);
  fscanf (fpp, "%lf\n", &mmp.n3);
  fscanf (fpp, "%lf\n", &mmp.d[0]);
  fclose (fpp);

  /*  read scanning parameters from scan parameter file  */
  fpp = fopen_r ("parameters/scan.par");

  fscanf (fpp, "%d\n", &scan_par.N);
  fscanf (fpp, "%lf\n", &scan_par.t);
  fscanf (fpp, "%lf\n", &scan_par.f_pr);
  fscanf (fpp, "%lf\n", &scan_par.f_cam);
  fscanf (fpp, "%d\n", &scan_par.N_slice);
  fscanf (fpp, "%lf\n", &scan_par.s);
  fscanf (fpp, "%lf\n", &scan_par.index);
  fclose (fpp);

 

  /* read illuminated layer data */
  fpp = fopen_r ("parameters/criteria.par");
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

  fpp = fopen_r ("parameters/sequence.par");

  for (i=0; i<4; i++)		fscanf (fpp, "%s\n", seq_name[i]);
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);

  /*  create file names  */
  for (i=0; i<n_img; i++)
    {
      strcpy (img_lp_name[i], img_name[i]); strcat (img_lp_name[i], "_lp");
      strcpy (img_hp_name[i], img_name[i]); strcat (img_hp_name[i], "_hp");
      strcpy (img_ori[i], img_cal[i]);  strcat (img_ori[i], ".ori");
      strcpy (img_addpar[i], img_cal[i]); strcat (img_addpar[i],".addpar");
    }

  /*  read orientation and additional parameters  */
  for (i=0; i<n_img; i++)
    {
      read_ori (&Ex[i], &I[i], img_ori[i]);
      rotation_matrix (Ex[i], Ex[i].dm);

      fp1 = fopen_r (img_addpar[i]);
      fscanf (fp1,"%lf %lf %lf %lf %lf %lf %lf",
	      &ap[i].k1, &ap[i].k2, &ap[i].k3, &ap[i].p1, &ap[i].p2,
	      &ap[i].scx, &ap[i].she);
      fclose (fp1);
    }

  /* read and display original images */
  for (i=0; i<n_img; i++)
    {
      /* reading */
      sprintf(val, "camcanvas %d", i+1);
      Tcl_Eval(interp, val);

      read_image (interp, img_name[i], img[i]);
      sprintf(val, "newimage %d", i+1);

      Tcl_Eval(interp, val);
      sprintf(val, "keepori %d", i+1);
      Tcl_Eval(interp, val);
    }

  if (!trackallocflag2)
    {
      for (i=0; i<4; i++)
	{
	  mega[i]=(P *) calloc(sizeof(P),M);
	  c4[i]=(corres *) calloc(sizeof(corres),M);
	  for (k=0; k<4; k++) { t4[i][k]=(target *) calloc(sizeof (target),M);}
	}
      trackallocflag2=1;
    }

  return TCL_OK;

}

int pre_processing_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i_img, sup;

  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;

  sprintf(val, "Filtering with Highpass");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");

  /* read support of unsharp mask */
  fpp = fopen ("parameters/unsharp_mask.par", "r");
  if ( fpp == 0) { sup = 12;}
  else	{ fscanf (fpp, "%d\n", &sup); fclose (fpp); }

  for (i_img=0; i_img<n_img; i_img++)
    {
      highpass (img_name[i_img], img[i_img], img[i_img], sup, 0, chfield, i_img);

      if (display) {
      img_handle = Tk_FindPhoto( interp, "temp");
      Tk_PhotoGetImage (img_handle, &img_block);
      tclimg2cimg (interp, img[i_img], &img_block);

      sprintf(val, "newimage %d", i_img+1);
      Tcl_GlobalEval(interp, val);
      }
    }

  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");

  return TCL_OK;

}


int detection_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int	       	i, i_img;
  int	       	xmin, pft_version=3;
  char val[256];

  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;

  /* process info */
  sprintf(val, "Detection of Particles");
  Tcl_Eval(interp, ".text delete 2");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text insert 2 $tbuf");

  if (display) {
    for (i_img=0; i_img<n_img; i_img++)
      {
	img_handle = Tk_FindPhoto( interp, "temp");
	Tk_PhotoGetImage (img_handle, &img_block);
	tclimg2cimg (interp, img[i_img], &img_block);
	sprintf(val, "newimage %d", i_img+1);
	Tcl_Eval(interp, val);
      }
  }

  strcpy(val, "");

  /* xmin set to 10 so v_line is not included in detection, in future xmin should
     be set to 0, peakfitting has to be changed too */
  xmin=0;

  /*  read pft version  */
  fpp = fopen ("parameters/pft_version", "r");
  if (fpp)
    {
      fscanf (fpp, "%d\n", &pft_version);
      fclose (fpp);
    }


  /* reset zoom values */
  for (i_img=0; i_img<n_img; i_img++)
    {
      zoom_x[i_img] = imx/2; zoom_y[i_img] = imy/2;  zoom_f[i_img] = 1;
    }

  /*copy images because the target recognition will set greyvalues to 0*/
  for (i_img=0; i_img<n_img; i_img++)
    {
      copy_images (img[i_img], img0[i_img]);
    }

  /* target recognition */
  for (i_img=0; i_img<n_img; i_img++)
    {
      switch (pft_version)
	{
	case 3:	/* pft with profile and distance check */
	  /* newest version */
	  xmin=0; /* vertical line restriction */
	  num[i_img] = peak_fit_new (interp, img[i_img],
				     "parameters/targ_rec.par",
				     xmin, imx, 1, imy, pix[i_img], i_img);
	  break;

	case 0:	/* without peak fitting technique */
	  simple_connectivity (interp, img[i_img], img0[i_img],
			       "parameters/targ_rec.par",
			       xmin, imx, 1, imy, pix[i_img], i_img, &num[i_img]);
	  break;

	case 1:	/* with old (but fast) peak fitting technique */
	  targ_rec (interp, img[i_img], img0[i_img],
		    "parameters/targ_rec.par",
		    xmin, imx, 1, imy, pix[i_img], i_img, &num[i_img]);
	  break;
	}

      sprintf (buf,"%d: %d,  ", i_img+1, num[i_img]);
      strcat(val, buf);

      /* proper sort of targets in y-direction for later binary search */
      /* and for dimitris' tracking */
      quicksort_target_y (pix[i_img], num[i_img]);

      /* reorganize target numbers */
      for (i=0; i<num[i_img]; i++)  pix[i_img][i].pnr = i;
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
  int	i, i_img;
  double x,y;

  puts ("\nTransformation to metric coordinates\n");

  /* rearrange point numbers after manual deletion of points */
  for (i_img=0; i_img<n_img; i_img++)
    for (i=0; i<num[i_img]; i++)  pix[i_img][i].pnr = i;
  /* transformations pixel coordinates -> metric coordinates */
  /* transformations metric coordinates -> corrected metric coordinates */
  for (i_img=0; i_img<n_img; i_img++)
    {
      for (i=0; i<num[i_img]; i++)
	{
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
    {
      quicksort_coord2d_x (geo[i_img], num[i_img]);
    }


  /* init multimedia radial displacement LUTs */
  /* ======================================== */

  if ( !mmp.lut && (mmp.n1 != 1 || mmp.n2[0] != 1 || mmp.n3 != 1))
    {
      puts ("Init multimedia displacement LUTs");
      for (i_img=0; i_img<n_img; i_img++) init_mmLUT(i_img);
      mmp.lut = 1;
    }

  correspondences_4 ( interp);

  /* --------------- */
  /* save pixel coords for tracking */
  for (i_img=0; i_img<n_img; i_img++)
    {
      sprintf (filename, "%s_targets", img_name[i_img]);
      if (openfile)
	  {
		  }
	  fp1 = fopen (filename, "w");
	  
      fprintf(fp1,"%d\n", num[i_img]);
      for (i=0; i<num[i_img]; i++)
	{
	  fprintf (fp1, "%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n", pix[i_img][i].pnr, pix[i_img][i].x,
		   pix[i_img][i].y, pix[i_img][i].n ,
		   pix[i_img][i].nx ,pix[i_img][i].ny,
		   pix[i_img][i].sumg, pix[i_img][i].tnr);
	}
		if (closefile)
		{
			}
		fclose (fp1);	
		
    }

  return TCL_OK;
}


int determination_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int  	i, j, n;
  int  	p[4];
  double  x[4], y[4], X,Y,Z;
  double  Zlo = 1e20, Zhi = -1e20;

  puts ("Determinate");

  sprintf (buf, "Point positioning (l.sq.)");
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");

  fp1 = fopen (res_name, "w");

  if ( ! fp1)
    {
      sprintf(res_name,"res/dt_lsq");
      fp1 = fopen (res_name, "w");
    }
  if ( ! fp1)
    {
      printf ("cannot find dir: res,  data written to dt_lsq in same dir\n");
      sprintf (res_name, "dt_lsq");
      fp1 = fopen (res_name, "w");
    }
  /* create dump file for rdb */
  if (examine == 4)
    {
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

  rmsX = 0; rmsY = 0; rmsZ = 0;	mean_sigma0 = 0;

  for (i=0; i<match; i++)
    {
      for (j=0; j<4; j++)
	if (con[i].p[j] >= 0)	p[j] = geo[j][con[i].p[j]].pnr;
	else		       	p[j] = -1;

      for (j=0, n=0; j<4; j++)
	{
	  if (p[j] > -1)
	    {
	      x[j] = crd[j][p[j]].x;	y[j] = crd[j][p[j]].y;
	      n++;
	    }
	  else
	    {
	      x[j] = -1e10;	y[j] = -1e10;
	      if (p[j] == -2)	n = -100;
	    }
	}

      /* take only points which are matched in all images */
      /* or triplets/quadruplets which result from object model */
      /* e.g.: quad -> n=4; model triplet -> n=3; model pair -> n=2;
	 unrestricted triplet -> n<0; unrestricted pair -> n<0 */
      /*     if (n_img > 2  &&  n < 3)	continue; */

      /* ################################# */
      /* take only points which are matched in all images */
      /* or triplets/quadruplets which result from object model */
      /* e.g.: quad -> n=4; model triplet -> n=3; model pair -> n=2;
	 unrestricted triplet -> n<0; unrestricted pair -> n<0 */
      if ((n_img > 2 && num[0]>64 && num[1]>64 && num[2]>64 && num[3]>64)
	  &&  n < 3)	continue;

/* hack due to problems with approx in det_lsq: */
X = 0.0; Y = 0.0; Z = (Zmin_lay[0]+Zmax_lay[0])/2.0;
for (j=0; j<n_img; j++) { X += Ex[j].x0; Y += Ex[j].y0; }
X /= n_img; Y /= n_img;
      /* ******************************** */


      det_lsq (Ex, I, ap, mmp,
	       x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], &X, &Y, &Z);


      /* write a sequential point number,
	 sumg, if the point was used, and the 3D coordinates */
      fprintf (fp1, "%4d", i+1);


      /*
      if (p[0] > -1)	fprintf (fp1, "  %4d", pix[0][p[0]].sumg);
      else			fprintf (fp1, "   %4d", -1);
      if (p[1] > -1)	fprintf (fp1, "  %4d", pix[1][p[1]].sumg);
      else			fprintf (fp1, "  %4d", -1);
      if (p[2] > -1)	fprintf (fp1, "  %4d", pix[2][p[2]].sumg);
      else			fprintf (fp1, "  %4d", -1);
      if (p[3] > -1)	fprintf (fp1, "  %4d", pix[3][p[3]].sumg);
      else			fprintf (fp1, "  %4d", -1);
      */

      fprintf (fp1, " %9.3f %9.3f %9.3f", X, Y, Z);
      if (p[0] > -1)	fprintf (fp1, " %4d", pix[0][p[0]].pnr);
      else			fprintf (fp1, " %4d", -1);
      if (p[1] > -1)	fprintf (fp1, " %4d", pix[1][p[1]].pnr);
      else			fprintf (fp1, " %4d", -1);
      if (p[2] > -1)	fprintf (fp1, " %4d", pix[2][p[2]].pnr);
      else			fprintf (fp1, " %4d", -1);
      if (p[3] > -1)	fprintf (fp1, " %4d\n", pix[3][p[3]].pnr);
      else			fprintf (fp1, " %4d\n", -1);

      /* write data as new points to dump for rdb */
      if (examine == 4)
	{
	  fprintf (fp2, "%d %10.3f %10.3f %10.3f   %d    ", i, X, Y, Z, 3);
	  for (j=0; j<n_img; j++)
	    if (x[j] != -1e10)
	      fprintf (fp2, "%4d %8.5f %8.5f    ", i, x[j], y[j]);
	    else
	      fprintf (fp2, "%4d %8.5f %8.5f    ", -999, x[j], y[j]);
	  fprintf (fp2, "\n");
	  fclose (fp2);
	}

      if (Z < Zlo)  Zlo = Z;   if (Z > Zhi)  Zhi = Z;
    }

  fclose (fp1);

  rmsX = sqrt(rmsX/match); rmsY = sqrt(rmsY/match); rmsZ = sqrt(rmsZ/match);
  mean_sigma0 = sqrt (mean_sigma0/match);

  sprintf (buf, "Match: %d, => sigma0 = %4.2f micron, RMS = %5.3f/%5.3f/%5.3f mm", match, mean_sigma0*1000, rmsX, rmsY, rmsZ);
  puts (buf);
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");

  /* sort coordinates for binary search in epi line segment drawing */
  for (i=0; i<n_img; i++)  quicksort_coord2d_x (geo[0], num[0]);

  puts ("Determinate done\n");

  return TCL_OK;

}


int sequence_proc_c  (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int     i, j, ok, k;
  int slicepos=0;
  char    seq_ch[10], seq_name[4][128];
  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;
  double slice_step;
  double slicethickness;
  double zdim, z_cen_slice[MAXSLICES+1];
  double zzzz[20];
  
  /* ************************ */
  
  double X1, Y1, Z1, X2, Y2, Z2;
  foundpix p16[16];
  int h, step, step2, count;
  int zaehler, zaehler2, quali, invol, okay;
  int philf[4][4];
  double x1[4], y1[4], x2[4], y2[4], xn[4], yn[4];
  double l_quader;
  double distance;
  double xr[4], xl[4], yd[4], yu[4];
  double Ymin=0, Ymax=0;
  char nn_name[128], img_filename[128];
  int reverse=1;
  double offs[50];
  double deltaoffset; /* due to velocity of lightsheet */

  /* ************************ */
 
  double x,y,z,xx,yy,zz,mx,my,mz;
  double X0,Y0,Z0;
  double V0x,V0y,V0z,acc,deltaT;
  double leftlimit,rightlimit;
  int Nparticles,time,initialtime,finaltime,counter;
  char filen[4][128],resname[128];
  double XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax;
  double a,b,c,X,Y,Z;
  int npixtot,npixx,npixy,sumgray;
  double left,right,maxdist,intersection[40];
  int zero=0,one=1,seqnumber;
  double margin;

  /* ************************ */

  int resolution,n,Value;
  int MaxIndex[MAXSLICES+1],MaxValue[MAXSLICES+1];
  double deltaz,zeta;
  int PointHist[MAXSLICES+1][MAXRESOLUTION+1]={0};

  int summ,arm,NumPoints[MAXSLICES+1];
  double slice_centre[MAXSLICES+1],variance[MAXSLICES+1],sheet_thickness_95[MAXSLICES+1];
  int iter_first,iter_last,iter=2,first,last,yes=0,counter2,dummy;
  /* ************************ */

  
  fpp = fopen_r ("parameters/sequence.par");
  for (i=0; i<4; i++)
  {
  fscanf (fpp, "%s\n", seq_name[i]);     /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  //fscanf (fpp,"%d\n", &iter_first);
  //fscanf (fpp,"%d\n", &iter_last);
  }
  fclose (fpp);

  display = atoi(argv[1]);



  /* scanning ptv ************** */
  printf("\nObject volume is scanned in %d slices!\n", scan_par.N_slice);

  /* read illuminated Volume */
  fpp = fopen_r ("parameters/criteria.par");
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


  Zmin0_lay[0] = Zmin_lay[0];
  Zmin0_lay[1] = Zmin_lay[1];
  Zmax0_lay[0] = Zmax_lay[0];
  Zmax0_lay[1] = Zmax_lay[1];

  mmp.nlay = 1;

if(1)
{
  printf("\n[1] Create slice_centre.par?\n[2] Do Sequence?\n");
  scanf("%d",&iter);


  if (iter==1)
  {
  printf("\nEnter Sequence_first and Sequence_last for iteration. e.g. 1002, 1021\n");
  scanf("%d %d",&iter_first,&iter_last);

	
//printf("\nzdim: %f, max: %f, min: %f, st: %f\n", zdim,Zmax_lay[0], Zmin_lay[0], scan_par.N_slice);

// z_cen_slice[0]=-79.;
/*
speedsheet=zdim*scan_par.f_pr
deltadicke=speedsheet/scan_par.f_cam
deltaoffset=deltadicke/2.
*/

/*   // \commented out on 15/02/2013
   zdim=Zmax_lay[0]-Zmin_lay[0];
   deltaoffset = zdim*scan_par.f_pr/scan_par.f_cam/2.;

   for (j=1; j<scan_par.N_slice+1; j++)
   {
		 sheet_offset(j);

		 offs[j]=scan_par.offset;
		
		 if (reverse==1)
		 {
		 z_cen_slice[j-1] = Zmin_lay[0] + zdim*0.5 - scan_par.offset + deltaoffset;
		 // number of black images: 1, hence j+1
		 }
		 else
		 {
         z_cen_slice[j] = Zmin_lay[0] + zdim*0.5 + scan_par.offset - deltaoffset;
         }			
   }
   //  /commented out on 15/02/2013   */

   //for scanner
   zdim=Zmax_lay[0]-Zmin_lay[0];
   for (j=0; j<scan_par.N_slice; j++)
   {
		 if (reverse==1)
		 {
			 z_cen_slice[j] = Zmax_lay[0] - j*zdim/scan_par.N_slice - scan_par.s/2.;
		 }
		 else
		 {
			 z_cen_slice[j] = Zmin_lay[0] + j*zdim/scan_par.N_slice + scan_par.s/2.;
		 }
   }
   // for scanner

   printf("slicenumber, slicecentre\n");
   for (i=0;i<scan_par.N_slice;i++)
   {
      printf("%d\t %f\n",i,z_cen_slice[i]);
   }
   printf("Slicethickness=\t%f",scan_par.s);

 first=iter_first;
 last=iter_last;
 printf("\nEnter any number to continue\n");
 scanf("%d",&dummy);
 printf("First Iteration...");
 counter2=1;
 }  // end of if iter==1


if (iter!=1)
{
   fpp = fopen ("parameters/slice_centre.par","r");
       for (i=0;i<scan_par.N_slice;i++)
	   {      
          fscanf (fpp, "%lf", &z_cen_slice[i],"\t");
	   }
   fclose (fpp);

   counter2=3;
   first=seq_first;
   last=seq_last;
}

backtosequence:



  for (i=first; i<=last; i++)
  {
	if ((i-first) % (scan_par.N_slice)== 0)
	{
		slicepos = 0;		
		Zmin_lay[0] = Zmin0_lay[0];
		Zmin_lay[1] = Zmin0_lay[1];
		Zmax_lay[0] = Zmax0_lay[0];
		Zmax_lay[1] = Zmax0_lay[1];
	}



	printf("\nstep: %d, zslice[j]: %f, slicepos: %d\n", i,z_cen_slice[slicepos], slicepos );

	Zmax_lay[0]= z_cen_slice[slicepos] - scan_par.s/2.0;
	Zmin_lay[0]= z_cen_slice[slicepos] + scan_par.s/2.0;
	Zmax_lay[1]= z_cen_slice[slicepos] - scan_par.s/2.0;
	Zmin_lay[1]= z_cen_slice[slicepos] + scan_par.s/2.0;
	printf("in sequence zslice[j]: %f, zmin0: %f, zmax0: %f\n",
	z_cen_slice[slicepos], Zmax_lay[0],Zmin_lay[0] );

	slicepos++;
	
	if (i < 10)             sprintf (seq_ch, "00%1d", i);
    else if (i < 100)       sprintf (seq_ch, "0%2d",  i);
    else       sprintf (seq_ch, "%3d",  i);

    for (j=0; j<n_img; j++)
	{
	  sprintf (img_name[j], "%s%s", seq_name[j], seq_ch);
	  sprintf (img_lp_name[j], "%s%s_lp", seq_name[j], seq_ch);
	  sprintf (img_hp_name[j], "%s%s_hp", seq_name[j], seq_ch);
	}

      if (chfield == 0)       sprintf (res_name, "res/rt_is.%s", seq_ch);
      else            sprintf (res_name, "res/rt_is.%s_%1d", seq_ch, chfield);

      sprintf (buf, "\nImages:");
      for (j=0; j<n_img; j++) sprintf (buf, "%s  %s", buf, img_name[j]);
      puts (buf);

      /* calling function for each sequence-n-tupel */
      /* read and display original images */

      for (k=0; k<n_img; k++)
	{
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

      if (hp_flag) {
	pre_processing_c (clientData, interp, argc, argv);
	puts("\nHighpass switched on\n");
      } else { puts("\nHighpass switched off\n"); }
      if (display) {Tcl_Eval(interp, "update idletasks");}
      detection_proc_c (clientData, interp, argc, argv);
      if (display) {Tcl_Eval(interp, "update idletasks");}
      correspondences_proc_c (clientData, interp, argc, argv);
      if (display) {Tcl_Eval(interp, "update idletasks");}
      determination_proc_c (clientData, interp, argc, argv);


      /* delete unneeded files */

      for (j=0; j<n_img; j++)
	{
	  ok = remove (img_lp_name[j]);
	  ok = remove (img_hp_name[j]);
	}
    }
  /* reset of display flag */
  display = 1;



/* ************************* 01 12 2004 ************************************ */


if (counter2==1 || counter2==2) // CALC CENTRE
{
resolution=200;

	for (i=0;i<scan_par.N_slice;i++)
	{
	   for (j=0;j<resolution;j++)
	   {
	   PointHist[i][j]=0;
	   }
	}

    counter2++;
	
    deltaz=(Zmax0_lay[0]-Zmin0_lay[0])/resolution;

    for (i=0;i<scan_par.N_slice;i++)
	{
		//puts(seq_ch);
        read_ascii_data(i+first);

        for (h=0; h<m[3]; h++)
        {
		   if(mega[3][h].x[2]<Zmax0_lay[0] && mega[3][h].x[2]>Zmin0_lay[0])
		   { 
		      n=abs((mega[3][h].x[2]-Zmax0_lay[0])/deltaz);
			  if (n==0) {n=1;}
			  PointHist[i][n]=PointHist[i][n]+1;
		   }
        }
        
        summ=0;
        arm=0;
        for (j=0;j<resolution;j++) //mean
		   {
			  summ=summ+PointHist[i][j];
			  arm=arm+j*PointHist[i][j];
		   }
	    
		n=(int)((double) arm/summ);  
		slice_centre[i]=Zmax0_lay[0]-n*deltaz;
        NumPoints[i]=summ;

		arm=0;
        for (j=0;j<resolution;j++) //variance
		   {			  
			  arm=arm+(j-n)*(j-n)*PointHist[i][j];
		   }

	    n=(int)(1.0/((double) summ -1 )*(double) arm);
        variance[i]=sqrt((double) n)*deltaz;
        sheet_thickness_95[i]=4*variance[i];


	}

    printf("slicenumber, slicecentre, 95%thickness, Number of Points in Slice\n");
    for (i=0;i<scan_par.N_slice;i++)
	{
	printf("%d\t %f\t %f\t %d\n",i,slice_centre[i],sheet_thickness_95[i],NumPoints[i]);
	//replace centres
	z_cen_slice[i]=slice_centre[i];
	}

   // printf("\nCurrent Sheet-thickness: %f",scan_par.s);
    printf("\nSheet-thickness =\t%f",scan_par.s);
	printf("\nReset Sheet-thickness [1]=yes [2]=no?");
	

	scanf("\n%d",&yes);
	
    if (yes==1)
	{
		printf("Type new thickness:\t");
        scanf("%f",&scan_par.s);
	}

    if (counter2==2) printf("Second Iteration...");

	if (counter2==3)
	{

	fpp = fopen ("parameters/slice_centre.par","w");
       for (i=0;i<scan_par.N_slice;i++)
	   {      
          fprintf (fpp, "%3.6f\t", slice_centre[i]);     /* name of sequence */
	   }
    fclose (fpp);

    first=seq_first;
    last=seq_last;

	fpp = fopen ("parameters/scan.par","w");

    fprintf (fpp, "%d\n", scan_par.N);
    fprintf (fpp, "%lf\n", scan_par.t);
    fprintf (fpp, "%lf\n", scan_par.f_pr);
    fprintf (fpp, "%lf\n", scan_par.f_cam);
    fprintf (fpp, "%d\n", scan_par.N_slice);
    fprintf (fpp, "%lf\n", scan_par.s);
    fprintf (fpp, "%lf\n", scan_par.index);
    fclose (fpp);

	printf("\nDo Sequence...");
	}
    display=0;
	goto backtosequence;
}






}   // end of if{0}



printf("\nCreate intersect.par [1]=yes [2]=no?");
scanf("\n%d",&yes);
if (yes==1)// FIND INTERSECTIONS
{     

	resolution=200;
    deltaz=(Zmax0_lay[0]-Zmin0_lay[0])/resolution;

	fp1 = fopen ("parameters/intersect.par","w");

    for (step2=seq_first; step2<seq_last; step2=step2+scan_par.N_slice)
	{

    fprintf(fp1,"%4d\t",step2);

       for (i=0;i<scan_par.N_slice;i++)
	   {
		   i=i+step2;

		   if (i < 10)             sprintf (seq_ch, "00%1d", i);
           else if (i < 100)       sprintf (seq_ch, "0%2d",  i);
           else       sprintf (seq_ch, "%3d",  i);
		   puts(seq_ch);
  
           read_ascii_data(i);

		   i=i-step2;

           for (h=0; h<m[3]; h++)
           {
			   if(mega[3][h].x[2]<Zmax0_lay[0] && mega[3][h].x[2]>Zmin0_lay[0])
			   { 
				   n=abs((mega[3][h].x[2]-Zmax0_lay[0])/deltaz);
				   if (n==0) {n=1;}
				   PointHist[i][n]=PointHist[i][n]+1;
		       }
           }

           MaxValue[i]=0;
		   MaxIndex[i]=0;

		   for (j=0;j<resolution;j++)
		   {
			  Value=PointHist[i][j];
			  if(Value>MaxValue[i])
			  {
			     MaxValue[i]=Value;
			     MaxIndex[i]=j;
			  }			
		   }

           if (MaxValue[i]==0) //no points at all in the slice
		   {
		      MaxIndex[i]=(int)(abs(z_cen_slice[i]-Zmax0_lay[0])/deltaz);
		   }

	   }


       for (i=0;i<scan_par.N_slice-1;i++)
	   {
          j=MaxIndex[i];
          while (PointHist[i][j]>PointHist[i+1][j])
		  {
              j=j+1;
			  if (j>MaxIndex[i+1]) {j=(int)((MaxIndex[i+1]-MaxIndex[i])/2.); break; }
		  }

          zeta=Zmax0_lay[0]-j*deltaz;
		  zeta=zeta/1000.;

          fprintf(fp1,"%3.6f\t",zeta);
        
	   }

    fprintf (fp1, "\n");


	}

fclose(fp1);


}



 //IF{0}

/* ************************* 01 12 2004 ************************************ */




  /* ************************* 07 06 2004 ************************************ */  
//}
printf("\nDo nearest neighbour stuff? [1]=yes [2]=no");
scanf("\n%d",&yes);
if (yes==1)   // nearest neighbour stuff
{

  l_quader=.3;
  maxdist=1000.0; // 0.2
  margin=100.0;//.25;

  fp4 = fopen_r ("parameters/intersect.par");
  //fpp = fopen_r ("parameters/slice_centre.par");

/* (sequence) nearest neighbour cycle */

for (step2=seq_first; step2<seq_last; step2++)
{
 
 if ((step2-seq_first) % (scan_par.N_slice)== 0)
 {
 Zmin_lay[0] = Zmin0_lay[0];
 Zmin_lay[1] = Zmin0_lay[1];
 Zmax_lay[0] = Zmax0_lay[0];
 Zmax_lay[1] = Zmax0_lay[1];


 read_ascii_data(step2);

 fscanf (fp4, "%ld\t", &seqnumber);

 for (i=0;i<scan_par.N_slice-1;i++)
 {
	  if (i==scan_par.N_slice-2)
	  {
	  fscanf (fp4, "%lf\t", &intersection[i]);
	  intersection[i]=intersection[i]*1000.;
	  }
	  else
	  {
	  fscanf (fp4, "%lf\n", &intersection[i]); 
	  intersection[i]=intersection[i]*1000.;
	  }
  }

  //fscanf (fpp, "%ld\t", &seqnumber);

   /*   for (i=0;i<scan_par.N_slice;i++)
	  {
		  if (i==scan_par.N_slice-1)
		  {
		  fscanf (fpp, "%lf\n", &z_cen_slice[i]);
		  z_cen_slice[i]=z_cen_slice[i]*1000.;
		  }
		  else
		  {
		  fscanf (fpp, "%lf\t", &z_cen_slice[i]);
          z_cen_slice[i]=z_cen_slice[i]*1000.;
		  }
	  }*/
      for (i=0;i<scan_par.N_slice;i++)
      {
		  if (i==0)
		  {
			  z_cen_slice[i]=Zmax0_lay[0]+(intersection[i]-Zmax0_lay[0])/2.;
		  }
          if (i==scan_par.N_slice-1)
		  {
			  z_cen_slice[i]=intersection[i-1]+(Zmin0_lay[0]-intersection[i-1])/2.;
		  }
		  if (i>0 && i<scan_par.N_slice-1)
		  {
              z_cen_slice[i]=intersection[i-1]+(intersection[i]-intersection[i-1])/2.;
		  }
	  }

/* sheet cycle */
  for (step=step2; step<(step2+scan_par.N_slice); step++)
  {
  
  count=0;
  
  if (step < 10)             sprintf (seq_ch, "00%1d", step);
  else if (step < 100)       sprintf (seq_ch, "0%2d",  step);
  else       sprintf (seq_ch, "%3d",  step);	  
    
  puts(seq_ch);

  rotate_dataset();
  if (step<seq_last)
  {
  read_ascii_data(step+1);
  }
  else
  {
  break;
  }

  volumedimension (&X_lay[1], &X_lay[0], &Ymax, &Ymin, &Zmax_lay[1], &Zmin_lay[0]);

/*  overlap_z = z_cen_slice[step-100+1] - z_cen_slice[step-100] + scan_par.s */



   /*
   if (step-step2 >0 && step-step2<scan_par.N_slice-1)
   {
     right  = z_cen_slice[step-step2]+(z_cen_slice[step-step2]-z_cen_slice[step-step2-1])/2.;
     left  = z_cen_slice[step-step2]-(z_cen_slice[step-step2+1]-z_cen_slice[step-step2])/2.;
   }
   if (step-step2==0)
   {
      right= z_cen_slice[step-step2]+(z_cen_slice[step-step2+1]-z_cen_slice[step-step2])/2.;
   }
   if (step-step2==scan_par.N_slice-1)
   {
      left = z_cen_slice[step-step2]-(z_cen_slice[step-step2]-z_cen_slice[step-step2-1])/2.;
   }
   */


   if (step-step2>0 && step-step2<scan_par.N_slice-1)
   {
      left = intersection[step-step2-1];
	  right = intersection[step-step2];
   }
   if (step-step2==0)
   {
      right=intersection[step-step2];
	  left=Zmax0_lay[0]+margin;
   }
   if (step-step2==scan_par.N_slice-1)
   {
      left=intersection[step-step2-1];
	  right=Zmin0_lay[0]-margin;
   }


	/* for each particle */  
   for (h=0; h<m[2]; h++)
   {   
   
   


    /* if particle in overlapped volume (assume 50% overlap!) */
    if (mega[2][h].x[2] < z_cen_slice[step-step2])
    {


    X1=Y1=Z1=X2=Y2=Z2=0;
    mega[2][h].inlist=0;

      for (i=0; i<16;i++)
	  {
	  p16[i].ftnr=-1;
	  p16[i].freq=0;
	    for(j=0;j<n_img;j++) p16[i].whichcam[j] =0;
	  }
	
    /* 3D-position */
	X1=mega[2][h].x[0];
	Y1=mega[2][h].x[1];
	Z1=mega[2][h].x[2];

      for (j=0; j<n_img; j++)
	  {
	  img_coord (X1, Y1, Z1, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
	  metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
	  x1[j]=xn[j];
	  y1[j]=yn[j];
	  }

 	  
	tpar.dvxmax=tpar.dvymax=tpar.dvzmax=l_quader;
	tpar.dvxmin=tpar.dvymin=tpar.dvzmin=-l_quader;
    
	searchquader(X1, Y1, Z1, &xr, &xl, &yd, &yu);

    /* search in pix for candidates in next time step */
	  for (j=0;j<n_img;j++)
	  {
	  zaehler = candsearch_in_pix (t4[3][j], nt4[3][j], x1[j], y1[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	    for(k=0; k<4; k++)
		{
		p16[j*4+k].ftnr=t4[3][j][philf[j][k]].tnr;
		if(philf[j][k] != -999) p16[j*4+k].whichcam[j]=1;
		if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
		}
	  }

	/* end of search in pix */

    /* fill and sort candidate struct */
	sortwhatfound(&p16, &zaehler);
    mega[2][h].inlist=zaehler;

	/* take distance as decision criteria */

	  for (j=0;j<zaehler;j++)
	  {
	  X2=mega[3][p16[j].ftnr].x[0];
      Y2=mega[3][p16[j].ftnr].x[1];
	  Z2=mega[3][p16[j].ftnr].x[2];
	  
	  distance=sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1));
	  
      mega[2][h].decis[j]=distance;
      mega[2][h].linkdecis[j]=p16[j].ftnr;
	  }

	 /* end if particle in overlapped volume  */
	 }
	
   
     /* 11.06.04 creating new particle position */
	 /* *************************************************************** */
	if(0)
	{
	   /* reset img coord because of n_img smaller 4 */
	 for (j=0;j<4;j++) { x2[j]=-1e10; y2[j]=-1e10; }
     
	 for (j=0;j<n_img;j++)
	  {
      zaehler2 = candsearch_in_pixrest (t4[3][j], nt4[3][j], x1[j], y1[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
   
        if(zaehler2>0 )
		{
	    x2[j]=t4[3][j][philf[j][0]].x; y2[j]= t4[3][j][philf[j][0]].y;
        }

	  }

      quali=0;

      for (j=0;j<n_img;j++)
	  {
	    if (x2[j] !=-1e10 && y2[j] != -1e10)
		{
	    pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
        }
	  }

      if ( quali >= 2)
	  {

	  X2 = X1; Y2 =Y1; Z2 = Z1;
	  invol=0; okay=0;

	  det_lsq (Ex, I, ap, mmp,
			 x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X2, &Y2, &Z2);

       /* volume check */
	     if ( X_lay[0] < X2 && X2 < X_lay[1] && Ymin < Y2 && Y2 < Ymax &&
		     Zmin_lay[0] < Z2 && Z2 < Zmax_lay[1]) {invol=1;}

		/* displacement check */
		 if ( invol==1 &&
		     tpar.dvxmin < (X1-X2) && (X1-X2) < tpar.dvxmax &&
		     tpar.dvymin < (Y1-Y2) && (Y1-Y2) < tpar.dvymax &&
		     tpar.dvzmin < (Z1-Z2) && (Z1-Z2) < tpar.dvzmax )
		 {   okay=1;
		     /* end displacement check */

		   if (okay == 1)
		   {
		   count++;
		   puts("uh");

		   distance=sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1));
	  
           mega[2][h].decis[mega[2][h].inlist]=distance;
           mega[2][h].linkdecis[mega[2][h].inlist]=m[3];
		   mega[2][h].inlist++;
		   
		   mega[3][m[3]].x[0]=X2;
		   mega[3][m[3]].x[1]=Y2;
		   mega[3][m[3]].x[2]=Z2;
		   mega[3][m[3]].prev= -1;
		   mega[3][m[3]].next= -2;
		   mega[3][m[3]].prio= 2;
		   
		     for (j=0;j<n_img;j++)
			 {
			 c4[3][m[3]].p[j]=-1;
			   if(philf[j][0]!=-999)
			   {
			   t4[3][j][philf[j][0]].tnr=m[3];
			   c4[3][m[3]].p[j]= philf[j][0];
			   c4[3][m[3]].nr=m[3];
			   }
			 }
		   m[3]++;
		   
		   
		   }
		 okay=0;
		 }
	  invol=0;
	  }
	  quali=0;
	}
     /* *************************************************************** */
     /* end 11.06.04 creating new particle position */
   
    }
    /* end for each particle */



    /* sort decis and give preliminary "finaldecis"  */
      for (h=0;h<m[2];h++)
	  {
	    if(mega[2][h].inlist > 0 )
		{
	    sort(mega[2][h].inlist, &mega[2][h].decis, &mega[2][h].linkdecis);
      	mega[2][h].finaldecis=mega[2][h].decis[0];
	    mega[2][h].next=mega[2][h].linkdecis[0];
	    }
	  }
      
	 

	  /* create links with decision check */
      for (h=0;h<m[2];h++)
	{
	  if(mega[2][h].inlist > 0 ) {
	    /* best choice wasn't used yet, so link is created */
	    if ( mega[3][mega[2][h].next].prev == -1) {	mega[3][mega[2][h].next].prev=h;
		mega[3][mega[2][h].next].finaldecis=mega[2][h].finaldecis; }

	    /* best choice was already used by mega[2][mega[1][h].next].prev */
	    else {
	      /* check which is the better choice */
	      if ( mega[2][mega[3][mega[2][h].next].prev].finaldecis > mega[2][h].finaldecis)
		{
		  /*
		    printf("h ist besser, h: %4d jetzt: %5.3f vorher: %d %5.3f\n",
		    h, mega[1][h].finaldecis, mega[1][mega[2][mega[1][h].next].prev].next,
		    mega[1][mega[2][mega[1][h].next].prev].finaldecis); */

		  if (mega[2][mega[3][mega[2][h].next].prev].inlist>1) {
		    /*
		      printf("zweite Wahl fuer %d waere: %5.3f\n",
		      mega[1][mega[2][mega[1][h].next].prev].linkdecis[1],
		      mega[1][mega[2][mega[1][h].next].prev].decis[1]); */ }

		  /* remove link with prev */
		  mega[2][mega[3][mega[2][h].next].prev].next= -2;
		  mega[3][mega[2][h].next].prev=h;
		}
	      else {
		/*
		  printf("h ist schlechter, h: %4d jetzt: %5.3f vorher: %d %5.3f\n",
		  h, mega[1][h].finaldecis, mega[1][mega[2][mega[1][h].next].prev].next,
		  mega[1][mega[2][mega[1][h].next].prev].finaldecis); */
		if (mega[2][h].inlist>1) {
		  /*
		    printf("zweite Wahl fuer %d waere: %5.3f\n",
		    mega[1][h].linkdecis[1],
		    mega[1][h].decis[1]);  */ }

		mega[2][h].next=-2;}
	    }
	  }
	 
	} /* end of creation of links with decision check */
  
	if (chfield == 0)       sprintf (nn_name, "res/nn_is.%s", seq_ch);
    else            sprintf (nn_name, "res/nn_is.%s_%1d", seq_ch, chfield);

    if(0)
	{
	fp1 = fopen (nn_name, "w");
    
	  for (h=0;h<m[2];h++)
	  {
	    if (mega[2][h].next != -2 && mega[2][h].finaldecis<maxdist)
		{
	    fprintf (fp1, "%d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", h, mega[2][h].x[0], mega[2][h].x[1],
		   mega[2][h].x[2], mega[3][mega[2][h].next].x[0], mega[3][mega[2][h].next].x[1],
		   mega[3][mega[2][h].next].x[2], mega[2][h].decis[0]);
		}  
	  }
	
	fclose (fp1);
	}



    
	
	if (chfield == 0)       sprintf (res_name, "res/rt_is.%s", seq_ch);
    else            sprintf (res_name, "res/rt_is.%s_%1d", seq_ch, chfield);

    fp2 = fopen (res_name, "w");
    
	if (! fp2) printf("Can't open ascii file for writing\n");

    fprintf(fp2, "%d\n", m[2]);

     for(i=0; i<m[2]; i++)
	 {
	 if (mega[2][i].x[2] < left+margin && mega[2][i].x[2] > right-margin)
     {
       if (mega[2][i].next < 0 && mega[2][i].prev < 0)
	   {
       fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	     i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	     c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
	   }
	   else
	   {
	      if (mega[2][i].next > -1)
	      {
			  //mega[2][i].prev=-1;
	         if (mega[2][i].finaldecis<maxdist)
		     {
		        if (mega[2][i].x[2]<right)
			    {
                   fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	               i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	               // c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], zero);	               
				   c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
		        }
			    else
			    {
                   fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	               i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	               c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
			    }
			 }
			 else
			 {
                fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	            i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	            c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
			 }
		   
			 
	      } 
          else
	      {
	         if (mega[2][i].finaldecis<maxdist)
		     {
		        if (mega[2][i].x[2]>left)
			    {
                fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	            i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	            // c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], zero);
	            c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
		        }
			    else
			    {
                fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	            i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	            c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
			    }
		     }
			 else
			 {
             fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	         i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	         c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
			 }
	      }
	   }
	 }
	 else
	 {
	 fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d %4d\n",
	 i+1, mega[2][i].x[0], mega[2][i].x[1], mega[2][i].x[2],
	 // c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], zero);
	 c4[2][i].p[0], c4[2][i].p[1], c4[2][i].p[2], c4[2][i].p[3], one);
	 }
     }

     fclose(fp2);


    /* create/update of new targets- and new rt_is-files */
    if(0)
	{
    if (step < 10)             sprintf (seq_ch, "00%1d", step+1);
    else if (step < 100)       sprintf (seq_ch, "0%2d",  step+1);
    else       sprintf (seq_ch, "%3d",  step+1);

    if (count !=0)
	{
     for (j=0; j<n_img; j++)
	 {
	 sprintf (img_name[j], "%s%s", seq_name[j], seq_ch);
	 }

     if (chfield == 0)       sprintf (res_name, "res/rt_is.%s", seq_ch);
     else            sprintf (res_name, "res/rt_is.%s_%1d", seq_ch, chfield);

     if (chfield == 0)       sprintf (nn_name, "res/nn_is.%s", seq_ch);
     else            sprintf (nn_name, "res/nn_is.%s_%1d", seq_ch, chfield);

     fp2 = fopen (res_name, "w");
       if (! fp2) printf("Can't open ascii file for writing\n");

     fprintf(fp2, "%d\n", m[3]);

     for(i=0; i<m[3]; i++)
     {
       fprintf(fp2, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d\n",
	     i+1, mega[3][i].x[0], mega[3][i].x[1], mega[3][i].x[2],
	     c4[3][i].p[0], c4[3][i].p[1], c4[3][i].p[2], c4[3][i].p[3]);
     }

     fclose(fp2);

     for (i=0; i<n_img; i++)
     {
       
	   sprintf (img_filename, "%s_targets", img_name[i]);	 
	   fp3= fopen (img_filename, "w");
       if (! fp3) printf("Can't open ascii file: %s\n", img_name[i]);

       fprintf (fp3, "%d\n", nt4[3][i]);
       for (j=0; j<nt4[3][i]; j++)
	 {
	   fprintf (fp3, "%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n",
	 	  t4[3][i][j].pnr, t4[3][i][j].x,
		  t4[3][i][j].y, t4[3][i][j].n ,
		  t4[3][i][j].nx ,t4[3][i][j].ny,
		  t4[3][i][j].sumg, t4[3][i][j].tnr);
	 }
     fclose (fp3);
    }

  }
	}

	}
 /* end sheet cycle */
 }
}


fclose(fp4);
//fclose(fpp);
}
/* end sequence nearest neighbour cycle */
  /* ************************* end 07 06 2004 *********************************** */


 return TCL_OK;
}


int calibration_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, j, sel, i_img, k, n, sup;
  int intx1, inty1, intx2, inty2;
  coord_2d    	apfig1[11][11];	/* regular grid for ap figures */
  coord_2d     	apfig2[11][11];	/* ap figures */
  coord_3d     	fix4[4];       	/* object points for preorientation */
  coord_2d     	crd0[4][4];    	/* image points for preorientation */
  char	       	filename[256], val[256];
  const char *valp;
  
  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;

  /* read support of unsharp mask */
  fp1 = fopen ("parameters/unsharp_mask.par", "r");
  if (! fp1)	sup = 12;
  else	{ fscanf (fp1, "%d\n", &sup); fclose (fp1); }

  /* Get Selection value from TclTk */

  valp = Tcl_GetVar(interp, "sel",  TCL_GLOBAL_ONLY);
  sel = atoi (valp);

  switch (sel)
    {
    case 1: /*  read calibration parameter file  */
      fp1 = fopen_r ("parameters/cal_ori.par");
      fscanf (fp1,"%s\n", fixp_name);
      for (i=0; i<4; i++)
	{
	  fscanf (fp1, "%s\n", img_name[i]);
	  fscanf (fp1, "%s\n", img_ori0[i]);
	}
      fscanf (fpp, "%d\n", &tiff_flag);
      fscanf (fp1, "%d\n", &chfield);
      fclose (fp1);

      /*  create file names  */
      for (i=0; i<n_img; i++)
	{
	  strcpy (img_ori[i], img_name[i]);
	  strcat (img_ori[i], ".ori");
	  strcpy (img_addpar0[i], img_name[i]);
	  strcat (img_addpar0[i], ".addpar0");
	  strcpy (img_addpar[i], img_name[i]);
	  strcat (img_addpar[i], ".addpar");
	  strcpy (img_hp_name[i], img_name[i]);
	  strcat (img_hp_name[i], "_hp");
	}

      for (i=0; i<n_img; i++)
	{

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
      for (i=0; i<n_img; i++)
	{
	  zoom_x[i] = imx/2; zoom_y[i] = imy/2; zoom_f[i] = 1;
	}

     /* copy images because the target recognition
	 will set greyvalues to zero */

     for (i=0; i<n_img; i++)
	{
	  copy_images (img[i], img0[i]);
	}


      /* target recognition */
      for (i=0; i<n_img; i++)
	{
	  targ_rec (interp, img[i], img0[i], "parameters/detect_plate.par",
		    0, imx, 1, imy, pix[i], i, &num[i]);

	  sprintf (buf,"image %d: %d,  ", i+1, num[i]);
	  strcat(val, buf);

	  if (num[i] > nmax)  exit (1);
	}

      /* save pixel coord as approx. for template matching */
      if (examine)	for (i=0; i<n_img; i++)
	{
	  sprintf (filename, "%s_pix", img_name[i]);
	  fp1 = fopen (filename, "w");
	  for (j=0; j<num[i]; j++)
	    fprintf (fp1, "%4d  %8.3f  %8.3f\n",
		     pix[i][j].pnr, pix[i][j].x, pix[i][j].y);

	  fclose (fp1);
	}

      sprintf(buf,"Number of detected targets, interaction enabled");
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");
      Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 3");
      Tcl_Eval(interp, ".text insert 3 $tbuf");
      break;


    case 3:	pp1=0;	pp2=0;	pp3=0;	pp4=0;

      for (i=0; i<n_img; i++)
	{
	  sprintf (buf, "%d targets remain", num[i]);
	  puts (buf);
	}
      fp1 = fopen_r ("parameters/man_ori.par");
      for (i=0; i<n_img; i++)
	{
	  fscanf (fp1, "%d %d %d %d\n", &nr[i][0], &nr[i][1], &nr[i][2], &nr[i][3]);
	}
      fclose (fp1);

      for (i=0; i<n_img; i++)
	{
	  sprintf(val, "measure %d %d %d %d %d", nr[i][0], nr[i][1], nr[i][2], nr[i][3], i+1);
	  Tcl_Eval(interp, val);
	  valp = Tcl_GetVar(interp, "px0",  TCL_GLOBAL_ONLY);
	  pix0[i][0].x = atoi (valp);
	  valp = Tcl_GetVar(interp, "py0",  TCL_GLOBAL_ONLY);
	  pix0[i][0].y = atoi (valp);
	  valp = Tcl_GetVar(interp, "px1",  TCL_GLOBAL_ONLY);
	  pix0[i][1].x = atoi (valp);
	  valp = Tcl_GetVar(interp, "py1",  TCL_GLOBAL_ONLY);
	  pix0[i][1].y = atoi (valp);
	  valp = Tcl_GetVar(interp, "px2",  TCL_GLOBAL_ONLY);
	  pix0[i][2].x = atoi (valp);
	  valp = Tcl_GetVar(interp, "py2",  TCL_GLOBAL_ONLY);
	  pix0[i][2].y = atoi (valp);
	  valp = Tcl_GetVar(interp, "px3",  TCL_GLOBAL_ONLY);
	  pix0[i][3].x = atoi (valp);
	  valp = Tcl_GetVar(interp, "py3",  TCL_GLOBAL_ONLY);
	  pix0[i][3].y = atoi (valp);
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
      fp1 = fopen_r ("parameters/man_ori.par");
      for (i=0; i<n_img; i++)
	{
	  fscanf (fp1, "%d %d %d %d\n",
		  &nr[i][0], &nr[i][1], &nr[i][2], &nr[i][3]);
	}
      fclose (fp1);

      /* read coordinates of pre-clicked points */
      fp1 = fopen ("man_ori.dat", "r");
      if (! fp1)	break;
      for (i_img=0; i_img<n_img; i_img++)	for (i=0; i<4; i++)
	{
	  fscanf (fp1, "%lf %lf\n",
		  &pix0[i_img][i].x, &pix0[i_img][i].y);
	  drawcross (interp,  (int) pix0[i_img][i].x,
		     (int) pix0[i_img][i].y, cr_sz+2, i_img, "red");
	  draw_pnr (interp, (int) pix0[i_img][i].x, (int) pix0[i_img][i].y,
		    nr[i_img][i], i_img, "red");

	}
      fclose (fp1);

      break;


    case 5: puts ("Sort grid points");
      for (i=0; i<n_img; i++)
	{
	  /* read control point coordinates for man_ori points */
	  fp1 = fopen_r (fixp_name);
	  k = 0;
	  while ( fscanf (fp1, "%d %lf %lf %lf", &fix[k].pnr,
			  &fix[k].x, &fix[k].y, &fix[k].z) != EOF) k++;
	  fclose (fp1);
	  nfix = k;

	  /* take clicked points from control point data set */
	  for (j=0; j<4; j++)	for (k=0; k<nfix; k++)
	    {
	      if (fix[k].pnr == nr[i][j])	fix4[j] = fix[k];
	    }

	  /* get approx for orientation and ap */
	  read_ori (&Ex[i], &I[i], img_ori0[i]);
	  fp1 = fopen (img_addpar0[i], "r");
	  if (! fp1)  fp1 = fopen ("addpar.raw", "r");

	  if (fp1) {
	    fscanf (fp1, "%lf %lf %lf %lf %lf %lf %lf",
		    &ap[i].k1,&ap[i].k2,&ap[i].k3,
		    &ap[i].p1,&ap[i].p2,
		    &ap[i].scx,&ap[i].she);
	    fclose (fp1);} else {
	      printf("no addpar.raw\n");
	      ap[i].k1=ap[i].k2=ap[i].k3=ap[i].p1=ap[i].p2=ap[i].she=0.0;
	      ap[i].scx=1.0;
	    }


	  /* transform clicked points */
	  for (j=0; j<4; j++)
	    {
	      pixel_to_metric (pix0[i][j].x, pix0[i][j].y,
			       imx,imy, pix_x, pix_y,
			       &crd0[i][j].x, &crd0[i][j].y,
			       chfield);
	      correct_brown_affin (crd0[i][j].x, crd0[i][j].y, ap[i],
				   &crd0[i][j].x, &crd0[i][j].y);
	    }

	  /* raw orientation with 4 points */
	  raw_orient (Ex[i], I[i], ap[i], mmp, 4, fix4, crd0[i], &Ex[i]);
	  sprintf (filename, "raw%d.ori", i);
	  write_ori (Ex[i], I[i], filename);
	 
	  /* sorting of detected points by back-projection */
	  sortgrid_man (interp, Ex[i], I[i], ap[i], mmp,
			imx,imy, pix_x,pix_y,
			nfix, fix, num[i], pix[i], chfield, i);

	  /* adapt # of detected points */
	  num[i] = nfix;

	  for (j=0; j<nfix; j++)
	    {
	      if (pix[i][j].pnr < 0)	continue;
	      intx1 = (int) pix[i][j].x ;
	      inty1 = (int) pix[i][j].y ;

	      drawcross (interp, intx1, inty1, cr_sz, i, "white");
	      draw_pnr (interp, intx1, inty1, fix[j].pnr, i, "white");
	    }
	}

      /* dump dataset for rdb */
      if (examine == 4)
	{
	  /* create filename for dumped dataset */
	  sprintf (filename, "dump_for_rdb");
	  fp1 = fopen (filename, "w");

	  /* write # of points to file */
	  fprintf (fp1, "%d\n", nfix);

	  /* write point and image coord to file */
	  for (i=0; i<nfix; i++)
	    {
	      fprintf (fp1, "%4d %10.3f %10.3f %10.3f   %d    ",
		       fix[i].pnr, fix[i].x, fix[i].y, fix[i].z, 0);
	      for (i_img=0; i_img<n_img; i_img++)
		{
		  if (pix[i_img][i].pnr >= 0)
		    {
		      /* transform pixel coord to metric */
		      pixel_to_metric (pix[i_img][i].x,
				       pix[i_img][i].y, imx,imy, pix_x, pix_y,
				       &crd[i_img][i].x, &crd[i_img][i].y,
				       chfield);
		      fprintf (fp1, "%4d %8.5f %8.5f    ",
			       pix[i_img][i].pnr,
			       crd[i_img][i].x, crd[i_img][i].y);
		    }
		  else
		    {
		      fprintf (fp1, "%4d %8.5f %8.5f    ",
			       pix[i_img][i].pnr, 0.0, 0.0);
		    }
		}
	      fprintf (fp1, "\n");
	    }
	  fclose (fp1);
	  printf ("dataset dumped into %s\n", filename);
	}
      break;




    case 6: puts ("Orientation"); strcpy(buf, "");

      for (i_img=0; i_img<n_img; i_img++)
	{
	  for (i=0; i<nfix ; i++)
	    {
	      pixel_to_metric (pix[i_img][i].x, pix[i_img][i].y,
			       imx,imy, pix_x, pix_y,
			       &crd[i_img][i].x, &crd[i_img][i].y,
			       chfield);
	      crd[i_img][i].pnr = pix[i_img][i].pnr;
	    }

	  /* save data for special use of resection routine */
	  if (examine == 4)
	    {
	      printf ("try write resection data to disk\n");
	      /* point coordinates */
	      sprintf (filename, "resect_%s.fix", img_name[i_img]);
	      write_ori (Ex[i_img], I[i_img], img_ori[i_img]);
	      fp1 = fopen (filename, "w");
	      for (i=0; i<nfix; i++)
		fprintf (fp1, "%3d  %10.5f  %10.5f  %10.5f\n",
			 fix[i].pnr, fix[i].x, fix[i].y, fix[i].z);
	      fclose (fp1);

	      /* metric image coordinates */
	      sprintf (filename, "resect_%s.crd", img_name[i_img]);
	      fp1 = fopen (filename, "w");
	      for (i=0; i<nfix; i++)
		fprintf (fp1,
			 "%3d  %9.5f  %9.5f\n", crd[i_img][i].pnr,
			 crd[i_img][i].x, crd[i_img][i].y);
	      fclose (fp1);

	      /* orientation and calibration approx data */
	      write_ori (Ex[i_img], I[i_img], "resect.ori0");
	      fp1 = fopen ("resect.ap0", "w");
	      fprintf (fp1, "%f %f %f %f %f %f %f",
		       ap[i_img].k1, ap[i_img].k2, ap[i_img].k3,
		       ap[i_img].p1, ap[i_img].p2,
		       ap[i_img].scx, ap[i_img].she);
	      fclose (fp1);
	      printf ("resection data written to disk\n");
	    }


	  /* resection routine */
	  /* ================= */

	  if (examine != 4)
	    orient (interp, Ex[i_img], I[i_img], ap[i_img], mmp,
		    nfix, fix, crd[i_img],
		    &Ex[i_img], &I[i_img], &ap[i_img], i_img);

	  /* ================= */


	  /* resection with dumped datasets */
	  if (examine == 4)
	    {

	      printf("Resection with dumped datasets? (y/n)");
	      scanf("%s",buf);
	      if (buf[0] != 'y')	continue;
	      strcpy (buf, "");

	      /* read calibration frame datasets */
	      for (n=0, nfix=0, dump_for_rdb=0; n<100; n++)
		{
		  sprintf (filename, "resect.fix%d", n);
		  fp1 = fopen (filename, "r");
		  if (! fp1)	continue;

		  printf("reading file: %s\n", filename);
		  printf ("reading dumped resect data #%d\n", n);
		  k = 0;
		  while ( fscanf (fp1, "%d %lf %lf %lf",
				  &fix[nfix+k].pnr, &fix[nfix+k].x,
				  &fix[nfix+k].y, &fix[nfix+k].z)
			  != EOF) k++;
		  fclose (fp1);
		  /* read metric image coordinates */
		  sprintf (filename, "resect_%d.crd%d", i_img, n);
		  printf("reading file: %s\n", filename);
		  fp1 = fopen (filename, "r");
		  for (i=nfix; i<nfix+k; i++)
		    fscanf (fp1, "%d %lf %lf",
			    &crd[i_img][i].pnr,
			    &crd[i_img][i].x, &crd[i_img][i].y);
		  nfix += k;
		}

	      /* resection */
	      orient (interp, Ex[i_img], I[i_img], ap[i_img], mmp,
		      nfix, fix, crd[i_img],
		      &Ex[i_img], &I[i_img], &ap[i_img], i_img);
	    }


	  /* save orientation and additional parameters */
	  write_ori (Ex[i_img], I[i_img], img_ori[i_img]);
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
      Tcl_Eval(interp, ".text insert 2 \"...done, sigma0 for each image -> \"");
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text insert 3 $tbuf");

      break;

    case 7: checkpoint_proc (interp);
      sprintf(val,"blue: planimetry,   yellow: height");
      Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 3");
      Tcl_Eval(interp, ".text insert 3 $tbuf");
      break;


    case 8: /* draw additional parameter figures */

      Tcl_Eval(interp, "clearcam");

      /*  read orientation and additional parameters  */
      for (i=0; i<n_img; i++)	read_ori (&Ex[i], &I[i], img_ori[i]);
      for (i=0; i<n_img; i++)
	{
	  fp1 = fopen_r (img_addpar[i]);
	  fscanf (fp1,"%lf %lf %lf %lf %lf %lf %lf",
		  &ap[i].k1, &ap[i].k2, &ap[i].k3,
		  &ap[i].p1, &ap[i].p2, &ap[i].scx, &ap[i].she);
	  fclose (fp1);
	}
      for (i_img=0; i_img<n_img; i_img++)
	{
	  /* create undistorted grid */
	  for (i=0; i<11; i++)	for (j=0; j<11; j++)
	    {
	      apfig1[i][j].x = i * imx/10;
	      apfig1[i][j].y = j * imy/10;
	    }
	  /* draw undistorted grid */
	  for (i=0; i<10; i++)	for (j=0; j<10; j++)
	    {
	      intx1 = (int) apfig1[i][j].x;
	      inty1 = (int) apfig1[i][j].y;
	      intx2 = (int) apfig1[i+1][j].x;
	      inty2 = (int) apfig1[i][j+1].y;
	      drawvector (interp, intx1, inty1, intx2, inty1, 1, i_img, "black");
	      drawvector (interp, intx1, inty1, intx1, inty2, 1, i_img, "black");
	    }
	  for (j=0; j<10; j++)
	    {
	      intx1 = (int) apfig1[10][j].x;
	      inty1 = (int) apfig1[10][j].y;
	      inty2 = (int) apfig1[10][j+1].y;
	      drawvector (interp, intx1, inty1, intx1, inty2, 1, i_img, "black");
	    }
	  for (i=0; i<10; i++)
	    {
	      intx1 = (int) apfig1[i][10].x;
	      inty1 = (int) apfig1[i][10].y;
	      intx2 = (int) apfig1[i+1][10].x;
	      drawvector (interp, intx1, inty1, intx2, inty1, 1, i_img, "black");
	    }
	  /* distort grid */
	  for (i=0; i<11; i++)	for (j=0; j<11; j++)
	    {
	      /* transform to metric, distort and re-transform */
	      pixel_to_metric (apfig1[i][j].x, apfig1[i][j].y,
			       imx,imy, pix_x,pix_y,
			       &apfig2[i][j].x, &apfig2[i][j].y, chfield);
	      distort_brown_affin (apfig2[i][j].x, apfig2[i][j].y,
				   ap[i_img], &apfig2[i][j].x, &apfig2[i][j].y);
	      metric_to_pixel (apfig2[i][j].x, apfig2[i][j].y,
			       imx,imy, pix_x,pix_y,
			       &apfig2[i][j].x, &apfig2[i][j].y, chfield);
	      /* exaggerate distortion by factor 5 */
	      apfig2[i][j].x = 5*apfig2[i][j].x - 4*apfig1[i][j].x;
	      apfig2[i][j].y = 5*apfig2[i][j].y - 4*apfig1[i][j].y;

	    }
	  /* draw distorted grid */
	  for (i=0; i<10; i++)	for (j=0; j<10; j++)
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
	  for (j=0; j<10; j++)
	    {
	      intx1 = (int) apfig2[10][j].x;
	      inty1 = (int) apfig2[10][j].y;
	      intx2 = (int) apfig2[10][j+1].x;
	      inty2 = (int) apfig2[10][j+1].y;
	      drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
	    }
	  for (i=0; i<10; i++)
	    {
	      intx1 = (int) apfig2[i][10].x;
	      inty1 = (int) apfig2[i][10].y;
	      intx2 = (int) apfig2[i+1][10].x;
	      inty2 = (int) apfig2[i+1][10].y ;
	      drawvector (interp, intx1, inty1, intx2, inty2, 3, i_img, "magenta");
	    }
	}

      break;
    }
  return TCL_OK;
}


int quit_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, k;

  for (i=0; i<n_img; i++) { free (img[i]); free (img0[i]); }
  free (zoomimg);

  /* delete unneeded files */
  for (i=0; i<n_img; i++)      	k = remove (img_lp_name[i]);
  return TCL_OK;
}

