/****************************************************************************

Routine:				segmentation.c

Author/Copyright:		Hans-Gerd Maas

Address:		       	Institute of Geodesy and Photogrammetry
				       	ETH - Hoenggerberg
				       	CH - 8093 Zurich

Creation Date:			1988/89

Description:	  target recognition with adaptive threshold
		  based on high pass filter and connectivity analysis
						with peak fitting technique

Routines contained:		unsharp_mask () ...

****************************************************************************/

#include "ptv.h"

void highpass (pic_name, img, img_hp, dim_lp, filter_hp, field, nr)
char	        pic_name[256];  /* image name */
unsigned char  *img;
unsigned char  *img_hp;       	/* highpass filtered image */
int             dim_lp;	       	/* dimension of subtracted lowpass image */
int	        filter_hp;     	/* flag for additional filtering of _hp */
int    		field;	       	/* field to be used */
int	      	nr;	       	/* image number for display */

{
  register int	i;
  FILE	       	*fp;
  unsigned char *img_lp;
  char	       	lp_name[256], hp_name[256];
  register unsigned char *ptr1, *ptr2, *ptr3;
  void 	unsharp_mask ();

  sprintf (lp_name, "%s_lp", pic_name);
  sprintf (hp_name, "%s_hp", pic_name);

  /* allocate memory */

  img_lp = (unsigned char *) calloc (imgsize, 1);
  if ( ! img_lp)
    {
      puts ("calloc for img_lp --> error");
      exit (1);
    }

  unsharp_mask (dim_lp, img, img_lp);

  if (examine == 3)
    {
      fp = fopen (lp_name, "w");  /*  save lowpass image */
      if (tiff_flag) { write_tiff (lp_name, img_lp, imx, imy); }
      else { fwrite (img_lp, 1, imgsize, fp); }
      printf("low pass: %s will be deleted when quitting ptv\n", lp_name);
      fclose (fp);
    }

  /*  subtract lowpass from original  (=> highpass)  */
  for (ptr1=img, ptr2=img_lp, ptr3=img_hp, i=0; i<imgsize;
       ptr1++, ptr2++, ptr3++, i++)
    {
      if (*ptr1 > *ptr2)  *ptr3 = *ptr1 - *ptr2;
      else  *ptr3 = 0;
    }


  /* consider field mode */
  if (field == 1 || field == 2)  split (img_hp, field);


  /* filter highpass image, if wanted */
  switch (filter_hp)
    {
    case 0: break;
    case 1: lowpass_3 (img_hp, img_hp);	break;
    case 2: filter_3 (img_hp, img_hp);	break;
    }

  /* save highpass image for later use */
  if (examine == 3)
    {  fp = fopen (hp_name, "w");
  if (tiff_flag) { write_tiff (hp_name, img_hp, imx, imy); }
  else { fwrite (img_hp, 1, imgsize, fp);}
  fclose (fp);
  }

  free (img_lp);
}

void targ_rec (interp, img0, img, par_file, xmin,xmax,ymin,ymax, pix, nr, num)

Tcl_Interp* interp;
unsigned char	*img, *img0;   	/* image data, image to be set to zero */
char	       	par_file[];    	/* name of parameter file */
int    		xmin,xmax,ymin,ymax;	/* search area */
target	       	pix[];		       	/* pixel coord array, global */
int	       	nr;		       	/* image number for display */
int	       	*num;	       		/* number of detections */

/*  thresholding and center of gravity with a peak fitting technique  */
/*  uses 4 neighbours for connectivity and 8 to find local maxima     */

{
  register int	i, j, m;
  int	      	n=0, n_wait=0, n_targets=0, sumg, sumg_min;
  int	        numpix;
  int	       	thres, gvthres[4], disco;
  int	       	cr_sz;			       /* size of crosses to be drawn */
  int	       	nnmin,nnmax, nx,ny, nxmin,nxmax, nymin,nymax;
  int	       	xa,ya,xb,yb, x4[4],y4[4], xn,yn;
  double	       	x, y;
  FILE	       	*fpp;

  register unsigned char	gv, gvref;

  targpix	       	waitlist[2048];

  /* read image name, filter dimension and threshold from parameter file */
  fpp = fopen_r (par_file);
  fscanf (fpp, "%d", &gvthres[0]);      /* threshold for binarization 1.image */
  fscanf (fpp, "%d", &gvthres[1]);      /* threshold for binarization 2.image */
  fscanf (fpp, "%d", &gvthres[2]);      /* threshold for binarization 3.image */
  fscanf (fpp, "%d", &gvthres[3]);      /* threshold for binarization 4.image */
  fscanf (fpp, "%d", &disco);    	/* threshold value for discontinuity */
  fscanf (fpp, "%d  %d", &nnmin, &nnmax);  /* min. and max. number of  */
  fscanf (fpp, "%d  %d", &nxmin, &nxmax);	 /* pixels per target, */
  fscanf (fpp, "%d  %d", &nymin, &nymax);  /* abs, in x, in y    */
  fscanf (fpp, "%d", &sumg_min);	       	/* min. sumg	       */
  fscanf (fpp, "%d", &cr_sz);    	    	/* size of crosses     */
  fclose (fpp);

  /* give thres value refering to image number */
  thres=gvthres[nr];

  /*  thresholding and connectivity analysis in image  */

  for (i=ymin; i<ymax; i++)  for (j=xmin; j<xmax; j++)
    {
      gv = *(img + i*imx + j);
      if ( gv > thres)
	if (	gv >= *(img + i*imx + j-1)
		&&	gv >= *(img + i*imx + j+1)
		&&	gv >= *(img + (i-1)*imx + j)
		&&	gv >= *(img + (i+1)*imx + j)
		&&	gv >= *(img + (i-1)*imx + j-1)
		&&	gv >= *(img + (i+1)*imx + j-1)
		&&	gv >= *(img + (i-1)*imx + j+1)
		&&	gv >= *(img + (i+1)*imx + j+1) )

	  /* => local maximum, 'peak' */
	  {
	    yn=i;  xn=j;
	    xn = xn;
	    sumg = gv;  *(img + i*imx + j) = 0;
	    xa = xn;  xb = xn;  ya = yn;  yb = yn;
	    gv -= thres;
	    x = (xn) * gv;
	    y = yn * gv;
	    numpix = 1;
	    waitlist[0].x = j;  waitlist[0].y = i;  n_wait = 1;

	    while (n_wait > 0)
	      {
		gvref = *(img0 + imx*(waitlist[0].y) + (waitlist[0].x));

		x4[0] = waitlist[0].x - 1;  y4[0] = waitlist[0].y;
		x4[1] = waitlist[0].x + 1;  y4[1] = waitlist[0].y;
		x4[2] = waitlist[0].x;  y4[2] = waitlist[0].y - 1;
		x4[3] = waitlist[0].x;  y4[3] = waitlist[0].y + 1;

		for (n=0; n<4; n++)
		  {
		    xn = x4[n];  yn = y4[n];
		    xn = xn;
		    gv = *(img + imx*yn + xn);

		    /* conditions for threshold, discontinuity, image borders */
		    /* and peak fitting */
		    if (   (gv > thres)
			   && (xn>=xmin)&&(xn<xmax) && (yn>=ymin)&&(yn<ymax)
			   && (gv <= gvref+disco)
			   && (gvref + disco >= *(img0 + imx*(yn-1) + xn))
			   && (gvref + disco >= *(img0 + imx*(yn+1) + xn))
			   && (gvref + disco >= *(img0 + imx*yn + (xn-1)))
			   && (gvref + disco >= *(img0 + imx*yn + (xn+1)))  )
		      {
			sumg += gv;  *(img + imx*yn + xn) = 0;
			if (xn < xa)	xa = xn;	if (xn > xb)	xb = xn;
			if (yn < ya)	ya = yn;	if (yn > yb)	yb = yn;
			waitlist[n_wait].x = xn;	waitlist[n_wait].y = yn;
			x = x + (xn) * (gv - thres);
			y = y + yn * (gv - thres);
			numpix++;	n_wait++;
		      }
		  }

		n_wait--;
		for (m=0; m<n_wait; m++)  waitlist[m] = waitlist[m+1];
		waitlist[n_wait].x = 0;  waitlist[n_wait].y = 0;

	      }	/*  end of while-loop  */


	    /* check whether target touches image borders */
	    if (xa==xmin || ya==ymin || xb==xmax-1 || yb==ymax-1)	continue;


	    /* get targets extensions in x and y */
	    nx = xb - xa + 1;  ny = yb - ya + 1;

	    if (   (numpix >= nnmin) && (numpix <= nnmax)
		   && (nx >= nxmin) && (nx <= nxmax)
		   && (ny >= nymin) && (ny <= nymax)
		   && (sumg > sumg_min)			 )
	      {
		pix[n_targets].n = numpix;
		pix[n_targets].nx = nx;
		pix[n_targets].ny = ny;
		pix[n_targets].sumg = sumg;
		sumg -= (numpix*thres);
		x /= sumg;	x += 0.5;	y /= sumg;	y += 0.5;
		pix[n_targets].x = x;
		pix[n_targets].y = y;
		pix[n_targets].pnr = n_targets++;

		xn = x;  yn = y;
		drawcross (interp, (int) xn, (int) yn, cr_sz, nr, "Blue");

	      }
	  }	/*  end of if-loop  */
    }
  *num = n_targets;
}





void simple_connectivity (interp, img0, img, par_file, xmin,xmax,ymin,ymax, pix, nr, num)
Tcl_Interp* interp;
unsigned char	*img, *img0;   	/* image data, image to be set to zero */
char	       	par_file[];    	/* name of parameter file */
int	       	xmin,xmax,ymin,ymax;	/* search area */
target	       	pix[];        	/* pixel coord array, global */
int    	       	nr;    	       	/* image number for display */
int	       	*num;	       	/* number of detections */

/*  thresholding and center of gravity with a peak fitting technique  */
/*  uses 4 neighbours for connectivity and 8 to find local maxima     */

{
  register int	i, j, m;
  int         	n=0, n_wait=0, n_targets=0, sumg, sumg_min;
  int	       	numpix;
  int	       	thres;
  int	      	cr_sz; 	      	/* size of crosses to be drawn */
  int          	nnmin,nnmax, nx,ny, nxmin,nxmax, nymin,nymax;
  int	       	xa,ya,xb,yb, x4[4],y4[4], xn,yn;
  double       	x, y;
  FILE	       	*fpp;

  register unsigned char  gv, gvref;

  targpix	waitlist[2048];

  /* read image name, threshold and shape limits from parameter file */
  fpp = fopen_r (par_file);
  fscanf (fpp, "%d", &thres);    	/* threshold for binarization  	 */
  fscanf (fpp, "%d", &n);		/* threshold value for discontinuity */
  fscanf (fpp, "%d  %d", &nnmin, &nnmax);	 /* min. and max. number of  */
  fscanf (fpp, "%d  %d", &nxmin, &nxmax);	 /* pixels per target, 	*/
  fscanf (fpp, "%d  %d", &nymin, &nymax);	 /* abs, in x, in y   	*/
  fscanf (fpp, "%d", &sumg_min);	       	 /* min. sumg 		*/
  fscanf (fpp, "%d", &cr_sz);		 /* size of crosses	*/
  fclose (fpp);


  /*  thresholding and connectivity analysis in image  */

  for (i=ymin; i<ymax; i++)  for (j=xmin; j<xmax; j++)
    {
      gv = *(img + i*imx + j);
      if ( gv > 2*thres)
	{
	  yn=i;  xn=j;
	  xn = xn;
	  sumg = gv;  *(img + i*imx + j) = 0;
	  xa = xn;  xb = xn;  ya = yn;  yb = yn;
	  gv -= thres;
	  x = (xn) * gv;
	  y = yn * gv;
	  numpix = 1;
	  waitlist[0].x = j;  waitlist[0].y = i;  n_wait = 1;

	  while (n_wait > 0)
	    {
	      gvref = *(img0 + imx*(waitlist[0].y) + (waitlist[0].x));

	      x4[0] = waitlist[0].x - 1;  y4[0] = waitlist[0].y;
	      x4[1] = waitlist[0].x + 1;  y4[1] = waitlist[0].y;
	      x4[2] = waitlist[0].x;  y4[2] = waitlist[0].y - 1;
	      x4[3] = waitlist[0].x;  y4[3] = waitlist[0].y + 1;

	      for (n=0; n<4; n++)
		{
		  xn = x4[n];  yn = y4[n];
		  xn = xn;
		  gv = *(img + imx*yn + xn);

		  /* conditions for threshold, discontinuity, image borders */
		  /* and peak fitting */
		  if (   (gv > thres)
			 && (xn>=xmin)&&(xn<xmax) && (yn>=ymin)&&(yn<ymax)  )
		    {
		      sumg += gv;  *(img + imx*yn + xn) = 0;
		      if (xn < xa)	xa = xn;	if (xn > xb)	xb = xn;
		      if (yn < ya)	ya = yn;	if (yn > yb)	yb = yn;
		      waitlist[n_wait].x = xn;	waitlist[n_wait].y = yn;
		      x = x + (xn ) * (gv - thres);
		      y = y + yn * (gv - thres);
		      numpix++;	n_wait++;
		    }
		}

	      n_wait--;
	      for (m=0; m<n_wait; m++)  waitlist[m] = waitlist[m+1];
	      waitlist[n_wait].x = 0;  waitlist[n_wait].y = 0;

	    }	/*  end of while-loop  */

	  nx = xb - xa + 1;  ny = yb - ya + 1;

	  if (   (numpix >= nnmin) && (numpix <= nnmax)
		 && (nx >= nxmin) && (nx <= nxmax)
		 && (ny >= nymin) && (ny <= nymax)
		 && (sumg > sumg_min)			 )
	    {
	      pix[n_targets].n = numpix;
	      pix[n_targets].nx = nx;
	      pix[n_targets].ny = ny;
	      pix[n_targets].sumg = sumg;
	      sumg -= (numpix*thres);
	      x /= sumg;	x += 0.5;	y /= sumg;	y += 0.5;
	      pix[n_targets].x = x;
	      pix[n_targets].y = y;
	      pix[n_targets].tnr = -1;
	      pix[n_targets].pnr = n_targets++;

	      xn =  x;  yn = y;
	      drawcross (xn, yn, cr_sz, 8);
	    }
	}	/*  end of if-loop  */
    }

  *num = n_targets;
}


