/*******************************************************************

Routine:      	orientation.c

Author/Copyright:      	Hans-Gerd Maas

Address:      	Institute of Geodesy and Photogrammetry
                ETH - Hoenggerberg
	       	CH - 8093 Zurich

Creation Date:	summer 1988

Description:   	computes 6 parameters of exterior orientation
	     	and different sets of additional parameters
	      	from given approximate values, fixpoint- and image
	    	coordinates in ASCII files

Functiones used:   	img_affin, rotation_matrix, adjlib.a

Routines contained:

Related routines:

******************************************************************/
#include "ptv.h"

void orient (interp, Ex0, I0, ap0, mm, nfix, fix, crd, Ex, I, ap, nr)
Tcl_Interp*     interp;
Exterior	Ex0, *Ex;	/* exterior orientation, approx and result */
Interior	I0, *I;		/* interior orientation, approx and result */
ap_52		ap0, *ap;	/* add. parameters, approx and result */
mm_np		mm;    		/* multimedia parameters */
int	       	nfix;		/* # of object points */
coord_3d	fix[];		/* object point data */
coord_2d	crd[];		/* image coordinates */
int	       	nr;  		/* image number for residual display */
{
  int  	i,j,n, itnum, stopflag, n_obs=0;
  int  	useflag, ccflag, scxflag, sheflag, xhflag, yhflag,
    k1flag, k2flag, k3flag, p1flag, p2flag;
  int  	intx1, intx2, inty1, inty2;
  double       	dm = 0.00001,  drad = 0.0000001;
  double       	X[1800][16], Xh[1800][16], y[1800], yh[1800], ident[10],
    XPX[16][16], XPy[16], beta[16], Xbeta[1800],
    resi[1800], omega=0, sigma0, sigmabeta[16],
    P[1800], p, sumP, pixnr[3600];
  double 	Xp, Yp, Zp, xp, yp, xpd, ypd, r, qq;
  FILE 	*fp1;
  int dummy, multi;


  /* read, which parameters shall be used */
  fp1 = fopen_r ("parameters/orient.par");
  fscanf (fp1,"%d", &useflag);
  fscanf (fp1,"%d", &ccflag);
  fscanf (fp1,"%d", &xhflag);
  fscanf (fp1,"%d", &yhflag);
  fscanf (fp1,"%d", &k1flag);
  fscanf (fp1,"%d", &k2flag);
  fscanf (fp1,"%d", &k3flag);
  fscanf (fp1,"%d", &p1flag);
  fscanf (fp1,"%d", &p2flag);
  fscanf (fp1,"%d", &scxflag);
  fscanf (fp1,"%d", &sheflag);
  fclose (fp1);


  /* init X, y (set to zero) */
  for (i=0; i<1800; i++)
    {
      for (j=0; j<16; j++)  X[i][j] = 0;
      y[i] = 0;  P[i] = 0;
    }

  /* init identities */
  ident[0] = I0.cc;  ident[1] = I0.xh;  ident[2] = I0.yh;
  ident[3]=ap0.k1; ident[4]=ap0.k2; ident[5]=ap0.k3;
  ident[6]=ap0.p1; ident[7]=ap0.p2;
  ident[8]=ap0.scx; ident[9]=ap0.she;

  /* main loop, program runs through it, until none of the beta values
     comes over a threshold and no more points are thrown out
     because of their residuals */

  puts ("\n\nbegin of iterations");
  itnum = 0;  stopflag = 0;
  while ((stopflag == 0) && (itnum < 20))
    {
      printf ("\n\n%2d. iteration\n", ++itnum);

      for (i=0, n=0; i<nfix; i++)  if (crd[i].pnr == fix[i].pnr)
	{
	  /* use only certain points as control points */
	  switch (useflag)
	    {
	    case 1: if ((fix[i].pnr % 2) == 0)  continue;  break;
	    case 2: if ((fix[i].pnr % 2) != 0)  continue;  break;
	    case 3: if ((fix[i].pnr % 3) == 0)  continue;  break;
	    }

	  /* check for correct correspondence */
	  if (crd[i].pnr != fix[i].pnr)	continue;

	  /* do not use the corner points of plate 85 */
	  /*if (nfix == 85  &&  fix[i].pnr == 1)	continue;
	  if (nfix == 85  &&  fix[i].pnr == 7)	continue;
	  if (nfix == 85  &&  fix[i].pnr == 43)	continue;
	  if (nfix == 85  &&  fix[i].pnr == 49)	continue;*/


	  pixnr[n/2] = i;		/* for drawing residuals */
	  Xp = fix[i].x;  Yp = fix[i].y;  Zp = fix[i].z;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xp,&yp);


	  /* derivatives of add. parameters */

	  r = sqrt (xp*xp + yp*yp);

	  X[n][7] = ap0.scx;
	  X[n+1][7] = sin(ap0.she);

	  X[n][8] = 0;
	  X[n+1][8] = 1;

	  X[n][9] = ap0.scx * xp * r*r;
	  X[n+1][9] = yp * r*r;

	  X[n][10] = ap0.scx * xp * pow(r,4.0);
	  X[n+1][10] = yp * pow(r,4.0);

	  X[n][11] = ap0.scx * xp * pow(r,6.0);
	  X[n+1][11] = yp * pow(r,6.0);

	  X[n][12] = ap0.scx * (2*xp*xp + r*r);
	  X[n+1][12] = 2 * xp * yp;

	  X[n][13] = 2 * ap0.scx * xp * yp;
	  X[n+1][13] = 2*yp*yp + r*r;

	  qq =  ap0.k1*r*r; qq += ap0.k2*pow(r,4.0);
	  qq += ap0.k3*pow(r,6.0);
	  qq += 1;
	  X[n][14] = xp * qq + ap0.p1 * (r*r + 2*xp*xp) + 2*ap0.p2*xp*yp;
	  X[n+1][14] = 0;

	  X[n][15] = -cos(ap0.she) * yp;
	  X[n+1][15] = -sin(ap0.she) * yp;



	  /* numeric derivatives */

	  Ex0.x0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][0]      = (xpd - xp) / dm;
	  X[n+1][0] = (ypd - yp) / dm;
	  Ex0.x0 -= dm;

	  Ex0.y0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][1]      = (xpd - xp) / dm;
	  X[n+1][1] = (ypd - yp) / dm;
	  Ex0.y0 -= dm;

	  Ex0.z0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][2]      = (xpd - xp) / dm;
	  X[n+1][2] = (ypd - yp) / dm;
	  Ex0.z0 -= dm;

	  Ex0.omega += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][3]      = (xpd - xp) / drad;
	  X[n+1][3] = (ypd - yp) / drad;
	  Ex0.omega -= drad;

	  Ex0.phi += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][4]      = (xpd - xp) / drad;
	  X[n+1][4] = (ypd - yp) / drad;
	  Ex0.phi -= drad;

	  Ex0.kappa += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][5]      = (xpd - xp) / drad;
	  X[n+1][5] = (ypd - yp) / drad;
	  Ex0.kappa -= drad;

	  I0.cc += dm;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I0, ap0, mm, &xpd,&ypd);
	  X[n][6]      = (xpd - xp) / dm;
	  X[n+1][6] = (ypd - yp) / dm;
	  I0.cc -= dm;

	  y[n]   = crd[i].x - xp;
	  y[n+1] = crd[i].y - yp;

	  n += 2;
	}
      n_obs = n;


      /* identities */

      for (i=0; i<10; i++)  X[n_obs+i][6+i] = 1;

      y[n_obs+0] = ident[0] - I0.cc;		y[n_obs+1] = ident[1] - I0.xh;
      y[n_obs+2] = ident[2] - I0.yh;		y[n_obs+3] = ident[3] - ap0.k1;
      y[n_obs+4] = ident[4] - ap0.k2;		y[n_obs+5] = ident[5] - ap0.k3;
      y[n_obs+6] = ident[6] - ap0.p1;		y[n_obs+7] = ident[7] - ap0.p2;
      y[n_obs+8] = ident[8] - ap0.scx;	y[n_obs+9] = ident[9] - ap0.she;



      /* weights */
      for (i=0; i<n_obs; i++)  P[i] = 1;
      if ( ! ccflag)  P[n_obs+0] = 1e20;
      if ( ! xhflag)  P[n_obs+1] = 1e20;
      if ( ! yhflag)  P[n_obs+2] = 1e20;
      if ( ! k1flag)  P[n_obs+3] = 1e20;
      if ( ! k2flag)  P[n_obs+4] = 1e20;
      if ( ! k3flag)  P[n_obs+5] = 1e20;
      if ( ! p1flag)  P[n_obs+6] = 1e20;
      if ( ! p2flag)  P[n_obs+7] = 1e20;
      if ( ! scxflag) P[n_obs+8] = 1e20;
      if ( ! sheflag) P[n_obs+9] = 1e20;


      n_obs += 10;  sumP = 0;
      for (i=0; i<n_obs; i++)	       	/* homogenize */
	{
	  p = sqrt (P[i]);
	  for (j=0; j<16; j++)  Xh[i][j] = p * X[i][j];
	  yh[i] = p * y[i];  sumP += P[i];
	}



      /* Gauss Markoff Model */

      ata (Xh, XPX, n_obs, 16);

      matinv (XPX, 16);

      atl (XPy, Xh, yh, n_obs, 16);

      matmul (beta, XPX, XPy, 16,16,1);

      stopflag = 1;
      puts ("\n==> beta :\n");
      for (i=0; i<16; i++)
	{
	  printf ("%10.6f  ",beta[i]);
	  if (fabs (beta[i]) > 0.01)  stopflag = 0;	/* more iterations */
	}
      printf ("\n\n");

      Ex0.x0 += beta[0];  Ex0.y0 += beta[1];  Ex0.z0 += beta[2];
      Ex0.omega += beta[3];  Ex0.phi += beta[4];  Ex0.kappa += beta[5];
      I0.cc += beta[6];  I0.xh += beta[7];  I0.yh += beta[8];
      ap0.k1 += beta[9];  ap0.k2 += beta[10];  ap0.k3 += beta[11];
      ap0.p1 += beta[12];  ap0.p2 += beta[13];
      ap0.scx += beta[14];  ap0.she += beta[15];
    }



  /* compute residuals etc. */

  matmul ( Xbeta, X, beta, n_obs, 16, 1);

  omega = 0;
  for (i=0; i<n_obs; i++)
    {
      resi[i] = Xbeta[i] - y[i];  omega += resi[i] * P[i] * resi[i];
    }
  sigma0 = sqrt (omega / (n_obs - 16));

  for (i=0; i<16; i++)  sigmabeta[i] = sigma0 * sqrt(XPX[i][i]);


  /* correlations between parameters */
  if (examine)	for (i=0; i<16; i++)
    {
      for (j=0; j<16; j++)
	printf ("%6.2f",
		XPX[i][j] / (sqrt(XPX[i][i]) * sqrt(XPX[j][j])));
      printf ("\n");
    }


  /* print results */
  printf ("\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||");
  printf ("\n\nResults after %d iterations:\n\n", itnum);
  printf ("sigma0 = %6.2f micron\n", sigma0*1000);
  printf ("X0 =    %8.3f   +/- %8.3f\n", Ex0.x0, sigmabeta[0]);
  printf ("Y0 =    %8.3f   +/- %8.3f\n", Ex0.y0, sigmabeta[1]);
  printf ("Z0 =    %8.3f   +/- %8.3f\n", Ex0.z0, sigmabeta[2]);
  printf ("omega = %8.4f   +/- %8.4f\n", Ex0.omega*ro, sigmabeta[3]*ro);
  printf ("phi   = %8.4f   +/- %8.4f\n", Ex0.phi*ro, sigmabeta[4]*ro);
  printf ("kappa = %8.4f   +/- %8.4f\n", Ex0.kappa*ro, sigmabeta[5]*ro);
  printf ("camera const  = %8.5f   +/- %8.5f\n", I0.cc, sigmabeta[6]);
  printf ("xh            = %8.5f   +/- %8.5f\n", I0.xh, sigmabeta[7]);
  printf ("yh            = %8.5f   +/- %8.5f\n", I0.yh, sigmabeta[8]);
  printf ("k1            = %8.5f   +/- %8.5f\n", ap0.k1, sigmabeta[9]);
  printf ("k2            = %8.5f   +/- %8.5f\n", ap0.k2, sigmabeta[10]);
  printf ("k3            = %8.5f   +/- %8.5f\n", ap0.k3, sigmabeta[11]);
  printf ("p1            = %8.5f   +/- %8.5f\n", ap0.p1, sigmabeta[12]);
  printf ("p2            = %8.5f   +/- %8.5f\n", ap0.p2, sigmabeta[13]);
  printf ("scale for x'  = %8.5f   +/- %8.5f\n", ap0.scx, sigmabeta[14]);
  printf ("shearing      = %8.5f   +/- %8.5f\n", ap0.she*ro, sigmabeta[15]*ro);


  fp1 = fopen_r ("parameters/examine.par");
  fscanf (fp1,"%d\n", &dummy);
  fscanf (fp1,"%d\n", &multi);
  fclose (fp1);
  if (dummy==1){
      examine=4;
  }
  else{
      examine=0;
  }
  

  /* show original images with residual vectors (requires globals) */
  sprintf (val, "%d: %5.2f micron, ", nr+1, sigma0*1000);
  strcat(buf,val);

//if(multi==0){
  read_image (interp, img_name[nr], img[nr]);
  sprintf(val, "newimage %d", nr+1);
  Tcl_Eval(interp, val);
//}

  puts (buf);

//if(multi==0){
  for (i=0; i<n_obs-10; i+=2)
    {
      n = pixnr[i/2];
      intx1 = (int) pix[nr][n].x;
      inty1 = (int) pix[nr][n].y;
      intx2 = intx1 + resi[i]*5000;
      inty2 = inty1 + resi[i+1]*5000;

      drawcross (interp, intx1, inty1, 3, nr , "orange");
      drawvector (interp, intx1, inty1, intx2, inty2, 1, nr , "red");
    }
//}



  if (stopflag)
    {
      rotation_matrix (Ex0, Ex0.dm);
      *Ex = Ex0;	*I = I0;	*ap = ap0;
    }
  else	puts ("orientation does not converge");
}




void raw_orient (Ex0, I, ap, mm, nfix, fix, crd, Ex,only_show)
Exterior  Ex0, *Ex;
Interior  I;
ap_52	  ap;
mm_np	  mm;
int	  nfix;
coord_3d  fix[];
coord_2d  crd[];
{
  double		X[8][6], y[8],
    XPX[6][6], XPy[6], beta[6];
  double 		Xp, Yp, Zp, xp, yp, xpd, ypd;
  int     	i,j,n, itnum, stopflag, n_obs=0;
  double		dm = 0.0001,  drad = 0.000001;
 

  /* init X, y (set to zero) */
  for (i=0; i<8; i++)
    {
      for (j=0; j<6; j++)  X[i][j] = 0;
      y[i] = 0;
    }

  ap.k1 = 0;	ap.k2 = 0;	ap.k3 = 0;	ap.p1 = 0;	ap.p2 = 0;
  ap.scx = 1; ap.she = 0;


  /* main loop, program runs through it, until none of the beta values
     comes over a threshold and no more points are thrown out
     because of their residuals */

  itnum = 0;  stopflag = 0;

///////////make a menu so one see the raw guess!!!!!
  if(only_show==1) stopflag=1;
/////// Beat Lüthi 9. Mai 2007

  while ((stopflag == 0) && (itnum < 20))
    {
      ++itnum;

      for (i=0, n=0; i<nfix; i++)  if (crd[i].x != -999)
	{
	  Xp = fix[i].x;  Yp = fix[i].y;  Zp = fix[i].z;
	  rotation_matrix (Ex0, Ex0.dm);

	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xp,&yp);

	  /* numeric derivatives */

	  Ex0.x0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][0]      = (xpd - xp) / dm;
	  X[n+1][0] = (ypd - yp) / dm;
	  Ex0.x0 -= dm;

	  Ex0.y0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][1]	  = (xpd - xp) / dm;
	  X[n+1][1] = (ypd - yp) / dm;
	  Ex0.y0 -= dm;

	  Ex0.z0 += dm;
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][2]	  = (xpd - xp) / dm;
	  X[n+1][2] = (ypd - yp) / dm;
	  Ex0.z0 -= dm;

	  Ex0.omega += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][3]	  = (xpd - xp) / drad;
	  X[n+1][3] = (ypd - yp) / drad;
	  Ex0.omega -= drad;

	  Ex0.phi += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][4]	  = (xpd - xp) / drad;
	  X[n+1][4] = (ypd - yp) / drad;
	  Ex0.phi -= drad;

	  Ex0.kappa += drad;
	  rotation_matrix (Ex0, Ex0.dm);
	  img_coord (Xp,Yp,Zp, Ex0,I, ap, mm, &xpd,&ypd);
	  X[n][5]	  = (xpd - xp) / drad;
	  X[n+1][5] = (ypd - yp) / drad;
	  Ex0.kappa -= drad;

	  y[n]   = crd[i].x - xp;
	  y[n+1] = crd[i].y - yp;

	  n += 2;
	}
      n_obs = n;

      /* Gauss Markoff Model */

      ata (X, XPX, n_obs, 6);
      matinv (XPX, 6);
      atl (XPy, X, y, n_obs, 6);
      matmul (beta, XPX, XPy, 6,6,1);

      stopflag = 1;
      for (i=0; i<6; i++)  if (fabs (beta[i]) > 0.1)  stopflag = 0;

      Ex0.x0 += beta[0];  Ex0.y0 += beta[1];  Ex0.z0 += beta[2];
      Ex0.omega += beta[3];  Ex0.phi += beta[4];
      Ex0.kappa += beta[5];
    }

  if (stopflag)
    {
      *Ex = Ex0;
      rotation_matrix (*Ex, Ex->dm);
    }
  else puts ("raw orientation impossible");
}
