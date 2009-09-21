#include "ptv.h"

double epi_line (xl, yl, Ex1, I1, G1, Ex2, I2, G2)

double    xl, yl;
Exterior  Ex1, Ex2;
Interior  I1, I2;
Glass     G1, G2;

{
  int i,j;
  double m2;
  double vect1[3], vect2[3], vect3[3], nk[3], n2[3],
    p1l[3], K2[3], k2[3], D2t[3][3];
  void crossprod ();

  /* base O1 -> O2 */
  vect1[0] = Ex2.x0 - Ex1.x0;
  vect1[1] = Ex2.y0 - Ex1.y0;
  vect1[2] = Ex2.z0 - Ex1.z0;

  /* coordinates of arbitrary point P1 in image1 */
  p1l[0] = xl;  p1l[1] = yl;	p1l[2] = - I1.cc;

  /* beam O1 -> P in space */
  matmul (vect2, Ex1.dm, p1l, 3,3,1);

  /* normale to epipolar plane */
  crossprod (vect1,vect2,nk);

  /* normale to image2 */
  vect3[0] = 0;	vect3[1] = 0;	vect3[2] = - I2.cc;

  /* normale to image 2, in space */
  matmul (n2, Ex2.dm, vect3, 3,3,1);

  /* epipolar line in image2, in space */
  crossprod (nk,n2,K2);

  /* epipolar line in image2 */
  for (i=0; i<3; i++)  for (j=0; j<3; j++)  D2t[i][j] = Ex2.dm[j][i];
  matmul (k2, D2t, K2, 3,3,1);
  m2 = k2[1] / k2[0];
  return (m2);
}



int epi_mm (x1, y1, Ex1, I1, G1, Ex2, I2, G2, mmp, xmin, ymin, xmax, ymax)

double     x1, y1;	  	    /* input coord */
Exterior   Ex1, Ex2;        /* orientation data */
Interior   I1, I2;	      	/* orientation data */
Glass      G1, G2;	      	/* glass data */
mm_np	   mmp;		        /* multimed param. (layers) */
double	   *xmin, *ymin, *xmax, *ymax;    /* output search window */

{
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a, b, c, xa,ya,xb,yb;
  double X1,Y1,Z1, X, Y, Z;
  double Zmin, Zmax;
  double Zmin_neu, Zmax_neu,ok;

  //ray_tracing    (x1,y1, Ex1, I1,     mmp, &X1, &Y1, &Z1, &a, &b, &c);
  ray_tracing_v2 (x1,y1, Ex1, I1, G1, mmp, &X1, &Y1, &Z1, &a, &b, &c);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = Zmin_lay[0]
    + (X1-X_lay[0]) * (Zmin_lay[1]-Zmin_lay[0]) / (X_lay[1]-X_lay[0]);
  Zmax = Zmax_lay[0]
    + (X1-X_lay[0]) * (Zmax_lay[1]-Zmax_lay[0]) / (X_lay[1]-X_lay[0]);


  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  //changed by Beat Lüthi 17. March 2008 to keep points inside illuminated area defined in main parameters
  ok=1.;
  if( X<X_lay[0]){
     Zmin_neu = (X_lay[0]-X1) * c/a + Z1;
	 if (Zmin < Zmin_neu && Zmin_neu < Zmax){
        Z = Zmin_neu;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
	 }
	 else{
        ok=0.;
	 }
  }
  else{
     if( X>X_lay[1]){
        Zmin_neu = (X_lay[1]-X1) * c/a + Z1;
	    if (Zmin < Zmin_neu && Zmin_neu < Zmax){
           Z = Zmin_neu;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
	    }
	    else{
           ok=0.;
	    }
     }
  }
  //img_xy_mm_geo_old (X,Y,Z, Ex2, I2,     mmp, &xa, &ya);
  img_xy_mm_geo     (X,Y,Z, Ex2, I2, G2, mmp, &xa, &ya);

  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  //changed by Beat Lüthi 17. March 2008
  if( X<X_lay[0]){
     Zmax_neu = (X_lay[0]-X1) * c/a + Z1;
	 if (Zmax > Zmax_neu && Zmax_neu > Zmin){
        Z = Zmax_neu;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
	 }
	 else{
        ok=0.;
	 }
  }
  else{
     if( X>X_lay[1]){
        Zmax_neu = (X_lay[1]-X1) * c/a + Z1;
	    if (Zmax > Zmax_neu && Zmax_neu > Zmin){
           Z = Zmax_neu;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
	    }
	    else{
           ok=0.;
	    }
     }
  }
  //img_xy_mm_geo_old (X,Y,Z, Ex2, I2,     mmp, &xb, &yb);
  img_xy_mm_geo     (X,Y,Z, Ex2, I2, G2, mmp, &xb, &yb);
   

  /*  ==> window given by xa,ya,xb,yb  */
  if(ok==1){
     *xmin = xa;  *ymin = ya;  *xmax = xb;  *ymax = yb;
  }
  else{
     *xmin = 0.;  *ymin = 0.;  *xmax = 0.;  *ymax = 0.;
  }

  return (0);
}

int epi_mm_2D (x1, y1, Ex1, I1, G1, mmp, xp,yp,zp)

double     x1, y1;	  	/* input coord */
Exterior   Ex1;           	/* orientation data */
Interior   I1;	      	/* orientation data */
Glass      G1;	      	/* glass data */
mm_np	   mmp;		        /* multimed param. (layers) */
double *xp, *yp, *zp;
//double	   *xmin, *ymin, *xmax, *ymax;    /* output search window */

{
  /*  ray tracing gives the point of exit and the direction
      cosines at the waterside of the glass;
      min. and max. depth give window in object space,
      which can be transformed into _2 image
      (use img_xy_mm because of comparison with img_geo)  */

  double a, b, c, xa,ya,xb,yb;
  double X1,Y1,Z1,X,Y,Z;
  
  double Zmin, Zmax;

  ray_tracing_v2 (x1,y1, Ex1, I1, G1, mmp, &X1, &Y1, &Z1, &a, &b, &c);

  /* calculate min and max depth for position (valid only for one setup) */
  Zmin = Zmin_lay[0]
    + (X1-X_lay[0]) * (Zmin_lay[1]-Zmin_lay[0]) / (X_lay[1]-X_lay[0]);
  Zmax = Zmax_lay[0]
    + (X1-X_lay[0]) * (Zmax_lay[1]-Zmax_lay[0]) / (X_lay[1]-X_lay[0]);
 

  Z = 0.5*(Zmin+Zmax);   
  X = X1 + (Z-Z1) * a/c;   
  Y = Y1 + (Z-Z1) * b/c;
  
  *xp=X; *yp=Y; *zp=Z;

  return (0);
}



void find_candidate (crd, pix, num, xa,ya,xb,yb,eps, n, nx, ny, sumg,
					 cand, count, nr)

/*  binarized search in a x-sorted coord-set, exploits shape information  */
/*  gives messages (in examination)  */

coord_2d	crd[];	       		/* metric coordinates */
target		pix[];		       	/* pixel data for correlation */
int	       	num, *count, nr;
double		xa, ya, xb, yb, eps;
int	       	n, nx, ny, sumg;
candidate	cand[];

{
  register int	j;
  int          	j0, dj, p2;
  double       	m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin = -4.40, xmax = 4.40, ymin = -2.94, ymax = 2.94;
  /* max. sensor format (HR 480 / Maxscan) */
  char	       	str[64], buf[32];



  /* define sensor format for search interrupt */
  xmin = (-1) * pix_x * imx/2;	xmax = pix_x * imx/2;
  ymin = (-1) * pix_y * imy/2;	ymax = pix_y * imy/2;
  xmin -= I[nr].xh;	ymin -= I[nr].yh;
  xmax -= I[nr].xh;	ymax -= I[nr].yh;
  correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);


  if (nr != 0)
    {
      strcpy (str, "");	strcpy (buf, "");   puts (str);
    }


  for (j=0; j<8; j++)		       	/* initialize */
    {
      cand[j].pnr = -999;  cand[j].tol = -999;  cand[j].corr = -999;
    }


  m = (yb-ya)/(xb-xa);  b = ya - m*xa; 	/* line equation: y = m*x + b */


  if (xa > xb)			       	/* sort search window */
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }


  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax))  /* sensor area */
    {

      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - eps))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to truncation */


      for (j=j0, *count=0; j<num; j++)         	/* candidate search */
	{
	  if (crd[j].x > xb+eps)  return;     	/* stop search */

	  if ((crd[j].y > ya-eps) && (crd[j].y < yb+eps))
	    {
	      if ((crd[j].x > xa-eps) && (crd[j].x < xb+eps))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
		  if (d < eps)
		    {
		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;

		      /* empirical correlation coefficient
			 from shape and brightness parameters */
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      /* create a tendency to prefer those matches
			 with brighter targets */
		      corr *= ((double) (sumg + pix[p2].sumg));

		      if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg)
			{
			  cand[*count].pnr = p2;
			  cand[*count].tol = d;
			  cand[*count].corr = corr;
			  (*count)++;
			  if (nr > 0)
			    {
			      sprintf (buf, "%3.0f/%3.1f + ", corr, d*1000);
			      strcat (str, buf);  puts (str);
			    }
			}
		    }
		}
	    }
	}
      if (*count == 0  &&  nr > 0)	puts ("- - -");
    }
  else  *count = -1;     		/* out of sensor area */
}






void find_candidate_plus (crd, pix, num, xa,ya,xb,yb, eps, n, nx, ny, sumg,
						  cand, count, nr)

/*  binarized search in a x-sorted coord-set, exploits shape information  */

coord_2d	crd[];
target		pix[];
int    		num, *count;
double		xa, ya, xb, yb, eps;
int    		n, nx, ny, sumg;
candidate	cand[];
int	       	nr;	       	/* image number for ap etc. */

{
  register int	j;
  int dummy;
  int	       	j0, dj, p2;
  double      	m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax;

  /* define sensor format for search interrupt */
  xmin = (-1) * pix_x * imx/2;	xmax = pix_x * imx/2;
  ymin = (-1) * pix_y * imy/2;	ymax = pix_y * imy/2;
  xmin -= I[nr].xh;	ymin -= I[nr].yh;
  xmax -= I[nr].xh;	ymax -= I[nr].yh;
  correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);


  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = -999;  cand[j].corr = -999;
    }


  /* line equation: y = m*x + b */
  if (xa == xb)	xa += 1e-10;
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax))  /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - eps))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;		       	/* due to trunc */

      for (j=j0, *count=0; j<num; j++)			/* candidate search */
	{
	  if (crd[j].x > xb+eps)  return;		/* finish search */

	  if ((crd[j].y > ya-eps) && (crd[j].y < yb+eps))
	    {
	      if ((crd[j].x > xa-eps) && (crd[j].x < xb+eps))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));

		  if ( d < eps )
		    {
		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;

		      /* empirical correlation coefficient
			 from shape and brightness parameters */
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      /* create a tendency to prefer those matches
			 with brighter targets */
		      corr *= ((double) (sumg + pix[p2].sumg));

		      if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg)
			{
				if ( *count < maxcand) {
			  cand[*count].pnr = j;
			  cand[*count].tol = d;
			  cand[*count].corr = corr;
			  (*count)++;
		  } else {
			  dummy=(int)maxcand;
			  printf("in find_candidate_plus: count > maxcand\n");}
			}
		    }
		}
	    }
	}
    }

  else  *count = -1;	 	/* out of sensor area */
}




void find_candidate_plus_msg (crd, pix, num, xa,ya,xb,yb,eps, n, nx, ny, sumg,
							  cand, count, i12)

/*  binarized search in a x-sorted coord-set, exploits shape information  */
/*  gives messages (in examination)  */

coord_2d	crd[];
target		pix[];
int    		num, *count, i12;
double		xa, ya, xb, yb, eps;
int    		n, nx, ny, sumg;
/*
candidate	cand[3];
*/
candidate	cand[];

{
  register int	j;
  int	       	j0, dj, p2;
  double        m, b, d, temp, qn, qnx, qny, qsumg, corr;
  double       	xmin, xmax, ymin, ymax;

  /* define sensor format for search interrupt */
  xmin = (-1) * pix_x * imx/2;	xmax = pix_x * imx/2;
  ymin = (-1) * pix_y * imy/2;	ymax = pix_y * imy/2;
  xmin -= I[i12].xh;	ymin -= I[i12].yh;
  xmax -= I[i12].xh;	ymax -= I[i12].yh;
  correct_brown_affin (xmin,ymin, ap[i12], &xmin,&ymin);
  correct_brown_affin (xmax,ymax, ap[i12], &xmax,&ymax);

  for (j=0; j<4; j++)
    {
      cand[j].pnr = -999;  cand[j].tol = 999;
    }
  m = (yb-ya)/(xb-xa);  b = ya - m*xa;   /* line equation: y = m*x + b */

  if (xa > xb)
    {
      temp = xa;  xa = xb;  xb = temp;
    }
  if (ya > yb)
    {
      temp = ya;  ya = yb;  yb = temp;
    }

  if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) /* sensor area */
    {
      /* binarized search for start point of candidate search */
      for (j0=num/2, dj=num/4; dj>1; dj/=2)
	{
	  if (crd[j0].x < (xa - eps))  j0 += dj;
	  else  j0 -= dj;
	}
      j0 -= 12;  if (j0 < 0)  j0 = 0;  	/* due to trunc */

      for (j=j0, *count=0; j<num; j++) 	/* candidate search */
	{
	  if (crd[j].x > xb+eps)  return;      	/* finish search */

	  if ((crd[j].y > ya-eps) && (crd[j].y < yb+eps))
	    {
	      if ((crd[j].x > xa-eps) && (crd[j].x < xb+eps))
		{
		  d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
		  if (d < eps)
		    {
		      p2 = crd[j].pnr;
		      if (n  < pix[p2].n)      	qn  = (double) n/pix[p2].n;
		      else		       	qn  = (double) pix[p2].n/n;
		      if (nx < pix[p2].nx)	qnx = (double) nx/pix[p2].nx;
		      else		       	qnx = (double) pix[p2].nx/nx;
		      if (ny < pix[p2].ny)	qny = (double) ny/pix[p2].ny;
		      else		       	qny = (double) pix[p2].ny/ny;
		      if (sumg < pix[p2].sumg)
			qsumg = (double) sumg/pix[p2].sumg;
		      else	qsumg = (double) pix[p2].sumg/sumg;


		      /* empirical correlation coefficient
			 from shape and brightness parameters */
		      corr = (4*qsumg + 2*qn + qnx + qny);
		      /* create a tendency to prefer those matches
			 with brighter targets */
		      corr *= ((double) (sumg + pix[p2].sumg));

		      if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg)
			{
			  if (*count>=maxcand)
			    { printf("More candidates than (maxcand): %d\n",*count); return; }
			  cand[*count].pnr = p2;
			  cand[*count].tol = d;
 			  cand[*count].corr = corr;
			  (*count)++;
			  printf ("%d %3.0f/%3.1f \n", p2, corr, d*1000);
			}
		    }
		}
	    }
	}
      if (*count == 0)  puts ("- - -");
    }
  else  *count = -1;		       	       /* out of sensor area */

}





void crossprod (a,b,c)

double  a[3], b[3], c[3];

{
	c[0] = a[1] * b[2]  -  a[2] * b[1];
	c[1] = a[2] * b[0]  -  a[0] * b[2];
	c[2] = a[0] * b[1]  -  a[1] * b[0];
}
