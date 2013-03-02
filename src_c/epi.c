/*
Copyright (c) 1990-2011 ETH Zurich

See the file license.txt for copying permission.
*/

#include "ptv.h"

#ifdef EVER_CALLED		// Unused function, ad holten 12-2012
double epi_line (double xl, double yl, Exterior Ex1, Interior I1, Glass G1,
				 Exterior Ex2, Interior I2, Glass G2)
{
	int i,j;
	double m2, vect1[3], vect2[3], vect3[3], nk[3], n2[3];
	double p1l[3], K2[3], k2[3], D2t[3][3];
	void crossprod (double a[3], double b[3], double c[3]);

	/* base O1 -> O2 */
	vect1[0] = Ex2.x0 - Ex1.x0;
	vect1[1] = Ex2.y0 - Ex1.y0;
	vect1[2] = Ex2.z0 - Ex1.z0;

	/* coordinates of arbitrary point P1 in image1 */
	p1l[0] = xl;  p1l[1] = yl;	  p1l[2] = -I1.cc;

	/* beam O1 -> P in space */
	matmul (vect2, Ex1.dm, p1l, 3,3,1);

	/* normale to epipolar plane */
	crossprod (vect1,vect2,nk);

	/* normale to image2 */
	vect3[0] = 0;  vect3[1] = 0;  vect3[2] = -I2.cc;

	/* normale to image 2, in space */
	matmul (n2, Ex2.dm, vect3, 3,3,1);

	/* epipolar line in image2, in space */
	crossprod (nk,n2,K2);

	/* epipolar line in image2 */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			D2t[i][j] = Ex2.dm[j][i];
	matmul (k2, D2t, K2, 3,3,1);
	m2 = k2[1] / k2[0];
	return (m2);
}
#endif


int epi_mm (double x1, double y1, 
			Exterior Ex1, Interior I1, Glass G1, 
			Exterior Ex2, Interior I2, Glass G2, mm_np mmp, 
			double* xmin, double* ymin, double* xmax, double* ymax)
	// double	x1, y1		input coord
	// Exterior Ex1, Ex2	orientation data
	// Interior I1, I2		orientation data
	// Glass	G1, G2; 	glass data
	// mm_np	mmp;		multimed param. (layers)
	// double	*xmin, *ymin, *xmax, *ymax		output search window
{
	/* ray tracing gives the point of exit and the direction cosines
	   at the waterside of the glass; 
	   min. and max. depth give window in object space, which can be 
	   transformed into _2 image
	  (use img_xy_mm because of comparison with img_geo)  */

	double a, b, c, xa,ya,xb,yb;
	double X1,Y1,Z1, X, Y, Z;
	double Zmin, Zmax;

	//ray_tracing	 (x1,y1, Ex1, I1,	  mmp, &X1, &Y1, &Z1, &a, &b, &c);
	ray_tracing_v2 (x1,y1, Ex1, I1, G1, mmp, &X1, &Y1, &Z1, &a, &b, &c);

	/* calculate min and max depth for position (valid only for one setup) */
	Zmin = Zmin_lay[0]
		   + (X1-X_lay[0]) * (Zmin_lay[1]-Zmin_lay[0]) / (X_lay[1]-X_lay[0]);
	Zmax = Zmax_lay[0]
		   + (X1-X_lay[0]) * (Zmax_lay[1]-Zmax_lay[0]) / (X_lay[1]-X_lay[0]);


	Z = Zmin;	X = X1 + (Z-Z1) * a/c;	 Y = Y1 + (Z-Z1) * b/c;
	//img_xy_mm_geo_old (X,Y,Z, Ex2, I2,	 mmp, &xa, &ya);
	img_xy_mm_geo (X,Y,Z, Ex2, I2, G2, mmp, &xa, &ya);

	Z = Zmax;	X = X1 + (Z-Z1) * a/c;	 Y = Y1 + (Z-Z1) * b/c;
	//img_xy_mm_geo_old (X,Y,Z, Ex2, I2,	 mmp, &xb, &yb);
	img_xy_mm_geo (X,Y,Z, Ex2, I2, G2, mmp, &xb, &yb);

	/*	==> window given by xa,ya,xb,yb  */

	*xmin = xa;  *ymin = ya;  *xmax = xb;  *ymax = yb;

	return (0);
}

int epi_mm_2D (double x1, double y1, Exterior Ex1, Interior I1, Glass G1, 
			   mm_np mmp, double *xp, double *yp, double *zp)
	// double	  x1, y1	input coord
	// Exterior   Ex1		orientation data
	// Interior   I1		orientation data
	// Glass	  G1		glass data
	// mm_np	  mmp		multimed param. (layers)
	// double *xp, *yp, *zp;
{
	/* ray tracing gives the point of exit and the direction cosines
	   at the waterside of the glass; 
	   min. and max. depth give window in object space, which can be 
	   transformed into _2 image */
	double a, b, c;
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



#ifdef EVER_CALLED		// Unused function, ad holten 12-2012
void find_candidate (coord_2d crd[], target pix[], int num, 
					 double xa, double ya, double xb, double yb, 
					 double eps, int n, int nx, int ny, int sumg,
					 candidate cand[], int *count, int nr)
	// binarized search in a x-sorted coord-set, exploits shape information
	// gives messages (in examination)
	// coord_2d crd[]	metric coordinates
	// target	pix[]	pixel data for correlation
{
	register int j;
	int 		 j0, dj, p2;
	double		 m, b, d, temp, qn, qnx, qny, qsumg, corr;
	double		 xmin = -4.40, xmax = 4.40, ymin = -2.94, ymax = 2.94;
	/* max. sensor format (HR 480 / Maxscan) */
	char		 str[64], buf[32];

	/* define sensor format for search interrupt */
	xmin  = -pix_x * imx/2;  xmax = pix_x * imx/2;
	ymin  = -pix_y * imy/2;  ymax = pix_y * imy/2;
	xmin -= I[nr].xh;	 ymin -= I[nr].yh;
	xmax -= I[nr].xh;	 ymax -= I[nr].yh;
	correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
	correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);

	if (nr != 0) {
		strcpy (str, "");
		strcpy (buf, "");
		puts (str);
	}
	for (j=0; j<8; j++) {					/* initialize */
		cand[j].pnr  = -999;  
		cand[j].tol  = -999;  
		cand[j].corr = -999;
	}

	m = (yb-ya)/(xb-xa);  b = ya - m*xa;	/* line equation: y = m*x + b */


	if (xa > xb) { temp = xa; xa = xb; xb = temp; }	/* sort search window */
	if (ya > yb) { temp = ya; ya = yb; yb = temp; }

	if ((xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) {  /* sensor area */
		/* binarized search for start point of candidate search */
		for (j0=num/2, dj=num/4; dj>1; dj/=2) {
			if (crd[j0].x < (xa - eps))  j0 += dj;
			else  j0 -= dj;
		}
		j0 -= 12;  if (j0 < 0)	j0 = 0; 		/* due to truncation */

		for (j=j0, *count=0; j<num; j++) {		/* candidate search */
			if (crd[j].x > xb+eps)	return; 	/* stop search */

			if ((crd[j].y > ya-eps) && (crd[j].y < yb+eps)) {
				if ((crd[j].x > xa-eps) && (crd[j].x < xb+eps)) {
					d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
					if (d < eps) {
						p2 = crd[j].pnr;
						if (n  < pix[p2].n) 	 qn  = (double) n/pix[p2].n;
						else					 qn  = (double) pix[p2].n/n;
						if (nx < pix[p2].nx)	 qnx = (double) nx/pix[p2].nx;
						else					 qnx = (double) pix[p2].nx/nx;
						if (ny < pix[p2].ny)	 qny = (double) ny/pix[p2].ny;
						else					 qny = (double) pix[p2].ny/ny;
						if (sumg < pix[p2].sumg) qsumg = (double) sumg/pix[p2].sumg;
						else					 qsumg = (double) pix[p2].sumg/sumg;

						/* empirical correlation coefficient from shape and 
						   brightness parameters */
						corr = (4*qsumg + 2*qn + qnx + qny);
						/* create a tendency to prefer those matches
						   with brighter targets */
						corr *= ((double) (sumg + pix[p2].sumg));

						if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg) {
							cand[*count].pnr = p2;
							cand[*count].tol = d;
							cand[*count].corr = corr;
							(*count)++;
							if (nr > 0) {
								sprintf (buf, "%3.0f/%3.1f + ", corr, d*1000);
								strcat (str, buf);	puts (str);
							}
						}
					}
				}
			}
		}
		if (*count == 0 && nr > 0)	puts ("- - -");
	}
	else  *count = -1;	/* out of sensor area */
}
#endif

void find_candidate_plus (coord_2d crd[], target pix[], int num, 
					 double xa, double ya, double xb, double yb, 
					 double eps, int n, int nx, int ny, int sumg,
					 candidate cand[], int *count, int nr, const char** argv)

	// binarized search in a x-sorted coord-set, exploits shape information
	// int nr	image number for ap etc.
{
	register int j;
	int dummy, j0, dj, p2;
	double m, b, d, temp, qn, qnx, qny, qsumg, corr;
	double xmin, xmax, ymin, ymax,particle_size;
	int dumbbell=0;
	double tol_band_width;

	//Beat Mai 2010 for dumbbell
	if (atoi(argv[1])==3)
		dumbbell=1;

	if (dumbbell==0) {
	  ///// here is new Beat version of April 2010
		if (nx>ny) particle_size = nx;
		else	   particle_size = ny;
		tol_band_width = eps*0.5*(pix_x+pix_y)*particle_size;
	}
	else
		tol_band_width=eps;
	
	if (tol_band_width < 0.06)
		tol_band_width = 0.06;

	/* define sensor format for search interrupt */
	xmin  = -pix_x * imx/2;    xmax = pix_x * imx/2;
	ymin  = -pix_y * imy/2;    ymax = pix_y * imy/2;
	xmin -= I[nr].xh;	 ymin -= I[nr].yh;
	xmax -= I[nr].xh;	 ymax -= I[nr].yh;
	correct_brown_affin (xmin,ymin, ap[nr], &xmin,&ymin);
	correct_brown_affin (xmax,ymax, ap[nr], &xmax,&ymax);

	for (j=0; j<4; j++) { /* initialize, why 8 in the prev version, adh ????*/
		cand[j].pnr  = -999;
		cand[j].tol  = -999;
		cand[j].corr = -999;
	}

	/* line equation: y = m*x + b */
	if (xa == xb)  xa += 1e-10;
	m = (yb-ya)/(xb-xa);  b = ya - m*xa;

	if (xa > xb) { temp = xa; xa = xb; xb = temp; }
	if (ya > yb) { temp = ya; ya = yb; yb = temp; }

	if ((xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) { /* sensor area */
		/* binarized search for start point of candidate search */
		for (j0=num/2, dj=num/4; dj>1; dj/=2) {
			if (crd[j0].x < (xa - tol_band_width))	j0 += dj;
			else  j0 -= dj;
		}
		j0 -= 12;  if (j0 < 0)	j0 = 0; 				/* due to trunc */

		for (j=j0, *count=0; j<num; j++) {				/* candidate search */
			if (crd[j].x > xb+tol_band_width) return;	/* finish search */

			//  ad holten, 12-2012 : merged the two if's
			//  if ((crd[j].y > ya-tol_band_width) && (crd[j].y < yb+tol_band_width)) {
			//	    if ((crd[j].x > xa-tol_band_width) && (crd[j].x < xb+tol_band_width)) {
			if (crd[j].y > ya-tol_band_width && crd[j].y < yb+tol_band_width &&
				crd[j].x > xa-tol_band_width && crd[j].x < xb+tol_band_width)
			{
				d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
	  
				// Beat: modified in April 2010 to allow for better treatment of 
				// different sized traced particles, in particular colloids and tracers
				// old : if (d < eps) {
				// new : // if (nx>ny) particle_size=nx;
				//		 // else	   particle_size=ny;
				//		 if (d < tol_band_width) {

				if (d < tol_band_width) {
					p2 = crd[j].pnr;
					if (n  < pix[p2].n) 	 qn  = (double) n/pix[p2].n;
					else					 qn  = (double) pix[p2].n/n;
					if (nx < pix[p2].nx)	 qnx = (double) nx/pix[p2].nx;
					else					 qnx = (double) pix[p2].nx/nx;
					if (ny < pix[p2].ny)	 qny = (double) ny/pix[p2].ny;
					else					 qny = (double) pix[p2].ny/ny;
					if (sumg < pix[p2].sumg) qsumg = (double) sumg/pix[p2].sumg;
					else					 qsumg = (double) pix[p2].sumg/sumg;

					// empirical correlation coefficient from shape and 
					// brightness parameters 
					corr = (4*qsumg + 2*qn + qnx + qny);
					// create a tendency to prefer those matches
					// with brighter targets 
					corr *= ((double) (sumg + pix[p2].sumg));

					if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg) {
						if (*count < maxcand) {
							cand[*count].pnr = j;
							cand[*count].tol = d;
							cand[*count].corr = corr;
							(*count)++;
						} else {
							dummy = (int)maxcand;
							printf("in find_candidate_plus: count > maxcand\n");
						}
					}
				}
			}
		}
	}
	
	else  *count = -1;		   /* out of sensor area */
}



void find_candidate_plus_msg (coord_2d crd[], target pix[], int num, 
					 double xa, double ya, double xb, double yb, 
					 double eps, int n, int nx, int ny, int sumg,
					 candidate cand[], int *count, int i12)
	//	binarized search in a x-sorted coord-set, exploits shape information  */
	//	gives messages (in examination)  */
{
	register int j;
	int    j0, dj, p2;
	double m, b, d, temp, qn, qnx, qny, qsumg, corr;
	double xmin, xmax, ymin, ymax, tol_band_width, particle_size;

	/* define sensor format for search interrupt */
	xmin  = -pix_x * imx/2;    xmax = pix_x * imx/2;
	ymin  = -pix_y * imy/2;    ymax = pix_y * imy/2;
	xmin -= I[i12].xh;	  ymin -= I[i12].yh;
	xmax -= I[i12].xh;	  ymax -= I[i12].yh;
	correct_brown_affin (xmin,ymin, ap[i12], &xmin,&ymin);
	correct_brown_affin (xmax,ymax, ap[i12], &xmax,&ymax);

	if (nx>ny) particle_size=nx;
	else	   particle_size=ny;
	tol_band_width = eps*0.5*(pix_x+pix_y)*particle_size;

	for (j=0; j<4; j++) {
		cand[j].pnr = -999;
		cand[j].tol = 999;
	}
	m = (yb-ya)/(xb-xa);  b = ya - m*xa;   /* line equation: y = m*x + b */

	if (xa > xb) { temp = xa; xa = xb; xb = temp; }
	if (ya > yb) { temp = ya; ya = yb; yb = temp; }

	if ( (xb>xmin) && (xa<xmax) && (yb>ymin) && (ya<ymax)) {	/* sensor area */
		/* binarized search for start point of candidate search */
		for (j0=num/2, dj=num/4; dj>1; dj/=2) {
			if (crd[j0].x < (xa - tol_band_width))	j0 += dj;
			else  j0 -= dj;
		}
		j0 -= 12;  if (j0 < 0)	j0 = 0; 				/* due to trunc */

		for (j=j0, *count=0; j<num; j++) {				/* candidate search */
			if (crd[j].x > xb+tol_band_width)  return;	/* finish search */

			// ad holten, 12-2012 : merged two if's to a single one 
			if (crd[j].y > ya-tol_band_width && crd[j].y < yb+tol_band_width &&
				crd[j].x > xa-tol_band_width && crd[j].x < xb+tol_band_width)
			{
				d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));
				if ( d < tol_band_width ){
					p2 = crd[j].pnr;
					if (n  < pix[p2].n) 	 qn  = (double) n/pix[p2].n;
					else					 qn  = (double) pix[p2].n/n;
					if (nx < pix[p2].nx)	 qnx = (double) nx/pix[p2].nx;
					else					 qnx = (double) pix[p2].nx/nx;
					if (ny < pix[p2].ny)	 qny = (double) ny/pix[p2].ny;
					else					 qny = (double) pix[p2].ny/ny;
					if (sumg < pix[p2].sumg) qsumg = (double) sumg/pix[p2].sumg;
					else					 qsumg = (double) pix[p2].sumg/sumg;

					// empirical correlation coefficient from shape and 
					// brightness parameters 
					corr = (4*qsumg + 2*qn + qnx + qny);
					// create a tendency to prefer those matches
					// with brighter targets 
					corr *= ((double) (sumg + pix[p2].sumg));

					if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg) {
						if (*count >= maxcand) {
							printf("More candidates than (maxcand): %d\n",*count);
							return;
						}
						cand[*count].pnr = p2;
						cand[*count].tol = d;
						cand[*count].corr = corr;
						(*count)++;
						printf ("%d %3.0f/%3.1f \n", p2, corr, d*1000);
					}
				}
			}
		}
		if (*count == 0)  puts ("- - -");
	}
	else  *count = -1;							/* out of sensor area */
}

#ifdef EVER_CALLED		// Unused function, ad holten 12-2012
void crossprod (double a[3], double b[3], double c[3])
{
	c[0] = a[1] * b[2]	-  a[2] * b[1];
	c[1] = a[2] * b[0]	-  a[0] * b[2];
	c[2] = a[0] * b[1]	-  a[1] * b[0];
}
#endif