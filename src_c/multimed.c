/****************************************************************************

Routine:	       	multimed.c

Author/Copyright:      	Hans-Gerd Maas

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
	       	       	CH - 8093 Zurich

Creation Date:          20.4.88
	
Description:			
1 point, exterior orientation, multimedia parameters   	==> X,Y=shifts
		  (special case 3 media, medium 2 = plane parallel)
							
Routines contained:		-

****************************************************************************/

#include "ptv.h"

double get_mmf_from_mmLUT ();


void  multimed (ex,mm,X,Y,Z,Xq,Yq)
Exterior  ex;
mm_3p	  mm;
double    X, Y, Z, *Xq,*Yq;
{
  int    i, it=0;
  double beta1, beta2, beta3, r, rbeta, rdiff, rq, mmf;
  double ocf=1.0;      	/* over compensation factor for faster convergence */
   
  if (mm.n1==1 && mm.n2==1 && mm.n3==1)		/* 1-medium case */
    {
      *Xq = X;	*Yq = Y;	return;
    }
    
  /* interpolation in mmLUT, if selected (requires some global variables) */
  if (mm.lut)
    {
      /* check, which is the correct image */
      for (i=0; i<n_img; i++)
	if (Ex[i].x0 == ex.x0  &&  Ex[i].y0 == ex.y0  &&  Ex[i].z0 == ex.z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)
	{
	  *Xq = ex.x0 + (X-ex.x0) * mmf;
	  *Yq = ex.y0 + (Y-ex.y0) * mmf;
	  return;
	}
    }
  
  /* iterative procedure */
  
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      beta2 = asin (sin(beta1) * mm.n1/mm.n2);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d) * tan(beta1) + mm.d * tan(beta2) - Z * tan(beta3);
      rdiff = r - rbeta;
      rdiff *= ocf;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      puts ("Multimed stopped after 40. Iteration");
      if (fabs(rdiff) < 0.1)	
	{
	  *Xq = ex.x0 + (X-ex.x0) * rq/r;
	  *Yq = ex.y0 + (Y-ex.y0) * rq/r;
	}
      else
	{
	  *Xq = X; *Yq = Y;
	}
      return;
    }
   
  if (r != 0)
    {
      *Xq = ex.x0 + (X-ex.x0) * rq/r;
      *Yq = ex.y0 + (Y-ex.y0) * rq/r;
    }
  else
    {
      *Xq = X;
      *Yq = Y;
    }
}



double multimed_r (ex,mm,X,Y,Z)
/* calculates and returns the radial shift */
Exterior  ex;
mm_3p	  mm;
double	  X, Y, Z;
{
  int  	 i, it=0;
  double beta1, beta2, beta3, r, rbeta, rdiff, rq, mmf;
  double ocf=1.0;      	/* over compensation factor for faster convergence */
  
  if (mm.n1==1 && mm.n2==1 && mm.n3==1)	return (1.0);	/* 1-medium case */
    
  /* interpolation in mmLUT, if selected (requires some global variables) */
  if (mm.lut)
    {
      /* check, which is the correct image */
      for (i=0; i<n_img; i++)
	if (Ex[i].x0 == ex.x0  &&  Ex[i].y0 == ex.y0  &&  Ex[i].z0 == ex.z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)	return (mmf);
    }
  
  /* iterative procedure */
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      beta2 = asin (sin(beta1) * mm.n1/mm.n2);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d) * tan(beta1) + mm.d * tan(beta2) - Z * tan(beta3);
      rdiff = r - rbeta;
      rdiff *= ocf;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      puts ("Multimed_r stopped after 40. Iteration");	return (1.0);
    }
  
  if (r != 0)	return (rq/r);	else return (1.0);
}


void  multimed_nlay_true (ex,mm,X,Y,Z,Xq,Yq)
Exterior  ex;
mm_np	  mm;
double    X, Y, Z, *Xq,*Yq;
{
  int		i, it=0;
  double	beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  double	ocf=1.00;		/* over compensation factor for faster convergence */
  
  /* interpolation in mmLUT, if selected (requires some global variables) */
  if (mm.lut)
    {
      /* check, which is the correct image */
      for (i=0; i<n_img; i++)
	if (Ex[i].x0 == ex.x0  &&  Ex[i].y0 == ex.y0  &&  Ex[i].z0 == ex.z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)
	{
	  *Xq = ex.x0 + (X-ex.x0) * mmf;
	  *Yq = ex.y0 + (Y-ex.y0) * mmf;
	  return;
	}
    }
   
  /* iterative procedure */
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rdiff *= ocf;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      *Xq = X; *Yq = Y;
      puts ("Multimed_nlay_true stopped after 40. Iteration");	return;
    }
    
  if (r != 0)
    {
      *Xq = ex.x0 + (X-ex.x0) * rq/r;
      *Yq = ex.y0 + (Y-ex.y0) * rq/r;
    }
  else
    {
      *Xq = X;
      *Yq = Y;
    }
}



void  multimed_nlay (ex,mm,X,Y,Z,Xq,Yq)
Exterior	ex;
mm_np		mm;
double  	X, Y, Z, *Xq,*Yq;
{
  int		i, it=0;
  double	beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  
  /* interpolation in mmLUT, if selected (requires some global variables) */
  if (mm.lut)
    {    
      /* check, which is the correct image */
      for (i=0; i<n_img; i++)
	if (Ex[i].x0 == ex.x0  &&  Ex[i].y0 == ex.y0  &&  Ex[i].z0 == ex.z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)
	{
	  *Xq = ex.x0 + (X-ex.x0) * mmf;
	  *Yq = ex.y0 + (Y-ex.y0) * mmf;
	  return;
	}
    }
  
  /* iterative procedure (if mmLUT does not exist or has no entry) */
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      *Xq = X; *Yq = Y;
      puts ("Multimed_nlay stopped after 40. Iteration");	return;
    }
    
  if (r != 0)
    {
      *Xq = ex.x0 + (X-ex.x0) * rq/r;
      *Yq = ex.y0 + (Y-ex.y0) * rq/r;
    }
  else
    {
      *Xq = X;
      *Yq = Y;
    }
}



double multimed_r_nlay (ex,mm,X,Y,Z)
/* calculates and returns the radial shift */
Exterior	ex;
mm_np		mm;
double		X, Y, Z;
{
  int  	i, it=0;
  double beta1, beta2[32], beta3, r, rbeta, rdiff, rq, mmf;
  double ocf=1.0; /* over compensation factor for faster convergence */
  
  
  /* 1-medium case */
  if (mm.n1==1 && mm.nlay == 1 && mm.n2[0]==1 && mm.n3==1) return (1.0);
    
  /* interpolation in mmLUT, if selected (requires some global variables) */
  if (mm.lut)
    {
      /* check, which is the correct image */
      for (i=0; i<n_img; i++)
	if (Ex[i].x0 == ex.x0  &&  Ex[i].y0 == ex.y0  &&  Ex[i].z0 == ex.z0)
	  break;
      
      mmf = get_mmf_from_mmLUT (i, X,Y,Z);
      
      if (mmf > 0)	return (mmf);
    }
 
  /* iterative procedure */
  r = sqrt ((X-ex.x0)*(X-ex.x0)+(Y-ex.y0)*(Y-ex.y0));
  rq = r;
  
  do
    {
      beta1 = atan (rq/(ex.z0-Z));
      for (i=0; i<mm.nlay; i++)	beta2[i] = asin (sin(beta1) * mm.n1/mm.n2[i]);
      beta3 = asin (sin(beta1) * mm.n1/mm.n3);
      
      rbeta = (ex.z0-mm.d[0]) * tan(beta1) - Z * tan(beta3);
      for (i=0; i<mm.nlay; i++)	rbeta += (mm.d[i] * tan(beta2[i]));
      rdiff = r - rbeta;
      rdiff *= ocf;
      rq += rdiff;
      it++;
    }
  while (((rdiff > 0.001) || (rdiff < -0.001))  &&  it < 40);
  
  if (it >= 40)
    {
      puts ("Multimed_r_nlay stopped after 40. Iteration");	return (1.0);
    }
  
  if (r != 0)	return (rq/r);	else return (1.0);
}



void init_mmLUT (i_cam)
int    	i_cam;
{
  register int	i,j, nr, nz;
  double       	X,Y,Z, R, X1,Y1,Z1, Zmin, Rmax=0,Zmax, a,b,c;
  double       	x,y, *Ri,*Zi;
  double       	rw = 2;
    
  /* find extrema of imaged object volume */
  /* ==================================== */
  
  /* find extrema in depth */
  
  fpp = fopen ("parameters/criteria.par", "r");
  fscanf (fpp, "%lf\n", &X);
  fscanf (fpp, "%lf\n", &Zmin);
  fscanf (fpp, "%lf\n", &Zmax);
  fscanf (fpp, "%lf\n", &X);
  fscanf (fpp, "%lf\n", &Z);	if (Z < Zmin)	Zmin = Z;
  fscanf (fpp, "%lf\n", &Z);	if (Z > Zmax)	Zmax = Z;
  fclose (fpp);
  
  Zmin -= fmod (Zmin, rw);
  Zmax += (rw - fmod (Zmax, rw));
  
  /* intersect with image vertices rays */
  
  pixel_to_metric (0., 0., imx,imy, pix_x,pix_y, &x,&y, chfield);
  x = x - I[i_cam].xh;
  y = y - I[i_cam].yh;
  correct_brown_affin (x, y, ap[i_cam], &x,&y);
  ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	
  
  if (R > Rmax)	Rmax = R;
  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));
  
  if (R > Rmax)	Rmax = R;
  pixel_to_metric (0., (double) imy, imx,imy, pix_x,pix_y, &x,&y, chfield);
  x = x - I[i_cam].xh;
  y = y - I[i_cam].yh;
  correct_brown_affin (x, y, ap[i_cam], &x,&y);
  ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));
  
  if (R > Rmax)	Rmax = R;
  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));
  
  if (R > Rmax)	Rmax = R;
  
  pixel_to_metric ((double) imx, 0., imx,imy, pix_x,pix_y, &x,&y, chfield);
  x = x - I[i_cam].xh;
  y = y - I[i_cam].yh;
  correct_brown_affin (x, y, ap[i_cam], &x,&y);
  ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));
  
  if (R > Rmax)	Rmax = R;
  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	

  if (R > Rmax)	Rmax = R;
  pixel_to_metric ((double) imx, (double) imy,
		   imx,imy, pix_x,pix_y, &x,&y, chfield);
  x = x - I[i_cam].xh;
  y = y - I[i_cam].yh;
  correct_brown_affin (x, y, ap[i_cam], &x,&y);
  ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
  Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	
  
  if (R > Rmax)	Rmax = R;
  Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
  R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
	      + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	
  
  if (R > Rmax)	Rmax = R;
  
  /* round values (-> enlarge) */
  Rmax += (rw - fmod (Rmax, rw));
    
  /* get # of rasterlines in r,z */
  nr = Rmax/rw + 1;
  nz = (Zmax-Zmin)/rw + 1;
 
  /* create twodimensional mmLUT structure */
  mmLUT[i_cam].origin.x = Ex[i_cam].x0;
  mmLUT[i_cam].origin.y = Ex[i_cam].y0;
  mmLUT[i_cam].origin.z = Zmin;
  mmLUT[i_cam].nr = nr;
  mmLUT[i_cam].nz = nz;
  mmLUT[i_cam].rw = rw;
  mmLUT[i_cam].data = (double *) malloc (nr*nz * sizeof (double));
   
  /* fill mmLUT structure */
  /* ==================== */
  
  Ri = (double *) malloc (nr * sizeof (double));
  for (i=0; i<nr; i++)	Ri[i] = i*rw;
  Zi = (double *) malloc (nz * sizeof (double));
  for (i=0; i<nz; i++)	Zi[i] = Zmin + i*rw;
  
  for (i=0; i<nr; i++)	for (j=0; j<nz; j++)
    {
      mmLUT[i_cam].data[i*nz + j]
	= multimed_r_nlay (Ex[i_cam], mmp,
			   Ri[i]+Ex[i_cam].x0, Ex[i_cam].y0, Zi[j]);
    }
}



double get_mmf_from_mmLUT (i_cam, X,Y,Z)
int		i_cam;
double	X,Y,Z;
{
  int		i, ir,iz, nr,nz, v4[4];
  double	R, sr,sz, rw, mmf=1;
  
  rw =  mmLUT[i_cam].rw;
  
  Z -= mmLUT[i_cam].origin.z; sz = Z/rw; iz = (int) sz;	sz -= iz;
  X -= mmLUT[i_cam].origin.x;
  Y -= mmLUT[i_cam].origin.y;
  R = sqrt (X*X + Y*Y);	sr = R/rw; ir = (int) sr; sr -= ir;
    
  nz =  mmLUT[i_cam].nz;
  nr =  mmLUT[i_cam].nr;
    
  /* check whether point is inside camera's object volume */
  if (ir > nr)				return (0);
  if (iz < 0  ||  iz > nz)	return (0);
  
  /* bilinear interpolation in r/z box */
  /* ================================= */
  
  /* get vertices of box */
  v4[0] = ir*nz + iz;
  v4[1] = ir*nz + (iz+1);
  v4[2] = (ir+1)*nz + iz;
  v4[3] = (ir+1)*nz + (iz+1);
  
  /* 2. check wther point is inside camera's object volume */
  /* important for epipolar line computation */
  for (i=0; i<4; i++)
    if (v4[i] < 0  ||  v4[i] > nr*nz)	return (0);
  
  /* interpolate */
  mmf = mmLUT[i_cam].data[v4[0]] * (1-sr)*(1-sz)
    + mmLUT[i_cam].data[v4[1]] * (1-sr)*sz
    + mmLUT[i_cam].data[v4[2]] * sr*(1-sz)
    + mmLUT[i_cam].data[v4[3]] * sr*sz;
  
  return (mmf);
}



void volumedimension (xmax, xmin, ymax, ymin, zmax, zmin)
double *xmax, *xmin, *ymax, *ymin, *zmax, *zmin;
{
  int	i_cam;
  double X,Y,Z, R, X1,Y1,Z1, Zmin, Rmax=0,Zmax, a,b,c;
  double x,y;
    
  /* find extrema of imaged object volume */
  /* ==================================== */
  
  fpp = fopen ("parameters/criteria.par", "r");
  fscanf (fpp, "%lf\n", &X);
  fscanf (fpp, "%lf\n", &Zmin);
  fscanf (fpp, "%lf\n", &Zmax);
  fscanf (fpp, "%lf\n", &X);
  fscanf (fpp, "%lf\n", &Z);	if (Z < Zmin)	Zmin = Z;
  fscanf (fpp, "%lf\n", &Z);	if (Z > Zmax)	Zmax = Z;
  fclose (fpp);

  *zmin=Zmin;
  *zmax=Zmax;

  for (i_cam=0;i_cam<n_img;i_cam++)
    {  
      /* intersect with image vertices rays */
      pixel_to_metric (0.0, 0.0, imx,imy, pix_x,pix_y, &x,&y, chfield);
      x = x - I[i_cam].xh;
      y = y - I[i_cam].yh;
      correct_brown_affin (x, y, ap[i_cam], &x,&y);
      ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
      Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;

      if (R > Rmax)	Rmax = R;
      Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      pixel_to_metric (0.0, (double) imy, imx,imy, pix_x,pix_y, &x,&y, chfield);
      x = x - I[i_cam].xh;
      y = y - I[i_cam].yh;
      correct_brown_affin (x, y, ap[i_cam], &x,&y);
      ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
      Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      
      pixel_to_metric ((double) imx, 0., imx,imy, pix_x,pix_y, &x,&y, chfield);
      x = x - I[i_cam].xh;
      y = y - I[i_cam].yh;
      correct_brown_affin (x, y, ap[i_cam], &x,&y);
      ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
      Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      pixel_to_metric ((double) imx, (double) imy,
		       imx,imy, pix_x,pix_y, &x,&y, chfield);
      x = x - I[i_cam].xh;
      y = y - I[i_cam].yh;
      correct_brown_affin (x, y, ap[i_cam], &x,&y);
      ray_tracing (x,y, Ex[i_cam], I[i_cam], mmp, &X1, &Y1, &Z1, &a, &b, &c);
      Z = Zmin;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));	

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y;
      
      if (R > Rmax)	Rmax = R;
      Z = Zmax;   X = X1 + (Z-Z1) * a/c;   Y = Y1 + (Z-Z1) * b/c;
      R = sqrt (  (X-Ex[i_cam].x0)*(X-Ex[i_cam].x0)
		  + (Y-Ex[i_cam].y0)*(Y-Ex[i_cam].y0));           

      if ( X > *xmax) *xmax=X;
      if ( X < *xmin) *xmin=X;
      if ( Y > *ymax) *ymax=Y;
      if ( Y < *ymin) *ymin=Y; 
    }
}
