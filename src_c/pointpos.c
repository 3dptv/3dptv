/****************************************************************************

Routine:	       	pointpos.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	July 1988, June 1989
	
Description:	       	point positioning with least squares adjustment
		       	(multimedia, 2 or 3 channels)
	
Routines contained:		

****************************************************************************/
#include "ptv.h"

void det_lsq_old (Ex, I, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp)
Exterior	Ex[4];
Interior	I[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp;
{
  int		i, j, itnum, stopflag, n_obs;
  double	X[8][3], XtX[3][3], y[8], Xty[3], beta[3],
    Xbeta[8], edach[8], omega=0., sigma0, sigmabeta[3];
  double	Zx, Zy, N, mmf;
  double	xp1, yp1;
  double	xp[4], yp[4];
  
  double	multimed_r_nlay ();
  
  /* static parameters for examination of Qxx, Qvv */
  static int   	Sn;
  static double	Ssigma0, SQxx[3][3], SQvv[8][8];
  double	Xt[3][8], buf[3][8], Qvv[8][8];
  
  xp[0] = x1;	xp[1] = x2;	xp[2] = x3;	xp[3] = x4;
  yp[0] = y1;	yp[1] = y2;	yp[2] = y3;	yp[3] = y4;

  /* approximate values: rt_is with the uncorrected image coordinates */
  /* of the first 2 rays		*/
  
  /*
    for (i=0; i<n_img; i++)	if (xp[i] > -1e10  &&  yp[i] > -1e10)
    {
    ray_tracing (xp[i], yp[i], Ex[i], I[i],
    mm, &X1, &Y1, &Z1, &a1, &b1, &c1);	break;
    }
    for (i=i+1; i<n_img; i++)	if (xp[i] > -1e10  &&  yp[i] > -1e10)
    {
    ray_tracing (xp[i], yp[i], Ex[i], I[i],
    mm, &X2, &Y2, &Z2, &a2, &b2, &c2);	break;
    }
    
    intersect_rt (X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, Xp, Yp, Zp);
    */  
  
  /* hack due to problems with approx in det_lsq: */
  
  /************************************************************************/  
  /* least squares adjustment */
  itnum = 0;	stopflag = 0;
  while ((stopflag == 0) && (itnum < 10))
    {
      /********************************************************************/
      /* X-matrix  +  observation vector from image coordinates */
      
      for (i=0; i<4; i++)
	{
	  if (xp[i] > -1e10  &&  yp[i] > -1e10)
	    {
	      /* derivatives d(x(i),y(i)) / d(X,Y,Z)   	(Kraus)	*/
	      
	      Zx =  Ex[i].dm[0][0] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][0] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][0] * (*Zp - Ex[i].z0);
	      Zy =  Ex[i].dm[0][1] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][1] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][1] * (*Zp - Ex[i].z0);
	      N  =  Ex[i].dm[0][2] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][2] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][2] * (*Zp - Ex[i].z0);
	      
	      for (j=0; j<3; j++)
		X[2*i][j] = ((-I[i].cc)/(N*N))
		  * (N*Ex[i].dm[j][0] - Zx*Ex[i].dm[j][2]);
	      for (j=0; j<3; j++)
		X[2*i+1][j] = ((-I[i].cc)/(N*N))
		  * (N*Ex[i].dm[j][1] - Zy*Ex[i].dm[j][2]);
	            
	      /* multimedia factor (radial shift) */

	      mmf = multimed_r_nlay (Ex[i], mm, *Xp, *Yp, *Zp);

	      X[2*i][0] *= mmf;	X[2*i][1] *= mmf;
	      X[2*i+1][0] *= mmf;	X[2*i+1][1] *= mmf;
	      
	      /* observation vector */
	      img_coord_old (*Xp,*Yp,*Zp, Ex[i], I[i], ap[i], mm, &xp1, &yp1);
	      
	      y[2*i] = xp[i] - xp1;	y[2*i+1] = yp[i] - yp1;
	    }
	  
	  else							/* set values to 0 */
	    {
	      for (j=0; j<3; j++)
		{
		  X[2*i][j] = 0; X[2*i+1][j] = 0;
		}
	      y[2*i] = 0;  y[2*i+1] = 0;
	    }
	}
      
      /********************************************************************/
      
      /* Gaussian algorithm */
      
      ata (X, XtX, 8, 3);
      matinv (XtX, 3);
      atl (Xty, X, y, 8, 3);
      matmul (beta, XtX, Xty, 3,3,1);
      
      stopflag = 1;
      for (i=0; i<3; i++)  
	if (fabs(beta[i]) > 0.005)   stopflag = 0;
      
      *Xp += beta[0];  *Yp += beta[1];  *Zp += beta[2];	itnum++;
    }
  
  /************************************************************************/
   
  /* residuals */
  
  /* number of observations for sigma0 */
  for (i=0, n_obs=0; i<n_img; i++)
    if (xp[i] > -1e10  &&  yp[i] > -1e10)  n_obs += 2;

  matmul (Xbeta, X, beta, 8, 3, 1);
  for (i=0; i<n_obs; i++)
    {
      edach[i] = Xbeta[i] - y[i];
      omega += edach[i] * edach[i];
    }
  sigma0 = sqrt(omega /(n_obs-3));
  for (i=0; i<3; i++)  sigmabeta[i] = sigma0 * sqrt(XtX[i][i]);
  

  mean_sigma0 += sigma0 * sigma0;
  rmsX += sigmabeta[0] * sigmabeta[0];
  rmsY += sigmabeta[1] * sigmabeta[1];
  rmsZ += sigmabeta[2] * sigmabeta[2];

  if (examine == 5)
    {
      /* statistical examinations of Qxx, Qvv */
      
      Sn++;
      Ssigma0 += sigma0;
      
      for (i=0; i<3; i++)  for (j=0; j<3; j++)	SQxx[i][j] += XtX[i][j];
      
      mat_transpose (X, Xt, 8, 3);
      matmul (buf, XtX, Xt, 3,3,8);
      matmul (Qvv, X, buf, 8,3,8);
      for (i=0; i<8; i++)  for (j=0; j<8; j++)	Qvv[i][j] *= -1;
      for (i=0; i<8; i++)  Qvv[i][i] = 1 + Qvv[i][i];
      
      for (i=0; i<8; i++)  for (j=0; j<8; j++)	SQvv[i][j] += Qvv[i][j];
      
      if (Sn == match)
	{
	  printf ("mean values over %d points", Sn);
	  printf ("mean sigma0 = %6.3f micron\n\n", 1000*Ssigma0/Sn);
	  puts ("mean Qxx:");
	  for (i=0; i<3; i++)
	    {
	      for (j=0; j<3; j++) printf ("%6.3f   ", SQxx[i][j]/Sn);
	      printf ("\n");
	    }
	  puts ("mean Qvv:");
	  for (i=0; i<8; i++)
	    {
	      for (j=0; j<8; j++) printf ("%6.3f   ", SQvv[i][j]/Sn);
	      printf ("\n");
	    }
	}
    }
}

void dist_to_ray(x, y, Ex, I, G, ap, mm, Xp,Yp,Zp, dist)

Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
double		x, y,Xp,Yp,Zp, *dist;
{
  double    gX[4],gY[4],gZ[4],a[4],b[4],c[4];
  double    x01,x02,x03,x12,x13,x23;
  double    y01,y02,y03,y12,y13,y23;
  double    z01,z02,z03,z12,z13,z23;
  
    
  	
  *dist=1;

}
 
void pos_from_ray(Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp, dist)

Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp, *dist;
{
  double    gX[4],gY[4],gZ[4],a[4],b[4],c[4];
  double    x01,x02,x03,x12,x13,x23;
  double    y01,y02,y03,y12,y13,y23;
  double    z01,z02,z03,z12,z13,z23;
  
    
  ray_tracing_v2 (x1, y1, Ex[0], I[0], G[0], mmp, &gX[0], &gY[0], &gZ[0], &a[0], &b[0], &c[0]);
  ray_tracing_v2 (x2, y2, Ex[1], I[1], G[1], mmp, &gX[1], &gY[1], &gZ[1], &a[1], &b[1], &c[1]);
  ray_tracing_v2 (x3, y3, Ex[2], I[2], G[2], mmp, &gX[2], &gY[2], &gZ[2], &a[2], &b[2], &c[2]);
  ray_tracing_v2 (x4, y4, Ex[3], I[3], G[3], mmp, &gX[3], &gY[3], &gZ[3], &a[3], &b[3], &c[3]);

  //something of my own to determine X,Y,Z, and rms dist
  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[1], I[1], G[1],      gX[1], gY[1], gZ[1], a[1], b[1], c[1], &x01,&y01,&z01);

  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[2], I[2], G[2],      gX[2], gY[2], gZ[2], a[2], b[2], c[2], &x02,&y02,&z02);

  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x03,&y03,&z03);

  point_line_line(Ex[1], I[1], G[1], mmp, gX[1], gY[1], gZ[1], a[1], b[1], c[1],
	              Ex[2], I[2], G[2],      gX[2], gY[2], gZ[2], a[2], b[2], c[2], &x12,&y12,&z12);

  point_line_line(Ex[1], I[1], G[1], mmp, gX[1], gY[1], gZ[1], a[1], b[1], c[1],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x13,&y13,&z13);

  point_line_line(Ex[2], I[2], G[2], mmp, gX[2], gY[2], gZ[2], a[2], b[2], c[2],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x23,&y23,&z23);

  *Xp=(1./6.)*(x01+x02+x03+x12+x13+x23);  
  *Yp=(1./6.)*(y01+y02+y03+y12+y13+y23);  
  *Zp=(1./6.)*(z01+z02+z03+z12+z13+z23);	
  *dist=1;

}

void det_lsq_3d (Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp)
Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp;
{

	    int     i,count_inner=0,n,m,flag[4];
	    double  d_inner=0.,x,y;
	    double X[4],Y[4],Z[4],a[4],b[4],c[4],dist,dist_error,X_pos[6],Y_pos[6],Z_pos[6],XX,YY,ZZ,si0,sqX,sqY,sqZ;
	    
        //new det_lsq function, bloody fast!
		flag[0]=0;flag[1]=0;flag[2]=0;flag[3]=0;
		if(x1>-999){
			flag[0]=1;
			x = x1 - I[0].xh;
	        y = y1 - I[0].yh;
	        //correct_brown_affin (x, y, ap[0], &x, &y);
		    ray_tracing_v2 (x,y, Ex[0], I[0], G[0], mmp, &X[0], &Y[0], &Z[0], &a[0], &b[0], &c[0]);
		}		
		if(x2>-999){
			flag[1]=1;
			x = x2 - I[1].xh;
	        y = y2 - I[1].yh;
	        //correct_brown_affin (x, y, ap[1], &x, &y);
		    ray_tracing_v2 (x,y, Ex[1], I[1], G[1], mmp, &X[1], &Y[1], &Z[1], &a[1], &b[1], &c[1]);
		}		
		if(x3>-999){
			flag[2]=1;
			x = x3 - I[2].xh;
	        y = y3 - I[2].yh;
	        //correct_brown_affin (x, y, ap[2], &x, &y);
		    ray_tracing_v2 (x,y, Ex[2], I[2], G[2], mmp, &X[2], &Y[2], &Z[2], &a[2], &b[2], &c[2]);
		}		
		if(x4>-999){
			flag[3]=1;
			x = x4 - I[3].xh;
	        y = y4 - I[3].yh;
	        //correct_brown_affin (x, y, ap[3], &x, &y);
		    ray_tracing_v2 (x,y, Ex[3], I[3], G[3], mmp, &X[3], &Y[3], &Z[3], &a[3], &b[3], &c[3]);
		}

		count_inner=0;
		for (n=0;n<n_img;n++){
			for(m=n+1;m<n_img;m++){
				if(flag[n]==1 && flag[m]==1){
                    mid_point(X[n],Y[n],Z[n],a[n],b[n],c[n],X[m],Y[m],Z[m],a[m],b[m],c[m],&dist,&XX,&YY,&ZZ);
                    d_inner += dist;
					X_pos[count_inner]=XX;Y_pos[count_inner]=YY;Z_pos[count_inner]=ZZ;
					count_inner++;
				}
			}
		}
        d_inner/=(double)count_inner;		
		XX=0.;YY=0.;ZZ=0.;
		for(i=0;i<count_inner;i++){
           XX+=X_pos[i]; 
		   YY+=Y_pos[i];
		   ZZ+=Z_pos[i];
		}
		XX/=(double)count_inner;YY/=(double)count_inner;ZZ/=(double)count_inner;
		//end of new det_lsq
		*Xp=XX;
		*Yp=YY;
		*Zp=ZZ;

		//statistics
		si0=0.;sqX=0.;sqY=0.;sqZ=0.;
		for(i=0;i<count_inner;i++){
           si0+=pow(X_pos[i]-XX,2.)+pow(Y_pos[i]-YY,2.)+pow(Z_pos[i]-ZZ,2.);
           sqX+=pow(X_pos[i]-XX,2.); 
		   sqY+=pow(Y_pos[i]-YY,2.);
		   sqZ+=pow(Z_pos[i]-ZZ,2.);		   
		}
		si0/=(double)count_inner;sqX/=(double)count_inner;sqY/=(double)count_inner;sqZ/=(double)count_inner;
		
		mean_sigma0 += pow(si0,0.5);
        rmsX += pow(sqX,0.5);
        rmsY += pow(sqY,0.5);
        rmsZ += pow(sqZ,0.5);
		//end of statistics

}

void det_lsq (Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp)
Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp;
{
  int		i, j, itnum, stopflag, n_obs;
  double	X[8][3], XtX[3][3], y[8], Xty[3], beta[3],
    Xbeta[8], edach[8], omega=0., sigma0, sigmabeta[3];
  double	Zx, Zy, N, mmf;
  double	xp1, yp1;
  double	xp[4], yp[4];
  Exterior	Ex_t[4];
  double Xp_t,Yp_t,Zp_t;
  double X_t[8][3],cross_p[3],cross_c[3];


  double    multimed_r_nlay ();
  double	multimed_r_nlay_v2 ();
  
  /* static parameters for examination of Qxx, Qvv */
  static int   	Sn;
  static double	Ssigma0, SQxx[3][3], SQvv[8][8];
  double	Xt[3][8], buf[3][8], Qvv[8][8];
  
  xp[0] = x1;	xp[1] = x2;	xp[2] = x3;	xp[3] = x4;
  yp[0] = y1;	yp[1] = y2;	yp[2] = y3;	yp[3] = y4;

  /* approximate values: rt_is with the uncorrected image coordinates */
  /* of the first 2 rays		*/
  
  /*
    for (i=0; i<n_img; i++)	if (xp[i] > -1e10  &&  yp[i] > -1e10)
    {
    ray_tracing (xp[i], yp[i], Ex[i], I[i],
    mm, &X1, &Y1, &Z1, &a1, &b1, &c1);	break;
    }
    for (i=i+1; i<n_img; i++)	if (xp[i] > -1e10  &&  yp[i] > -1e10)
    {
    ray_tracing (xp[i], yp[i], Ex[i], I[i],
    mm, &X2, &Y2, &Z2, &a2, &b2, &c2);	break;
    }
    
    intersect_rt (X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, Xp, Yp, Zp);
    */  
  
  /* hack due to problems with approx in det_lsq: */
  
  /************************************************************************/  
  /* least squares adjustment */
  itnum = 0;	stopflag = 0;
  while ((stopflag == 0) && (itnum < 10))
    {
      /********************************************************************/
      /* X-matrix  +  observation vector from image coordinates */
      
      for (i=0; i<4; i++)
	{
	  if (xp[i] > -1e10  &&  yp[i] > -1e10)
	    {
	      /* derivatives d(x(i),y(i)) / d(X,Y,Z)   	(Kraus)	*/
	      
	      Zx =  Ex[i].dm[0][0] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][0] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][0] * (*Zp - Ex[i].z0);
	      Zy =  Ex[i].dm[0][1] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][1] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][1] * (*Zp - Ex[i].z0);
	      N  =  Ex[i].dm[0][2] * (*Xp - Ex[i].x0)
		+ Ex[i].dm[1][2] * (*Yp - Ex[i].y0)
		+ Ex[i].dm[2][2] * (*Zp - Ex[i].z0);
	      
	      for (j=0; j<3; j++)
		X[2*i][j] = ((-I[i].cc)/(N*N))
		  * (N*Ex[i].dm[j][0] - Zx*Ex[i].dm[j][2]);
	      for (j=0; j<3; j++)
		X[2*i+1][j] = ((-I[i].cc)/(N*N))
		  * (N*Ex[i].dm[j][1] - Zy*Ex[i].dm[j][2]);
	            
	      /* multimedia factor (radial shift) */
          //trans

		  trans_Cam_Point(Ex[i],mm,G[i],*Xp,*Yp,*Zp,&Ex_t[i],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
	      mmf = multimed_r_nlay_v2 (Ex_t[i], Ex[i], mm, Xp_t, Yp_t, Zp_t);

          X[2*i][0] *= mmf;	X[2*i][1] *= mmf;///
	      X[2*i+1][0] *= mmf;	X[2*i+1][1] *= mmf;///
	      
	      /* observation vector */
	      img_coord (*Xp,*Yp,*Zp, Ex[i], I[i], G[i], ap[i], mm, &xp1, &yp1);
	      
	      y[2*i] = xp[i] - xp1;	y[2*i+1] = yp[i] - yp1;
	    }
	  
	  else							/* set values to 0 */
	    {
	      for (j=0; j<3; j++)
		{
		  X[2*i][j] = 0; X[2*i+1][j] = 0;
		}
	      y[2*i] = 0;  y[2*i+1] = 0;
	    }
	}
      
      /********************************************************************/
      
      /* Gaussian algorithm */
      
      ata (X, XtX, 8, 3);
      matinv (XtX, 3);
      atl (Xty, X, y, 8, 3);
      matmul (beta, XtX, Xty, 3,3,1);
      
      stopflag = 1;
      for (i=0; i<3; i++)  
	if (fabs(beta[i]) > 0.002)   stopflag = 0;
      
      *Xp += beta[0];  *Yp += beta[1];  *Zp += beta[2];	itnum++;
    }
  
  /************************************************************************/
   
  /* residuals */
  
  /* number of observations for sigma0 */
  for (i=0, n_obs=0; i<n_img; i++)
    if (xp[i] > -1e10  &&  yp[i] > -1e10)  n_obs += 2;

  matmul (Xbeta, X, beta, 8, 3, 1);
  for (i=0; i<n_obs; i++)
    {
      edach[i] = Xbeta[i] - y[i];
      omega += edach[i] * edach[i];
    }
  sigma0 = sqrt(omega /(n_obs-3));
  for (i=0; i<3; i++)  sigmabeta[i] = sigma0 * sqrt(XtX[i][i]);
  

  mean_sigma0 += sigma0 * sigma0;
  rmsX += sigmabeta[0] * sigmabeta[0];
  rmsY += sigmabeta[1] * sigmabeta[1];
  rmsZ += sigmabeta[2] * sigmabeta[2];

  if (examine == 5)
    {
      /* statistical examinations of Qxx, Qvv */
      
      Sn++;
      Ssigma0 += sigma0;
      
      for (i=0; i<3; i++)  for (j=0; j<3; j++)	SQxx[i][j] += XtX[i][j];
      
      mat_transpose (X, Xt, 8, 3);
      matmul (buf, XtX, Xt, 3,3,8);
      matmul (Qvv, X, buf, 8,3,8);
      for (i=0; i<8; i++)  for (j=0; j<8; j++)	Qvv[i][j] *= -1;
      for (i=0; i<8; i++)  Qvv[i][i] = 1 + Qvv[i][i];
      
      for (i=0; i<8; i++)  for (j=0; j<8; j++)	SQvv[i][j] += Qvv[i][j];
      
      if (Sn == match)
	{
	  printf ("mean values over %d points", Sn);
	  printf ("mean sigma0 = %6.3f micron\n\n", 1000*Ssigma0/Sn);
	  puts ("mean Qxx:");
	  for (i=0; i<3; i++)
	    {
	      for (j=0; j<3; j++) printf ("%6.3f   ", SQxx[i][j]/Sn);
	      printf ("\n");
	    }
	  puts ("mean Qvv:");
	  for (i=0; i<8; i++)
	    {
	      for (j=0; j<8; j++) printf ("%6.3f   ", SQvv[i][j]/Sn);
	      printf ("\n");
	    }
	}
    }
}


void det_lsq_3 (Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, Xp, Yp, Zp)
Exterior	Ex[3];
Interior	I[3];
Glass   	G[3];
ap_52		ap[3];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, *Xp,*Yp,*Zp;
{
  int		i, itnum, stopflag;
  double	X[6][3], XtX[3][3], y[6], Xty[3], beta[3],
    Xbeta[6], edach[6], omega=0., sigma0, sigmabeta[3];
  double	Zx, Zy, N, mmf;
  double	xp1, yp1, xp2, yp2, xp3, yp3;
  double	X1,Y1,Z1, a1,b1,c1,X2,Y2,Z2, a2,b2,c2;
  
  double	multimed_r_nlay_v2 ();
  Exterior	Ex_t[3];
  double *Xp_t,*Yp_t,*Zp_t,cross_p[3],cross_c[3];
  double	X_t[6][3];
  
  
  /* approximate values (rt_is with the uncorrected image coordinates) */
  ray_tracing_v2 (x1, y1, Ex[0], I[0],G[0],
	       mm, &X1, &Y1, &Z1, &a1, &b1, &c1);
  ray_tracing_v2 (x2, y2, Ex[1],I[1],G[1],
	       mm, &X2, &Y2, &Z2, &a2, &b2, &c2);
  intersect_rt (X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, Xp, Yp, Zp);
  
  /************************************************************************/
  
  itnum = 0;	stopflag = 0;
  while ((stopflag == 0) && (itnum < 10))
    {
      
      /********************************************************************/
      /* derivatives d(x',y') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[0].dm[0][0] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][0] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][0] * (*Zp - Ex[0].z0);
      Zy =  Ex[0].dm[0][1] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][1] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][1] * (*Zp - Ex[0].z0);
      N  =  Ex[0].dm[0][2] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][2] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][2] * (*Zp - Ex[0].z0);
      
      X[0][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][0] - Zx*Ex[0].dm[0][2]);
      X[0][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][0] - Zx*Ex[0].dm[1][2]);
      X[0][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][0] - Zx*Ex[0].dm[2][2]);
      X[1][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][1] - Zy*Ex[0].dm[0][2]);
      X[1][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][1] - Zy*Ex[0].dm[1][2]);
      X[1][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][1] - Zy*Ex[0].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[0],mm,G[0],*Xp,*Yp,*Zp,&Ex_t[0],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[0], Ex[0], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  
      X[0][0] *= mmf;	X[0][1] *= mmf;///
      X[0][0] *= mmf;	X[0][1] *= mmf;	
      X[1][0] *= mmf;	X[1][1] *= mmf;///
	 
      
      
      /********************************************************************/
      /* derivatives d(x'',y'') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[1].dm[0][0] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][0] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][0] * (*Zp - Ex[1].z0);
      Zy =  Ex[1].dm[0][1] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][1] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][1] * (*Zp - Ex[1].z0);
      N  =  Ex[1].dm[0][2] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][2] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][2] * (*Zp - Ex[1].z0);
      
      X[2][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][0] - Zx*Ex[1].dm[0][2]);
      X[2][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][0] - Zx*Ex[1].dm[1][2]);
      X[2][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][0] - Zx*Ex[1].dm[2][2]);
      X[3][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][1] - Zy*Ex[1].dm[0][2]);
      X[3][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][1] - Zy*Ex[1].dm[1][2]);
      X[3][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][1] - Zy*Ex[1].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[1],mm,G[1],*Xp,*Yp,*Zp,&Ex_t[1],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[1], Ex[1], mm, *Xp_t, *Yp_t, *Zp_t);
	  
      X[2][0] *= mmf;	X[2][1] *= mmf;///
      X[3][0] *= mmf;	X[3][1] *= mmf;///
      
      
      /********************************************************************/
      /* derivatives d(x''',y''') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[2].dm[0][0] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][0] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][0] * (*Zp - Ex[2].z0);
      Zy =  Ex[2].dm[0][1] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][1] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][1] * (*Zp - Ex[2].z0);
      N  =  Ex[2].dm[0][2] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][2] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][2] * (*Zp - Ex[2].z0);
      
      X[4][0] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[0][0] - Zx*Ex[2].dm[0][2]);
      X[4][1] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[1][0] - Zx*Ex[2].dm[1][2]);
      X[4][2] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[2][0] - Zx*Ex[2].dm[2][2]);
      X[5][0] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[0][1] - Zy*Ex[2].dm[0][2]);
      X[5][1] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[1][1] - Zy*Ex[2].dm[1][2]);
      X[5][2] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[2][1] - Zy*Ex[2].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[2],mm,G[2],*Xp,*Yp,*Zp,&Ex_t[2],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[2], Ex[2], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  
      X[4][0] *= mmf;	X[4][1] *= mmf;///
	  
      X[5][0] *= mmf;	X[5][1] *= mmf;///
	      
      
      /********************************************************************/
      
      img_coord (*Xp,*Yp,*Zp, Ex[0],I[0], G[0], ap[0], mm, &xp1, &yp1);
      img_coord (*Xp,*Yp,*Zp, Ex[1],I[1], G[1], ap[1], mm, &xp2, &yp2);
      img_coord (*Xp,*Yp,*Zp, Ex[2],I[2], G[2], ap[2], mm, &xp3, &yp3);
      
      y[0] = x1 - xp1;	y[2] = x2 - xp2;	y[4] = x3 - xp3;
      y[1] = y1 - yp1;	y[3] = y2 - yp2;	y[5] = y3 - yp3;
      
      
      /********************************************************************/
      /* Gaussian algorithm */
      
      ata (X, XtX, 6, 3);
      matinv (XtX, 3);
      atl (Xty, X, y, 6, 3);
      matmul (beta, XtX, Xty, 3,3,1);
      
      stopflag = 1;
      for (i=0; i<3; i++)  
	if (fabs(beta[i]) > 0.005)  stopflag = 0;
      
      *Xp += beta[0];  *Yp += beta[1];  *Zp += beta[2];	itnum++;
    }
  
  
  /************************************************************************/
  
  
  /* residuals */
  
  matmul (Xbeta, X, beta, 6, 3, 1);
  for (i=0; i<6; i++)
    {
      edach[i] = Xbeta[i] - y[i];
      omega += edach[i] * edach[i];
    }
  sigma0 = sqrt(omega /(6-3));
  for (i=0; i<3; i++)  sigmabeta[i] = sigma0 * sqrt(XtX[i][i]);
  
  
  mean_sigma0 += sigma0 * sigma0;
  rmsX += sigmabeta[0] * sigmabeta[0];
  rmsY += sigmabeta[1] * sigmabeta[1];
  rmsZ += sigmabeta[2] * sigmabeta[2];
}


void det_lsq_4 (Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp)
Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp;
{
  int		i, itnum, stopflag;
  double	X[8][3], XtX[3][3], y[8], Xty[3], beta[3],
    Xbeta[8], edach[8], omega=0., sigma0, sigmabeta[3];
  double	Zx, Zy, N, mmf;
  double	xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4;
  double	X1,Y1,Z1, a1,b1,c1,X2,Y2,Z2, a2,b2,c2;
  
  double	multimed_r_nlay_v2 ();
  Exterior	Ex_t[4];
  double	X_t[8][3];
  double *Xp_t,*Yp_t,*Zp_t,cross_p[3],cross_c[3];
  
  
  /* approximate values (rt_is with the uncorrected image coordinates) */
  ray_tracing_v2 (x1, y1, Ex[0], I[0],G[0],
	       mm, &X1, &Y1, &Z1, &a1, &b1, &c1);
  ray_tracing_v2 (x2, y2, Ex[1], I[1],G[1],
	       mm, &X2, &Y2, &Z2, &a2, &b2, &c2);
  intersect_rt (X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, Xp, Yp, Zp);
  
  /************************************************************************/
  
  itnum = 0;	stopflag = 0;
  while ((stopflag == 0) && (itnum < 10))
    {
      
      /********************************************************************/
      /* derivatives d(x',y') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[0].dm[0][0] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][0] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][0] * (*Zp - Ex[0].z0);
      Zy =  Ex[0].dm[0][1] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][1] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][1] * (*Zp - Ex[0].z0);
      N  =  Ex[0].dm[0][2] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][2] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][2] * (*Zp - Ex[0].z0);
      
      X[0][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][0] - Zx*Ex[0].dm[0][2]);
      X[0][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][0] - Zx*Ex[0].dm[1][2]);
      X[0][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][0] - Zx*Ex[0].dm[2][2]);
      X[1][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][1] - Zy*Ex[0].dm[0][2]);
      X[1][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][1] - Zy*Ex[0].dm[1][2]);
      X[1][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][1] - Zy*Ex[0].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[0],mm,G[0],*Xp,*Yp,*Zp,&Ex_t[0],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[0], Ex[0], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  
      X[0][0] *= mmf;	X[0][1] *= mmf;///
      X[1][0] *= mmf;	X[1][1] *= mmf;///
	        
      
      /********************************************************************/
      /* derivatives d(x'',y'') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[1].dm[0][0] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][0] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][0] * (*Zp - Ex[1].z0);
      Zy =  Ex[1].dm[0][1] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][1] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][1] * (*Zp - Ex[1].z0);
      N  =  Ex[1].dm[0][2] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][2] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][2] * (*Zp - Ex[1].z0);
      
      X[2][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][0] - Zx*Ex[1].dm[0][2]);
      X[2][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][0] - Zx*Ex[1].dm[1][2]);
      X[2][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][0] - Zx*Ex[1].dm[2][2]);
      X[3][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][1] - Zy*Ex[1].dm[0][2]);
      X[3][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][1] - Zy*Ex[1].dm[1][2]);
      X[3][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][1] - Zy*Ex[1].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[1],mm,G[1],*Xp,*Yp,*Zp,&Ex_t[1],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[1], Ex[1], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  
      X[2][0] *= mmf;	X[2][1] *= mmf;///
      X[3][0] *= mmf;	X[3][1] *= mmf;///
	 
      
      /********************************************************************/
      /* derivatives d(x''',y''') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[2].dm[0][0] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][0] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][0] * (*Zp - Ex[2].z0);
      Zy =  Ex[2].dm[0][1] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][1] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][1] * (*Zp - Ex[2].z0);
      N  =  Ex[2].dm[0][2] * (*Xp - Ex[2].x0)
	+ Ex[2].dm[1][2] * (*Yp - Ex[2].y0)
	+ Ex[2].dm[2][2] * (*Zp - Ex[2].z0);
      
      X[4][0] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[0][0] - Zx*Ex[2].dm[0][2]);
      X[4][1] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[1][0] - Zx*Ex[2].dm[1][2]);
      X[4][2] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[2][0] - Zx*Ex[2].dm[2][2]);
      X[5][0] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[0][1] - Zy*Ex[2].dm[0][2]);
      X[5][1] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[1][1] - Zy*Ex[2].dm[1][2]);
      X[5][2] = ((-I[2].cc)/(N*N)) * (N*Ex[2].dm[2][1] - Zy*Ex[2].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[2],mm,G[2],*Xp,*Yp,*Zp,&Ex_t[2],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[2], Ex[2], mm, *Xp_t, *Yp_t, *Zp_t);
	  
      X[4][0] *= mmf;	X[4][1] *= mmf;///
      X[5][0] *= mmf;	X[5][1] *= mmf;///
      
      
      /********************************************************************/
      /* derivatives d(x'''',y'''') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[3].dm[0][0] * (*Xp - Ex[3].x0)
	+ Ex[3].dm[1][0] * (*Yp - Ex[3].y0)
	+ Ex[3].dm[2][0] * (*Zp - Ex[3].z0);
      Zy =  Ex[3].dm[0][1] * (*Xp - Ex[3].x0)
	+ Ex[3].dm[1][1] * (*Yp - Ex[3].y0)
	+ Ex[3].dm[2][1] * (*Zp - Ex[3].z0);
      N  =  Ex[3].dm[0][2] * (*Xp - Ex[3].x0)
	+ Ex[3].dm[1][2] * (*Yp - Ex[3].y0)
	+ Ex[3].dm[2][2] * (*Zp - Ex[3].z0);
      
      X[6][0] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[0][0] - Zx*Ex[3].dm[0][2]);
      X[6][1] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[1][0] - Zx*Ex[3].dm[1][2]);
      X[6][2] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[2][0] - Zx*Ex[3].dm[2][2]);
      X[7][0] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[0][1] - Zy*Ex[3].dm[0][2]);
      X[7][1] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[1][1] - Zy*Ex[3].dm[1][2]);
      X[7][2] = ((-I[3].cc)/(N*N)) * (N*Ex[3].dm[2][1] - Zy*Ex[3].dm[2][2]);
      //trans
	  trans_Cam_Point(Ex[3],mm,G[3],*Xp,*Yp,*Zp,&Ex_t[3],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[3], Ex[3], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  
      X[6][0] *= mmf;	X[6][1] *= mmf;///
      X[7][0] *= mmf;	X[7][1] *= mmf;///
	
      
      /********************************************************************/
      
      img_coord (*Xp,*Yp,*Zp, Ex[0],I[0], G[0], ap[0], mm, &xp1, &yp1);
      img_coord (*Xp,*Yp,*Zp, Ex[1],I[1], G[1], ap[1], mm, &xp2, &yp2);
      img_coord (*Xp,*Yp,*Zp, Ex[2],I[2], G[2], ap[2], mm, &xp3, &yp3);
      img_coord (*Xp,*Yp,*Zp, Ex[3],I[3], G[3], ap[3], mm, &xp4, &yp4);
      
      y[0] = x1 - xp1;	y[1] = y1 - yp1;
      y[2] = x2 - xp2;	y[3] = y2 - yp2;
      y[4] = x3 - xp3;	y[5] = y3 - yp3;
      y[6] = x4 - xp4;	y[7] = y4 - yp4;
      
      
      /********************************************************************/
      /* Gaussian algorithm */
      
      ata (X, XtX, 8, 3);
      matinv (XtX, 3);
      atl (Xty, X, y, 8, 3);
      matmul (beta, XtX, Xty, 3,3,1);
      
      stopflag = 1;
      for (i=0; i<3; i++)  
	if (fabs(beta[i]) > 0.005)  stopflag = 0;
      
      *Xp += beta[0];  *Yp += beta[1];  *Zp += beta[2];	itnum++;
    }
  
  
  /************************************************************************/
  
  
  /* residuals */
  
  matmul (Xbeta, X, beta, 8, 3, 1);
  for (i=0; i<8; i++)
    {
      edach[i] = Xbeta[i] - y[i];
      omega += edach[i] * edach[i];
    }
  sigma0 = sqrt(omega /(8-3));
  for (i=0; i<3; i++)  sigmabeta[i] = sigma0 * sqrt(XtX[i][i]);
  
  
  mean_sigma0 += sigma0 * sigma0;
  rmsX += sigmabeta[0] * sigmabeta[0];
  rmsY += sigmabeta[1] * sigmabeta[1];
  rmsZ += sigmabeta[2] * sigmabeta[2];
}

void det_lsq_2 (Ex, I, G, ap, mm, x1, y1, x2, y2, Xp, Yp, Zp)
Exterior	Ex[2];
Interior	I[2];
Glass   	G[2];
ap_52		ap[2];
mm_np		mm;
double		x1, y1, x2, y2, *Xp,*Yp,*Zp;
{
  int		i, itnum, stopflag;
  double	X[4][3], XtX[3][3], y[4], Xty[3], beta[3],
    Xbeta[4], edach[4], omega=0., sigma0, sigmabeta[3];
  double	Zx, Zy, N, mmf;
  double	xp1, yp1, xp2, yp2;
  double	X1,Y1,Z1, a1,b1,c1, X2,Y2,Z2, a2,b2,c2;
  
  double	multimed_r_nlay_v2 ();
  Exterior	Ex_t[2];
  double	X_t[4][3];
  double *Xp_t,*Yp_t,*Zp_t,cross_p[3],cross_c[3];
  
  
  /* approximate values (rt_is with the uncorrected image coordinates) */
  ray_tracing_v2 (x1, y1, Ex[0], I[0],G[0],
	       mm, &X1, &Y1, &Z1, &a1, &b1, &c1);
  ray_tracing_v2 (x2, y2, Ex[1],I[1],G[1],
	       mm, &X2, &Y2, &Z2, &a2, &b2, &c2);
  intersect_rt (X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, Xp, Yp, Zp);
  /************************************************************************/
  
  itnum = 0;	stopflag = 0;
  while ((stopflag == 0) && (itnum < 10))
    {
      
      /********************************************************************/
      /* derivatives d(x',y') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[0].dm[0][0] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][0] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][0] * (*Zp - Ex[0].z0);
      Zy =  Ex[0].dm[0][1] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][1] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][1] * (*Zp - Ex[0].z0);
      N  =  Ex[0].dm[0][2] * (*Xp - Ex[0].x0)
	+ Ex[0].dm[1][2] * (*Yp - Ex[0].y0)
	+ Ex[0].dm[2][2] * (*Zp - Ex[0].z0);
      
      X[0][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][0] - Zx*Ex[0].dm[0][2]);
      X[0][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][0] - Zx*Ex[0].dm[1][2]);
      X[0][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][0] - Zx*Ex[0].dm[2][2]);
      X[1][0] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[0][1] - Zy*Ex[0].dm[0][2]);
      X[1][1] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[1][1] - Zy*Ex[0].dm[1][2]);
      X[1][2] = ((-I[0].cc)/(N*N)) * (N*Ex[0].dm[2][1] - Zy*Ex[0].dm[2][2]);
      
	  //trans
	  trans_Cam_Point(Ex[0],mm,G[0],*Xp,*Yp,*Zp,&Ex_t[0],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[0], Ex[0], mm, *Xp_t, *Yp_t, *Zp_t);
	  // do stuff
      X[0][0] *= mmf;	X[0][1] *= mmf;	X[1][0] *= mmf;	X[1][1] *= mmf;
      
      
      /********************************************************************/
      /* derivatives d(x'',y'') / d(X,Y,Z)  (Kraus) */
      
      Zx =  Ex[1].dm[0][0] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][0] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][0] * (*Zp - Ex[1].z0);
      Zy =  Ex[1].dm[0][1] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][1] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][1] * (*Zp - Ex[1].z0);
      N  =  Ex[1].dm[0][2] * (*Xp - Ex[1].x0)
	+ Ex[1].dm[1][2] * (*Yp - Ex[1].y0)
	+ Ex[1].dm[2][2] * (*Zp - Ex[1].z0);
      
      X[2][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][0] - Zx*Ex[1].dm[0][2]);
      X[2][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][0] - Zx*Ex[1].dm[1][2]);
      X[2][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][0] - Zx*Ex[1].dm[2][2]);
      X[3][0] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[0][1] - Zy*Ex[1].dm[0][2]);
      X[3][1] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[1][1] - Zy*Ex[1].dm[1][2]);
      X[3][2] = ((-I[1].cc)/(N*N)) * (N*Ex[1].dm[2][1] - Zy*Ex[1].dm[2][2]);
      
	  //trans
	  trans_Cam_Point(Ex[1],mm,G[1],*Xp,*Yp,*Zp,&Ex_t[1],&Xp_t,&Yp_t,&Zp_t,&cross_p,&cross_c);
      mmf = multimed_r_nlay_v2 (Ex_t[1], Ex[1], mm, *Xp_t, *Yp_t, *Zp_t);
	  
	  trans_Cam_Point(Ex[i],mm,G[i],X[2][0],X[2][1],X[2][2],&Ex_t[i],
			          &X_t[2][0],&X_t[2][1],&X_t[2][2],&cross_p,&cross_c);
      X_t[2][0] *= mmf;	X_t[2][1] *= mmf;///
	  back_trans_Point(X_t[2][0],X_t[2][1],X_t[2][2],mm, G[i],cross_p,cross_c,
			          &X[2][0],&X[2][1],&X[2][2]);

	  trans_Cam_Point(Ex[i],mm,G[i],X[3][0],X[3][1],X[3][2],&Ex_t[i],
			          &X_t[3][0],&X_t[3][1],&X_t[3][2],&cross_p,&cross_c);
      X_t[3][0] *= mmf;	X_t[3][1] *= mmf;///
	  back_trans_Point(X_t[3][0],X_t[3][1],X_t[3][2],mm, G[i],cross_p,cross_c,
			          &X[3][0],&X[3][1],&X[3][2]);
      
      
      /********************************************************************/
      
      img_coord (*Xp,*Yp,*Zp, Ex[0],I[0], G[0], ap[0], mm, &xp1, &yp1);
      img_coord (*Xp,*Yp,*Zp, Ex[1],I[1], G[1], ap[1], mm, &xp2, &yp2);
      
      y[0] = x1 - xp1;	y[2] = x2 - xp2;
      y[1] = y1 - yp1;	y[3] = y2 - yp2;
      
      
      /********************************************************************/
      /* Gaussian algorithm */
      
      ata (X, XtX, 4, 3);
      matinv (XtX, 3);
      atl (Xty, X, y, 4, 3);
      matmul (beta, XtX, Xty, 3,3,1);
      
      stopflag = 1;
      for (i=0; i<3; i++)  
	if (fabs(beta[i]) > 0.005)  stopflag = 0;
      
      *Xp += beta[0];  *Yp += beta[1];  *Zp += beta[2];	itnum++;
    }
  
  
  /************************************************************************/  
  /* residuals */
  
  matmul (Xbeta, X, beta, 4, 3, 1);
  for (i=0; i<4; i++)
    {
      edach[i] = Xbeta[i] - y[i];
      omega += edach[i] * edach[i];
    }
  sigma0 = sqrt(omega /(4-3));
  for (i=0; i<3; i++)  sigmabeta[i] = sigma0 * sqrt(XtX[i][i]);
  
  
  mean_sigma0 += sigma0 * sigma0;
  rmsX += sigmabeta[0] * sigmabeta[0];
  rmsY += sigmabeta[1] * sigmabeta[1];
  rmsZ += sigmabeta[2] * sigmabeta[2];
}
