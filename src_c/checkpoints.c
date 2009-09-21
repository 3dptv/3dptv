/****************************************************************************

Routine:	       	checkpoints.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	1988
	
Description:	       	calibration quality control from check points
		       	on calibration plate
	
Routines contained:		

****************************************************************************/

#include "ptv.h"
	
void checkpoint_proc (Tcl_Interp* interp)
{
  int       i, count1 = 0, useflag;
  int       intx1, inty1, intx2, inty2;
  double    X, Y, Z;
  double    sigmax1=0, sigmay1=0, sigmaz1=0;
  FILE      *fp1, *fp3;
  coord_3d  fix[256];
    
  puts ("check points");
  
  for (i=0; i<n_img; i++)
    {
      read_ori (&Ex[i], &I[i], img_ori[i]);
      
      fp1 = fopen_r (img_addpar[i]);
      fscanf (fp1,"%lf %lf %lf %lf %lf %lf %lf",
	      &ap[i].k1, &ap[i].k2, &ap[i].k3, &ap[i].p1, &ap[i].p2,
	      &ap[i].scx, &ap[i].she);
      fclose (fp1);
    }
  
  /* read calibration plate coordinates */
  fp3 = fopen_r (fixp_name);
  for (i=0; i<nfix; i++)  fscanf (fp3,"%d %lf %lf %lf",
			       	  &fix[i].pnr,&fix[i].x,&fix[i].y,&fix[i].z);
  fclose (fp3);
  
  /* read, which points shall be used  (those not used for orientation) */
  fp1 = fopen_r ("parameters/orient.par");
  fscanf (fp1,"%d", &useflag);
  fclose (fp1);  
 
  rmsX = 0;	rmsY = 0;	rmsZ = 0;	mean_sigma0 = 0;
    
  for (i=0; i<nfix; i++)
    {    
      /* use (all) (even) (odd) point numbers as check points  */
      switch (useflag)
	{
	case 0:	break;
	case 1:	if ((fix[i].pnr % 2) != 0)	continue;	break;
	case 2:	if ((fix[i].pnr % 2) == 0)	continue;	break;
	case 3:	if ((fix[i].pnr % 3) != 0)	continue;	break;
	}
      
      /* safety check */
      if (n_img > 0  &&  crd[0][i].pnr != fix[i].pnr)	continue;
      if (n_img > 1  &&  crd[1][i].pnr != fix[i].pnr)	continue;
      if (n_img > 2  &&  crd[2][i].pnr != fix[i].pnr)	continue;
      if (n_img > 3  &&  crd[3][i].pnr != fix[i].pnr)	continue;
      
      /* do not use the corner points of plate 85 */
      if (nfix == 85  &&  fix[i].pnr == 1)	continue;
      if (nfix == 85  &&  fix[i].pnr == 7)	continue;
      if (nfix == 85  &&  fix[i].pnr == 43)	continue;
      if (nfix == 85  &&  fix[i].pnr == 49)	continue;
           
      switch (n_img)
	{
	case 2: det_lsq_2 (Ex, I, G, ap, mmp,
			   crd[0][i].x,crd[0][i].y,
			   crd[1][i].x,crd[1][i].y,
			   &X,&Y,&Z);			break;
	case 3: det_lsq_3 (Ex, I, G, ap, mmp,
			   crd[0][i].x,crd[0][i].y,
			   crd[1][i].x,crd[1][i].y,
			   crd[2][i].x,crd[2][i].y,
			   &X,&Y,&Z);			break;
	case 4: det_lsq_4 (Ex, I, G, ap, mmp,
			   crd[0][i].x,crd[0][i].y,
			   crd[1][i].x,crd[1][i].y,
			   crd[2][i].x,crd[2][i].y,
			   crd[3][i].x,crd[3][i].y,
			   &X,&Y,&Z);			break;
	}
           
      /*  compute deviations  */
      sigmax1 += (X - fix[i].x) * (X - fix[i].x);
      sigmay1 += (Y - fix[i].y) * (Y - fix[i].y);
      sigmaz1 += (Z - fix[i].z) * (Z - fix[i].z);
      printf ("%3d  -   %7.3f    %7.3f   %7.3f   |   %6.3f   %6.3f  %6.3f\n",
	      fix[i].pnr, X,Y,Z, X-fix[i].x, Y-fix[i].y, Z-fix[i].z);
      
      
      /* draw residual vectors into img0 */
      intx1 = (int) pix[0][i].x;
      inty1 = (int) pix[0][i].y;
      intx2 = intx1 + (fix[i].x - X)*5*Ex[0].z0/I[0].cc;
      inty2 = inty1 + (fix[i].y - Y)*5*Ex[0].z0/I[0].cc;
      drawvector (interp,intx1, inty1, intx2, inty2, 1, 0, "yellow");
      inty2 = inty1 + (fix[i].z - Z)*5*Ex[0].z0/I[0].cc;
      drawvector (interp, intx1, inty1, intx1, inty2, 1, 0, "blue");
      
      count1++;
    }
  
  
  /*  compute a priori RMS  */
  rmsX = sqrt(rmsX/count1);
  rmsY = sqrt(rmsY/count1);
  rmsZ = sqrt(rmsZ/count1);
  mean_sigma0 = sqrt (mean_sigma0/count1);	
  printf ("RMS from %d checkpoints:\n", count1);
  printf ("a priori => sigma0 = %4.2f micron, RMS = %6.3f/%6.3f/%6.3f mm",
	  mean_sigma0*1000, rmsX, rmsY, rmsZ);
  
  
  /*  compute a posteriori RMS from check points  */
  puts ("\n accuracies from check points:\n");
  sigmax1 = sqrt (sigmax1/count1);
  sigmay1 = sqrt (sigmay1/count1);
  sigmaz1 = sqrt (sigmaz1/count1);
  printf ("a posteriori    =>                RMS = %6.3f/%6.3f/%6.3f mm\n",
	  sigmax1, sigmay1, sigmaz1);
  
  
  sprintf (buf, "%d check points => mx=%6.3f mm, my=%6.3f mm, mz=%6.3f mm",
	   count1, sigmax1, sigmay1, sigmaz1);    puts (buf);
}
