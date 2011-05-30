/****************************************************************************

Routine:	       	sortgrid.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	22.6.88

Description:	       	reads objects, detected by detection etc.,
		       	sorts them with respect to the 13 wires of the grid,
		       	detects missing points
		       	and writes the result to  a new file
			
		       	does not work in each imaginable case !
****************************************************************************/

#include "ptv.h"

void just_plot (interp, Ex, I, G, ap, mm, imx,imy, pix_x,pix_y,
				   nfix,fix, field, n_img)

Tcl_Interp* interp;
Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
int			nfix, imx, imy, field, n_img;
double		pix_x, pix_y;
coord_3d	fix[];


{
  int	       	i, j;
  int	       	intx, inty;
  double       	xp, yp, eps=10.0;
  // target       	old[512]; Alex, 17.09.09, working on Wesleyan data
  target       	old[1024]; 
  
  
  
  
  
  
  
  
  
  /* reproject all calibration plate points into pixel space
     and search a detected target nearby */
  
  for (i=0; i<nfix; i++)
    {
      img_coord (fix[i].x, fix[i].y, fix[i].z,  Ex, I, G, ap, mm, &xp,&yp);
      metric_to_pixel (xp, yp, imx,imy, pix_x,pix_y, &xp, &yp, field);
      
      /* draw projected points for check purpuses */
      
      intx = (int) xp;
      inty = (int) yp;

	  printf ("coord of point %d: %d, %d\n", i,intx,inty);
      
      drawcross (interp, intx, inty, cr_sz+1, n_img, "yellow");
     draw_pnr (interp, intx, inty, fix[i].pnr, n_img, "yellow");
      
      
    }
}



void sortgrid_man (interp, Ex, I, G, ap, mm, imx,imy, pix_x,pix_y,
				   nfix,fix, num,pix, field, n_img)

Tcl_Interp* interp;
Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
int			nfix, num, imx, imy, field, n_img;
double		pix_x, pix_y;
coord_3d	fix[];
target		pix[];

{
  int	       	i, j;
  int	       	intx, inty;
  double       	xp, yp, eps=10.0;
//  target       	old[512];
  target       	old[1024];
  
  
  
  /* copy and re-initialize pixel data before sorting */
  for (i=0; i<num; i++)	old[i] = pix[i];
  for (i=0; i<nfix; i++)
    {
      pix[i].pnr = -999;  pix[i].x = -999;  pix[i].y = -999;
      pix[i].n = 0; pix[i].nx = 0; pix[i].ny = 0;
      pix[i].sumg = 0;
    }
  
  
  fpp = fopen ("parameters/sortgrid.par", "r");
  if (fpp) {
    fscanf (fpp, "%lf", &eps);
    printf ("Sortgrid search radius: %.1f pixel (from sortgrid.par)\n",eps);
    fclose (fpp);
  }
  else {
    printf ("parameters/sortgrid.par does not exist, ");
    printf ("using default search radius 10 pixel\n");
  }
  
  
  /* reproject all calibration plate points into pixel space
     and search a detected target nearby */
  
  for (i=0; i<nfix; i++)
    {
      img_coord (fix[i].x, fix[i].y, fix[i].z,  Ex, I, G, ap, mm, &xp,&yp);
      metric_to_pixel (xp, yp, imx,imy, pix_x,pix_y, &xp, &yp, field);
      
      /* draw projected points for check purpuses */
      
      intx = (int) xp;
      inty = (int) yp;

	  printf ("coord of point %d: %d, %d\n", i,intx,inty);
      
      drawcross (interp, intx, inty, cr_sz+1, n_img, "cyan");
      
      if (xp > -eps  &&  yp > -eps  &&  xp < imx+eps  &&  yp < imy+eps)
	{
	  j = nearest_neighbour_pix (old, num, xp, yp, eps);
	  
	  if (j != -999)
	    {
	      pix[i] = old[j];  pix[i].pnr = fix[i].pnr;
	    }
	}
    }
}

void sortgrid_file (interp, Ex, I, G, ap, mm, imx,imy, pix_x,pix_y,
				   nfix,fix, num,pix, field, n_img)


///Beat and Debashish January 2011

Tcl_Interp* interp;
Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
int			nfix, num, imx, imy, field, n_img;
double		pix_x, pix_y;
coord_3d	fix[];
target		pix[];



{
  int	       	i, j,k, n_sel;
  int	       	intx, inty;
  double       	xp, yp, eps=10.0;
//  target       	old[512];
  target       	old[1024];
  char	file_sort[256];
  int   dummy,detection_pnr[1000];
  
 
  sprintf (file_sort, "for_sortgrid.%1d", n_img);
  
  fpp = fopen_r (file_sort);
  k = 0;
  while ( fscanf (fpp, "%d %lf %lf %lf %d", &fix[k].pnr,
	  &fix[k].x, &fix[k].y, &fix[k].z, &dummy) != EOF){
          detection_pnr[k]=dummy;
		  k++;
  }
  fclose (fpp);

  for (i=0; i<num; i++)	old[i] = pix[i];
  for (i=0; i<nfix; i++)
    {
      pix[i].pnr = -999;  pix[i].x = -999;  pix[i].y = -999;
      pix[i].n = 0; pix[i].nx = 0; pix[i].ny = 0;
      pix[i].sumg = 0;
  }

  for (i=0; i<nfix; i++){ //...all yellow numbers selected in the sortgrtid.x files      
	  //j = nearest_neighbour_pix in detected points (blue number)
	  // 
	  // if (j != -999)
	  //    pix[i] = old[j];                 pix[i].pnr = fix[i].pnr;
	  if(detection_pnr[i] !=-999){
	        pix[i] = old[detection_pnr[i]];  pix[i].pnr = fix[i].pnr;
	  }
  }
}
