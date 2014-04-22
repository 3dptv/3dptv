#include "ptv.h"

int mouse_proc_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)

{
  int     i, j, click_x, click_y, n, zf, kind;
  double  x, y, xa, ya;
  double  xa12, xb12, ya12, yb12;
  int     k, pt1, intx1, inty1, count, intx2, inty2, pt2;
  candidate cand[maxcand];
  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;
 
  if (zoom_f[0] == 1) {zf = 2;} else { zf = zoom_f[0];}

  click_x = atoi(argv[1]);
  click_y = atoi(argv[2]);

  n = atoi(argv[3]);
  kind = atoi(argv[4]);
  if (examine)	zf *= 2;
  if (argc == 6 ) zf = atoi(argv[5]);
  
  switch (kind) 
    {

/* -------------------------- MIDDLE MOUSE BUTTON ---------------------------------- */

    case 1:
      
      zoom_x[n] = (click_x - imx/2)/zoom_f[n] + zoom_x[n];
      zoom_y[n] = (click_y - imy/2)/zoom_f[n] + zoom_y[n];
      zoom_f[n] *= 2;
      
      if (zoom_x[n] < imx/(2*zoom_f[n]))	  zoom_x[n] = imx/(2*zoom_f[n]);
      if (zoom_x[n] > imx-imx/(2*zoom_f[n]))	  zoom_x[n] = imx-imx/(2*zoom_f[n]);
      if (zoom_y[n] < imy/(2*zoom_f[n]))	  zoom_y[n] = imy/(2*zoom_f[n]);
      if (zoom_y[n] > imy-imy/(2*zoom_f[n]))      zoom_y[n] = imy-imy/(2*zoom_f[n]);
      
      zoom (img[n], zoomimg, zoom_x[n], zoom_y[n], zoom_f[n]);
          
      img_handle = Tk_FindPhoto( interp, "temp");
      Tk_PhotoGetImage (img_handle, &img_block);
      tclimg2cimg (interp, zoomimg, &img_block);
  
      sprintf(val, "newimage %d", n+1);
      Tcl_Eval(interp, val);
 
      break;
          
      
    case 2:	/* zoom images at corresponding position by factor zf */
      
      if (n == 0)
	{
	  /* zoom and show img1 */
	  zoom_x[0] = (click_x-imx/2)/zoom_f[0] + zoom_x[0];
	  zoom_y[0] = (click_y-imy/2)/zoom_f[0] + zoom_y[0];
	  zoom_f[0] = zf;

	  if (zoom_x[0] < imx/(2*zoom_f[0])) zoom_x[0] = imx/(2*zoom_f[0]);
	  if (zoom_x[0] > imx-imx/(2*zoom_f[0])) zoom_x[0] = imx-imx/(2*zoom_f[0]);
	  if (zoom_y[0] < imy/(2*zoom_f[0]))  zoom_y[0] = imy/(2*zoom_f[0]);
	  if (zoom_y[0] > imy-imy/(2*zoom_f[0]))  zoom_y[0] = imy-imy/(2*zoom_f[0]);

	  zoom (img[0], zoomimg, zoom_x[0], zoom_y[0], zoom_f[0]);
 
	  img_handle = Tk_FindPhoto( interp, "temp");
	  Tk_PhotoGetImage (img_handle, &img_block);
	  tclimg2cimg (interp, zoomimg, &img_block);
     
	  sprintf(val, "newimage %d", 1);
	  Tcl_Eval(interp, val);

	  
	  if (argc != 6 ) {
	    mark_detections (interp, 0);  
	    mark_correspondences (interp, 0);
	    mark_corr (interp, 0);
	  }
	  
	  /* transform center of zoomed img1 */
	  x = (double) zoom_x[0];  y = (double) zoom_y[0];
	  pixel_to_metric (x,y, imx,imy, pix_x,pix_y, &x,&y, chfield);
	  x -= I[n].xh;	y -= I[n].yh;
	  correct_brown_affin (x, y, ap[0], &x, &y);
	  
	  for (i=1; i<n_img; i++)
	    {
	      /* calculate middle of epipolar band */
	      /* for zoom positions in the other images */
	      epi_mm (x,y, Ex[0], I[0], Ex[i], I[i], mmp,
		      &xa12, &ya12, &xb12, &yb12);
	      xa = (xa12 + xb12) / 2;  ya = (ya12 + yb12) / 2;
	      distort_brown_affin (xa, ya, ap[i], &xa, &ya);
	      xa += I[i].xh;	ya += I[i].yh;
	      metric_to_pixel (xa,ya, imx,imy, pix_x,pix_y,
			       &xa,&ya, chfield);
	      zoom_x[i] = (int) xa;  zoom_y[i] = (int) ya;
	      zoom_f[i] = zf;
	      if (zoom_x[i] < imx/(2*zoom_f[i]))
		zoom_x[i] = imx/(2*zoom_f[i]);
	      if (zoom_x[i] > imx - imx/(2*zoom_f[i]))
		zoom_x[i] = imx - imx/(2*zoom_f[i]);
	      if (zoom_y[i] < imy/(2*zoom_f[i]))
		zoom_y[i] = imy/(2*zoom_f[i]);
	      if (zoom_y[i] > imy - imy/(2*zoom_f[i]))
		zoom_y[i] = imy - imy/(2*zoom_f[i]);
	      zoom (img[i], zoomimg, zoom_x[i],zoom_y[i],zoom_f[i]);

	      img_handle = Tk_FindPhoto( interp, "temp");
	      Tk_PhotoGetImage (img_handle, &img_block);
	      tclimg2cimg (interp, zoomimg, &img_block);
	      
	      sprintf(val, "newimage %d", i+1);
	      Tcl_Eval(interp, val);
	      
	      if (argc != 6 ) {
		mark_detections (interp, i);  
		mark_correspondences (interp, i);
		mark_corr (interp, i);
	      }
	      
	    }
	} else { 
	  sprintf(val,"Zooming of corresp. areas is only possible by clicking in image 1");
	  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
	  Tcl_Eval(interp, ".text delete 2");
	  Tcl_Eval(interp, ".text insert 2 $tbuf");
	}
      break;


    case 6:
      
      /* Zoom with high zf for tracking */
      zoom_x[n] = (click_x - imx/2)/zoom_f[n] + zoom_x[n];
      zoom_y[n] = (click_y - imy/2)/zoom_f[n] + zoom_y[n];
      zoom_f[n] = atoi(argv[5]);

      if (zoom_x[n] < imx/(2*zoom_f[n]))	  zoom_x[n] = imx/(2*zoom_f[n]);
      if (zoom_x[n] > imx-imx/(2*zoom_f[n]))	  zoom_x[n] = imx-imx/(2*zoom_f[n]);
      if (zoom_y[n] < imy/(2*zoom_f[n]))	  zoom_y[n] = imy/(2*zoom_f[n]);
      if (zoom_y[n] > imy-imy/(2*zoom_f[n]))      zoom_y[n] = imy-imy/(2*zoom_f[n]);
      
      zoom (img[n], zoomimg, zoom_x[n], zoom_y[n], zoom_f[n]);
          
      img_handle = Tk_FindPhoto( interp, "temp");
      Tk_PhotoGetImage (img_handle, &img_block);
      tclimg2cimg (interp, zoomimg, &img_block);
     
      sprintf(val, "newimage %d", n+1);
      Tcl_Eval(interp, val);

      break;
/*----------------------- RIGHT MOUSE BUTTON ------------------------------*/

    case 3: /* generate epipolar line segments */
      
      /* get geometric coordinates of nearest point in img[n] */
      x = (float) (click_x - imx/2)/zoom_f[n] + zoom_x[n];
      y = (float) (click_y - imy/2)/zoom_f[n] + zoom_y[n];

      pixel_to_metric (x,y, imx,imy, pix_x,pix_y, &x,&y, chfield);
      x -= I[n].xh;	y -= I[n].yh;
      correct_brown_affin (x, y, ap[n], &x, &y);
      k = nearest_neighbour_geo (geo[n], num[n], x, y, 0.05);
      if (k == -999)
	{	  
	  sprintf (buf, "no point near click coord ! Click again!");  puts (buf);
	  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	  Tcl_Eval(interp, ".text delete 2");
	  Tcl_Eval(interp, ".text insert 2 $tbuf");
	  return TCL_OK;
	}
      pt1 = geo[n][k].pnr;

      intx1 = (int) ( imx/2 + zoom_f[n] * (pix[n][pt1].x-zoom_x[n]));
      inty1 = (int) ( imy/2 + zoom_f[n] * (pix[n][pt1].y-zoom_y[n]));
      drawcross (interp, intx1, inty1, cr_sz+2, n, "BlueViolet");

      sprintf (buf, "%d %d %d %d %d\n", pt1, pix[n][pt1].nx, pix[n][pt1].ny,
	       pix[n][pt1].n, pix[n][pt1].sumg);  puts (buf);
	       	       	       
	       for (i=0; i<n_img; i++)	 if (i != n)
		 {
		   /* calculate epipolar band in img[i] */
		   epi_mm (geo[n][k].x,geo[n][k].y,
			   Ex[n],I[n], Ex[i],I[i], mmp,
			   &xa12, &ya12, &xb12, &yb12);
		   
		   /* search candidate in img[i] */
		   printf("\ncandidates in img: %d\n", i);
		   find_candidate_plus_msg (geo[i], pix[i], num[i],
					    xa12, ya12, xb12, yb12, eps0,
					    pix[n][pt1].n, pix[n][pt1].nx, pix[n][pt1].ny,
					    pix[n][pt1].sumg, cand, &count, i);

		   distort_brown_affin (xa12,ya12, ap[i], &xa12,&ya12);
		   distort_brown_affin (xb12,yb12, ap[i], &xb12,&yb12);
		   xa12 += I[i].xh;	ya12 += I[i].yh;
		   xb12 += I[i].xh;	yb12 += I[i].yh;
		   metric_to_pixel (xa12, ya12, imx,imy, pix_x,pix_y,
				    &xa12, &ya12, chfield);
		   metric_to_pixel (xb12, yb12, imx,imy, pix_x,pix_y,
				    &xb12, &yb12, chfield);
		   intx1 = (int) ( imx/2 + zoom_f[i] * (xa12 - zoom_x[i]));
		   inty1 = (int) ( imy/2 + zoom_f[i] * (ya12 - zoom_y[i]));
		   intx2 = (int) ( imx/2 + zoom_f[i] * (xb12 - zoom_x[i]));
		   inty2 = (int) ( imy/2 + zoom_f[i] * (yb12 - zoom_y[i]));

		   if ( n == 0 ) sprintf( val,"yellow");
		   if ( n == 1 ) sprintf( val,"green");
		   if ( n == 2 ) sprintf( val,"red");
		   if ( n == 3 ) sprintf( val,"blue");

		   drawvector ( interp, intx1, inty1, intx2, inty2, 1, i, val);

                   for (j=0; j<count; j++)
                     {
                       pt2 = cand[j].pnr;
                       intx2 = (int) ( imx/2 + zoom_f[i] * (pix[i][pt2].x - zoom_x[i]));
                       inty2 = (int) ( imy/2 + zoom_f[i] * (pix[i][pt2].y - zoom_y[i]));
                       drawcross (interp, intx2, inty2, cr_sz+2, i, "orange");
                     }
   
		   
		 }

	       break;

	       	       
    case 4: /* delete points, which should not be used for orientation */

      
      j = kill_in_list (interp, n, num[n], click_x, click_y);
      if (j != -1)
	{
	  num[n] -= 1;
	  sprintf (buf, "point %d deleted", j);  puts (buf);
	  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	  Tcl_Eval(interp, ".text delete 2 2");
	  Tcl_Eval(interp, ".text insert 2 $tbuf");
	}
      else {
	  sprintf (buf, "no point near click coord !");  puts (buf);
	  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
	  Tcl_Eval(interp, ".text delete 2 2");
	  Tcl_Eval(interp, ".text insert 2 $tbuf");
      }
      break;

/*------------------------ LEFT MOUSE BUTTON ------------------------------*/

    case 5: /* measure coordinates and grey value */

      x = (float) (click_x - imx/2)/zoom_f[n] + zoom_x[n];
      y = (float) (click_y - imy/2)/zoom_f[n] + zoom_y[n];

      sprintf (buf, "   %6.2f    %6.2f    %s", x, y, argv[5]); puts (buf);
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");
      break;	       
    }
  return TCL_OK;
}
