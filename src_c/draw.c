/****************************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	end of 97
	
Description:	       	different drawing function to interact
                        with Tcl/Tk
	
Routines contained:     drawcross, drawvector, draw_pnr, mark_detections
                        mark_correspondences, mark_corr, mark_track_c

****************************************************************************/
#include "ptv.h"


int drawcross (Tcl_Interp* interp, int x0, int y0, int size, int imgnr, char color[256])
{
  char val[256];
  sprintf(val, "markparticle %d %d %d %d %s", x0, y0, size, imgnr+1, color);
  Tcl_Eval(interp, val );
  return TCL_OK;
}

int drawvector (Tcl_Interp* interp, int x0, int y0, int x1, int y1, int width, int imgnr, char color[256])
{
  char val[256];
  sprintf(val, "drawline %d %d %d %d %d %d %s", x0, y0, x1, y1, width, imgnr+1, color);
  Tcl_Eval(interp, val );
  return TCL_OK;
}

int draw_pnr (Tcl_Interp* interp, int x, int y, int pnr, int imgnr, char color[256])
{
  char val[256];
  sprintf (val, "%d", pnr);
  sprintf(buf, "drawtext %d %d %s %d %s", x, y, val, imgnr+1, color); 
  Tcl_Eval(interp, buf );
  return TCL_OK; 
}

int draw_value (Tcl_Interp* interp, int x, int y, double pnr, int imgnr, char color[256])
{
  char val[256];
  sprintf (val, "%5.3f", pnr);
  sprintf(buf, "drawtext %d %d %s %d %s", x, y, val, imgnr+1, color); 
  Tcl_Eval(interp, buf );
  return TCL_OK; 
}


void mark_detections (Tcl_Interp* interp, int nr)
/* draws crosses for detected points in a displayed image */
{
  int  	i,limx, limy, intx, inty;
  
  if (num[nr] == 0)
    {
      printf ("No points detected");  return;
    }
  limy = imy/(2*zoom_f[nr]);
  limx = imx/(2*zoom_f[nr]);
  for (i=0; i<num[nr]; i++)
    {
      if (   (fabs(pix[nr][i].x-zoom_x[nr]) < limx)
	     && (fabs(pix[nr][i].y-zoom_y[nr]) < limy))
	{
	  intx = (int)(imx/2+zoom_f[nr]*(pix[nr][i].x-zoom_x[nr]));
	  inty = (int)(imy/2+zoom_f[nr]*(pix[nr][i].y-zoom_y[nr]));
	  drawcross (interp, intx, inty, cr_sz , nr, "blue");
	}
    }
}

void mark_correspondences (Tcl_Interp* interp, int nr)
/* draws crosses and numbers for corresponding points in a displayed window */
{
  int  	i,j, pnr, lim, sum, intx, inty;
  double  x, y;
  
  if (match == 0) return;
  
  lim = imx/(2*zoom_f[nr]);
  
  for (i=0; i<match; i++)
    {
      pnr = geo[nr][con[i].p[nr]].pnr;
      if (pnr < 0 || con[i].p[nr] < 0)	continue;
      
      x = pix[nr][pnr].x;  y = pix[nr][pnr].y;
      if ((fabs (x-zoom_x[nr]) < lim) && (fabs (y-zoom_y[nr]) < lim))
	{
	  intx = (int) ( imx/2 + zoom_f[nr] * (x-zoom_x[nr]));
	  inty = (int) ( imy/2 + zoom_f[nr] * (y-zoom_y[nr]));
	  
	  /* check whether quadruplet, triplet or pair -> select color */
	  for (j=0, sum=0; j<4; j++)	if (con[i].p[j] > 0) sum++; 
	  if ( sum == 2 ) sprintf(buf ,"yellow");
	  if ( sum == 3 ) sprintf(buf ,"green");
	  if ( sum == 4 ) sprintf(buf ,"red");

	  drawcross (interp, intx, inty, cr_sz, nr, buf);

	  /* draw point number */
	  if ( examine  && zoom_f[nr] > 2 )
	    {
	      /* number of established correspondence */
	      /*
	      draw_pnr (interp, intx+3 , inty+3, i, nr, "white"); 
	      */
	       } 
	}
    }
}




void mark_corr (Tcl_Interp* interp, int nr)
/* draws crosses and numbers for corresponding points in a displayed window */
{
  int  	i,j, pnr, sum, intx, inty;
  double  x, y;
  
  if (match == 0)	return;
  
  for (i=0; i<match; i++)
    {
      pnr = geo[nr][con[i].p[nr]].pnr;
      if (pnr < 0|| con[i].p[nr] < 0)	continue;
      
      x = pix[nr][pnr].x;  y = pix[nr][pnr].y;
      
      intx = (int) ( imx/2 + zoom_f[nr] * (x-zoom_x[nr]));
      inty = (int) ( imy/2 + zoom_f[nr] * (y-zoom_y[nr]));
      
      /* check whether quadruplet, triplet or pair -> select color */
      for (j=0, sum=0; j<4; j++)	if (con[i].p[j] > 0) { sum++; }
      if ( sum == 2 ) sprintf(buf ,"yellow");
      if ( sum == 3 ) sprintf(buf ,"green");
      if ( sum == 4 ) sprintf(buf ,"red");

      /* i is the number of the established correspondence */
      draw_pnr (interp, intx+5 , inty, i, nr, "white"); 
      /*
      hilfi = geo[nr][con[i].p[nr]].pnr;
      drawcross (interp, intx, inty, cr_sz, nr, buf);
      */
      /* hilfi is the number of the detected point, see in targetlist */
      /*
      draw_pnr (interp, intx+25 , inty-15, hilfi, nr, "red"); 
      */
    }
}

int mark_track_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
/* draws crosses for detected points in a displayed image */
{
  char  seq_name[4][128];
  int   i_img, i_seq, h, intx, inty;

  cr_sz = atoi(Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY));

  fpp = fopen_r ("parameters/sequence.par");
  for (i_img=0; i_img<4; i_img++) 
    { 
      fscanf (fpp, "%s\n", seq_name[i_img]); 
    }
  /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);
  
  sprintf (buf, "Show detected particles "); puts (buf);
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");
  
  /* track sequence */
  for (i_seq=seq_first; i_seq<=seq_last; i_seq++)
    {
      read_ascii_data(i_seq);
      /* treat the cameras one after the other */
      for (i_img=0; i_img<n_img; i_img++)
	{
	  for (h=0; h<nt4[3][i_img]; h++)
	    {
	      if ( ( fabs(t4[3][i_img][h].x-zoom_x[i_img]) < imx/(2*zoom_f[i_img]))
		   && ( fabs(t4[3][i_img][h].y-zoom_y[i_img]) < imy/(2*zoom_f[i_img])) )
		{		    
		  intx = (int)(imx/2+zoom_f[i_img]*(t4[3][i_img][h].x-zoom_x[i_img]));
		  inty = (int)(imy/2+zoom_f[i_img]*(t4[3][i_img][h].y-zoom_y[i_img]));		  
		  if (t4[3][i_img][h].tnr>-1)
		    { 
		      drawcross ( interp, intx, inty, cr_sz+1, i_img, "green");
		      if (zoom_f[i_img] >= 6) {
		      draw_pnr ( interp, intx, inty+10, i_seq, i_img, "orange");
		      draw_pnr ( interp, intx, inty, t4[3][i_img][h].tnr, i_img, "green");
		      }
		    } else { drawcross ( interp, intx, inty, cr_sz, i_img, "blue"); }
		}
	    }
	  Tcl_Eval(interp, "update idletasks");	      
	}
    }

  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");

  return TCL_OK;  
}


int trajectories_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
/* draws crosses for detected points in a displayed image */
{
  int   k, intx1, inty1, intx2, inty2;
  int i, anz1, anz2, m, j;
  FILE *fp1;
  char val[256];
  vector *line1, *line2;
  double color;
  coord_2d p1[4], p2[4];

  cr_sz = atoi(Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY));

  fpp = fopen_r ("parameters/sequence.par");

  /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);
  
  sprintf (buf, "Show trajectories "); puts (buf);
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");


  for (i=seq_first; i<seq_last;i++)
    {
      if (i < 10)             sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/ptv_is.%2d",  i);
      else       sprintf (val, "res/ptv_is.%3d",  i);
 
      fp1 = fopen (val, "r");
      
      color = ((double)(i-seq_first))/((double)(seq_last-2-seq_first));
      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) malloc (anz1 * sizeof (vector));
      for (j=0;j<anz1;j++) {
	fscanf (fp1, "%d\n", &line1[j].p);
	fscanf (fp1, "%d\n", &line1[j].n);
	fscanf (fp1, "%lf\n", &line1[j].x1);
	fscanf (fp1, "%lf\n", &line1[j].y1);
	fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)             sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)       sprintf (val, "res/ptv_is.%2d",  i+1);
      else       sprintf (val, "res/ptv_is.%3d",  i+1);
      
      fp1 = fopen (val, "r");      
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) calloc (anz2, sizeof (vector));
      
      for (j=0;j<anz2;j++) {
	fscanf (fp1, "%d\n", &line2[j].p);
	fscanf (fp1, "%d\n", &line2[j].n);
	fscanf (fp1, "%lf\n", &line2[j].x1);
	fscanf (fp1, "%lf\n", &line2[j].y1);
	fscanf (fp1, "%lf\n", &line2[j].z1);
      }
      fclose (fp1);
      
      for(j=0;j<anz1;j++) { 	
	m = line1[j].n;

	if (m >= 0)  {	  
	  for (k=0; k<n_img; k++)
	    {
	      img_coord (line1[j].x1, line1[j].y1, line1[j].z1, Ex[k],I[k], G[k], ap[k], mmp, &p1[k].x, &p1[k].y);
	      metric_to_pixel (p1[k].x, p1[k].y, imx,imy, pix_x,pix_y, &p1[k].x, &p1[k].y, chfield);
	      
	      img_coord (line2[m].x1, line2[m].y1, line2[m].z1, Ex[k],I[k], G[k], ap[k], mmp, &p2[k].x, &p2[k].y);
	      metric_to_pixel (p2[k].x, p2[k].y, imx,imy, pix_x,pix_y, &p2[k].x, &p2[k].y, chfield); 
	      
	      if ( fabs( p2[k].x-zoom_x[k]) < imx/(2*zoom_f[k])
		   && ( fabs(p2[k].y-zoom_y[k]) < imy/(2*zoom_f[k])) )
		{	
		  intx1 = (int)(imx/2+zoom_f[k]*(p1[k].x-zoom_x[k]));
		  inty1 = (int)(imy/2+zoom_f[k]*(p1[k].y-zoom_y[k]));
		  intx2 = (int)(imx/2+zoom_f[k]*(p2[k].x-zoom_x[k]));
		  inty2 = (int)(imy/2+zoom_f[k]*(p2[k].y-zoom_y[k]));

		  drawcross ( interp, intx1, inty1, cr_sz+1, k, "blue");	 
		  drawcross ( interp, intx2, inty2, cr_sz+1, k, "red");
		  drawvector (interp, intx1, inty1, intx2, inty2, 2, k, "green");
		}	
	    }
	}
	}
	  Tcl_Eval(interp, "update idletasks");	      


      strcpy(val, "");
      free(line1); free(line2);
    }  /* end of sequence loop */
  
  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf"); 
  
  return TCL_OK;
}
