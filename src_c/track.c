/*******************************************************************

Routine:	      	track.c

Author/Copyright:     	Jochen Willneff

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
	               	CH - 8093 Zurich

Creation Date:		Beginning: February '98
                        End: far away

Description:   	        Tracking of particles in image- and objectspace

Routines contained:    	trackcorr_c

*******************************************************************/
#include "ptv.h"

void write_added();
void write_addedback();

int trackcorr_c (ClientData clientData, Tcl_Interp* interp,
		 int argc, const char** argv)
{
  char  val[256], buf[256];
  int i, j, h, k, mm, kk, step, okay=0, invol=0;
  int zaehler1, zaehler2,zaehler1_2, zaehler2_2, philf[4][4];
  int count1=0, count2=0, count3=0, count5=0, count6=0, lost =0, zusatz=0;
  int intx0, intx1, inty0, inty1;
  int intx2, inty2,intx3, inty3;
  int quali=0, help=-1;
  int sl;
  double x1[4], y1[4], x2[4], y2[4], angle, acc, angle0, acc0, lmax, dl;
  double xr[4], xl[4], yd[4], yu[4], angle1, acc1;
  double X1, Y1, Z1, X0, Y0, Z0, X2, Y2, Z2;
  double X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6;
  double xp[4], yp[4], xc[4], yc[4], xn[4], yn[4];
  double rr, Ymin=0, Ymax=0;
  double npart=0, nlinks=0;
  
  double sliceposition=0, leftlimit=0, rightlimit=0;
  double zdimens, z_centerslice[41], deltaoffset, deltazeta;
  int reverse=1;
  int prevsl, nextsl;
  int seqnumber;
  double intersection[40];



  foundpix *w, *wn, p16[16], p16_bis[16];

  display = atoi(argv[1]);

  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 \"Track established correspondences\"");

  /* read data */
  readseqtrackcrit ();
 


  zdimens=Zmax_lay[0]-Zmin_lay[0];
  
	
// z_centerslice[0]=-111.882;




/*
speedsheet=zdim*scan_par.f_pr
deltadicke=speedsheet/scan_par.f_cam
deltaoffset=deltadicke/2.
*/

deltaoffset = zdimens*scan_par.f_pr/scan_par.f_cam/2.;

for (j=1; j<scan_par.N_slice+1; j++) ///doesn't use them; centres are (intersect_(i+1)-intersect_(i))/*2
	{
		sheet_offset(j);

		
		if (reverse==1)
		{
		z_centerslice[j-1] = Zmin_lay[0] + zdimens*0.5 - scan_par.offset + deltaoffset;
		}
		else
		{
        z_centerslice[j-1] = Zmin_lay[0] + zdimens*0.5 + scan_par.offset - deltaoffset;
        }			
	}




  /*Alloc space, if checkflag for mega, c4, t4 is zero */
  if (!trackallocflag)
    {
      for (i=0; i<4; i++)
	{
		for (j=0; j<scan_par.N_slice; j++)
	  {
	  mega_slice[i][j]=(P *) calloc(sizeof(P),M);	
	  
	  c4_slice[i][j]=(corres *) calloc(sizeof(corres),M);
	  
	  for (k=0; k<4; k++) { t4_slice[i][k][j]=(target *) calloc(sizeof (target),M);}
	  }
	  }
	trackallocflag=1;
    }
  
    

  /*load again first data sets*/
  step = seq_first;
  read_ascii_data_first(step);
  rotate_dataset_first();
  read_ascii_data_first(step+scan_par.N_slice);
  rotate_dataset_first();
  read_ascii_data_first(step+2*scan_par.N_slice);

  lmax=sqrt((tpar.dvxmin-tpar.dvxmax)*(tpar.dvxmin-tpar.dvxmax)
	    +(tpar.dvymin-tpar.dvymax)*(tpar.dvymin-tpar.dvymax)
	    +(tpar.dvzmin-tpar.dvzmax)*(tpar.dvzmin-tpar.dvzmax));

  volumedimension (&X_lay[1], &X_lay[0], &Ymax, &Ymin, &Zmax_lay[1], &Zmin_lay[0]);



  fpp = fopen_r ("parameters/intersect.par");
  //fp1 = fopen_r ("parameters/slice_centre.par");
  



  /* sequence loop */
  for (step = seq_first; step < seq_last; step=step+scan_par.N_slice)
    {
      sprintf (buf, "Time step: %d, seqnr: %d, Particle info:", (step- seq_first), step);
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");

      count1=0; lost=0; zusatz=0;

      count5=count6=0;  

      /* read intersections */

	  fscanf (fpp, "%ld\t", &seqnumber);

      for (i=0;i<scan_par.N_slice-1;i++)
	  {
		  if (i==scan_par.N_slice-2)
		  {
		  fscanf (fpp, "%lf\n", &intersection[i]); 
		  intersection[i]=intersection[i]*1000.;
		  }
		  else
		  {
		  fscanf (fpp, "%lf\t", &intersection[i]); 
		  intersection[i]=intersection[i]*1000.;
		  }
	  }
         
     // fscanf (fp1, "%ld\t", &seqnumber);

     /* for (i=0;i<scan_par.N_slice;i++)
	  {
		  if (i==scan_par.N_slice-1)
		  {
		  fscanf (fp1, "%lf\n", &z_centerslice[i]);
		  z_centerslice[i]=z_centerslice[i]*1000.;
		  }
		  else
		  {
		  fscanf (fp1, "%lf\t", &z_centerslice[i]);
		  z_centerslice[i]=z_centerslice[i]*1000.;
		  }
	  }*/
	  for (i=0;i<scan_par.N_slice;i++)
      {
		  if (i==0)
		  {
			  z_centerslice[i]=Zmax_lay[0]+(intersection[i]-Zmax_lay[0])/2.;
		  }
          if (i==scan_par.N_slice-1)
		  {
			  z_centerslice[i]=intersection[i-1]+(Zmin_lay[0]-intersection[i-1])/2.;
		  }
		  if (i>0 && i<scan_par.N_slice-1)
		  {
              z_centerslice[i]=intersection[i-1]+(intersection[i]-intersection[i-1])/2.;
		  }
	  }
      /*for (i=0;i<scan_par.N_slice-1;i++)
	  {
	  intersection[i]=z_centerslice[i]+(z_centerslice[i+1]-z_centerslice[i])/2.;
	  }*/

	  /* slice loop */
	for (sl=0;sl<scan_par.N_slice;sl++)
	{

	  if (reverse==0)
	  {
	     prevsl=sl-1;
	     nextsl=sl+1;
	     leftlimit = intersection[sl];
	     rightlimit = intersection[sl+1];	
	  }
	  else
	  {
	     prevsl=sl+1;
	     nextsl=sl-1;
	     if (sl>0 && sl<scan_par.N_slice-1)
	     {
	        leftlimit = intersection[sl-1];
	        rightlimit = intersection[sl];
	     }
	     if (sl==0)
	     {
	        rightlimit=intersection[sl];
	     }
         if (sl==scan_par.N_slice-1)
	     {
	        leftlimit=intersection[sl-1];
	     }
	  }
	    
	   

      /* try to track correspondences from previous 0 - corp, variable h */
      for (h=0; h<m_slice[1][sl]; h++)
	{


	  
	  if (mega_slice[1][sl][h].prio==1)
	  {
	  
      count5++;
		  
	  X1=Y1=Z1=X0=Y0=Z0=X2=Y2=Z2=X5=Y5=Z5=X3=Y3=Z3=X4=Y4=Z4=X6=Y6=Z6=-999;
	  mega_slice[1][sl][h].inlist=0;

	  for (i=0; i<16;i++)
	    {
	      p16[i].ftnr=-1;
	      p16[i].freq=0;
	      p16[i].slice=-1;
          for(j=0;j<n_img;j++) p16[i].whichcam[j] =0;
	    }

		for (i=0; i<16;i++)
	    {
	      p16_bis[i].ftnr=-1;
	      p16_bis[i].freq=0;
          p16_bis[i].slice=-1;
		  for(j=0;j<n_img;j++) p16_bis[i].whichcam[j] =0;
		 }

	  /* 3D-position */
	  X1=mega_slice[1][sl][h].x[0];
	  Y1=mega_slice[1][sl][h].x[1];
	  Z1=mega_slice[1][sl][h].x[2];

	  /* use information from previous to locate new search position
	     and to calculate values for search area */
	  if (mega_slice[1][sl][h].prev>=0) {
	    X0=mega_slice[0][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].x[0];
	    Y0=mega_slice[0][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].x[1];
		Z0=mega_slice[0][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].x[2];
	    X2=2*X1-X0;
	    Y2=2*Y1-Y0;
	    Z2=2*Z1-Z0;

	    for (j=0; j<n_img; j++)
	      {
		img_coord (X2, Y2, Z2, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
		metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
		x1[j]=xn[j];
		y1[j]=yn[j];
	      }
	  } else {  X2=X1; Y2=Y1; Z2=Z1;
	  for (j=0;j<n_img;j++) {
	    if (c4_slice[1][sl][h].p[j] == -1) {
	      img_coord (X2, Y2, Z2, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
	      metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
	      x1[j]=xn[j];
	      y1[j]=yn[j];
	    } else { x1[j]=t4_slice[1][j][sl][c4_slice[1][sl][h].p[j]].x; 
		         y1[j]=t4_slice[1][j][sl][c4_slice[1][sl][h].p[j]].y; } }
	  }

	  /* calculate searchquader and reprojection in image space */
	  searchquader(X2, Y2, Z2, &xr, &xl, &yd, &yu);

	  /* mark search quader in image */
	  /*
	    for (j=0;j<n_img;j++)
	    {
	    intx0 = (int)(imx/2+zoom_f[j]*(x1[j]-xl[j]-zoom_x[j]));
	    inty0 = (int)(imy/2+zoom_f[j]*(y1[j]+yd[j]-zoom_y[j]));
	    intx1 = (int)(imx/2+zoom_f[j]*(x1[j]-xl[j]-zoom_x[j]));
	    inty1 = (int)(imy/2+zoom_f[j]*(y1[j]-yu[j]-zoom_y[j]));
	    intx2 = (int)(imx/2+zoom_f[j]*(x1[j]+xr[j]-zoom_x[j]));
	    inty2 = (int)(imy/2+zoom_f[j]*(y1[j]-yu[j]-zoom_y[j]));
	    intx3 = (int)(imx/2+zoom_f[j]*(x1[j]+xr[j]-zoom_x[j]));
	    inty3 = (int)(imy/2+zoom_f[j]*(y1[j]+yd[j]-zoom_y[j]));

	    drawvector (interp, intx0, inty0, intx1, inty1, 1, j, "white");
	    drawvector (interp, intx1, inty1, intx2, inty2, 1, j, "white");
	    drawvector (interp, intx3, inty3, intx2, inty2, 1, j, "white");
	    drawvector (interp, intx0, inty0, intx3, inty3, 1, j, "white");
	    }
	  */
	  /* search in pix for candidates in next time step */
	  
	  zaehler1_2=0;
	  
	  	
	  /* if searchquader in 2 and 1 */
	  if (Z2 >= z_centerslice[sl])
		 
	  { 
	  for (j=0;j<n_img;j++)
	    {
	    zaehler1 = candsearch_in_pix (t4_slice[2][j][sl], nt4_slice[2][j][sl], x1[j], y1[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][sl][p16[j*4+k].ftnr].prio!=1)
		    {
            p16[j*4+k].ftnr = -1;
			
		    }
		  }
		  }
		  }
		if (sl>0)
		
		{

		zaehler1_2 = candsearch_in_pix (t4_slice[2][j][nextsl], nt4_slice[2][j][nextsl], x1[j], y1[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16_bis[j*4+k].whichcam[j]=1;
		  p16_bis[j*4+k].ftnr=t4_slice[2][j][nextsl][philf[j][k]].tnr;
		  
		 
		  p16_bis[j*4+k].slice=nextsl;
		 
		
		  if(p16_bis[j*4+k].ftnr!=-1)
		  { 
		 
		  if (mega_slice[2][nextsl][p16_bis[j*4+k].ftnr].prio!=1)
		    {
			p16_bis[j*4+k].ftnr=-1;
			
			}
		 
  
		  }
		  }
		  }
	    }
		}
	  }    
	  else       /* if searchquader in 1 and 0 */
	  {
		
	  for (j=0;j<n_img;j++)
	    {
	    zaehler1 = candsearch_in_pix (t4_slice[2][j][sl], nt4_slice[2][j][sl], x1[j], y1[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][sl][p16[j*4+k].ftnr].prio!=1)
		    {
            p16[j*4+k].ftnr = -1;
			
		    }
		  }
		  }
		  }
		if (sl<scan_par.N_slice-1)
		  
		{
		zaehler1_2 = candsearch_in_pix (t4_slice[2][j][prevsl], nt4_slice[2][j][prevsl], x1[j], y1[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16_bis[j*4+k].whichcam[j]=1;
		  p16_bis[j*4+k].ftnr=t4_slice[2][j][prevsl][philf[j][k]].tnr;
		  
		 
		  p16_bis[j*4+k].slice=prevsl;
		 
		  
		  if (p16_bis[j*4+k].ftnr != -1)
		  {
		
		  if (mega_slice[2][prevsl][p16_bis[j*4+k].ftnr].prio!=1)
		    {
			p16_bis[j*4+k].ftnr=-1;
			
			}
		  
		 
		  
		  }
		  }
		  }
	    }
		}
	  }


	  /* end of search in pix */

	  /* fill and sort candidate struct */
	  sortwhatfound(&p16, &zaehler1);
      sortwhatfound(&p16_bis, &zaehler1_2);	  

    

	  w = (foundpix *) calloc ((zaehler1+zaehler1_2), sizeof (foundpix));

	  if ((zaehler1+zaehler1_2) > 0) count2++; 
	  for (i=0; i<zaehler1;i++)
	    {
	      w[i].ftnr = p16[i].ftnr;
	      w[i].freq = p16[i].freq;
	      w[i].slice = p16[i].slice;
		  for (j=0; j<n_img; j++) w[i].whichcam[j] = p16[i].whichcam[j];
	    }
		for (i=0; i<zaehler1_2;i++)
	    {
	      w[i+zaehler1].ftnr = p16_bis[i].ftnr;
	      w[i+zaehler1].freq = p16_bis[i].freq;
	      w[i+zaehler1].slice = p16_bis[i].slice;
		  for (j=0; j<n_img; j++) w[i+zaehler1].whichcam[j] = p16_bis[i].whichcam[j];
	    }
	  /*end of candidate struct */

      zaehler1=zaehler1+zaehler1_2;


  /* ******************************************************* */
	  /* check for what was found */
	  for (mm=0; mm<zaehler1;mm++) /* zaehler1-loop */
	    {
	      /* search for found corr of current the corr in next
		 with predicted location */

	      /*reset p16 value for new search */
	      for (i=0; i<16;i++)
		{
		  p16[i].ftnr=-1;
		  p16[i].freq=0;
		  for(j=0;j<n_img;j++) p16[i].whichcam[j] =0;
		  p16[i].slice =-1;
          }
          for (i=0; i<16;i++)
	    {
	      p16_bis[i].ftnr=-1;
	      p16_bis[i].freq=0;
	      for(j=0;j<n_img;j++) p16_bis[i].whichcam[j] =0;
	   	  p16_bis[i].slice=-1;
          }

	      /* found 3D-position */
	      X3=mega_slice[2][w[mm].slice][w[mm].ftnr].x[0];
	      Y3=mega_slice[2][w[mm].slice][w[mm].ftnr].x[1];
	      Z3=mega_slice[2][w[mm].slice][w[mm].ftnr].x[2];

	      if (mega_slice[1][sl][h].prev>=0) {
		X5=0.5*(5.0*X3-4.0*X1+X0);
		Y5=0.5*(5.0*Y3-4.0*Y1+Y0);
		Z5=0.5*(5.0*Z3-4.0*Z1+Z0);
	      } else {
		X5=2*X3-X1;
		Y5=2*Y3-Y1;
		Z5=2*Z3-Z1; }

	      searchquader(X5, Y5, Z5, &xr, &xl, &yd, &yu);

	      for (j=0;j<n_img;j++) {
		img_coord (X5, Y5, Z5, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
		metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
		x2[j]=xn[j];
		y2[j]=yn[j];
	      }
    
	  /* search for candidates in next time step */
	  zaehler2_2=0;
	  

      /* if searchquader in 2 and 1 */
	  if (Z5>=z_centerslice[sl])
	  {
	  for (j=0;j<n_img;j++)
	    {
	    zaehler2 = candsearch_in_pix (t4_slice[3][j][sl], nt4_slice[3][j][sl], x2[j], y2[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
		  if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  if( t4_slice[3][j][sl][philf[j][k]].tnr != -1)
		  {
          p16[j*4+k].ftnr=t4_slice[3][j][sl][philf[j][k]].tnr;
		  
		 
		  p16[j*4+k].slice=sl;
         
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[3][sl][p16[j*4+k].ftnr].prio!=1)
		    {
            p16[j*4+k].ftnr = -1;
		    }
		  }
		  }
		  }
		  }
	    if (sl>0)
		 
		{
		zaehler2_2 = candsearch_in_pix (t4_slice[3][j][nextsl], nt4_slice[3][j][nextsl], x2[j], y2[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
			  p16_bis[j*4+k].whichcam[j]=1;
		  if( t4_slice[3][j][nextsl][philf[j][k]].tnr != -1)
		  {
		  
		  p16_bis[j*4+k].ftnr=t4_slice[3][j][nextsl][philf[j][k]].tnr;
		  
		 
		  p16_bis[j*4+k].slice=nextsl;
		 
		  
		  if(p16_bis[j*4+k].ftnr!=-1)
		  {
			
		  if (mega_slice[3][nextsl][p16_bis[j*4+k].ftnr].prio!=1)
		    {
			p16_bis[j*4+k].ftnr=-1;
			}
		 
			
		  
		  }
		  }
		  }
		  }
	    }
		}
	  }
	  else    /* if searchquader in 1 and 0 */
	  {
	  for (j=0;j<n_img;j++)
	    {
	    zaehler2 = candsearch_in_pix (t4_slice[3][j][sl], nt4_slice[3][j][sl], x2[j], y2[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
		  if(philf[j][k] != -999)
		  {p16[j*4+k].whichcam[j]=1;
		  if( t4_slice[3][j][sl][philf[j][k]].tnr != -1)
		  {
		  p16[j*4+k].ftnr=t4_slice[3][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  
		   if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[3][sl][p16[j*4+k].ftnr].prio!=1)
		    {
            p16[j*4+k].ftnr = -1;
		    }
		  }
		  }
		  }
		  }
		if (sl<scan_par.N_slice-1)
		 
		{
		zaehler2_2 = candsearch_in_pix (t4_slice[3][j][prevsl], nt4_slice[3][j][prevsl], x2[j], y2[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
		  if(philf[j][k] != -999)
		  {p16_bis[j*4+k].whichcam[j]=1;
		  if( t4_slice[3][j][prevsl][philf[j][k]].tnr != -1)
		  {
		  p16_bis[j*4+k].ftnr=t4_slice[3][j][prevsl][philf[j][k]].tnr;
		
		 
		  p16_bis[j*4+k].slice=prevsl;
		  
		  
		  if (p16_bis[j*4+k].ftnr!=-1)
		  {
		 
		  if (mega_slice[3][prevsl][p16_bis[j*4+k].ftnr].prio)
		    {
			p16_bis[j*4+k].ftnr=-1;
			}
		 
			
		  }
		  }
		  }
		  }
	    }
		}
	  }



	    
	  
	  /* end of search in pix */
	  
	  /* fill and sort candidate struct */

	      sortwhatfound(&p16, &zaehler2);
		  sortwhatfound(&p16_bis, &zaehler2_2);

		 

	      wn = (foundpix *) calloc ((zaehler2+zaehler2_2), sizeof (foundpix));
	      if ((zaehler2+zaehler2_2) > 0) count3++;

	      for (i=0; i<zaehler2;i++)
		{
		  wn[i].ftnr = p16[i].ftnr;
		  wn[i].slice = p16[i].slice;
		  wn[i].freq = p16[i].freq;
		  for (j=0; j<n_img; j++) wn[i].whichcam[j] = p16[i].whichcam[j];
		}
        for (i=0; i<zaehler2_2;i++)
		{
		  wn[zaehler2+i].ftnr = p16_bis[i].ftnr;
		  wn[zaehler2+i].freq = p16_bis[i].freq;
		  wn[zaehler2+i].slice = p16_bis[i].slice;
		  for (j=0; j<n_img; j++) wn[zaehler2+i].whichcam[j] = p16_bis[i].whichcam[j];
		}
		zaehler2=zaehler2+zaehler2_2;

	      /*end of candidate struct */
	      /* ************************************************ */
	      for (kk=0; kk < zaehler2; kk++)  { /* zaehler2-loop */
                 
		X4=mega_slice[3][wn[kk].slice][wn[kk].ftnr].x[0];
		Y4=mega_slice[3][wn[kk].slice][wn[kk].ftnr].x[1];
		Z4=mega_slice[3][wn[kk].slice][wn[kk].ftnr].x[2];

		okay=0; rr=1000000; quali=0; dl=0;
		acc=2*tpar.dacc;angle=2*tpar.dangle;
		acc0=2*tpar.dacc;angle0=2*tpar.dangle;
		acc1=2*tpar.dacc;angle1=2*tpar.dangle;

		/* displacement check */
		if ( tpar.dvxmin < (X4-X3) && (X4-X3) < tpar.dvxmax &&
		     tpar.dvymin < (Y4-Y3) && (Y4-Y3) < tpar.dvymax &&
		     tpar.dvzmin < (Z4-Z3) && (Z4-Z3) < tpar.dvzmax ) { okay=1;
		     /* end displacement check */

		     if ( okay ==1 ) {
		       dl=(sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3))
			   +sqrt((X4-X3)*(X4-X3)+(Y4-Y3)*(Y4-Y3)+(Z4-Z3)*(Z4-Z3)))/2;

		       angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle1, &acc1);

		       if (mega_slice[1][sl][h].prev>=0) {
			 angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle0, &acc0);
		       } else { acc0=acc1; angle0=angle1;}

		       acc=(acc0+acc1)/2; angle=(angle0+angle1)/2;
		       quali=wn[kk].freq+w[mm].freq;
		       rr=1000000;

		       if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
			 {
				
			   rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
			   mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
			   mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=w[mm].ftnr;
			   mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=w[mm].slice;
			   mega_slice[1][sl][h].inlist++;
			   /*
printf("alt mit kk h: %d, X3: %6.3f %6.3f %6.3f, dl: %5.3f, ac:  %5.3f, an: %5.3f, quali: %d, rr: %5.3f\n",
h, X3, Y3, Z3, dl, acc, angle, quali, rr);
			   */
			 }
		       okay=0;
		     }
		}
	      }/* end of zaehler2-loop */
	      okay=0;

	      /* creating new particle position */
	      /* *************************************************************** */

if(1)
{



	      for (j=0;j<n_img;j++) {
		img_coord (X5, Y5, Z5, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
		metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
	      }

	      /* reset img coord because of n_img smaller 4 */
	      for (j=0;j<4;j++) { x2[j]=-1e10; y2[j]=-1e10; }
          
		  help=0;
	      
		  /* search for unused candidates in next time step */
	     
		
		  /* use fix distance to define xl, xr, yu, yd instead of searchquader */
		  xl[j]= xr[j]= yu[j]= yd[j] = 3.0;

		  
		  
		  
		  
		  /* if Z5 in 1 */	 
     	  if (Z5 < leftlimit && Z5 > rightlimit)
		  { 
			  for (j=0;j<n_img;j++)
			  {
              help=0;
			  zaehler2 = candsearch_in_pixrest (t4_slice[3][j][sl], nt4_slice[3][j][sl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler2>0 && philf[j][0]!=-999)
   			    { 
		        x2[j]=t4_slice[3][j][sl][philf[j][0]].x; y2[j]= t4_slice[3][j][sl][philf[j][0]].y;
		        }
			  }
		  }
		 
		  
		  /* if Z5 in 0 */
	     if (Z5<=rightlimit && sl<scan_par.N_slice-1)
	     {
			for (j=0;j<n_img;j++)
			{
		    if (reverse==0)
			{
			help=-1;
			}
			else
			{
			help=1;
			}
			zaehler2 = candsearch_in_pixrest (t4_slice[3][j][prevsl], nt4_slice[3][j][prevsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler2>0  && philf[j][0]!=-999)
			  {
		      x2[j]=t4_slice[3][j][prevsl][philf[j][0]].x; y2[j]= t4_slice[3][j][prevsl][philf[j][0]].y;
			  }
			}
	      }
	  
         /* if Z5 in 2 */
	     if (Z5>=leftlimit && sl>0)
         {
			for (j=0;j<n_img;j++)
			{

		    if (reverse==0)
			{
			help=1;
			}
			else
			{
			help=-1;
			}

			zaehler2 = candsearch_in_pixrest (t4_slice[3][j][nextsl], nt4_slice[3][j][nextsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler2>0  && philf[j][0]!=-999)
			  {
		      x2[j]=t4_slice[3][j][nextsl][philf[j][0]].x; y2[j]= t4_slice[3][j][nextsl][philf[j][0]].y;
			  }
			}
	      }



		
	      quali=0;

	      for (j=0;j<n_img;j++)
		{
		  if (x2[j] !=-1e10 && y2[j] != -1e10) {
		    pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
		  }
		}

	      if ( quali >= 2) {

		X4 = X5; Y4 =Y5; Z4 = Z5;
		invol=0; okay=0;

		det_lsq (Ex, I, ap, mmp,
			 x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X4, &Y4, &Z4);

		/* volume check */
		if ( X_lay[0] < X4 && X4 < X_lay[1] &&
		     Ymin < Y4 && Y4 < Ymax &&
		     Zmin_lay[0] < Z4 && Z4 < Zmax_lay[1]) {invol=1;}

		/* displacement check */
		if ( invol==1 &&
		     tpar.dvxmin < (X3-X4) && (X3-X4) < tpar.dvxmax &&
		     tpar.dvymin < (Y3-Y4) && (Y3-Y4) < tpar.dvymax &&
		     tpar.dvzmin < (Z3-Z4) && (Z3-Z4) < tpar.dvzmax ) { okay=1;
		     /* end displacement check */
		     if (okay == 1) {
		       rr=1000000; dl=0;
		       acc=2*tpar.dacc;angle=2*tpar.dangle;
		       angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle, &acc);
		       dl=(sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3))
			   +sqrt((X4-X3)*(X4-X3)+(Y4-Y3)*(Y4-Y3)+(Z4-Z3)*(Z4-Z3)))/2;

		       if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
			 {
			   rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali+w[mm].freq);
			   /*
printf("neu ohne prev h: %d, X3: %6.3f %6.3f %6.3f, dl: %5.3f, ac:  %5.3f, an: %5.3f, quali: %d, rr: %5.3f\n",
h, X4, Y4, Z4, dl, acc, angle, quali+w[mm].freq, rr);
			   */

			   mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
			   mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=w[mm].ftnr;
			   mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=w[mm].slice;
			   mega_slice[1][sl][h].inlist++;

			   if (tpar.add) {
			     mega_slice[3][sl+help][m[3]].x[0]=X4;
			     mega_slice[3][sl+help][m[3]].x[1]=Y4;
			     mega_slice[3][sl+help][m[3]].x[2]=Z4;
			     mega_slice[3][sl+help][m[3]].prev= -1;
			     mega_slice[3][sl+help][m[3]].next= -2;
			     mega_slice[3][sl+help][m[3]].prio= 1;
			    



			     for (j=0;j<n_img;j++)
			       {
				 c4_slice[3][sl+help][m[3]].p[j]=-1;
				 if(philf[j][0]!=-999) {
				   t4_slice[3][j][sl+help][philf[j][0]].tnr=m_slice[3][sl+help];
				   c4_slice[3][sl+help][m_slice[3][sl+help]].p[j]= philf[j][0];
				   c4_slice[3][sl+help][m_slice[3][sl+help]].nr=m_slice[3][sl+help];
				 }
			       }
			     m_slice[3][sl+help]++; zusatz++; 

             }
			 }
		     } okay=0; }
		invol=0;
	      }

	      quali=0;
}
	      /* end of creating new particle position */
	      /* *************************************************************** */


	      /* try to link if kk is not found/good enough and prev exist */
	      if ( mega_slice[1][sl][h].inlist == 0 && mega_slice[1][sl][h].prev>=0 ) {
		acc=2*tpar.dacc;angle=2*tpar.dangle;
		if ( tpar.dvxmin < (X3-X1) && (X3-X1) < tpar.dvxmax &&
		     tpar.dvymin < (Y3-Y1) && (Y3-Y1) < tpar.dvymax &&
		     tpar.dvzmin < (Z3-Z1) && (Z3-Z1) < tpar.dvzmax ) { okay=1;

		     if ( okay ==1 )
		       {
			 rr=1000000; quali=0;
			 quali=w[mm].freq;
			 angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);
			 dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			     +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

			 if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
			   {
			     rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
			     mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
			     mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=w[mm].ftnr;
			     mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=w[mm].slice;
				 mega_slice[1][sl][h].inlist++;
			     /*		    
printf("alt ohne kk h: %d, X3: %6.3f %6.3f %6.3f, dl: %5.3f, ac:  %5.3f, an: %5.3f, quali: %d, rr: %5.3f\n",
h, X3, Y3, Z3, dl, acc, angle, quali, rr);
			     */     
			   }
		       }
		}
	      }
	      okay=0;

	      free(wn);
	    } /* end of zaehler1-loop */

	  /* ******************************************************************************/
	  /* begin of inlist still zero */

/*
if(0)
{
*/

       if (tpar.add) {
	    if ( mega_slice[1][sl][h].inlist == 0 && mega_slice[1][sl][h].prev>=0 )
	      {
		/*
		  printf("h: %d\n", h);
		  printf("X0: %6.3f  %6.3f  %6.3f, prev: %d\n", X0, Y0, Z0, mega[1][h].prev);
		  printf("X1: %6.3f  %6.3f  %6.3f, h: %d\n", X1, Y1, Z1, h);
		  printf("X2: %6.3f  %6.3f  %6.3f\n", X2, Y2, Z2);
		  printf("X3: %6.3f  %6.3f  %6.3f, ftnr: %d\n", X3, Y3, Z3, w[mm].ftnr);
		*/
		for (j=0;j<n_img;j++) {
		  img_coord (X2, Y2, Z2, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
		  metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
		}
        
        help=0;

		/* reset img coord because of n_img smaller 4 */
		for (j=0;j<4;j++) { x2[j]=-1e10; y2[j]=-1e10; }

		/* search for unused candidates in next time step */
		
		    /* use fix distance to define xl, xr, yu, yd instead of searchquader */
		    xl[j]= xr[j]= yu[j]= yd[j] = 3.0;







           /* if Z2 in 1 */	 
     	  if (Z2< leftlimit && Z2 > rightlimit)
		  {
			for (j=0;j<n_img;j++)
		    {
            zaehler2 = candsearch_in_pixrest (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler2>0)
   			    { 
		        x2[j]=t4_slice[2][j][sl][philf[j][0]].x; y2[j]= t4_slice[2][j][sl][philf[j][0]].y;
		        }
			 help=0;
			}
		  }
		 
		  
		  /* if Z2 in 0 */
	     if (Z2<=rightlimit && sl<scan_par.N_slice-1)
	     {
		    for (j=0;j<n_img;j++)
		    {
			 zaehler2 = candsearch_in_pixrest (t4_slice[2][j][prevsl], nt4_slice[2][j][prevsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler2>0 )
			  {
		      x2[j]=t4_slice[2][j][prevsl][philf[j][0]].x; y2[j]= t4_slice[2][j][prevsl][philf[j][0]].y;
		      }
			  if (reverse==0)
			{
			help=-1;
			}
			else
			{
			help=1;
			}
			
			}
	      }
	  
         /* if Z2 in 2 */
	     if (Z2>=leftlimit && sl > 0)
         {
		    for (j=0;j<n_img;j++)
		    {
			zaehler2 = candsearch_in_pixrest (t4_slice[2][j][nextsl], nt4_slice[2][j][nextsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler2>0 )
			  {
		      x2[j]=t4_slice[2][j][nextsl][philf[j][0]].x; y2[j]= t4_slice[2][j][nextsl][philf[j][0]].y;
		      }

		    if (reverse==0)
			{
			help=1;
			}
			else
			{
			help=-1;
			}
			
			}
	      }



	  
		quali=0;

		for (j=0;j<n_img;j++)
		  {
		    if (x2[j] !=-1e10 && y2[j] != -1e10) {
		      pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
		    }
		  }

		if (quali>=2) {
		  X3 = X2; Y3 =Y2; Z3 = Z2;
		  invol=0; okay=0;

		  det_lsq (Ex, I, ap, mmp,
			   x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X3, &Y3, &Z3);

		  /* in volume check */
		  if ( X_lay[0] < X3 && X3 < X_lay[1] &&
		       Ymin < Y3 && Y3 < Ymax &&
		       Zmin_lay[0] < Z3 && Z3 < Zmax_lay[1]) {invol=1;}

         // if (help=-1 && Z3>rightlimit) {invol=0;}
		 // if (help=1 && Z3<leftlimit) {invol=0;}
         // if (help=0 && (Z3>leftlimit || Z3<rightlimit)) {invol=0;}


		  /* displacement check */
		  if ( invol==1 &&
		       tpar.dvxmin < (X2-X3) && (X2-X3) < tpar.dvxmax &&
		       tpar.dvymin < (Y2-Y3) && (Y2-Y3) < tpar.dvymax &&
		       tpar.dvzmin < (Z2-Z3) && (Z2-Z3) < tpar.dvzmax ) { okay=1;
		       /* end displacement check */
		       if (okay == 1) {
			 rr=1000000; dl=0;
			 acc=2*tpar.dacc;angle=2*tpar.dangle;
			 angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);
			 dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			     +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

			 if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
			   {
			     rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);

			     mega_slice[2][sl+help][m[2]].x[0]=X3;
			     mega_slice[2][sl+help][m[2]].x[1]=Y3;
			     mega_slice[2][sl+help][m[2]].x[2]=Z3;
			     mega_slice[2][sl+help][m[2]].prev= -1;
			     mega_slice[2][sl+help][m[2]].next= -2;
			     mega_slice[2][sl+help][m[2]].prio= 1;

			     mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
			     mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=m_slice[2][sl+help];
			     mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=sl+help;
				 mega_slice[1][sl][h].inlist++;



			     for (j=0;j<n_img;j++)
			       {
				 c4_slice[2][sl+help][m_slice[2][sl+help]].p[j]=-1;
				 if(philf[j][0]!=-999) {
				   t4_slice[2][j][sl+help][philf[j][0]].tnr=m_slice[2][sl+help];
				   c4_slice[2][sl+help][m_slice[2][sl+help]].p[j]= philf[j][0];
				   c4_slice[2][sl+help][m_slice[2][sl+help]].nr=m_slice[2][sl+help];
				 }
			       }
			     m_slice[2][sl+help]++; zusatz++;
			   }
		       } okay=0; }
		  invol=0;
		}
	      }
	  }

	  /* end of inlist still zero */
	  /***********************************/

	  free(w);
	}
	} /* end of h-loop */

      /* sort decis and give preliminary "finaldecis"  */
      for (h=0;h<m_slice[1][sl];h++)
	{
	  if(mega_slice[1][sl][h].inlist > 0 ) {
	    sort_slice(mega_slice[1][sl][h].inlist, &mega_slice[1][sl][h].decis, &mega_slice[1][sl][h].linkdecis,
			&mega_slice[1][sl][h].slices);
      	    mega_slice[1][sl][h].finaldecis=mega_slice[1][sl][h].decis[0];
	    mega_slice[1][sl][h].next=mega_slice[1][sl][h].linkdecis[0];
		mega_slice[1][sl][h].next_slice=mega_slice[1][sl][h].slices[0];
	  }
	}

      /* create links with decision check */
      for (h=0;h<m_slice[1][sl];h++)
	{
	  if(mega_slice[1][sl][h].inlist > 0 ) {
	    /* best choice wasn't used yet, so link is created */
	    if ( mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev == -1)
		{	mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev=h;
		    mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice=sl;

			if(fabs(mega_slice[1][sl][h].x[2]-mega_slice[2][mega_slice[1][sl][h].next_slice]
			[mega_slice[1][sl][h].next].x[2])>1)
			{
				deltazeta=fabs(mega_slice[1][sl][h].x[2]-mega_slice[2][mega_slice[1][sl][h].next_slice]
			[mega_slice[1][sl][h].next].x[2]);
				printf("%.2f",deltazeta);///////
			}

		}

	    /* best choice was already used by mega[2][mega[1][h].next].prev */
	    else {
	      /* check which is the better choice */
	      if ( mega_slice[1][mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice]
		  [mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev].finaldecis > 
			  mega_slice[1][sl][h].finaldecis)
		{
		  /*
		    printf("h ist besser, h: %4d jetzt: %5.3f vorher: %d %5.3f\n",
		    h, mega[1][h].finaldecis, mega[1][mega[2][mega[1][h].next].prev].next,
		    mega[1][mega[2][mega[1][h].next].prev].finaldecis); */

		  if ( mega_slice[1][mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice]
		  [mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev].inlist>1) {
		    /*
		      printf("zweite Wahl fuer %d waere: %5.3f\n",
		      mega[1][mega[2][mega[1][h].next].prev].linkdecis[1],
		      mega[1][mega[2][mega[1][h].next].prev].decis[1]); */ }

		  /* remove link with prev */
		  mega_slice[1][mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice]
		    [mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev].next= -2;
		  mega_slice[1][mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice]
		    [mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev].next_slice=-2;
		  mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev=h;
		  mega_slice[2][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].prev_slice=sl;
		}
	      else {
		/*
		  printf("h ist schlechter, h: %4d jetzt: %5.3f vorher: %d %5.3f\n",
		  h, mega[1][h].finaldecis, mega[1][mega[2][mega[1][h].next].prev].next,
		  mega[1][mega[2][mega[1][h].next].prev].finaldecis); */
		if (mega_slice[1][sl][h].inlist>1) {
		  /*
		    printf("zweite Wahl fuer %d waere: %5.3f\n",
		    mega[1][h].linkdecis[1],
		    mega[1][h].decis[1]);  */ }

		mega_slice[1][sl][h].next=-2;
		mega_slice[1][sl][h].next_slice=-2;  
		}
	    }
	  }
	  if (mega_slice[1][sl][h].next != -2 ) count1++;
	} /* end of creation of links with decision check */

      /* ******** Draw links now ******** */
      if (display)
	for (h=0;h<m_slice[1][sl];h++)
	  {
	    if(mega_slice[1][sl][h].next != -2 ) {

	      strcpy(buf,"");
	      sprintf(buf ,"green");

	      for (j=0;j<n_img;j++)
		{
		  if(c4_slice[1][sl][h].p[j]>0 && c4_slice[2][mega_slice[1][sl][h].next_slice]
		    [mega_slice[1][sl][h].next].p[j]>0) {

		    xp[j]=t4_slice[1][j][sl][c4_slice[1][sl][h].p[j]].x;
		    yp[j]=t4_slice[1][j][sl][c4_slice[1][sl][h].p[j]].y;
		    xc[j]=t4_slice[2][j][mega_slice[1][sl][h].next_slice][c4_slice[2][mega_slice[1][sl][h].next_slice]
			  [mega_slice[1][sl][h].next].p[j]].x;
		    yc[j]=t4_slice[2][j][mega_slice[1][sl][h].next_slice][c4_slice[2][mega_slice[1][sl][h].next_slice]
			  [mega_slice[1][sl][h].next].p[j]].y;
		    predict (xp[j], yp[j], xc[j], yc[j], &xn[j], &yn[j]);

		    if ( ( fabs(xp[j]-zoom_x[j]) < imx/(2*zoom_f[j]))
			 && (fabs(yp[j]-zoom_y[j]) < imy/(2*zoom_f[j])))
		      {
			strcpy(val,"");
			sprintf(val ,"orange");

			intx0 = (int)(imx/2+zoom_f[j]*(xp[j]-zoom_x[j]));
			inty0 = (int)(imy/2+zoom_f[j]*(yp[j]-zoom_y[j]));
			intx1 = (int)(imx/2+zoom_f[j]*(xc[j]-zoom_x[j]));
			inty1 = (int)(imy/2+zoom_f[j]*(yc[j]-zoom_y[j]));
			intx2 = (int)(imx/2+zoom_f[j]*(xn[j]-zoom_x[j]));
			inty2 = (int)(imy/2+zoom_f[j]*(yn[j]-zoom_y[j]));

			drawcross(interp,intx0,inty0,cr_sz,j,"green");
			drawcross(interp,intx1,inty1,cr_sz+1,j,"yellow");
			drawcross(interp,intx2,inty2,cr_sz+1,j,"white");
			drawvector (interp, intx0, inty0, intx1, inty1, 2, j, buf);
			drawvector (interp, intx1, inty1, intx2, inty2, 1, j, "white");

			if (mega_slice[1][sl][h].finaldecis> 0.2) {
			  draw_pnr ( interp, intx0, inty0, h, j, "white");
			  draw_pnr ( interp, intx0, inty0+10, mega_slice[1][sl][h].next, j, val);
			  draw_value (interp, intx0, inty0 + 20,mega_slice[1][sl][h].finaldecis, j, val); }

		      }
		  }
		}
	    }
	  }

      count6=+.4*m_slice[2][sl];

      /* ******** End of Draw links now ******** */
      sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
	       step+sl, count5, count6, count1, count5-count1, zusatz);

      /* for the average of particles and links */
      npart = npart + count5;
      nlinks = nlinks + count1;
      /*
	printf("%s\n", buf);
      */

      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 3");
      Tcl_Eval(interp, ".text insert 3 $tbuf");

      Tcl_Eval(interp, "update idletasks");


    } /* end of slice loop*/

      rotate_dataset_first();
      write_ascii_data_slice(step);
      /*write_added(step);*/

	  if (step<=seq_last-4*scan_par.N_slice) read_ascii_data_first(step+3*scan_par.N_slice); 
	  

	  if (step>=seq_last-3*scan_par.N_slice) break;

    } /* end of sequence loop */

  /* average over master timesteps */
  npart /= ((seq_last-seq_first)/scan_par.N_slice);
  nlinks /= ((seq_last-seq_first)/scan_par.N_slice);
  printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
	  npart, nlinks, npart-nlinks);

  fclose(fpp);
  //fclose(fp1);

  /*rotate_dataset();
  write_ascii_data(step);
  write_added(step);*/

  /* reset of display flag */
  display = 1;

  return TCL_OK;
}








/*     track backwards */

int trackback_c (ClientData clientData, Tcl_Interp* interp,
		 int argc, const char** argv)
{
  char  buf[256];
  int i, j, h, k, step, okay=0, invol=0;
  int zaehler1, philf[4][4];
  int count1=0, count2=0, zusatz=0;
  int quali=0;
  double x2[4], y2[4], angle, acc, lmax, dl;
  double xr[4], xl[4], yd[4], yu[4];
  double X1, Y1, Z1, X0, Y0, Z0, X2, Y2, Z2, X5, Y5, Z5;
  double X3, Y3, Z3, X4, Y4, Z4, X6, Y6, Z6;
  double xn[4], yn[4];
  double rr, Ymin=0, Ymax=0;
  double npart=0, nlinks=0;
  foundpix *w, p16[16], p16_bis[16];

/************/
 
  int sl, zaehler1_2, help, stepshift;
  double sliceposition=0, leftlimit=0, rightlimit=0;
  double zdimens, z_centerslice[41];
  int count5, count6;
  int reverse=1;
  int prevsl, nextsl;

/************/

  display = atoi(argv[1]);

  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 \"Track established correspondences\"");

  /* read data */
  readseqtrackcrit ();


    zdimens=Zmax_lay[0]-Zmin_lay[0];
  
	
//	z_centerslice[0]=-77.299;  


for (j=0; j<scan_par.N_slice; j++)
	{
		sheet_offset(j);

		
		if (reverse==1)
		{
		z_centerslice[j+0] = Zmin_lay[0] + zdimens*0.5 - scan_par.offset;
		}
		else
		{
        z_centerslice[j] = Zmin_lay[0] + zdimens*0.5 + scan_par.offset;
        }			
	}


 



  /*Alloc space, if checkflag for mega, c4, t4 is zero */
  if (!trackallocflag)
    {
      for (i=0; i<4; i++)
	{
		for (j=0; j<scan_par.N_slice; j++)
	  {
	  mega_slice[i][j]=(P *) calloc(sizeof(P),M);	
	  
	  c4_slice[i][j]=(corres *) calloc(sizeof(corres),M);
	  
	  for (k=0; k<4; k++) { t4_slice[i][k][j]=(target *) calloc(sizeof (target),M);}
	  }
	  }
	trackallocflag=1;
    }


  stepshift=(seq_last-seq_first) % (scan_par.N_slice);
  seq_last =seq_last-stepshift;


  /*load again first data sets*/
  step = seq_last;
  read_ascii_data_back(step-scan_par.N_slice);
  rotate_dataset();
  read_ascii_data_back(step-2*scan_par.N_slice);
  rotate_dataset();
  read_ascii_data_back(step-3*scan_par.N_slice);
  rotate_dataset();
  read_ascii_data_back(step-4*scan_par.N_slice);

  lmax=sqrt((tpar.dvxmin-tpar.dvxmax)*(tpar.dvxmin-tpar.dvxmax)
	    +(tpar.dvymin-tpar.dvymax)*(tpar.dvymin-tpar.dvymax)
	    +(tpar.dvzmin-tpar.dvzmax)*(tpar.dvzmin-tpar.dvzmax));

  volumedimension (&X_lay[1], &X_lay[0], &Ymax, &Ymin, &Zmax_lay[1], &Zmin_lay[0]);

  /* sequence loop */
  for (step = seq_last; step >= seq_first+2*scan_par.N_slice; step=step-scan_par.N_slice)
    {
      sprintf (buf, "Time step: %d, seqnr: %d, Particle info:", step- seq_first, step);
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");

      count1=count5=count6=zusatz=0;




   /* slice loop */
   for (sl=1;sl<scan_par.N_slice-1;sl++)
   {


    
	 if (reverse==0)
	  {
	  prevsl=sl-1;
	  nextsl=sl+1;
	  leftlimit = z_centerslice[sl]-(z_centerslice[sl]-z_centerslice[sl+1])/2.;
	  rightlimit = z_centerslice[sl]+(z_centerslice[sl-1]-z_centerslice[sl])/2.;
	  }
	  else
	  {
	  prevsl=sl+1;
	  nextsl=sl-1;
	  leftlimit = z_centerslice[sl]-(z_centerslice[sl]-z_centerslice[sl-1])/2.;
      rightlimit = z_centerslice[sl]+(z_centerslice[sl+1]-z_centerslice[sl])/2.;
	  }	

      for (h=0; h<m_slice[1][sl]; h++)
	{


	



	  if ((mega_slice[1][sl][h].x[2]<leftlimit)&&(mega_slice[1][sl][h].x[2]>rightlimit))
	  {
      count5++;

	  if (mega_slice[1][sl][h].next>=0 && mega_slice[1][sl][h].prev==-1) {
	    X1=Y1=Z1=X0=Y0=Z0=X2=Y2=Z2=X5=Y5=Z5=X3=Y3=Z3=X4=Y4=Z4=X6=Y6=Z6=-999;

	    mega_slice[1][sl][h].inlist=0;
	    for (i=0; i<16;i++)
	      {
		p16[i].ftnr=-1;
		p16[i].freq=0;
		for(j=0;j<n_img;j++) p16[i].whichcam[j] =0;
		p16[i].slice=-1;
	      }

		  for (i=0; i<16;i++)
	    {
	      p16_bis[i].ftnr=-1;
	      p16_bis[i].freq=0;
	      for(j=0;j<n_img;j++) p16_bis[i].whichcam[j] =0;
		  p16_bis[i].slice=-1;
		}

	    /* 3D-position */
	    X1=mega_slice[1][sl][h].x[0];
	    Y1=mega_slice[1][sl][h].x[1];
	    Z1=mega_slice[1][sl][h].x[2];

	    /* use information from previous to locate new search position
	       and to calculate values for search area */
	    X0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[0];
	    Y0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[1];
	    Z0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[2];
	    X2=2*X1-X0;
	    Y2=2*Y1-Y0;
	    Z2=2*Z1-Z0;

	    for (j=0; j<n_img; j++)
	      {
		img_coord (X2, Y2, Z2, Ex[j],I[j], ap[j], mmp, &xn[j], &yn[j]);
		metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
	      }

	    /* calculate searchquader and reprojection in image space */
	    searchquader(X2, Y2, Z2, &xr, &xl, &yd, &yu);


	  zaehler1_2=0;
	  
	  /* if searchquader entirely in 1 */
	  if ((Z2+tpar.dvzmax) < leftlimit && (Z2+tpar.dvzmin) > rightlimit)
	  {
        for (j=0;j<n_img;j++)
	    {
	      zaehler1 = candsearch_in_pix (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		{
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][sl][p16[j*4+k].ftnr].x[2]>leftlimit || mega_slice[2][sl][w[i].ftnr].x[2]<rightlimit)
		    {
            p16[j*4+k].ftnr = -1;
		    }
		  }
		}
		}
	    }
	  }
      /* if searchquader entirely in 0 */
	  if ((Z2+tpar.dvzmax)<rightlimit)
	  {
	  for (j=0;j<n_img;j++)
	    {
	      zaehler1 = candsearch_in_pix (t4_slice[2][j][prevsl], nt4_slice[2][j][prevsl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		{
	      if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][prevsl][philf[j][k]].tnr;
		  p16[j*4+k].slice=prevsl;
		  
		  
		  if(p16[j*4+k].ftnr!=-1)
		  {
		  if (mega_slice[2][prevsl][p16[j*4+k].ftnr].x[2]>rightlimit)
		    {
			p16[j*4+k].ftnr=-1;
			}
		  }
		 }
		}
	    }
	  }
	  /* if searchquader entirely in 2 */
	  if ((Z2+tpar.dvzmin)>leftlimit)
	  {
	  for (j=0;j<n_img;j++)
	    {
	      zaehler1 = candsearch_in_pix (t4_slice[2][j][nextsl], nt4_slice[2][j][nextsl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);

	      for(k=0; k<4; k++)
		{
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][nextsl][philf[j][k]].tnr;
		  p16[j*4+k].slice=nextsl;
		  
		  
		  if(p16[j*4+k].ftnr!=-1)
		  {
		  if (mega_slice[2][nextsl][p16[j*4+k].ftnr].x[2]<leftlimit)
		    {
			p16[j*4+k].ftnr=-1;
			}
		  }
		}
		}
	    }
	  }	
	  /* if searchquader in 2 and 1 */
	  if ((Z2+tpar.dvzmax) > leftlimit && (Z2+tpar.dvzmin) < leftlimit)
	  {
	  for (j=0;j<n_img;j++)
	    {
	    zaehler1 = candsearch_in_pix (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][sl][p16[j*4+k].ftnr].x[2]>leftlimit || mega_slice[2][sl][p16[j*4+k].ftnr].x[2]<rightlimit)
		    {
            p16[j*4+k].ftnr = -1;
		    }
		  }
		  }
		  }
		
		zaehler1_2 = candsearch_in_pix (t4_slice[2][j][nextsl], nt4_slice[2][j][nextsl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16_bis[j*4+k].whichcam[j]=1;
		  p16_bis[j*4+k].ftnr=t4_slice[2][j][nextsl][philf[j][k]].tnr;
		  p16_bis[j*4+k].slice=nextsl;
		  
		  
		  if(p16_bis[j*4+k].ftnr!=-1)
		  {
		  if (mega_slice[2][nextsl][p16_bis[j*4+k].ftnr].x[2]<leftlimit)
		    {
			p16_bis[j*4+k].ftnr=-1;
			}
		  }
		  }
		  }
	    }
	  }


      /* if searchquader in 1 and 0 */
	  if ((Z2+tpar.dvzmax) > rightlimit && (Z2+tpar.dvzmin) < rightlimit)
	  {
	  for (j=0;j<n_img;j++)
	    {
	    zaehler1 = candsearch_in_pix (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		  

	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16[j*4+k].whichcam[j]=1;
		  p16[j*4+k].ftnr=t4_slice[2][j][sl][philf[j][k]].tnr;
		  p16[j*4+k].slice=sl;
		  
		  
		  if (p16[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][sl][p16[j*4+k].ftnr].x[2]>leftlimit || mega_slice[2][sl][p16[j*4+k].ftnr].x[2]<rightlimit)
		    {
            p16[j*4+k].ftnr = -1;
		    }
		  }
		  }
		  }
		
		zaehler1_2 = candsearch_in_pix (t4_slice[2][j][prevsl], nt4_slice[2][j][prevsl], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
	      for(k=0; k<4; k++)
		  {
		  if(philf[j][k] == -999) p16_bis[j*4+k].ftnr=-1;
          if(philf[j][k] != -999)
		  {
		  p16_bis[j*4+k].whichcam[j]=1;
		  p16_bis[j*4+k].ftnr=t4_slice[2][j][prevsl][philf[j][k]].tnr;
		  p16_bis[j*4+k].slice=prevsl;
		  
		  
		  if (p16_bis[j*4+k].ftnr != -1)
		  {
		  if (mega_slice[2][prevsl][p16_bis[j*4+k].ftnr].x[2]>rightlimit)
		    {
			p16_bis[j*4+k].ftnr=-1;
			}
		  }
		  }
		  }
	    }
	  }


	  /* end of search in pix */

	  /* fill and sort candidate struct */
	  sortwhatfound(&p16, &zaehler1);
      sortwhatfound(&p16_bis, &zaehler1_2);	  

	  w = (foundpix *) calloc ((zaehler1+zaehler1_2), sizeof (foundpix));

	  if ((zaehler1+zaehler1_2) > 0) count2++; 
	  for (i=0; i<zaehler1;i++)
	    {
	      w[i].ftnr = p16[i].ftnr;
	      w[i].freq = p16[i].freq;
	      w[i].slice = p16[i].slice;
		  for (j=0; j<n_img; j++) w[i].whichcam[j] = p16[i].whichcam[j];
	    }
		for (i=0; i<zaehler1_2;i++)
	    {
	      w[i+zaehler1].ftnr = p16_bis[i].ftnr;
	      w[i+zaehler1].freq = p16_bis[i].freq;
	      w[i+zaehler1].slice = p16_bis[i].slice;
		  for (j=0; j<n_img; j++) w[i+zaehler1].whichcam[j] = p16_bis[i].whichcam[j];
	    }
	  /*end of candidate struct */

      zaehler1=zaehler1+zaehler1_2;




	    if (zaehler1 > 0) for (i=0; i<zaehler1;i++) {
			X3=mega_slice[2][w[i].slice][w[i].ftnr].x[0];
	        Y3=mega_slice[2][w[i].slice][w[i].ftnr].x[1];
	        Z3=mega_slice[2][w[i].slice][w[i].ftnr].x[2];

	      okay=0; acc=2*tpar.dacc;angle=2*tpar.dangle;rr=1000000;quali=0; dl=0;

	      /* displacement check */
	      if ( tpar.dvxmin < (X1-X3) && (X1-X3) < tpar.dvxmax &&
		   tpar.dvymin < (Y1-Y3) && (Y1-Y3) < tpar.dvymax &&
		   tpar.dvzmin < (Z1-Z3) && (Z1-Z3) < tpar.dvzmax ) { okay=1;
		   /* end displacement check */

		   if ( okay ==1 ) {
		     dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			 +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

		     quali=w[i].freq;
		     angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);

		     /* *********************check link *****************************/
		     if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
		       {
			 rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/quali;

			 mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
             mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=w[i].slice;
			 mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=w[i].ftnr;
			 mega_slice[1][sl][h].inlist++;
			 /*
printf("h: %d, old and ok: X3: %6.3f %6.3f %6.3f, ftnr: %d, prev: %d, next: %d, angle: %6.3f, acc: %6.3f, dl: %6.3f, rr: %6.3f, quali: %d\n", h, X3, Y3, Z3, w[i].ftnr, mega[2][w[i].ftnr].prev, mega[2][w[i].ftnr].next,angle, acc, dl, rr, quali);
			 */
		       }
		   }
	      } okay=0;
	    }

	    free(w);
	    /******************/
	    quali=0;
	    /* reset img coord because of n_img smaller 4 */
	    for (j=0;j<4;j++) { x2[j]=-1e10; y2[j]=-1e10;}

	    /* if old wasn't found try to create new particle position from rest */
	    if (tpar.add) {

		zaehler1=0;

	      if ( mega_slice[1][sl][h].inlist == 0)
		{
		  for (j=0;j<n_img;j++)
		    {
		      /* use fix distance to define xl, xr, yu, yd instead of searchquader */
		      xl[j]= xr[j]= yu[j]= yd[j] = 3.0;


        help=0;

        /* if searchquader entirely or partly in 1 */	 
     	  if ((Z2+tpar.dvzmax) < leftlimit && (Z2+tpar.dvzmin) > rightlimit)
		  {
			for (j=0;j<n_img;j++)
		    {
            zaehler1 = candsearch_in_pixrest (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler1>0)
   			    { 
		        x2[j]=t4_slice[2][j][sl][philf[j][0]].x; y2[j]= t4_slice[2][j][sl][philf[j][0]].y;
		        }
			 help=0;
			}
		  }
		  else if ((Z2+tpar.dvzmax) > leftlimit && (Z2+tpar.dvzmin) < leftlimit)
		  {
			  for (j=0;j<n_img;j++)
		      {
			  zaehler1 = candsearch_in_pixrest (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler1>0)
   			    { 
		        x2[j]=t4_slice[2][j][sl][philf[j][0]].x; y2[j]= t4_slice[2][j][sl][philf[j][0]].y;
		        }
			  help=0;
			  }
		  }
		  else if ((Z2+tpar.dvzmax) > rightlimit && (Z2+tpar.dvzmin) < rightlimit)
		  {
			  for (j=0;j<n_img;j++)
		      {
			  zaehler1 = candsearch_in_pixrest (t4_slice[2][j][sl], nt4_slice[2][j][sl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler1>0)
   			    { 
		        x2[j]=t4_slice[2][j][sl][philf[j][0]].x; y2[j]= t4_slice[2][j][sl][philf[j][0]].y;
		        }
			  help=0;
			  }
		  }
		  
		  /* if searchquader entirely in 0 */
	     if ((Z2+tpar.dvzmax)<rightlimit)
	     {
		    for (j=0;j<n_img;j++)
		    {
			 zaehler1 = candsearch_in_pixrest (t4_slice[2][j][prevsl], nt4_slice[2][j][prevsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler1>0 )
			  {
		      x2[j]=t4_slice[2][j][prevsl][philf[j][0]].x; y2[j]= t4_slice[2][j][prevsl][philf[j][0]].y;
		      }
			help=-1;
			}
	      }
	  
         /* if searchquader entirely in 2 */
	     if ((Z2+tpar.dvzmin)>leftlimit)
         {
		    for (j=0;j<n_img;j++)
		    {
			zaehler1 = candsearch_in_pixrest (t4_slice[2][j][nextsl], nt4_slice[2][j][nextsl], xn[j], yn[j],
						    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		      if(zaehler1>0 )
			  {
		      x2[j]=t4_slice[2][j][nextsl][philf[j][0]].x; y2[j]= t4_slice[2][j][nextsl][philf[j][0]].y;
		      }
			help=1;
			}
	      }



	  
		quali=0;
		
		
		  for (j=0;j<n_img;j++)
		    {
		      if (x2[j] !=-1e10 && y2[j] != -1e10) {
			pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
		      }
		    }

		  if (quali>=2) {
		    X3 = X2; Y3 =Y2; Z3 = Z2;
		    invol=0; okay=0;

		    det_lsq (Ex, I, ap, mmp,
			     x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X3, &Y3, &Z3);

		    /* volume check */
		    if ( X_lay[0] < X3 && X3 < X_lay[1] &&
			 Ymin < Y3 && Y3 < Ymax &&
			 Zmin_lay[0] < Z3 && Z3 < Zmax_lay[1]) {invol=1;}
            
			if (help=-1 && Z3>rightlimit) {invol=0;}
		    if (help=1 && Z3<leftlimit) {invol=0;}
            if (help=0 && (Z3>leftlimit || Z3<rightlimit)) {invol=0;}

		    okay=0; acc=2*tpar.dacc;angle=2*tpar.dangle;rr=1000000; dl=0;

		    /* displacement check */
		    if ( invol==1 &&
			 tpar.dvxmin < (X1-X3) && (X1-X3) < tpar.dvxmax &&
			 tpar.dvymin < (Y1-Y3) && (Y1-Y3) < tpar.dvymax &&
			 tpar.dvzmin < (Z1-Z3) && (Z1-Z3) < tpar.dvzmax ) { okay=1;
			 /* end displacement check */

			 if ( okay ==1 ) {
			   dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			       +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

			   angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);

			   if ( (acc<tpar.dacc && angle<tpar.dangle) || (acc<tpar.dacc/10) )
			     {
			       rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
			  
				   mega_slice[2][sl+help][m[2]].x[0]=X3;
			       mega_slice[2][sl+help][m[2]].x[1]=Y3;
			       mega_slice[2][sl+help][m[2]].x[2]=Z3;
			       mega_slice[2][sl+help][m[2]].prev= -1;
			       mega_slice[2][sl+help][m[2]].next= -2;
			       mega_slice[2][sl+help][m[2]].prio= 1;

			       mega_slice[1][sl][h].decis[mega_slice[1][sl][h].inlist]=rr;
			       mega_slice[1][sl][h].linkdecis[mega_slice[1][sl][h].inlist]=m_slice[2][sl+help];
			       mega_slice[1][sl][h].slices[mega_slice[1][sl][h].inlist]=sl+help;
				   mega_slice[1][sl][h].inlist++;

			   
			       for (j=0;j<n_img;j++)
				 {
				   c4_slice[2][sl+help][m_slice[2][sl+help]].p[j]=-1;
				   if(philf[j][0]!=-999) {
				     t4_slice[2][j][sl+help][philf[j][0]].tnr=m_slice[2][sl+help];
				     c4_slice[2][sl+help][m_slice[2][sl+help]].p[j]= philf[j][0];
				     c4_slice[2][sl+help][m_slice[2][sl+help]].nr=m_slice[2][sl+help];
				   }
				 }
			       mega_slice[1][sl][h].inlist++;
			       m_slice[2][sl+help]++;
			     }
			 } okay =0; }
		    invol=0;
		  }
		}
	    } /* end of if old wasn't found try to create new particle position from rest */
	  }
	 } 
	}/* end if fabs...*/
	}/* end of h-loop */

      /* sort decis  */
      for (h=0;h<m_slice[1][sl];h++)
	{
	  if(mega_slice[1][sl][h].inlist > 0 )
	  { sort_slice(mega_slice[1][sl][h].inlist, &mega_slice[1][sl][h].decis, 
		  &mega_slice[1][sl][h].linkdecis, &mega_slice[1][sl][h].slices); }
	}

      /* create links with decision check */
      count1=0; zusatz=0;
      for (h=0;h<m_slice[1][sl];h++)
	{
	  if (mega_slice[1][sl][h].inlist > 0 ) {
	    /* if old/new and unused prev == -1 and next == -2 link is created */
	    if ( mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev == -1 &&
			 mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].next == -2 )
	      {
		mega_slice[1][sl][h].finaldecis=mega_slice[1][sl][h].decis[0];
		mega_slice[1][sl][h].prev=mega_slice[1][sl][h].linkdecis[0];
        mega_slice[1][sl][h].prev_slice=mega_slice[1][sl][h].slices[0];
		mega_slice[2][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].next=h;
		mega_slice[2][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].next_slice=sl;
		zusatz++;
	      }

	    /* old which link to prev has to be checked */
	    if ( mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev != -1 &&
			 mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].next == -2 )
	      {
		X0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[0];
		Y0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[1];
		Z0=mega_slice[0][mega_slice[1][sl][h].next_slice][mega_slice[1][sl][h].next].x[2];

		X1=mega_slice[1][sl][h].x[0];
		Y1=mega_slice[1][sl][h].x[1];
		Z1=mega_slice[1][sl][h].x[2];

		X3=mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].x[0];
		Y3=mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].x[1];
		Z3=mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].x[2];

		X4=mega_slice[3][mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev_slice]
		       [mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev].x[0];
		Y4=mega_slice[3][mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev_slice]
		       [mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev].x[1];
		Z4=mega_slice[3][mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev_slice]
		       [mega_slice[2][mega_slice[1][sl][h].slices[0]][mega_slice[1][sl][h].linkdecis[0]].prev].x[2];

		X5=0.5*(5.0*X3-4.0*X1+X0);
		Y5=0.5*(5.0*Y3-4.0*Y1+Y0);
		Z5=0.5*(5.0*Z3-4.0*Z1+Z0);

		acc=2*tpar.dacc;angle=2*tpar.dangle;
		angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle, &acc);
		/*
printf("check X0: %6.3f %6.3f %6.3f\n", X0, Y0, Z0);
printf("check X1: %6.3f %6.3f %6.3f\n", X1, Y1, Z1);
printf("check X3: %6.3f %6.3f %6.3f\n", X3, Y3, Z3);
printf("check X4: %6.3f %6.3f %6.3f\n", X4, Y4, Z4);
printf("check X5: %6.3f %6.3f %6.3f\n", X5, Y5, Z5);
printf("check old, new: X3: %6.3f, Y3: %6.3f, Z3: %6.3f, angle: %6.3f, acc: %6.3f\n", X3, Y3, Z3, angle, acc);
*/
		if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) )
		  {
		    mega_slice[1][sl][h].finaldecis=mega_slice[1][sl][h].decis[0];
		    mega_slice[1][sl][h].prev=mega_slice[1][sl][h].linkdecis[0];
			mega_slice[1][sl][h].prev_slice=mega_slice[1][sl][h].slices[0];
		    mega_slice[2][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].next=h;
			mega_slice[2][mega_slice[1][sl][h].prev_slice][mega_slice[1][sl][h].prev].next_slice=sl;
		    zusatz++;
		  }
	      }
	  }

	  if (mega_slice[1][sl][h].prev != -1 ) count1++;
	} /* end of creation of links with decision check */


      count6=+.4*m_slice[2][sl];
      
	  sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
	       step+sl, count5, count6, count1, count5-count1, zusatz);

      /* for the average of particles and links */
      npart = npart + m_slice[1][sl];
      nlinks = nlinks + count1;


}  /* end of slice loop */




     
      /*
	printf("%s\n", buf);
      */

      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 3");
      Tcl_Eval(interp, ".text insert 3 $tbuf");
      Tcl_Eval(interp, "update idletasks");

      rotate_dataset_first();
	  write_ascii_data_slice(step-scan_par.N_slice);

	  /*write_addedback(step);*/

	  if(step>= seq_first+5*scan_par.N_slice) { read_ascii_data_back(step-5*scan_par.N_slice); }
	  else {break;}
  
	  



} /* end of sequence loop */




/* average over master timesteps */
  npart /= ((seq_last-seq_first)/scan_par.N_slice);
  nlinks /= ((seq_last-seq_first)/scan_par.N_slice);
  printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
	  npart, nlinks, npart-nlinks);

 /* rotate_dataset();
  write_ascii_data_slice(step);*/

  /*write_addedback(step);*/

  /* reset of display flag */
  display = 1;
  return TCL_OK;
}
