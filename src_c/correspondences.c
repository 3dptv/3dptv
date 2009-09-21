/****************************************************************************

Routine:	       	correspondences.c

Author/Copyright:      	Hans-Gerd Maas

Address:	      	Institute of Geodesy and Photogrammetry
	      		ETH - Hoenggerberg
	      		CH - 8093 Zurich

Creation Date:	       	1988/89

Description:	       	establishment of correspondences for 2/3/4 cameras

****************************************************************************/

#include "ptv.h"

#define maxcand 100

/****************************************************************************/
/*--------------- 4 camera model: consistent quadruplets -------------------*/
/****************************************************************************/

void correspondences_4 (Tcl_Interp* interp)
{
  int 	i,j,k,l,m,n,o,  i1,i2,i3;
  int   count, match0=0, match4=0, match3=0, match2=0, match1=0;
  int 	p1,p2,p3,p4, p31, p41, p42;
  int  	pt1;
  int 	tim[4][nmax];
  int  	intx, inty;
  double       	xa12,ya12,xb12,yb12,X,Y,Z;
  double       	corr;
  candidate   	cand[maxcand];
  n_tupel     	*con0;
  correspond  	*list[4][4];
/* ----------------------------------------------------------------------- */


  /* allocate memory for lists of correspondences */
  for (i1=0; i1<n_img-1; i1++)	for (i2=i1+1; i2<n_img; i2++)
    list[i1][i2] = (correspond *) malloc (num[i1] * sizeof (correspond));


  con0 = (n_tupel *) malloc (4*nmax * sizeof (n_tupel));

  /* ----------------------------------------------------------------------- */


printf("in corres zmin0: %f, zmax0: %f\n", Zmin_lay[0],Zmax_lay[0] );

  /*  initialize ...  */
  sprintf (buf,"Establishing correspondences");
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");

  match=0; match0=0; match2=0;

  for (i1=0; i1<n_img-1; i1++)
    for (i2=i1+1; i2<n_img; i2++)
      for (i=0; i<num[i1]; i++)
	{
	  list[i1][i2][i].p1 = 0;
	  list[i1][i2][i].n = 0;
	}

  for (i=0; i<nmax; i++)
    {
      for (j=0; j<4; j++) tim[j][i] = 0;
      for (j=0; j<4; j++) con0[i].p[j] = -1; con0[i].corr = 0;
    }


  /* -------------if only one cam and 2D--------- */ //by Beat Lüthi June 2007
  if(n_img==1){
	  
	  fp1 = fopen (res_name, "w");
	  fprintf (fp1, "%4d\n", num[0]);
	  for (i=0; i<num[0]; i++){
          o = epi_mm_2D (geo[0][i].x,geo[0][i].y,
		      Ex[0], I[0],  G[0], mmp,
		      &X,&Y,&Z);
          pix[0][geo[0][i].pnr].tnr=i;
		  fprintf (fp1, "%4d", i+1);
		  fprintf (fp1, " %9.3f %9.3f %9.3f", X, Y, Z);
          fprintf (fp1, " %4d", geo[0][i].pnr);
          fprintf (fp1, " %4d", -1);
          fprintf (fp1, " %4d", -1);
          fprintf (fp1, " %4d\n", -1);
	  }
	  fclose (fp1);
	  match1=num[0];
  }
  /* -------------end of only one cam and 2D ------------ */

  /* matching  1 -> 2,3,4  +  2 -> 3,4  +  3 -> 4 */

  for (i1=0; i1<n_img-1; i1++)	for (i2=i1+1; i2<n_img; i2++)
    {
      sprintf (buf, "Establishing correspondences  %d - %d", i1, i2);
      puts (buf);

      /* establish correspondences from num[i1] points of img[i1] to img[i2] */
      for (i=0; i<num[i1]; i++)	if (geo[i1][i].x != -999)
	{
	  /*o = epi_mm (geo[i1][i].x,geo[i1][i].y,
		      Ex[i1], I[i1], G[i1], Ex[i2], I[i2], G[i2], mmp,
		      &xa12, &ya12, &xb12, &yb12);

	  o = epi_mm (xa12, ya12,
		      Ex[i2], I[i2], G[i2], Ex[i1], I[i1], G[i1], mmp,
		      &xa12, &ya12, &xb12, &yb12);*/

      o = epi_mm (geo[i1][i].x,geo[i1][i].y,
		      Ex[i1], I[i1], G[i1], Ex[i2], I[i2], G[i2], mmp,
		      &xa12, &ya12, &xb12, &yb12);
	  

	  /* origin point in the list */
	  p1 = i;  list[i1][i2][p1].p1 = p1;	pt1 = geo[i1][p1].pnr;

	  /* search for a conjugate point in geo[i2] */
	  find_candidate_plus (geo[i2], pix[i2], num[i2],
			       xa12, ya12, xb12, yb12, eps0,
			       pix[i1][pt1].n,pix[i1][pt1].nx,pix[i1][pt1].ny,
			       pix[i1][pt1].sumg, cand, &count, i2);


	  /* write all corresponding candidates to the preliminary list */
	  /* of correspondences */
	  if (count > maxcand)	{ count = maxcand; }
	  for (j=0; j<count; j++)
	    {
	      list[i1][i2][p1].p2[j] = cand[j].pnr;
	      list[i1][i2][p1].corr[j] = cand[j].corr;
	      list[i1][i2][p1].dist[j] = cand[j].tol;
	    }
	  list[i1][i2][p1].n = count;
	}
    }

  /* repair memory fault (!?) */
  for (j=0; j<4; j++) for (i=0; i<nmax; i++) tim[j][i] = 0;


  /* ------------------------------------------------------------------ */
  /* ------------------------------------------------------------------ */

  /* search consistent quadruplets in the list */
  if (n_img == 4)
    {
      puts ("Search consistent quadruplets");
      for (i=0, match0=0; i<num[0]; i++)
	{
	  p1 = list[0][1][i].p1;
	  for (j=0; j<list[0][1][i].n; j++)
	    for (k=0; k<list[0][2][i].n; k++)
	      for (l=0; l<list[0][3][i].n; l++)
		{
		  p2 = list[0][1][i].p2[j];
		  p3 = list[0][2][i].p2[k];
		  p4 = list[0][3][i].p2[l];
		  for (m=0; m<list[1][2][p2].n; m++)
		    for (n=0; n<list[1][3][p2].n; n++)
		      {
			p31 = list[1][2][p2].p2[m];
			p41 = list[1][3][p2].p2[n];
			if (p3 == p31  &&  p4 == p41)
			  for (o=0; o<list[2][3][p3].n; o++)
			    {
			      p42 = list[2][3][p3].p2[o];
			      if (p4 == p42)
				{
				  corr = (list[0][1][i].corr[j]
					  + list[0][2][i].corr[k]
					  + list[0][3][i].corr[l]
					  + list[1][2][p2].corr[m]
					  + list[1][3][p2].corr[n]
					  + list[2][3][p3].corr[o])
				    / (list[0][1][i].dist[j]
				       + list[0][2][i].dist[k]
				       + list[0][3][i].dist[l]
				       + list[1][2][p2].dist[m]
				       + list[1][3][p2].dist[n]
				       + list[2][3][p3].dist[o]);
				  if (corr > corrmin)
				    {
				      /* accept as preliminary match */
				      con0[match0].p[0] = p1;
				      con0[match0].p[1] = p2;
				      con0[match0].p[2] = p3;
				      con0[match0].p[3] = p4;
				      con0[match0++].corr = corr;
				      if (match0 == 4*nmax)	/* security */
					{
					  printf ("Overflow in correspondences:");
					  printf (" > %d matches\n", match0);
					  i = num[0];
					}
				    }
				}
			    }
		      }
		}
	}


      /* -------------------------------------------------------------------- */

      /* sort quadruplets for match quality (.corr) */
      quicksort_con (con0, match0);

      /* -------------------------------------------------------------------- */

      /* take quadruplets from the top to the bottom of the sorted list */
      /* only if none of the points has already been used */
      for (i=0, match=0; i<match0; i++)
	{
	  p1 = con0[i].p[0];	if (p1 > -1)	if (++tim[0][p1] > 1)	continue;
	  p2 = con0[i].p[1];	if (p2 > -1)	if (++tim[1][p2] > 1)	continue;
	  p3 = con0[i].p[2];	if (p3 > -1)	if (++tim[2][p3] > 1)	continue;
	  p4 = con0[i].p[3];	if (p4 > -1)	if (++tim[3][p4] > 1)	continue;
	  con[match++] = con0[i];
	}

      match4 = match;
      sprintf (buf, "%d consistent quadruplets, ", match4);	puts (buf);
    }

  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */

  /* search consistent triplets :  123, 124, 134, 234 */
  if (n_img >= 3)
    {
      puts ("Search consistent triplets");
      match0=0;
      for (i1=0; i1<n_img-2; i1++)
	for (i2=i1+1; i2<n_img-1; i2++)
	  for (i3=i2+1; i3<n_img; i3++)
	    for (i=0; i<num[i1]; i++)
	      {
		p1 = list[i1][i2][i].p1;
		if (p1 > nmax  ||  tim[i1][p1] > 0)	continue;

		for (j=0; j<list[i1][i2][i].n; j++)
		  for (k=0; k<list[i1][i3][i].n; k++)
		    {
		      p2 = list[i1][i2][i].p2[j];
		      if (p2 > nmax  ||  tim[i2][p2] > 0)	continue;
		      p3 = list[i1][i3][i].p2[k];
		      if (p3 > nmax  ||  tim[i3][p3] > 0)	continue;

		      for (m=0; m<list[i2][i3][p2].n; m++)
			{
			  p31 = list[i2][i3][p2].p2[m];
			  if (p3 == p31)
			    {
			      corr = (list[i1][i2][i].corr[j]
				      + list[i1][i3][i].corr[k]
				      + list[i2][i3][p2].corr[m])
				/ (list[i1][i2][i].dist[j]
				   + list[i1][i3][i].dist[k]
				   + list[i2][i3][p2].dist[m]);
 			      if (corr > corrmin)
				{ for (n=0; n<n_img; n++) con0[match0].p[n] = -2;
				  con0[match0].p[i1] = p1;
				  con0[match0].p[i2] = p2;
				  con0[match0].p[i3] = p3;
				  con0[match0++].corr = corr;
				}
			    }
			}
		    }
	      }

      /* ----------------------------------------------------------------------- */

      /* sort triplets for match quality (.corr) */
      quicksort_con (con0, match0);

      /* ----------------------------------------------------------------------- */

      /* pragmatic version: */
      /* take triplets from the top to the bottom of the sorted list */
      /* only if none of the points has already been used */
      for (i=0; i<match0; i++)
	{
	  p1 = con0[i].p[0];	if (p1 > -1)	if (++tim[0][p1] > 1)	continue;
	  p2 = con0[i].p[1];	if (p2 > -1)	if (++tim[1][p2] > 1)	continue;
	  p3 = con0[i].p[2];	if (p3 > -1)	if (++tim[2][p3] > 1)	continue;
	  p4 = con0[i].p[3];	if (p4 > -1  && n_img > 3) if (++tim[3][p4] > 1) continue;

	  con[match++] = con0[i];
	}

      match3 = match - match4;
      sprintf (buf, "%d consistent quadruplets, %d triplets ", match4, match3);
      puts (buf);

      /* repair artifact (?) */
      if (n_img == 3)	for (i=0; i<match; i++)	con[i].p[3] = -1;
    }

  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */

  /* search consistent pairs :  12, 13, 14, 23, 24, 34 */
  /* only if an object model is available or if only 2 images are used */
  if(n_img>1){puts ("Search pairs");}


  match0 = 0;
  for (i1=0; i1<n_img-1; i1++)
    if ( n_img == 2 || (num[0] < 64 && num[1] < 64 && num[2] < 64 && num[3] < 64))
      for (i2=i1+1; i2<n_img; i2++)
	for (i=0; i<num[i1]; i++)
	  {
	    p1 = list[i1][i2][i].p1;
	    if (p1 > nmax  ||  tim[i1][p1] > 0)	continue;

	    /* take only unambigous pairs */
	    if (list[i1][i2][i].n != 1)	continue;

	    p2 = list[i1][i2][i].p2[0];
	    if (p2 > nmax  ||  tim[i2][p2] > 0)	continue;

	    corr = list[i1][i2][i].corr[0] / list[i1][i2][i].dist[0];

	    if (corr > corrmin)
	      {
		con0[match0].p[i1] = p1;
		con0[match0].p[i2] = p2;
		con0[match0++].corr = corr;
	      }
	  }

  /* ----------------------------------------------------------------------- */


  /* sort pairs for match quality (.corr) */
  quicksort_con (con0, match0);

  /* ----------------------------------------------------------------------- */


  /* take pairs from the top to the bottom of the sorted list */
  /* only if none of the points has already been used */
  for (i=0; i<match0; i++)
    {
      p1 = con0[i].p[0];	if (p1 > -1)	if (++tim[0][p1] > 1)	continue;
      p2 = con0[i].p[1];	if (p2 > -1)	if (++tim[1][p2] > 1)	continue;
      p3 = con0[i].p[2];	if (p3 > -1  && n_img > 2)
	if (++tim[2][p3] > 1)	continue;
      p4 = con0[i].p[3];	if (p4 > -1  && n_img > 3)
	if (++tim[3][p4] > 1)	continue;

      con[match++] = con0[i];
    }

  match2 = match-match4-match3;
  if(n_img==1){
     sprintf (buf, "determined %d points from 2D", match1);
     puts (buf);
  }
  else{
     sprintf (buf, "%d consistent quadruplets, %d triplets and %d unambigous pairs",
	      match4, match3, match2);
     puts (buf);
  }
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");

  /* ----------------------------------------------------------------------- */

  /* give each used pix the correspondence number */
  for (i=0; i<match; i++)
    {
      for (j=0; j<n_img; j++)
	{
	  p1 = geo[j][con[i].p[j]].pnr;
	  if (p1 > -1 && p1 < 1202590843)
	    {
	      pix[j][p1].tnr= i;
	    }
	}
    }

  /* draw crosses on canvas */
  if (display) {
    for (i=0; i<match4; i++)	       	/* red crosses for quadruplets */
      {
	for (j=0; j<n_img; j++)
	  {
	    p1 = geo[j][con[i].p[j]].pnr;
	    if (p1 > -1)
	      {
		if (   (fabs(pix[j][p1].x-zoom_x[j]) < imx/(2*zoom_f[j]))
		       && (fabs(pix[j][p1].y-zoom_y[j]) < imy/(2*zoom_f[j])))
		  {
		    intx = (int) ( imx/2 + zoom_f[j] * (pix[j][p1].x-zoom_x[j]));
		    inty = (int) ( imy/2 + zoom_f[j] * (pix[j][p1].y-zoom_y[j]));
		    drawcross (interp, intx, inty, cr_sz, j, "red");
		    if (zoom_f[j]>=2) draw_pnr (interp, intx+5 , inty+0, i, j, "white");
		  }
	      }
	  }
      }

    for (i=match4; i<match4+match3; i++)	/* green crosses for triplets */
      {
	for (j=0; j<n_img; j++)
	  {
	    p1 = geo[j][con[i].p[j]].pnr;
	    if (p1 > -1 && con[i].p[j] > -1)
	      {
		if (   (fabs(pix[j][p1].x-zoom_x[j]) < imx/(2*zoom_f[j]))
		       && (fabs(pix[j][p1].y-zoom_y[j]) < imy/(2*zoom_f[j])))
		  {
		    intx = (int) ( imx/2 + zoom_f[j] * (pix[j][p1].x-zoom_x[j]));
		    inty = (int) ( imy/2 + zoom_f[j] * (pix[j][p1].y-zoom_y[j]));
		    drawcross ( interp, intx, inty, cr_sz, j, "green" );
		    if (zoom_f[j]>=2) draw_pnr (interp, intx+5 , inty+0, i, j, "white");/* number of triplet */
		  }
		
	      }
	  }
      }
    for (i=match4+match3; i<match4+match3+match2; i++)
      {			      	/* yellow crosses for pairs */
	for (j=0; j<n_img; j++)
	  {
	    p1 = geo[j][con[i].p[j]].pnr;
	    if (p1 > -1 && con[i].p[j] > -1)
	      {
		if (   (fabs(pix[j][p1].x-zoom_x[j]) < imx/(2*zoom_f[j]))
		       && (fabs(pix[j][p1].y-zoom_y[j]) < imy/(2*zoom_f[j])))
		  {
		    intx = (int) ( imx/2 + zoom_f[j] * (pix[j][p1].x-zoom_x[j]));
		    inty = (int) ( imy/2 + zoom_f[j] * (pix[j][p1].y-zoom_y[j]));
		    drawcross (interp, intx, inty, cr_sz, j, "yellow");
		    if (zoom_f[j]>=2) draw_pnr (interp, intx+5 , inty+0, i, j, "white"); /* number of triplet */
		  }
	      }
	  }
      }
    
    for (j=0; j<n_img; j++)
      {
	
	for (i=0; i<num[j]; i++)
	  {			      	/* blue crosses for unused detections */
	    p1 = pix[j][i].tnr;
	    if (p1 == -1 )
	      {
		if (   (fabs(pix[j][i].x-zoom_x[j]) < imx/(2*zoom_f[j]))
		       && (fabs(pix[j][i].y-zoom_y[j]) < imy/(2*zoom_f[j])))
		  {
		    intx = (int) ( imx/2 + zoom_f[j] * (pix[j][i].x-zoom_x[j]));
		    inty = (int) ( imy/2 + zoom_f[j] * (pix[j][i].y-zoom_y[j]));
		    drawcross (interp, intx, inty, cr_sz, j, "blue");
		    
		  }
	      }
	  }
      }
  }
  /* ----------------------------------------------------------------------- */
  /* free memory for lists of correspondences */
  for (i1=0; i1<n_img-1; i1++)
    {
      for (i2=i1+1; i2<n_img; i2++) free (list[i1][i2]);
    }

  free (con0);

  sprintf (buf,"Correspondences done");
  Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 $tbuf");
}

