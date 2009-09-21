/****************************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	end of 99 ...going on

Description:	       	different search routines and track tools

Routines contained:     pix_in_next, candsearch_in_pix, searchposition,
                        predict, readseqtrackcrit, readtrackdata,
                        searchquader

****************************************************************************/
#include "ptv.h"

int pix_in_next (next, num, x, y, dl, dr, du, dd, found)
target  next[];
int     num;
double  x, y, dl, dr, du, dd;
int found[POSI];
{
  /* search of POSI near candidates in targetlist */

  int  j, j0, dj;
  int  zaehler=0;
  double  d, xmin, xmax, ymin, ymax;
  double dcand[POSI];
  int cand[POSI];
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;


  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].y > ymax)  break;	       	/* finish search */
      if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax && next[j].tnr>0)
	{
	  d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));
	  cand[zaehler]=next[j].tnr; dcand[zaehler]=d;
	  zaehler++;
	}
    }

  return (zaehler);
}


int candsearch_in_pix (next, num, x, y, dl, dr, du, dd, p)
target  next[];
int     num;
double  x, y, dl, dr, du, dd;
int p[4];
{
  /* search of four near candidates in targetlist */

  int  	  j, j0, dj, pnr = -999;
  int  zaehler=0, p1, p2, p3, p4;
  double  d, dmin=1e20, xmin, xmax, ymin, ymax;
  double d1, d2, d3, d4;
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;

  if(xmin<0.0) xmin=0.0;
  if(xmax>imx) xmax=imx;
  if(ymin<0.0) ymin=0.0;
  if(ymax>imy) ymax=imy;

  if(x<0.0) x=0.0;
  if(x>imx) x=imx;
  if(y<0.0) y=0.0;
  if(y>imy) y=imy;

  p1 = p2 = p3 = p4 = -999;
  d1 = d2 = d3 = d4 = dmin;

  if (x>=0.0 && x<=imx ) { if (y>=0.0 && y<=imy ) {

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].tnr != -1 ) {
	if (next[j].y > ymax )  break;	       	/* finish search */
	if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax )
	  {
	    d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));

	    if (d < dmin) { dmin = d; pnr = j;
	    }
	    if ( d < d1 )
	      {
		p4=p3; p3=p2; p2=p1; p1=j;
		d4=d3; d3=d2; d2=d1; d1=d;
	      }
	    else if ( d1 < d &&  d < d2 )
	      {
		p4=p3; p3=p2; p2=j;
		d4=d3; d3=d2; d2=d;
	      }
	    else if ( d2 < d && d < d3 )
	      {
		p4=p3; p3=j;
		d4=d3; d3=d;
	      }
	    else if ( d3 < d && d < d4 )
	      {
		p4=j;
		d4=d;
	      }
	  }
      }
    }

  p[0]=p1;
  p[1]=p2;
  p[2]=p3;
  p[3]=p4;
  for (j=0; j<4; j++) if ( p[j] != -999 ) {zaehler++; }
  } }
  return (zaehler);
}


int candsearch_in_pixrest (next, num, x, y, dl, dr, du, dd, p)
target  next[];
int     num;
double  x, y, dl, dr, du, dd;
int p[4];
{
  /* search of four near candidates in targetlist */

  int  	  j, j0, dj;
  int  zaehler=0, p1, p2, p3, p4;
  double  d, dmin=1e20, xmin, xmax, ymin, ymax;
  xmin = x - dl;  xmax = x + dr;  ymin = y - du;  ymax = y + dd;

  if(xmin<0.0) xmin=0.0;
  if(xmax>imx) xmax=imx;
  if(ymin<0.0) ymin=0.0;
  if(ymax>imy) ymax=imy;

  if(x<0.0) x=0.0;
  if(x>imx) x=imx;
  if(y<0.0) y=0.0;
  if(y>imy) y=imy;

  p1 = p2 = p3 = p4 = -999;

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (next[j0].y < ymin) j0 += dj;
      else j0 -= dj;
    }

  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */
  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (next[j].tnr == -1 ) {
	if (next[j].y > ymax )  break;	       	/* finish search */
	if (next[j].x > xmin  &&  next[j].x < xmax && next[j].y > ymin  &&  next[j].y < ymax )
	  {
	    d = sqrt ((x-next[j].x)*(x-next[j].x) + (y-next[j].y)*(y-next[j].y));
	    if (d < dmin) { dmin = d; p1 = j; }
	  }
      }
    }

  p[0]=p1;
  p[1]=p2;
  p[2]=p3;
  p[3]=p4;
  for (j=0; j<4; j++) if ( p[j] != -999 ) {zaehler++; }
  return (zaehler);
}



void predict (x1, y1, x2, y2, x3, y3)
double x1, y1, x2, y2, *x3, *y3;
{
  *x3 = 2*x2 - x1;
  *y3 = 2*y2 - y1;
  return;
}


void readseqtrackcrit ()
{
  int i_img;
  /* reads pixfiles and try to track particles in imagespace
     over the sequence */
  fpp = fopen_r ("parameters/sequence.par");
  for (i_img=0; i_img<4; i_img++) {
    fscanf (fpp, "%s\n", seq_name[i_img]);
  }
  /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);

  fpp = fopen_r ("parameters/track.par");
  fscanf (fpp, "%lf\n", &tpar.dvxmin);
  fscanf (fpp, "%lf\n", &tpar.dvxmax);
  fscanf (fpp, "%lf\n", &tpar.dvymin);
  fscanf (fpp, "%lf\n", &tpar.dvymax);
  fscanf (fpp, "%lf\n", &tpar.dvzmin);
  fscanf (fpp, "%lf\n", &tpar.dvzmax);
  fscanf (fpp, "%lf\n", &tpar.dangle);
  fscanf (fpp, "%lf\n", &tpar.dacc);
  fscanf (fpp,"%d\n", &tpar.add);
  /*
    fscanf (fpp,"%d\n", &tpar.dsumg);
    fscanf (fpp,"%d\n", &tpar.dn);
    fscanf (fpp,"%d\n", &tpar.dnx);
    fscanf (fpp,"%d\n", &tpar.dny);
  */
  fclose (fpp);

  /* read illuminated layer data */
  fpp = fopen_r ("parameters/criteria.par");
  fscanf (fpp, "%lf\n", &X_lay[0]);
  fscanf (fpp, "%lf\n", &Zmin_lay[0]);
  fscanf (fpp, "%lf\n", &Zmax_lay[0]);
  fscanf (fpp, "%lf\n", &X_lay[1]);
  fscanf (fpp, "%lf\n", &Zmin_lay[1]);
  fscanf (fpp, "%lf\n", &Zmax_lay[1]);
  fclose (fpp);
}



void searchquader(X, Y, Z, xr, xl, yd, yu)
double X, Y, Z, xr[4], xl[4], yd[4], yu[4];
{
  int k, i;
  double x, y, xz, yz;
  coord_3d quader[8], point;

  /* project quader in image space to define search area */
  for (k=0; k<8; k++)
    {
      quader[k].pnr=k;
    }
  /* calculation of quader points */
  point.pnr=0; point.x=X; point.y=Y; point.z=Z;

  quader[0].x=X+tpar.dvxmin; quader[0].y=Y+tpar.dvymin; quader[0].z=Z+tpar.dvzmin; /* --- */
  quader[1].x=X+tpar.dvxmax; quader[1].y=Y+tpar.dvymin; quader[1].z=Z+tpar.dvzmin; /* +-- */
  quader[2].x=X+tpar.dvxmin; quader[2].y=Y+tpar.dvymax; quader[2].z=Z+tpar.dvzmin; /* -+- */
  quader[3].x=X+tpar.dvxmax; quader[3].y=Y+tpar.dvymin; quader[3].z=Z+tpar.dvzmax; /* --+ */
  quader[4].x=X+tpar.dvxmax; quader[4].y=Y+tpar.dvymax; quader[4].z=Z+tpar.dvzmin; /* ++- */
  quader[5].x=X+tpar.dvxmax; quader[5].y=Y+tpar.dvymin; quader[5].z=Z+tpar.dvzmax; /* +-+ */
  quader[6].x=X+tpar.dvxmin; quader[6].y=Y+tpar.dvymax; quader[6].z=Z+tpar.dvzmax; /* -++ */
  quader[7].x=X+tpar.dvxmax; quader[7].y=Y+tpar.dvymax; quader[7].z=Z+tpar.dvzmax; /* +++ */

  /* calculation of search area */
  for (i=0; i<n_img; i++)
    {
      xr[i]=0;
      xl[i]=imx;
      yd[i]=0;
      yu[i]=imy;
      img_coord (point.x, point.y, point.z, Ex[i], I[i], G[i], ap[i], mmp, &xz,&yz);
      metric_to_pixel (xz,yz, imx,imy, pix_x,pix_y, &xz,&yz, chfield);

      for (k=0; k<8; k++)
	{
	  img_coord (quader[k].x, quader[k].y, quader[k].z, Ex[i], I[i], G[i], ap[i], mmp, &x,&y);
	  metric_to_pixel (x,y, imx,imy, pix_x,pix_y, &x,&y, chfield);

	  if (x <xl[i] ) xl[i]=x;
	  if (y <yu[i] ) yu[i]=y;
	  if (x >xr[i] ) xr[i]=x;
	  if (y >yd[i] ) yd[i]=y;
	}
      if (xl[i] < 0 ) xl[i]=0;
      if (yu[i] < 0 ) yu[i]=0;
      if (xr[i] > imx) xr[i]=imx;
      if (yd[i] > imy) yd[i]=imy;

      xr[i]=xr[i]-xz;
      xl[i]=xz-xl[i];
      yd[i]=yd[i]-yz;
      yu[i]=yz-yu[i];
    }
  return;

}



void sortwhatfound (foundpix item[16], int *zaehler)
{
  int i,j,m, different;
  foundpix temp;

  different=0;

  /* where what was found */
  for (i=0; i<16; i++)
    for (j=0; j<4; j++)
      for (m=0; m<4; m++)
	if(item[i].ftnr == item[4*j+m].ftnr)
	  {
	    item[i].whichcam[j]=1;
	  }

  /* how often was ftnr found */
  for (i=0; i<16; i++)
    for (j=0; j<n_img; j++)
      if (item[i].whichcam[j] == 1 && item[i].ftnr !=-1) item[i].freq++;

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }

  for (i=0; i<16; i++)
    for (j=i+1; j<16; j++)
      {
	if (item[i].ftnr == item[j].ftnr || item[j].freq <2)
	  {
	    item[j].freq=0;
	    item[j].ftnr=-1;
	  }
      }

  /* sort freq */
  for (i=1; i<16; ++i)  for (j=16-1; j>=i; --j)
    {
      if ( item[j-1].freq < item[j].freq )
	{
	  temp = *(item+j-1); *(item+j-1) = *(item+j); *(item+j) = temp;
	}
    }
  for (i=0; i<16; ++i) if(item[i].freq != 0) different++;
  *zaehler=different;

}

void angle_acc( X0, Y0, Z0,X1, Y1, Z1, X2, Y2, Z2, angle, acc )
double X0, Y0, Z0,X1, Y1, Z1, X2, Y2, Z2, *angle, *acc;
{
  double ds, da;
  coord_3d v0, v1, v2;

  /* calculation of angle and acceleration */
  v0.x=X1-X0; v0.y=Y1-Y0; v0.z=Z1-Z0;
  v1.x=X2-X0; v1.y=Y2-Y0; v1.z=Z2-Z0;
  v2.x=v1.x-v0.x; v2.y=v1.y-v0.y; v2.z=v1.z-v0.z;

  ds= sqrt( v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
  /* special case 200 gon */
  if (v1.x ==-v0.x && v1.y ==-v0.y && v1.z ==-v0.z)
    { *angle =200; } else {
  da=acos ((v0.x*v1.x+v0.y*v1.y+v0.z*v1.z)/
	   (sqrt( v1.x*v1.x+v1.y*v1.y+v1.z*v1.z)*
	    sqrt( v0.x*v0.x+v0.y*v0.y+v0.z*v0.z)));
  *angle=da*ro; }
  *acc=ds;

  return;
}
