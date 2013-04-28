#define sqr(x) x*x
#define maxcand 200 //Beat changed it on 090325
#define maxtrack 64 //Beat changed it on 090325

#define M 20000 //Beat changed it on 090325

#ifndef  M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define POSI 80 //Beat changed it on 090325

typedef	double	Dmatrix[3][3];	/* 3 x 3 rotation matrix */

typedef struct
{
  int pnr;
  double x, y, z;
    }
coord_3d;

typedef struct
{
  int pnr;
  double x, y;
}
coord_2d;

typedef struct
{
  double  x0, y0, z0;
  double  omega, phi, kappa;
  Dmatrix dm;
}
Exterior;

typedef struct
{
  double  dacc, dangle, dvxmax, dvxmin;
  double dvymax, dvymin, dvzmax, dvzmin;
  int  dsumg, dn, dnx, dny, add;
}
trackparameters;

typedef struct
{
  double xh, yh;
  double cc;
}
Interior;

typedef struct
{
  double vec_x,vec_y,vec_z;
}
Glass;

typedef struct
{
  double k1,k2,k3,p1,p2,scx,she;
}
ap_52;

typedef struct
{
  double n1,n2,n3,d;
  int    lut;
}
mm_3p;

typedef struct
{
  int  	  nlay;
  double  n1;
  double  n2[3];
  double  d[3];
  double  n3;
  int     lut;
}
mm_np;

typedef struct
{
  int     pnr;
  double  x, y;
  int     n, nx, ny, sumg;
  int     tnr;
}
target;

typedef struct
{
  int 	pos, status;
  short	xmin, xmax, ymin, ymax;
  int   n, sumg;
  double  x, y;
  int   unr, touch[4], n_touch;	/* unified with target unr, touching ... */
}
peak;

typedef struct
{
  short	       	x,y;
  unsigned char	g;
  short	       	tnr;
}
targpix;

typedef struct
{
  int  	left, right;
  double  tol, corr;
}
conjugate;

typedef struct WAIST {	// added, ad holten 04-2013
    double x;
    double y;
    double z;
    double dist;
} Waist;

typedef struct
{
  int  	pnr;
  double  tol, corr;
  Waist  waist;	/* distance perpendicular to epipolar line and coordiantes midpoint, ad holten 04-2013 */
}
candidate;

typedef struct
{
  int     p[4];
  double  corr;
}
n_tupel;

typedef struct
{
  int nr;
  int p[4];
}
corres;

typedef struct
{
  int    	p1;	       	/* point number of master point */
  int    	n;	       	/* # of candidates */
  int    	p2[maxcand];	/* point numbers of candidates */
  double	corr[maxcand];	/* feature based correlation coefficient */
  double	dist[maxcand];	/* distance perpendicular to epipolar line */
  Waist     waist[maxcand];	/* distance perpendicular to epipolar line and coordiantes midpoint, ad holten 04-2013 */
}
correspond;	       	/* correspondence candidates */

typedef struct
{
  coord_3d	origin;
  int          	nr,nz;
  double       	rw;
  double       	*data;
}
mm_LUT;


typedef struct /* struct for what was found to corres */
{
 int ftnr, freq, whichcam[4];
}
foundpix;

typedef struct
{
  int multi, h[maxtrack], freq[maxtrack];
  double quali[maxtrack];
}
currlist;

typedef struct
{
  int p,n;
  double x1, y1, z1;
  int type;
}
vector;

typedef struct
{
  int z1, c[maxcand], n[maxcand];
  double quali[maxcand];
}
prevlist;

typedef struct Pstruct
{
  float x[3]; /*coordinates*/
  int prev, next; /*pointer to prev or next link*/
  int prio; /*Prority of link is used for differen levels*/
  float decis[POSI]; /*Bin for decision critera of possible links to next dataset*/
  float finaldecis; /*final decision critera by which the link was established*/
  int linkdecis[POSI]; /* pointer of possible links to next data set*/
  int inlist; /* Counter of number of possible links to next data set*/
} P;

typedef struct
{
	int type;
	int x;
	int y;
	int x1;
	int y1;
	int size;
	int minzoom;
	double value;
	char color[32];
} Marker;

typedef struct
{
	double xc;			// normalized position of the view
	double yc;
	int fac;
	int vwx;
	int vwy;
	int fixed;
	int alpha;
} Zoompar;

typedef struct PTV_IS {    // point in the ptv_is files
    int prev;
    int next;
    double X;
    double Y;
    double Z;
} ptv_is;

typedef struct RT_IS {    // not existing in the ETH-code
    int    id;
    double X;
    double Y;
    double Z;
    int    tnr[4];
} rt_is;

typedef struct
{
  int pnr;
  double x, y, z;
  double r;
} sort_point;

typedef struct
{
	double xlay[2], ylay[2], zminlay[2], zmaxlay[2];
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double cnx, cny, cn, csumg, corrmin, eps0;
} critparameters;

/* ==== new typedefs needed for the polynomial mapping method === */

typedef struct {
	double pixx; 
	double pixy;
	double geox;
	double geoy;
} PNT;

typedef struct CORR {
    double X;
    double Y;
    double Z;
    int    pnr[4];
} CorrP;

typedef struct FITDATA {
    double x, y, z, u, v, w;
    int    pnr;
    double xp[4];
    double yp[4];
} FitData;

typedef struct FITGLOB {
    int order;
    int crossorder;
    int ndata;
    FitData  *data;
    coord_3d *pnt3;    
    coord_2d *pnt2;    
    coord_2d *in2D;
    coord_2d *out2D;
} FITGLOB;

typedef struct POLYFIT {
	int     order;
	int     crossorder;
	int     ncameras;
	int     nxparms1;
	int     nxparms2;
	int     nyparms1;
	int     nyparms2;
	double* Xparms1;
	double* Xparms2;
	double* Yparms1;
	double* Yparms2;
} POLYFIT;

typedef struct {
	double left;
	double right;
	double top;
	double bottom;
} Rect;

typedef struct Line3D {
	double xoff;
	double yoff;
	double zoff;
	double xrc;
	double yrc;
	// int    pnr;
} Line3D;

// --- types used in the shaking algorithm ---
typedef struct FITINFO {
	int     order;
	int     crossorder;
	int     nparms1;
	int     nparms2;
	double* parms1;
	double* parms2;
	double  result;
}	FitInfo;


POLYFIT* fitpixepi;
POLYFIT* fit3dpix;

#define MAX_FILENAME_LEN 1024
#define FILENAME_IN "res/rt_is"
#define FILENAME_OUT "res/ptv_is"
#define PARAMETERFILE "res/track.par"
#define COORDPARAFILE "res/coord.par"
#define STATUSFILE "res/track.out"
