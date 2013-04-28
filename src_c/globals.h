/*  global declarations for ptv  */

#define nmax 20240
#define MaxOri 6

enum MAPMETHODS {NONE, ETHZ, POLYN};	// mapping method in use
// zooming
enum MARKTYPES  {OBJ_MARK, OBJ_CROSS, OBJ_VECTOR, OBJ_PNR, OBJ_VALUE, OBJ_OVAL};

extern  Marker  *markers[4];
extern  int     nmarkers[4];
extern  int     nori;
extern  int     map_method;
extern  int     usingZplanes;
extern  int     examine3;

extern  int     n_img;                          /* no of images */
extern  int     hp_flag;                        /* flag for highpass */
extern  int     allCam_flag;                    /* flag for using all cams for points */
extern  int     tiff_flag;                      /* flag for tiff header */
extern  int     chfield;                        /* flag for field mode */
extern  int     nfix;                           /* no. of control points */
extern  int     num[];                          /* no. of particles per image */
extern  int     nump[4];                        /* no. of particles in previous image */
extern  int     numc[4];                        /* no. of particles in current image */
extern  int     numn[4];                        /* no. of particles in next image */
extern  int     n_trac[];                       /* no. of tracks */
extern  int     match;                          /* no. of matches */
extern  int     match2;                         /* no. of matches in 2nd pass */
extern  int     corp, corc, corn;               /* no. of correspondences in p,c,n */

extern  int     nr[4][4];                       /* point numbers for man. ori */
extern  int     imx, imy, imgsize;              /* image size */
// extern  int     zoom_x[], zoom_y[], zoom_f[];   /* zoom parameters */
extern  int     pp1, pp2, pp3, pp4,pp5;         /* for man. orientation */
extern  int     seq_first, seq_last;            /* 1. and last img of seq */
extern  int     demo_nr;                        /* for demo purposes */
extern  int     examine;                        /* extra output */
extern  int     dump_for_rdb;                   /* # of dumpfiles for rdb */
extern  int     cr_sz;                          /* size of crosses */
extern  int     display;                        /* display flag */
extern  int     m[4];
extern  int     trackallocflag;  /* checkflag if mega, c4, t4 already allocated */

extern  double  pix_x, pix_y;                   /* pixel size */
extern  double  pi, ro;                         /* pi, ro */
extern  double  cn, cnx, cny, csumg, eps0, corrmin; /* correspond. par */
extern  double  rmsX, rmsY, rmsZ, mean_sigma0;      /* a priori rms */
extern  double  X_lay[], Y_lay[], Zmin_lay[], Zmax_lay[];    /* illu. layer current slice */
extern  double  db_scale;                       /*dumbbell length, Beat Mai 2010*/ 

extern  FILE    *fp1, *fp2, *fp3, *fp4, *fpp;   /* file pointers */

extern  char  img_name[4][256];                 /* original image names */
extern  char  img_lp_name[4][256];              /* lowpass image names */
extern  char  img_hp_name[4][256];              /* highpass image names */
extern  char  img_cal[4][128];                  /* calibration image names */
extern  char  img_ori[4][128];                  /* image orientation data */
extern  char  img_ori0[4][128];                 /* orientation approx. values */
extern  char  img_addpar0[4][128];              /* approx. image additional parameters */
extern  char  img_addpar[4][128];               /* image additional parameters */
extern  char  fixp_name[128];                   /* testfield fixpoints */
extern  char  seq_name[4][128];                 /* sequence names */
extern  char  track_dir[128];                   /* directory with dap track data */
extern  char  res_name[128];                    /* result destination */
extern  char  buf[], val[];                     /* buffer */
extern  char  filename[128];

extern  unsigned char *img[];                   /* image data */
extern  unsigned char *zoomimg;                 /* zomm image data */

extern  Exterior  Ex[];                         /* exterior orientation */
extern  Interior  I[];                          /* interior orientation */
extern  Glass     G[];                          /* glass orientation */
extern  ap_52     ap[];                         /* add. parameters */
extern  mm_np     mmp;                          /* 3-media parameters */
extern  target    pix[4][nmax];                 /* target pixel data */
extern  target    pix0[4][12];                  /* pixel data for man_ori points */
extern  coord_2d  crd[4][nmax];                 /* (distorted) metric coordinates */
extern  coord_2d  geo[4][nmax];                 /* corrected metric coordinates */
extern  coord_3d  fix[];                        /* testfield points coordinates */
extern  n_tupel   con[];                        /* list of correspondences */
extern  mm_LUT    mmLUT[4];                     /* LUT for mm radial displacement */
extern  coord_3d  *p_c3d;
extern  target    *p[4];
extern  target    *c[4];
extern  target    *n[4];
extern  target    *t4[4][4];
extern  int       nt4[4][4];
extern  corres    *corrp;
extern  corres    *corrc;
extern  corres    *corrn;
extern  corres    *c4[4];
extern  trackparameters tpar;

extern  P         *mega[4];

extern  FILE      *fopen_r ();
extern  double    multimed_r ();

// --- change_parameter.c ---
int  parameter_panel_init(Tcl_Interp* interp);
int  done_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void save_method(Tcl_Interp* interp, int method, int usingplanes);

// --- checkpoints.c ---
void checkpoint_proc (Tcl_Interp* interp);

// --- checkpoints_poly.c ---
void checkpoint_poly (Tcl_Interp* interp);


// --- correspondences.c ---
void correspondences_4 (Tcl_Interp* interp, const char** argv);

// --- draw.c ---
int  drawcross (Tcl_Interp* interp, int x0, int y0, int size, int imgnr, char color[]);
int  drawvector (Tcl_Interp* interp, int x0, int y0, int x1, int y1, int width, int imgnr, char color[]);
int  draw_pnr (Tcl_Interp* interp, int x, int y, int pnr, int imgnr, char color[]);
int  draw_value (Tcl_Interp* interp, int x, int y, double pnr, int imgnr, char color[]);
void mark_detections (Tcl_Interp* interp, int nr);
void mark_correspondences (Tcl_Interp* interp, int nr);
void mark_corr (Tcl_Interp* interp, int nr);
int  mark_track_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
int  trajectories_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void add_mark(int x0, int y0, int imgnr);
void clear_drawnobjectslist();
int  clearmarkers_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void clicked_to_imgcoordinates(int *px, int *py, double x, double y, int n);
void draw_objects(Tcl_Interp* interp, int imgnr);
int  drawmarker(Tcl_Interp* interp, Marker mark, int imgnr);
int  drawoval (Tcl_Interp* interp, int x0, int y0, int size, int imgnr, char color[]);
int  draw_pnr_autohide (Tcl_Interp* interp, int x, int y, int pnr, int imgnr, char color[], int minzoomlevel, int zoomlevel);
int  draw_value_autohide (Tcl_Interp* interp, int x, int y, double pnr,	int imgnr, char color[], int minzoomlevel, int zoomlevel);
BOOL get_tclzoomparms(Tcl_Interp* interp, Zoompar *zoompar, int i_img);
void img_to_view_coordinates(int *px, int *py, double x, double y, int n);
BOOL isinview(double x, double y, int i_img);

// --- epi.c ---
int  epi_mm (double x1, double y1, Exterior Ex1, Interior I1, Glass G1, Exterior Ex2, Interior I2, Glass G2, mm_np mmp,
			double* xmin, double* ymin, double* xmax, double* ymax);
int  epi_mm_2D (double x1, double y1, Exterior Ex1, Interior I1, Glass G1, mm_np mmp, double *xp, double *yp, double *zp);
void find_candidate_plus (coord_2d crd[], target pix[], int num, double xa, double ya, double xb, double yb, 
			double eps, int n, int nx, int ny, int sumg, candidate cand[], int *count, int nr, const char** argv);
void find_candidate_plus_msg (coord_2d crd[], target pix[], int num, double xa, double ya, double xb, double yb, 
			double eps, int n, int nx, int ny, int sumg, candidate cand[], int *count, int i12);

// --- image_processing.c ---
void filter_3(unsigned char *img, unsigned char *img_lp);
void lowpass_3(unsigned char *img, unsigned char *img_lp);
void unsharp_mask(int n, unsigned char *img0, unsigned char *img_lp);
void zoom (unsigned char *img, unsigned char *zoomimg, int xm, int ym, int zf);
void split (unsigned char *img, int field);
void copy_images(unsigned char *img1, unsigned char *img2);
void subtract_mask(unsigned char *img, unsigned char *img_mask, unsigned char *img_new) ;

// --- imgcoord.c ---
void img_coord (double X, double Y, double Z, Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, double *x, double *y);
void img_xy_mm_geo (double X, double Y, double Z, Exterior Ex, Interior I, Glass G, mm_np mm, double *x, double *y);

// --- intersect.c ---
void intersect_rt (double X1, double Y1, double Z1, double a1, double b1, double c1, double X2, double Y2, double Z2,
				   double a2, double b2, double c2, double *X, double *Y, double *Z);

// --- jw_main.c ---
int  main(int argc, char **argv);
int  Tcl_AppInit(Tcl_Interp *interp);
int  jw_Init(Tcl_Interp *interp);
int  change_method_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);


// --- lsqadj.c ---
void  ata (void *a, void *ata, int m, int n);
void  ata_v2 (void *a, void *ata, int m, int n, int n_large);
void  atl (double *u, void *a, double *l, int m, int n);
void  atl_v2 (double *u, void *a, double *l, int m, int n, int n_large);
void  matinv (void *a, int n);
void  matinv_v2 (void *a, int n, int n_large);
void  matmul (double *a, void *b, double *c, int m, int n, int k);
void  matmul_v2 (double *a, void *b, double *c, int m, int n, int k, int m_large, int n_large);
void  transp (double a[], int m, int n);
void  mat_transpose (void *mat1, void *mat2, int m, int n);

// --- multimed.c ---
void   multimed_nlay (Exterior ex, mm_np mm, double X, double Y, double Z, double *Xq, double *Yq);
void   multimed_nlay_v2 (Exterior ex, Exterior ex_o, mm_np mm, double X, double Y, double Z, double  *Xq, double *Yq);
double multimed_r_nlay (Exterior ex, mm_np mm, double X, double Y, double Z);
void   trans_Cam_Point(Exterior ex, mm_np mm, Glass gl, double X, double Y, double Z, Exterior *ex_t, 
					   double *X_t, double *Y_t, double *Z_t, double *cross_p, double *cross_c);
void   trans_Cam_Point_back(Exterior ex, mm_np mm, Glass gl, double X, double Y, double Z, Exterior *ex_t, 
					   double *X_t, double *Y_t, double *Z_t, double *cross_p, double *cross_c);
void   back_trans_Point(double X_t, double Y_t, double Z_t, mm_np mm, Glass G,
					   double cross_p[], double cross_c[], double *X, double *Y, double *Z);
double multimed_r_nlay_v2 (Exterior ex, Exterior ex_o, mm_np mm, double X, double Y, double Z);
void   init_mmLUT (int i_cam);
double get_mmf_from_mmLUT (int i_cam, double X, double Y, double Z);
void   volumedimension (double *xmax, double *xmin, double *ymax, double *ymin, double *zmax, double *zmin);

// --- orientation.c ---
void prepare_eval (int n_img, int *n_fix);
void mid_point(double A1x, double A1y, double A1z, double Ux, double Uy, double Uz, double B1x, double B1y, 
			   double B1z, double Vx, double Vy, double Vz, double *dist, double *XX, double *YY, double *ZZ);
void eval_ori_v2 (double db_scale, double weight_scale, int n_img, int nfix, 
				  double *d_outer, double *av_dist_error, double *residual);
void orient_v5 (int n_img, int nfix, Exterior *Ex, Interior *I, Glass *G, ap_52 *ap);
void orient_v3 (Tcl_Interp* interp, Exterior Ex0, Interior I0, Glass G0, ap_52 ap0, mm_np mm, int nfix, 
				coord_3d fix[], coord_2d crd[], Exterior *Ex, Interior *I, Glass *G, ap_52 *ap, int nr);
void orient_v6 (Tcl_Interp* interp, Exterior Ex0, Interior I0, Glass G0, ap_52 ap0, mm_np mm, int nfix, coord_3d fix[],
				int good[], coord_2d crd[], Exterior *Ex, Interior *I, Glass *G, ap_52 *ap, int nr);
void raw_orient_v3 (Exterior Ex0, Interior I, Glass G0, ap_52 ap, mm_np mm, int nfix,
					coord_3d fix[], coord_2d crd[], Exterior *Ex, Glass *G, int only_show);

// --- parse.c ---
int read_strings(char* fname, int nparms, char val[][256]);
int read_allstrings(char* fname, char val[][256], int maxcnt);
BOOL parse_ptv_par(char* fname);			// x
BOOL parse_polyptv_par(char* fname);		// x
BOOL parse_criteria_par(char* fname);		// x
BOOL parse_polycriteria_par(char* fname);	// x
BOOL parse_polycriteria_par2(char* fname, critparameters *cpar);
BOOL parse_sequence_par(char* fname, char seq_name[][128], int *pseq_first, int *pseq_last);	// x
BOOL parse_shaking_par(char* fname, int *seq_first, int *seq_last, int *max_points, int *max_frames);
BOOL parse_rt_is(char* pathname, rt_is **plist, int *pcount);
BOOL parse_targetfile(char* pathname, target **tlist, int *pcount);
BOOL createfolder(char* pathname);

BOOL parse_body_calori_par(char* fname, char imgname[][256], char* fixpname, int deg3d2d[2], int deg2depi[2]);
BOOL parse_mult_calori_par(char* fname, int iplane, int *pncam, int *pnplanes, char imgname[][256], char* fixpname, int deg3d2d[2], int deg2depi[2]);

BOOL parse_body_manori_par(char* fname, int list[][MaxOri], int* pncam, int* pnori);

BOOL parse_multi_manori_par(char* fname, int iplane,  int list[][MaxOri], int *pncam, int *pnori);
BOOL read_2Dcoordinates(char *pixp_name, coord_2d **plist, int *pcount);
BOOL read_3Dcoordinates(char *fixp_name, coord_3d **plist, int *pcount);
BOOL read_target_crds(char *pathname, coord_3d **pfixp, int *pcount);
double get_sortradius(char *pathname, double defaultR);
void get_processingmethod(int* pMapping, int* pUsingplanes);
BOOL parse_track_par(char* fname, trackparameters *tpar);

// --- peakfitting.c ---
int  peak_fit_new (Tcl_Interp* interp, unsigned char *img, char par_file[], 
				  int xmin, int xmax, int ymin, int ymax, target pix[], int nr);

// --- pointpos.c ---
void  det_lsq_old (Exterior Ex[4], Interior I[4], ap_52 ap[4], mm_np mm, double x1, double y1, double x2, double y2, 
				  double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp);
void  dist_to_ray(double x, double y, Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, 
				 double Xp, double Yp, double Zp, double *dist);

void  pos_from_ray(Exterior Ex[4], Interior I[4], Glass G[4], ap_52 ap[4], mm_np mm, double x1, double y1, double x2, 
				 double y2, double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp, double *dist);
void  det_lsq_3d(Exterior Ex[4], Interior I[4], Glass G[4], ap_52 ap[4], mm_np mm, double x1, double y1, double x2,
				double y2, double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp);
void  det_lsq (Exterior Ex[4], Interior I[4], Glass G[4], ap_52 ap[4], mm_np mm, double x1, double y1, double x2, double y2, 
				double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp);
void  det_lsq_3 (Exterior Ex[3], Interior I[3], Glass G[3], ap_52 ap[3], mm_np mm, double x1, double y1,
				double x2, double y2, double x3, double y3, double *Xp, double *Yp, double *Zp);
void  det_lsq_4 (Exterior Ex[4], Interior I[4], Glass G[4], ap_52 ap[4], mm_np mm, double x1, double y1, double x2,
				double y2, double x3, double y3, double x4, double y4, double *Xp, double *Yp, double *Zp);
void  det_lsq_2 (Exterior Ex[2], Interior I[2], Glass G[2], ap_52 ap[2], mm_np mm, double x1, double y1, 
				double x2, double y2, double *Xp, double *Yp, double *Zp);

// --- poly_fit.c ---
void    boundingrect(Rect* prect, coord_2d* points, int npoints, int flag);
int     calibration_poly_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void    CreateZplaneFitdata(coord_2d *pixcrd, coord_2d *planecrd, FitInfo *pfit, double z, double x1, double x2, int nx, double y1, double y2, int ny);
void    destroy_POLYFIT_list(POLYFIT *fitlist);
void    det_lsq_3d_poly(double xp[], double yp[], int np, double *pX, double *pY, double *pZ);
void    draw_fitted_fixpoints2d(Tcl_Interp* interp, int iimg, FitInfo fit, sort_point* fix, int nfix);
int     epi_mm_poly(double xp1, double yp1, int ic1, int ic2, double* pxmin, double* pymin, double* pxmax, double* pymax);
BOOL    Fit2D_to_2D(FitInfo* pfi, int degree[2], coord_2d* Pin, coord_2d* Pout, int ndata);
void    Fit3D_to_2D(FitInfo* pfi, int degree[2], coord_3d *pnt3, coord_2d *pnt2, int ndata);
void    Fit3D_to_2D_body(FitInfo* pfi, coord_3d *pnt3, coord_2d *pnt2, int ndata);
void    find_candidate_poly_epiline(coord_2d crd[], target pix[], int num, double xa, double ya, double xb, double yb, 
                         double eps, target pt1, candidate cand[], int* pcount, const char** argv);
double  fitdata_Map2D_x(double X, double Y, FitInfo* pfit);
double  fitdata_Map2D_y(double X, double Y, FitInfo* pfit);
double  fitdata_Map3D_x(double X, double Y, double Z, FitInfo* pfit);
double  fitdata_Map3D_y(double X, double Y, double Z, FitInfo* pfit);
void    freefitparms(FitInfo* pfi);
void    GetFit_Pixels_Lines(FitInfo* pxFitinfo, FitInfo* pyFitinfo, FitInfo* fit, double* Zlev, int nz, Rect* prect);
void    GetFitPixels_Planes(FitInfo* fit, double* Zlev, int nz, Rect* prect, FitInfo* pfit3d2d);
void    getIOdata2D(double* iodata, double x, double y, int order, int crossorder);
void    getIOdata3D(double* iodata, double x, double y, double z, int order, int crossorder);
void    Line3D_mindist(Line3D line1, Line3D line2, Waist *pwaist);
double  Map2D_X(double x, double y, FitInfo* pfit);
double  Map2D_Y(double x, double y, FitInfo* pfit);
double  Map3D_x(double X, double Y, double Z, POLYFIT* pfit);
double  Map3D_y(double X, double Y, double Z, POLYFIT* pfit);
double  Mapper2D(double x, double y, int order, int crossorder, double* parms, int nparms);
void    parsepolycalibfile(char* fname, POLYFIT *pfit3dpix, POLYFIT *pfitpixepi);
void    saveparms_3D_to_pix(char* pathname, char* mode, FitInfo *pfit);
int     volumedimension_poly(double* xmax, double* xmin, double* ymax, double* ymin, double* zmax, double* zmin);


// --- poly_ptv.c ---
void ray_tracing_poly(double x, double y, int icam, Line3D* pline);


// --- poly_correspondences.c ---
void correspondences_poly(Tcl_Interp* interp, const char** argv);

// --- poly_shaking1.c ---
int poly_shaking_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);

// --- fit2D.c ---
int getNumParms2D(int order, int crossorder);
int getNumParms3D(int order, int crossorder);
void InvFitFunction3D(double index, double *iodata, int iosize);
void saveparms_2D_to_rays(char* pathname, char* mode, FitInfo *pXfit, FitInfo *pYfit);
double dotmult(double* v1, double *v2, int nelem);
void FitFunction2D(double index, double *iodata, int iosize);


// --- ptv.c ---
void  read_ascii_data(int filenumber);
void  read_targets(int i_img, int filenumber,  int *num);
void  write_ascii_data(int filenumber);
void  write_added(int filenumber);
void  write_addedback(int filenumber);
void  read_ascii_datanew(int filenumber);
void  write_ascii_datanew(int filenumber);
int   tracking(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void  level2(void);
void  level3(void);
void  sort(int n, float a[], int b[]);
void  rotate_dataset(void);
void  neighbours(float seekx[], float radi[], int nliste[], int *innliste, int set);

// --- nruril.c ---
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[], int ma, double **covar, 
	double *chisq, void (__cdecl *funcs)(double index, double *iodata, int iosize) );

// --- ray_tracing.c ---
void  ray_tracing (double x, double y, Exterior Ex, Interior I, mm_np mm,
				double *Xb2, double *Yb2, double *Zb2, double *a3, double *b3, double *c3);
void  point_line_line(Exterior Ex0, Interior I0, Glass G0, mm_np mm, double gX0, double gY0, double gZ0, 
				double a0, double b0, double c0, Exterior Ex1, Interior I1, Glass G1, double gX1, double gY1, 
				double gZ1, double a1, double b1, double c1, double *x, double *y, double *z);
void  norm_cross(double a[3], double b[3], double *n1, double *n2, double *n3);
void  dot(double a[3], double b[3], double *d);
void  modu(double a[3], double *m);
void  ray_tracing_v2 (double x, double y, Exterior Ex, Interior I, Glass G, mm_np mm, double *Xb2,
				double *Yb2, double* Zb2, double *a3, double *b3, double *c3);

// --- rotation.c ---
void rotation_matrix (Exterior	Ex, Dmatrix   dm);

// --- segmentation.c ---
void highpass (char* pic_name, unsigned char *img, unsigned char *img_hp, int dim_lp, int filter_hp, int field, int nr);
void simple_connectivity (Tcl_Interp* interp, unsigned char *img0, unsigned char *img, 
				char par_file[], int xmin, int xmax, int ymin, int ymax, target pix[], int nr, int *num);
void targ_rec (Tcl_Interp* interp, unsigned char *img0, unsigned char *img, char par_file[],
				int xmin, int xmax, int ymin, int ymax, target pix[], int nr, int *num);

// --- sortgrid.c ---
void just_plot (Tcl_Interp* interp, Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, int imx, int imy, 
				double pix_x, double pix_y, int nfix, coord_3d fix[], int field, int n_img);
void sortgrid_man (Tcl_Interp* interp, Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, int imx, int imy, 
				double pix_x, double pix_y, int nfix, coord_3d fix[], int num, target pix[], int field, int n_img);
void sortgrid_file (Tcl_Interp* interp, Exterior Ex, Interior I, Glass G, ap_52 ap, mm_np mm, int imx, int imy, 
				double pix_x, double pix_y, int nfix, coord_3d fix[], int num, target pix[], int field, int n_img);


// --- tools.c ---
void  write_ori (Exterior Ex, Interior I, Glass G, char filename[64]);
int   read_ori (Exterior *Ex, Interior *I, Glass *G, char filename[64]);
FILE  *fopen_r (CHAR filename[256]);
FILE  *fopen_rp (char *filename);
FILE *fopen_wp (char *filename);
int   read_image (Tcl_Interp* interp, char path[128], unsigned char *img);
int   write_tiff (const char path[256], unsigned char *data, int nx, int ny);
void  compose_name_plus_nr (char basename[256], char str[256], int nr, char filename[256]);
void  compose_name_plus_nr_str (char basename[256], char str[256], int nr, char filename[256]);
int   kill_in_list (Tcl_Interp* interp, int nr, int num, int ms_x, int ms_y);
int   nearest_neighbour_geo (coord_2d crd[], int num, double x, double y, double eps);
int   nearest_neighbour_pix (target pix[], int num, double x, double y, double eps);
void  quicksort_coord2d_x (coord_2d *crd, int num);
void  qs_coord2d_x (coord_2d *crd, int left, int right);
void  quicksort_target_y (target *pix, int num);
void  quicksort_con (n_tupel *con, int num);
void  tclimg2cimg (Tcl_Interp* interp, unsigned char *c_img, Tk_PhotoImageBlock *tcl_img);
void  cimg2tclimg (Tcl_Interp* interp, unsigned char *c_img, Tk_PhotoImageBlock *tcl_img, int alpha);
double distance(double x1, double y1, double x2, double y2);


// --- track.c ---
int  trackcorr_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
int  trackback_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv);
void pixelcoord_from_3Dpnt(double *px, double *py, int ic, double X, double Y, double Z);

// --- trafo.c ---
void  pixel_to_metric(double xp, double yp, int imx, int imy, 
					 double pix_x, double pix_y, double *xc, double *yc, int field);
void  metric_to_pixel(double xc, double yc, int imx, int imy, 
					 double pix_x, double pix_y, double *xp, double *yp, int field);
void  distort_brown_affin (double x, double y, ap_52 ap, double *x1, double *y1);
void  correct_brown_affin (double x, double y, ap_52 ap, double *x1, double *y1);


// --- ttools.c ---
int  pix_in_next (target next[], int num, double x, double y, 
				 double dl, double dr, double du, double dd, int found[POSI]);
int  candsearch_in_pix (target next[], int num, double x, double y, 
				 double dl, double dr, double du, double dd, int p[4]);
int  candsearch_in_pixrest (target  next[], int num, double x, double y, 
						   double dl, double dr, double du, double dd, int p[4]);
void  predict (double x1, double y1, double x2, double y2, double *x3, double *y3);
int   readseqtrackcrit();
void  searchquader(double X, double Y, double Z, 
				  double xr[4], double xl[4], double yd[4], double yu[4]);
void  sortwhatfound (foundpix item[16], int *zaehler);
void  angle_acc(double X0, double Y0, double Z0, double X1, double Y1, double Z1, 
			   double X2, double Y2, double Z2, double *angle, double *acc);
