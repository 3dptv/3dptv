/*  global declarations for ptv  */

#define nmax 20000

extern	int    	n_img;			       /* no of images */
extern  int		hp_flag;           	       /* flag for highpass */
extern	int    	tiff_flag;	               /* flag for tiff header */
extern	int    	chfield;	               /* flag for field mode */
extern	int    	nfix;	       		       /* no. of control points */
extern	int    	num[]; 			       /* no. of particles per image */
extern	int    	nump[4]; 	       	       /* no. of particles in previous image */
extern	int    	numc[4]; 	       	       /* no. of particles in current image */
extern	int    	numn[4]; 	       	       /* no. of particles in next image */
extern	int    	n_trac[];      		       /* no. of tracks */
extern	int    	match; 			       /* no. of matches */
extern	int    	match2;		    	       /* no. of matches in 2nd pass */
extern  int     corp, corc, corn;              /* no. of correspondences in p,c,n */

extern	int    	nr[4][4];		       /* point numbers for man. ori */
extern	int    	imx, imy, imgsize;     	       /* image size */
extern	int    	zoom_x[], zoom_y[], zoom_f[];  /* zoom parameters */
extern	int    	pp1, pp2, pp3, pp4,pp5;	       /* for man. orientation */
extern	int    	seq_first, seq_last;	       /* 1. and last img of seq */
extern	int    	demo_nr;		       /* for demo purposes */
extern	int    	examine;       		       /* extra output */
extern	int    	dump_for_rdb;		       /* # of dumpfiles for rdb */
extern  int     cr_sz;                         /* size of crosses */
extern  int     display;                       /* display flag */
extern  int     m[4];
extern  int     trackallocflag;  /* checkflag if mega, c4, t4 already allocated */

extern	double 	pix_x, pix_y;		     	/* pixel size */
extern	double 	pi, ro;				/* pi, ro */
extern	double 	cn, cnx, cny, csumg, eps0, corrmin;	 /* correspond. par */
extern	double 	rmsX, rmsY, rmsZ, mean_sigma0;		 /* a priori rms */
extern	double  X_lay[],   Zmin_lay[],   Zmax_lay[];         /* illu. layer current slice */

extern	FILE   	*fp1, *fp2, *fp3, *fp4, *fpp;	/* file pointers */

extern	char  img_name[4][256];	        /* original image names */
extern	char  img_lp_name[4][256];    	/* lowpass image names */
extern	char  img_hp_name[4][256];      /* highpass image names */
extern	char  img_cal[4][128];		/* calibration image names */
extern	char  img_ori[4][128];		/* image orientation data */
extern	char  img_ori0[4][128];		/* orientation approx. values */
extern	char  img_addpar0[4][128];	/* approx. image additional parameters */
extern	char  img_addpar[4][128];	/* image additional parameters */
extern	char  fixp_name[128];	/* testfield fixpoints */
extern	char  seq_name[4][128];      	/* sequence names */
extern	char  track_dir[128];	       	/* directory with dap track data */
extern	char  res_name[128];	       	/* result destination */
extern	char  buf[], val[];		       	/* buffer */
extern	char  filename[128];

extern	unsigned char	*img[];			/* image data */
extern	unsigned char	*zoomimg;		/* zomm image data */

extern	Exterior        Ex[];	       	/* exterior orientation */
extern	Interior       	I[];	        /* interior orientation */
extern	Glass       	G[];	        /* glass orientation */
extern	ap_52	       	ap[];	       	/* add. parameters */
extern	mm_np	       	mmp;	       	/* 3-media parameters */
extern	target	       	pix[4][nmax];  	/* target pixel data */
extern	target	       	pix0[4][4];    	/* pixel data for man_ori points */
extern	coord_2d       	crd[4][nmax];  	/* (distorted) metric coordinates */
extern	coord_2d       	geo[4][nmax];  	/* corrected metric coordinates */
extern	coord_3d       	fix[];	        /* testfield points coordinates */
extern	n_tupel	       	con[];	      	/* list of correspondences */
extern	mm_LUT	       	mmLUT[4];     	/* LUT for mm radial displacement */
extern  coord_3d        *p_c3d;
extern  target          *p[4];
extern  target          *c[4];
extern  target          *n[4];
extern  target          *t4[4][4];
extern  int             nt4[4][4];
extern  corres          *corrp;
extern  corres	       	*corrc;
extern  corres	       	*corrn;
extern  corres	       	*c4[4];
extern  trackparameters tpar;

extern P *mega[4];

extern	FILE	*fopen_r ();
extern	double	multimed_r ();


int  jw_Init();
int  init_proc_c();
int  start_proc_c();
int  pre_processing_c ();
int  detection_proc_c();
int  correspondences_proc_c ();
int  determination_proc_c ();
int  sequence_proc_c();
void checkpoint_proc ();
int  quit_proc_c ();
int  parameter_panel_init();
int  done_proc_c();

void read_ori();
void write_ori ();
void compose_name_plus_nr ();
int  mouse_proc_c ();
void read_image();
int  read_tiff ();
int  write_tiff ();
void copy_images();

int  peak_fit_new();
void highpass();
void zoom ();
void lowpass_n ();
void lowpass_3 ();
void filter_3 ();
void unsharp_mask ();
void split ();
void multimed_nlay ();
void multimed_nlay_v2 ();
void trans_Cam_Point ();
void back_trans_Point ();
void getG();

void quicksort_con ();
int  flow_demo_c();
void simple_connectivity();
void targ_rec();
void quicksort_target_y();
void transp();
void quicksort_coord2d_x();
void rotation_matrix();
void img_xy ();
void img_xy_mm ();
void img_xy_mm_geo ();
void pixel_to_metric();
void metric_to_pixel ();
void correct_brown_affin();
void distort_brown_affin ();
int  pix_in_next();
int  kill_in_list ();
void compose_name_plus_nr ();
void compose_name_plus_nr_str ();
int  nearest_neighbour_geo ();
int  nearest_neighbour_pix ();
void bubble_y ();
void init_mmLUT ();


void correspondences_4();
void det_lsq();
void ray_tracing();
void ray_tracing_v2();
void intersect_rt_3m();
void intersect_rt();
int  nearest_neighbour_pix();
void img_coord();
void raw_orient();
void raw_orient_v3();
int mod();
void getabcFromRot();
void sortgrid_man();
void just_plot();
void det_lsq ();
void det_lsq_2 ();
void det_lsq_3 ();
void det_lsq_4 ();
void orient();
void orient_v3();
void volumedimension();

double epi_line ();
int  epi_mm ();
int  epi_mm_2D ();
void find_candidate_plus_msg ();
void find_candidate_plus ();

int  drawcross();
int  drawvector ();
int  draw_pnr ();
int  draw_value ();
void mark_detections ();
void mark_correspondences ();
void mark_corr ();

void ata ();
void ata_v2 ();
void matinv ();
void atl ();
void matmul ();
void mat_transpose ();

void tracking_proc ();
int  find_track ();

int  candsearch_in_pix ();
int  candsearch_in_pixrest ();
int  corrtest();
void predict();
void readseqtrackcrit();
void searchquader();
void angle_acc();
void sortwhatfound();
void bubble_foundpix1 ();
void bubble_foundpix2 ();
void tclimg2cimg();

void Pix_crd ();
void Back_trafo ();
void crd_pix();
int  vrmltracks_c();
int  trajectories_c();
int  vrmldetections_c();
int  vrmldettracks_c();


int  tracking();
void read_parameter();
void level1();
void level2();
void level3();
void read_ascii_data();
void rotate_dataset();
void write_ascii_data();
void readcoord_parameter();
void change_coordinates();
void sort();
void neighbours();
int  seq_track_proc_c();

void read_ascii_datanew();
void write_ascii_datanew();


