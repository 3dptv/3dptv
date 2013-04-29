/********************************************************************

  Author/Copyright: 	Jochen Willneff
  
  Address:				Institute of Geodesy and Photogrammetry
						ETH - Hoenggerberg
						CH - 8093 Zurich
  
  Creation Date:		October'97	   
  
  Description:			editing of parameter files
  
  Routines contained:	parameter_panel_init, done_proc_c
  
*********************************************************************/

/*
Copyright (c) 1990-2011 ETH Zurich

See the file license.txt for copying permission.
*/

/* -------------------------------------------------------------------------
	Added code for using polynomials for position mapping.
	Ad Holten, 04-2013
-------------------------------------------------------------------------- */
#include "ptv.h"

int init_polycalib_panelparms(Tcl_Interp* interp);
int init_polycalib_panelparms(Tcl_Interp* interp);
int init_polymanori_panelparms(Tcl_Interp* interp);



int init_ethz_panelparms(Tcl_Interp* interp)
{
    char varname[20], val[100][256];
	int  i, ic, ip, rv;
    
	// The function read_strings() will print a message if not all the parameters are found.
	
	/* read 21 parameters from ptv.par */
	rv = read_strings("parameters/ptv.par", 21, val);

	Tcl_SetVar2(interp, "mp", "ncam",     val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg1",    val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal1",  val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg2",    val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal2",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg3",    val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal3",  val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg4",    val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal4",  val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "highpass", val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "allCam",   val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "tiff",     val[11], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "imx",      val[12], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "imy",      val[13], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pix_x",    val[14], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pix_y",    val[15], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "type",     val[16], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "air",      val[17], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "glass",    val[18], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "water",    val[19], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "thicknessglass", val[20], TCL_GLOBAL_ONLY);

	/* read 12 parameters from criteria.par */
	if (!read_strings("parameters/criteria.par", 12, val))
		return FALSE;

	Tcl_SetVar2(interp, "mp", "xmin",     val[0], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xminzmin", val[1], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xminzmax", val[2], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xmax",     val[3], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xmaxzmin", val[4], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xmaxzmax", val[5], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "nx",       val[6], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "ny",       val[7], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "npix",     val[8], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "sgv",      val[9], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "mincorr",  val[10], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "tolepi",   val[11], TCL_GLOBAL_ONLY);

	/* read 12 parameters from cal_ori.par */
	if (!read_strings("parameters/cal_ori.par", 12, val))
		return FALSE;

    Tcl_SetVar2(interp, "cp", "platecoord", val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "calp1",      val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ori1",       val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "calp2",      val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ori2",       val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "calp3",      val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ori3",       val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "calp4",      val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ori4",       val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "tiff",       val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pair",       val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "type",       val[11], TCL_GLOBAL_ONLY);

	/* read 16 parameters from man_ori.par */
	if (!read_strings("parameters/man_ori.par", 16, val))
		return FALSE;
	
	for (i=0; i<16; i++) {
		ic = i/4 + 1;
		ip = i%4 + 1;
		sprintf(varname, "p%d%d", ic, ip);		// ie p11, p12, ..
		Tcl_SetVar2(interp, "cp", varname, val[i], TCL_GLOBAL_ONLY);
	}

    /* read 12 parameters from orient.par */
	if (!read_strings("parameters/orient.par", 12, val))
		return FALSE;

    Tcl_SetVar2(interp, "cp", "pnrori", val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pdist",  val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "xp",     val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "yp",     val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "k1",     val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "k2",     val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "k3",     val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "p1",     val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "p2",     val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "scx",    val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "she",    val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "interf", val[11], TCL_GLOBAL_ONLY);

    /* changed parameters are read out all together not before "done_proc_c" */
    return TRUE;
}


int init_poly_panelparms(Tcl_Interp* interp)
{
    char val[100][256];
    int  rv;

	/* read 16 parameters from poly_ptv.par */
	rv = read_strings("parameters/poly_ptv.par", 16, val);
    Tcl_SetVar2(interp, "mp", "ncam",     val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "nimg",     val[0], TCL_GLOBAL_ONLY);        // at leats as many images as camera's, adh
    Tcl_SetVar2(interp, "mp", "fimg1",    val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal1",  val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg2",    val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal2",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg3",    val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal3",  val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "fimg4",    val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "camcal4",  val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "highpass", val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "allCam",   val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "tiff",     val[11], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "imx",      val[12], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "imy",      val[13], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pix_x",    val[14], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "pix_y",    val[15], TCL_GLOBAL_ONLY);


	/* read 12 parameters from poly_criteria.par */
	rv = read_strings("parameters/poly_criteria.par", 12, val);
	Tcl_SetVar2(interp, "mp", "xmin",     val[0], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "ymin",     val[1], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "zmin",     val[2], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "xmax",     val[3], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "ymax",     val[4], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "zmax",     val[5], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "nx",       val[6], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "ny",       val[7], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "npix",     val[8], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "sgv",      val[9], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "mincorr",  val[10], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "mp", "tolepi",   val[11], TCL_GLOBAL_ONLY);

	/* --- read calibration paramters from body_calori.par or mult_calori.par --- */
	rv = init_polycalib_panelparms(interp);

	/* read the point indices for manual orientations */
	rv = init_polymanori_panelparms(interp);

	return TRUE;
}

int init_polycalib_panelparms(Tcl_Interp* interp)
{
    char *fname, varname[32], val[100][256], *pval;
    int  i, io, ic, ne, ncam, nplanes, plane;

	/* --- read all the calibration paramters --- */
	fname = usingZplanes ? "parameters/mult_calori.par"
						 : "parameters/body_calori.par";
	ne = read_allstrings(fname, val, 100);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);

	// count per plane: 1+ncam names, i.e.: file with coordinates plus the image files
	if (nplanes<1 || ne < 6 + nplanes*(1+ncam)) {
		printf("Not enough parameters in %s\n", fname);
		if (nplanes == 0) nplanes = 1;	// define at least Tcl-globals for one plane.
	}

	if (ncam < n_img) {
		ncam = n_img;		// define at least for n_img cameras, Tcl-globals
		sprintf(val[0],"%d", n_img);
	}

	Tcl_SetVar2(interp, "cp", "ncam",    val[0], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "cp", "nplanes", val[1], TCL_GLOBAL_ONLY);
	i = 2;
	if (!usingZplanes) {
		for (ic=1; ic<=4; ic++) {
			pval = (ic<=ncam) ? val[i++] : "";
			sprintf(varname, "calp%d", ic);
			Tcl_SetVar2(interp, "cp", varname, pval, TCL_GLOBAL_ONLY);
		}
		Tcl_SetVar2(interp, "cp", "platecoord", val[i++], TCL_GLOBAL_ONLY);
	}
	else {
		for (plane=1; plane<=nplanes; plane++) {
			sprintf(varname, "platecoord%d", plane);
			Tcl_SetVar2(interp, "cp", varname, val[i++], TCL_GLOBAL_ONLY);
			for (ic=1; ic<=4; ic++) {
				pval = (ic<=ncam) ? val[i++] : "-";
				sprintf(varname, "cal%d%d", plane, ic);
				Tcl_SetVar2(interp, "cp", varname, pval, TCL_GLOBAL_ONLY);
			}
		}
		Tcl_SetVar2(interp, "cp", "plane_id", "1", TCL_GLOBAL_ONLY);	// start with the first plane
	}
	for (io=1; io<=4; io++) {				// degree of the polynomals to use
		sprintf(varname, "order%d", io);
		Tcl_SetVar2(interp, "cp", varname, val[i++], TCL_GLOBAL_ONLY);
	}
	return TRUE;
}

int init_polymanori_panelparms(Tcl_Interp* interp)
{
    char *fname, varname[32], val[100][256], *pval;
    int  i, ic, ip, ne, ncam, nori, nplanes, plane;

	/* read the point indices for manual orientations */
	fname = !usingZplanes ? "parameters/body_manori.par" 
						  : "parameters/mult_manori.par";
	ne   = read_allstrings(fname, val, 100);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);
	nori    = atoi(val[2]);

	if (nplanes<1 || ne < 3 + nplanes*ncam*nori) {
		printf("Not enough parameters in %s\n", fname);
		if (nori<3 || nplanes<1) {
			nplanes = 1;	// define Tcl-globals for one plane at least
			nori    = 4;
			sprintf(val[1],"%d", nplanes);
			sprintf(val[2],"%d", nplanes);
		}
	}
	if (ncam < n_img) {
		ncam = n_img;		// define Tcl-globals for n_img cameras
		sprintf(val[0],"%d", n_img);
	}
	Tcl_SetVar2(interp, "cp", "nplanes", val[1], TCL_GLOBAL_ONLY);
	Tcl_SetVar2(interp, "cp", "nori",    val[2], TCL_GLOBAL_ONLY);

	i = 3;
	for (plane=1; plane<=nplanes; plane++)
		for (ic=1; ic<=4 ; ic++)
			for (ip=1; ip<=nori; ip++) {
				if (!usingZplanes)
					sprintf(varname, "p%d%d", ic, ip);
				else
					sprintf(varname, "p%d%d%d", plane, ic, ip);
				pval = (ic<=ncam) ? val[i++] : "";
				Tcl_SetVar2(interp, "cp", varname, pval, TCL_GLOBAL_ONLY);
			}
    return TRUE;
}

int init_common_panelparms(Tcl_Interp* interp)
{
	/* set the rest of the Tcl panel variables */
	int ne;
    char val[50][256];

	/* read 15 parameters from targ_rec.par */
	if (!read_strings("parameters/targ_rec.par", 15, val))
		return FALSE;

    Tcl_SetVar2(interp, "mp", "partgv1",   val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "partgv2",   val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "partgv3",   val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "partgv4",   val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "partdisc",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pminnpix",  val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pmaxnpix",  val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pminnpixx", val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pmaxnpixx", val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pminnpixy", val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pmaxnpixy", val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "psumgv",    val[11], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "pcrossize", val[12], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "mask",      val[13], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "maskname",  val[14], TCL_GLOBAL_ONLY);


	if (!read_strings("parameters/pft_version.par", 1, val))
		return FALSE;
	Tcl_SetVar2(interp, "mp", "target", val[0], TCL_GLOBAL_ONLY);


	if (!read_strings("parameters/sequence.par", 6, val))
		return FALSE;

    Tcl_SetVar2(interp, "mp", "basename1", val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "basename2", val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "basename3", val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "basename4", val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "seqfirst",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "mp", "seqlast",   val[5], TCL_GLOBAL_ONLY);

    /* read 13 parameters from detect_plate.par */
	if (!read_strings("parameters/detect_plate.par", 13, val))
		return FALSE;

    Tcl_SetVar2(interp, "cp", "partgv1",  val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "partgv2",  val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "partgv3",  val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "partgv4",  val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "toldisc",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "minnpix",  val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "maxnpix",  val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "minnpixx", val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "maxnpixx", val[8], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "minnpixy", val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "maxnpixy", val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ppsumgv",  val[11], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "ppcrossize", val[12], TCL_GLOBAL_ONLY);

    /* read 4 parameters from shaking.par */
	if (!read_strings("parameters/shaking.par", 4, val))
		return FALSE;

    Tcl_SetVar2(interp, "cp", "first_shake",     val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "last_shake",      val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "maxPoints_shake", val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "maxFrames_shake", val[3], TCL_GLOBAL_ONLY);

    /* read 3 parameters from examine.par */
	ne = read_allstrings("parameters/examine.par", val, 3);
	if (ne>=2) {
		Tcl_SetVar2(interp, "cp", "examineFlag", val[0], TCL_GLOBAL_ONLY);
		Tcl_SetVar2(interp, "cp", "combineFlag", val[1], TCL_GLOBAL_ONLY);
		if (ne>2 )Tcl_SetVar2(interp, "mp", "examine3", val[2], TCL_GLOBAL_ONLY);
	}
	else {
		printf("parameters/examine.par is missing parameters, must be >=2\n");
		return FALSE;
	}

    /* read 22 parameters from track.par */
	if (!read_strings("parameters/track.par", 22, val))
		return FALSE;
    
    Tcl_SetVar2(interp, "tp", "dvxmin",  val[0], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dvxmax",  val[1], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dvymin",  val[2], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dvymax",  val[3], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dvzmin",  val[4], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dvzmax",  val[5], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dangle",  val[6], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "dacc",    val[7], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "add",     val[8], TCL_GLOBAL_ONLY);
    /* 5 textboxes */
    Tcl_SetVar2(interp, "tp", "maxnum",  val[9], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "deltat",  val[10], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "vmax",    val[11], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "linktol", val[12], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "jumptol", val[13], TCL_GLOBAL_ONLY);
    /* 8 checkboxes */
    Tcl_SetVar2(interp, "tp", "gluing",  val[14], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "track",   val[15], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "iterate", val[16], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "remkin",  val[17], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "wbin",    val[18], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "wtraj",   val[19], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "output",  val[20], TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "tp", "camconfig", val[21], TCL_GLOBAL_ONLY);

    return TRUE;
}


/********************************************************************/
/* parameter_panel_init : initializes the parameter panel */
/********************************************************************/
int parameter_panel_init(Tcl_Interp* interp)
{
	int rv; 

    // set the Tcl variables: cp(method) and cp(multi)
    Tcl_SetVar2(interp, "cp", "method", map_method == ETHZ ? "ETHZ" : "POLY", TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "cp", "multi",  usingZplanes? "1" : "0", TCL_GLOBAL_ONLY);

	/* initialize Tcl variables for the dialog boxes */
	if (map_method == ETHZ)
		rv  = init_ethz_panelparms(interp);
	else 
		rv  = init_poly_panelparms(interp);

	rv &= init_common_panelparms(interp);
	return rv;	// returns 1 on succes, 0 on error.

	/* changed parameters are read out all together not before "done_proc_c" */
}



/* -- Saving data to the parameter files ----------------------------------------------------- */

void save_method(Tcl_Interp* interp, int method, int usingplanes)
{
  	FILE* fp = fopen("parameters/method.par", "w");

    if (fp) {
	    fprintf(fp, "%s\n", method == ETHZ ? "ETHZ" : "POLY");
        fprintf(fp, "%s\n", usingplanes == 1 ? "multi" : "body");
        fclose(fp);
    }
}

BOOL save_ptv_par(Tcl_Interp* interp, char *fname)
{
	fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "ncam",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "fimg1",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "camcal1",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "fimg2",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "camcal2",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "fimg3",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "camcal3",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "fimg4",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "camcal4",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "highpass", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "allCam",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "tiff",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "imx",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "imy",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pix_x",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pix_y",    TCL_GLOBAL_ONLY));
	if (map_method == ETHZ) {
		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "type",     TCL_GLOBAL_ONLY));
		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "air",      TCL_GLOBAL_ONLY));
		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "glass",    TCL_GLOBAL_ONLY));
		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "water",    TCL_GLOBAL_ONLY));
		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "thicknessglass", TCL_GLOBAL_ONLY));
	}
    fclose (fp1); 
	return TRUE;
}

//BOOL save_poly_ptv_par(Tcl_Interp* interp, char *fname)
//{
//	const char *valp;
//	char varname[32];
//	int i, ncam;
//	fp1 = fopen_wp(fname);
//	if (!fp1) return FALSE;
//
//
//	valp = Tcl_GetVar2(interp, "mp", "ncam", TCL_GLOBAL_ONLY);
//    fprintf (fp1, "%s\n", valp);
//	ncam = atoi(valp);
//	for (i=1; i<=ncam; i++) {
//		sprintf(varname, "fimg%d", i);
//		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", varname, TCL_GLOBAL_ONLY));
//		sprintf(varname, "camcal%d", i);
//		fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", varname, TCL_GLOBAL_ONLY));
//	}
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "highpass", TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "allCam",   TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "tiff",     TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "imx",      TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "imy",      TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pix_x",    TCL_GLOBAL_ONLY));
//    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pix_y",    TCL_GLOBAL_ONLY));
//    fclose (fp1); 
//	return TRUE;
//}

BOOL save_targrec_par(Tcl_Interp* interp, char *fname)
{
	fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "partgv1",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "partgv2",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "partgv3",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "partgv4",   TCL_GLOBAL_ONLY));
	fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "partdisc",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pminnpix",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pmaxnpix",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pminnpixx", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pmaxnpixx", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pminnpixy", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pmaxnpixy", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "psumgv",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "pcrossize", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "mask",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "maskname",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "rel_disc",  TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

BOOL save_criteria_par(Tcl_Interp* interp, char* fname)
{
	fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    if (map_method == ETHZ) {
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmin",     TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xminzmin", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xminzmax", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmax",     TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmaxzmin", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmaxzmax", TCL_GLOBAL_ONLY));
    }
    else {
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmin", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "ymin", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "zmin", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "xmax", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "ymax", TCL_GLOBAL_ONLY));
        fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "zmax", TCL_GLOBAL_ONLY));
    }
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "nx",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "ny",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "npix",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "sgv",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "mincorr",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "tolepi",   TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

BOOL save_calori_par(Tcl_Interp* interp, char* fname)
{
	fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "platecoord", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "calp1",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ori1",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "calp2",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ori2",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "calp3",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ori3",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "calp4",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ori4",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "tiff",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pair",       TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "type",       TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

BOOL save_sequence_par(Tcl_Interp* interp, char* fname)
{
    if (! (fp1 = fopen_wp(fname)) ) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "basename1", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "basename2", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "basename3", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "basename4", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "seqfirst",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "seqlast",   TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

BOOL save_detectplate_par(Tcl_Interp* interp, char *fname)
{
    if (! (fp1 = fopen_wp(fname)) ) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "partgv1",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "partgv2",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "partgv3",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "partgv4",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "toldisc",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "minnpix",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "maxnpix",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "minnpixx", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "maxnpixx", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "minnpixy", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "maxnpixy", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ppsumgv",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "ppcrossize", TCL_GLOBAL_ONLY));

    fclose(fp1);
	return TRUE;
}

BOOL save_pft_version_par(Tcl_Interp* interp, char *fname)
{
    fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "mp", "target",  TCL_GLOBAL_ONLY));
    fclose(fp1);
	return TRUE;
}

BOOL save_manori_par(Tcl_Interp* interp, char *fname)
{
	char pname[10];
	int  ic, ip;

	fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

	for (ic=1; ic<=4; ic++) {
		for (ip=1; ip<=4; ip++) {
			sprintf(pname, "p%d%d", ic, ip);	// ie. p11, p12, ...
			fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", pname, TCL_GLOBAL_ONLY));
		}
	}
    fclose(fp1);
	return TRUE;
}


BOOL save_poly_calori_par(Tcl_Interp* interp, char *fname)
{
	char varname[32];
	const char *valp;
	int  i, ic, ncam, nplanes;

	ncam    = atoi(Tcl_GetVar2(interp, "cp", "ncam", TCL_GLOBAL_ONLY));
	nplanes = atoi(Tcl_GetVar2(interp, "cp", "nplanes", TCL_GLOBAL_ONLY));

	fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

	fprintf(fp1, "%d\n", ncam);
	fprintf(fp1, "%d\n", nplanes);
	for (ic=1; ic<=ncam; ic++) {
		sprintf(varname, "calp%d", ic);
		fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
	}
	fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", "platecoord", TCL_GLOBAL_ONLY));

	for (i=1; i<=4; i++) {
		sprintf(varname, "order%d", i);
		valp = Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY);
		fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
	}
    fclose (fp1);
	return TRUE;
}

BOOL save_poly_manori_par(Tcl_Interp* interp, char *fname)
{
	char varname[32];
	int  ic, ip, ncam, nplanes, nori;

	ncam    = atoi(Tcl_GetVar2(interp, "cp", "ncam", TCL_GLOBAL_ONLY));
	nplanes = atoi(Tcl_GetVar2(interp, "cp", "nplanes", TCL_GLOBAL_ONLY));
	nori    = atoi(Tcl_GetVar2(interp, "cp", "nori", TCL_GLOBAL_ONLY));

	fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

	fprintf(fp1, "%d\n", ncam);
	fprintf(fp1, "%d\n", nplanes);
	fprintf(fp1, "%d\n", nori);
	for (ic=1; ic<=ncam; ic++) {
		for (ip=1; ip<=nori; ip++) {
			sprintf(varname, "p%d%d", ic,ip);		// ie. p11, p12, ...
			fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
		}	
	}
    fclose (fp1);
	return TRUE;
}

BOOL save_mult_calori_par(Tcl_Interp* interp, char *fname, int iplane)
{
	char varname[32];
	int  i, plane, ic, ncam, nplanes;

	ncam    = atoi(Tcl_GetVar2(interp, "cp", "ncam", TCL_GLOBAL_ONLY));
	nplanes = atoi(Tcl_GetVar2(interp, "cp", "nplanes", TCL_GLOBAL_ONLY));

	fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

	fprintf(fp1, "%d\n", ncam);
	fprintf(fp1, "%d\n", nplanes);
	for (plane=1; plane<=nplanes; plane++) {
		sprintf(varname, "platecoord%d", plane);
		fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
		for (ic=1; ic<=ncam; ic++) {
			sprintf(varname, "cal%d%d", plane, ic);
			fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
		}
	}
	for (i=1; i<=4; i++) {
		sprintf(varname, "order%d", i);
		fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
	}

    fclose (fp1);
	return TRUE;
}

BOOL save_mult_manori_par(Tcl_Interp* interp, char *fname, int iplane)
{
	char varname[32];
	int  plane, ic, ip, ncam, nplanes, nori;

	ncam    = atoi(Tcl_GetVar2(interp, "cp", "ncam", TCL_GLOBAL_ONLY));
	nori    = atoi(Tcl_GetVar2(interp, "cp", "nori", TCL_GLOBAL_ONLY));
	nplanes = atoi(Tcl_GetVar2(interp, "cp", "nplanes", TCL_GLOBAL_ONLY));

	fp1 = fopen_wp(fname);
    if (!fp1) return FALSE;

	fprintf(fp1, "%d\n", ncam);
	fprintf(fp1, "%d\n", nplanes);
	fprintf(fp1, "%d\n", nori);
	for (plane=1; plane<=nplanes; plane++) {
		for (ic=1; ic<=ncam; ic++)
			for (ip=1; ip<=4; ip++)	{
				sprintf(varname, "p%d%d%d", plane,ic,ip);		// ie. p111, p112, ...
				fprintf(fp1, "%s\n", Tcl_GetVar2(interp, "cp", varname, TCL_GLOBAL_ONLY));
			}	
	}
    fclose (fp1);
	return TRUE;
}

int save_shaking_par(Tcl_Interp* interp, char *fname)
{
    fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "first_shake",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "last_shake",      TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "maxPoints_shake", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "maxFrames_shake", TCL_GLOBAL_ONLY));
    
    fclose (fp1);
	return TRUE;
}

int save_orient_par(Tcl_Interp* interp, char *fname)
{
    fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pnrori", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "pdist",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "xp",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "yp",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "k1",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "k2",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "k3",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "p1",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "p2",     TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "scx",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "she",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "interf", TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

int save_examine_par(Tcl_Interp* interp, char* fname)
{
    fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "examineFlag", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "combineFlag", TCL_GLOBAL_ONLY));
    // extra examine parameter for defining debugging output: examine3
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "cp", "examine3", TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}

int save_track_par(Tcl_Interp* interp, char *fname)
{
    fp1 = fopen_wp(fname);
	if (!fp1) return FALSE;

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvxmin", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvxmax", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvymin", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvymax", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvzmin", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dvzmax", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dangle", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "dacc",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "add",    TCL_GLOBAL_ONLY));

    /*
    valp = Tcl_GetVar2(interp, "tp", "dsumg", TCL_GLOBAL_ONLY);
    fprintf (fp1, "%s\n", valp);
    valp = Tcl_GetVar2(interp, "tp", "dn", TCL_GLOBAL_ONLY);
    fprintf (fp1, "%s\n", valp);
    valp = Tcl_GetVar2(interp, "tp", "dnx", TCL_GLOBAL_ONLY);
    fprintf (fp1, "%s\n", valp);
    valp = Tcl_GetVar2(interp, "tp", "dny", TCL_GLOBAL_ONLY);
    fprintf (fp1, "%s\n", valp);
    */

    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "maxnum",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "deltat",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "vmax",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "linktol", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "jumptol", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "gluing",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "track",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "iterate", TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "remkin",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "wbin",    TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "wtraj",   TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "output",  TCL_GLOBAL_ONLY));
    fprintf (fp1, "%s\n", Tcl_GetVar2(interp, "tp", "camconfig", TCL_GLOBAL_ONLY));

    fclose (fp1);
	return TRUE;
}


/****************************************************************************/
/* done procedure:
   reads the item contents
   and prints them (changed or unchanged) back to the files */
/****************************************************************************/

int done_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
{
    const char *valp;
	int rv, plane, method, usingplanes;

    // update the flags: map_method and usingZplanes
	valp = Tcl_GetVar2(interp, "cp", "method", TCL_GLOBAL_ONLY);
    map_method = (valp && valp[0] == 'P') ? POLYN : ETHZ;

	valp = Tcl_GetVar2(interp, "cp", "multi", TCL_GLOBAL_ONLY);
    usingZplanes = valp && valp[0] == '1';
    save_method(interp, map_method, usingZplanes);


    /* update all the parameter files with the curent values */
	if (map_method == ETHZ) {
		rv  = save_ptv_par      (interp, "parameters/ptv.par");
        rv &= save_criteria_par (interp, "parameters/criteria.par");
        rv &= save_calori_par   (interp, "parameters/cal_ori.par");
		rv &= save_manori_par   (interp, "parameters/man_ori.par");
        rv &= save_orient_par   (interp, "parameters/orient.par");

	}
	else {		// map_method = POLYN
		rv  = save_ptv_par (interp, "parameters/poly_ptv.par");
        rv &= save_criteria_par (interp, "parameters/poly_criteria.par");
		if (!usingZplanes) {
			rv &= save_poly_calori_par(interp, "parameters/cal_oribody");
			rv &= save_poly_manori_par(interp, "parameters/body_manori.par");
		}
		else {
			// valp = TclGetString(interp, "cp", "nlevels");
			// nlevels = atoi(valp);
			plane = atoi(Tcl_GetVar2(interp, "cp", "plane_id", TCL_GLOBAL_ONLY));
			rv &= save_mult_calori_par(interp, "parameters/mult_calori.par", plane);
			rv &= save_mult_manori_par(interp, "parameters/mult_manori.par.par", plane);
		}
	}
    
	rv  = save_targrec_par		(interp, "parameters/targ_rec.par");
	rv &= save_pft_version_par	(interp, "parameters/pft_version.par");
	rv &= save_sequence_par		(interp, "parameters/sequence.par");
	rv &= save_detectplate_par  (interp, "parameters/detect_plate.par");
    rv &= save_shaking_par      (interp, "parameters/shaking.par");
    rv &= save_examine_par      (interp, "parameters/examine.par");
	rv &= save_track_par        (interp, "parameters/track.par");
        
	get_processingmethod(&method, &usingplanes);
	if (method != map_method || usingplanes != usingZplanes)
		return TCL_ERROR;
    return TCL_OK;
}

