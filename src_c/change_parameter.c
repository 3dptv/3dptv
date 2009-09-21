/********************************************************************

  Author/Copyright:     Jochen Willneff
  
  Address:	      	Institute of Geodesy and Photogrammetry
                        ETH - Hoenggerberg
                        CH - 8093 Zurich
  
  Creation Date:        October'97     
  
  Description:          editing of parameter files
  
  Routines contained:   parameter_panel_init, done_proc_c
  
*********************************************************************/
#include "ptv.h"

/********************************************************************/
/* parameter_panel_init : initializes the parameter panel */
/********************************************************************/

int parameter_panel_init(Tcl_Interp* interp)
{
   char val[256];

	/* read 20 parameters from ptvppar */

	fp1 = fopen_r ("parameters/ptv.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "ncam", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "fimg1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "camcal1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "fimg2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "camcal2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "fimg3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "camcal3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "fimg4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "camcal4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "highpass", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "tiff", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "imx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "imy", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "pix_x", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "pix_y", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "type", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "air", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "glass", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "water", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "thicknessglass", val, TCL_GLOBAL_ONLY);
	fclose (fp1);

	/* read 13 parameters from targ_rec.par */

	fp1 = fopen_r ("parameters/targ_rec.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "partgv1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "partgv2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "partgv3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "partgv4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "partdisc", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pminnpix", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pmaxnpix", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pminnpixx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pmaxnpixx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pminnpixy", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pmaxnpixy", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "psumgv", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "pcrossize", val, TCL_GLOBAL_ONLY);
	fclose (fp1);

	/* read 12 parameters from criteria.par */

	fp1 = fopen_r ("parameters/criteria.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xmin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xminzmin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xminzmax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xmax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xmaxzmin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "xmaxzmax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "nx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "ny", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "npix", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "sgv", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "mincorr", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "tolepi", val, TCL_GLOBAL_ONLY);
	fclose (fp1);


	/* read 6 parameters from sequence.par */

	fp1 = fopen_r ("parameters/sequence.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "basename1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "basename2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "basename3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "basename4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "seqfirst", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "mp", "seqlast", val, TCL_GLOBAL_ONLY);
	
	fclose (fp1);

	/* read 11 parameters from cal_ori.par */

	fp1 = fopen_r ("parameters/cal_ori.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "platecoord", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "calp1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ori1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "calp2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ori2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "calp3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ori3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "calp4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ori4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "tiff", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "type", val, TCL_GLOBAL_ONLY);

	fclose (fp1);

	/* read 13 parameters from detect_plate.par */

	fp1 = fopen_r ("parameters/detect_plate.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "partgv1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "partgv2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "partgv3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "partgv4", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "toldisc", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "minnpix", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "maxnpix", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "minnpixx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "maxnpixx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "minnpixy", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "maxnpixy", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ppsumgv", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "ppcrossize", val, TCL_GLOBAL_ONLY);

	fclose (fp1);


	/* read 16 parameters from man_ori.par */

	fp1 = fopen_r ("parameters/man_ori.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p11", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p12", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p13", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p14", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p21", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p22", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p23", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p24", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p31", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p32", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p33", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p34", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p41", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p42", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p43", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p44", val, TCL_GLOBAL_ONLY);

	fclose (fp1);

	fp1 = fopen_r ("parameters/shaking.par");
    
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "first_shake", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "last_shake", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "maxPoints_shake", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "maxFrames_shake", val, TCL_GLOBAL_ONLY);
	
	fclose (fp1);

	/* read 11 parameters from orient.par */

	fp1 = fopen_r ("parameters/orient.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "pnrori", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "pdist", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "xp", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "yp", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "k1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "k2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "k3", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p1", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "p2", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "scx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "she", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "interf", val, TCL_GLOBAL_ONLY);

	fclose (fp1);


   
	fp1 = fopen_r ("parameters/examine.par");
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "examineFlag", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "cp", "combineFlag", val, TCL_GLOBAL_ONLY);
	fclose (fp1);


	/* read 5 parameters from track.par */

	fp1 = fopen_r ("parameters/track.par");

	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvxmin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvxmax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvymin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvymax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvzmin", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dvzmax", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dangle", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dacc", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "add", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	/*
	Tcl_SetVar2(interp, "tp", "dsumg", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dn", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dnx", val, TCL_GLOBAL_ONLY);
	fscanf (fp1, "%s", val);
	Tcl_SetVar2(interp, "tp", "dny", val, TCL_GLOBAL_ONLY);	
	fscanf (fp1, "%s", val);
	*/
	
	fclose (fp1);

	/* changed parameters are read out all together not before "done_proc_c" */

        return TCL_OK;
}


/****************************************************************************/
/* done procedure:
   reads the item contents
   and prints them (changed or unchanged) back to the files */
/****************************************************************************/


int done_proc_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv) 
{
	const char *valp;


	/* rewrite all parameter files */

	fp1 = fopen ("parameters/ptv.par", "w");
	
	valp = Tcl_GetVar2(interp, "mp", "ncam",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "fimg1",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "camcal1",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "fimg2",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "camcal2",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "fimg3",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "camcal3",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "fimg4",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "camcal4",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "highpass",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "tiff",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "imx",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "imy",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "pix_x",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "pix_y",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "type",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "air",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "glass",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "water",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "thicknessglass",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);


	fclose (fp1); 


	fp1 = fopen ("parameters/targ_rec.par", "w");

	valp = Tcl_GetVar2(interp, "mp", "partgv1",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "partgv2",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "partgv3",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "partgv4",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "partdisc", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pminnpix",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pmaxnpix",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pminnpixx",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pmaxnpixx",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pminnpixy",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pmaxnpixy",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "psumgv",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "pcrossize",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);


	fp1 = fopen ("parameters/criteria.par", "w");


	valp = Tcl_GetVar2(interp, "mp", "xmin",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "xminzmin",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "xminzmax",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "xmax",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "xmaxzmin",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "xmaxzmax",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "nx",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "ny",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "npix",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "sgv",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "mincorr",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "tolepi",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);


	fp1 = fopen ("parameters/cal_ori.par", "w");

	valp = Tcl_GetVar2(interp, "cp", "platecoord", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "calp1", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ori1", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "calp2", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ori2", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "calp3", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ori3", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "calp4", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ori4", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "tiff", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "type", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);


	fp1 = fopen ("parameters/sequence.par", "w");

	valp = Tcl_GetVar2(interp, "mp", "basename1", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "basename2", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "basename3", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "basename4", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "seqfirst", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "mp", "seqlast", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);


	fp1 = fopen ("parameters/detect_plate.par", "w");

	valp = Tcl_GetVar2(interp, "cp", "partgv1",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "partgv2",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "partgv3",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "partgv4",  TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "toldisc", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "minnpix", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "maxnpix", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "minnpixx", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "maxnpixx", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "minnpixy", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "maxnpixy", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ppsumgv", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "ppcrossize", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);

       
	fp1 = fopen ("parameters/man_ori.par", "w");

	valp = Tcl_GetVar2(interp, "cp", "p11", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p12", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p13", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p14", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p21", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p22", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p23", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p24", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p31", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p32", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p33", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p34", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p41", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p42", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p43", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p44", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);

	fp1 = fopen ("parameters/shaking.par", "w");

	valp = Tcl_GetVar2(interp, "cp", "first_shake", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "last_shake", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "maxPoints_shake", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "maxFrames_shake", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	
	fclose (fp1);


	fp1 = fopen ("parameters/orient.par", "w");

	valp = Tcl_GetVar2(interp, "cp", "pnrori", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "pdist", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "xp", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "yp", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "k1", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "k2", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "k3", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p1", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "p2", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "scx", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "she", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "interf", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);

	fclose (fp1);

	fp1 = fopen ("parameters/examine.par", "w");
	valp = Tcl_GetVar2(interp, "cp", "examineFlag", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "cp", "combineFlag", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	fclose (fp1);



	fp1 = fopen ("parameters/track.par", "w");

	valp = Tcl_GetVar2(interp, "tp", "dvxmin", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dvxmax", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dvymin", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dvymax", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dvzmin", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dvzmax", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dangle", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "dacc", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
	valp = Tcl_GetVar2(interp, "tp", "add", TCL_GLOBAL_ONLY);
	fprintf (fp1, "%s\n", valp);
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


	fclose (fp1);
        return TCL_OK;
}
