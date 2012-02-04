#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <tcl.h>
#include <tk.h>

#include "tiff.h"
#include "tiffio.h"

#include "typedefs.h"
#include "globals.h"

extern Tk_PhotoImageFormat tkImgFmtTIF;

extern Tcl_CmdProc init_proc_c;
extern Tcl_CmdProc start_proc_c;
extern Tcl_CmdProc done_proc_c;
extern Tcl_CmdProc detection_proc_c;
extern Tcl_CmdProc pre_processing_c;
extern Tcl_CmdProc correspondences_proc_c;
extern Tcl_CmdProc determination_proc_c;
extern Tcl_CmdProc sequence_proc_c;
extern Tcl_CmdProc calibration_proc_c;
extern Tcl_CmdProc restore_proc_c;
extern Tcl_CmdProc quit_proc_c;
extern Tcl_CmdProc mouse_proc_c;
extern Tcl_CmdProc flow_demo_c;
extern Tcl_CmdProc mark_track_c;
extern Tcl_CmdProc trackcorr_c;
extern Tcl_CmdProc trackback_c;
extern Tcl_CmdProc nearestinnext_c;
extern Tcl_CmdProc tracksequence_c;
extern Tcl_CmdProc tracking;
extern Tcl_CmdProc trajectories_c;
extern Tcl_CmdProc vrmltracks_c;
extern Tcl_CmdProc trajectories_c;
extern Tcl_CmdProc vrmldetections_c;
extern Tcl_CmdProc vrmldettracks_c;
extern Tcl_CmdProc seq_track_proc_c;
