

Installation and Handling of the PTV-Software
=============================================

The ptvmanual.pdf contains information to use the
software and gives a description of the input data
file which have to be provided.

The test data set contains the following:

Cam1*                   Cam3.addpar*            man_ori.dat*
Cam1.addpar*            Cam3.ori*               parameters/
Cam1.ori*               Cam4*                   ptvmanual.pdf*
Cam2*                   Cam4.addpar*            res/
Cam2.addpar*            Cam4.ori*               start.bat*
Cam2.ori*               calFieldApril.txt*
Cam3*                   img/


images for calibration, camera orientation data, files for
additional parameters as well as the coordinate file of the
points on the reference.
In man_ori.dat the manually measured image coordinates for
the pre-orientation (for calibration purpose) are stored.

- subdirectory /img  contains the image sequences 
- subdirectory /res  for storage of results
- subdirectory /parameters  contains the parameter files

The data mentioned above are the data for the experiment itself.
To avoid confusion this data should be kept separated to the
software data!

The software for PTV is stored under /tk84ptv.
To start the software >> click double on start.bat, which establishes
the link to the software. Project and software data should not be
confused. To start the software it is sufficient that a start.bat-file,
thus the software (and code) can be stored independently.


Tcl/Tk-Installation
===================

The install executable for Tcl/Tk is ActiceTcl8.4.2-win32-ix86.exe
in this directory. Or can be downloaded from a webpage.

Download under http://downloads.activestate.com/ActiveTcl/Windows/8.4.2/

After installing Tcl/Tk 8.4.2 all files with extension *.tcl should
appear with the Tcl/Tk-symbol (feather). Otherwise repeat installation.
Make sure that the flag for the extensions (*.tcl) is included in the
path. Otherwise not all needed dll-files can be found from arbitrary
directories on the PC.


Compilation of the source code of PTV
======================================

In the /tk84ptv directory You will find the following data:

index		script to generate tclIndex (/ might be missing in the generated tclIndex!)
ptv.tcl		main script to start graphical user interface (Windows)
ptvunix.tcl	dito for Unix
start		start file for Unix
tclIndex	Index with relative paths to Tcl functions
/src_c		source code directory
/src_tcl	tcl script directory

The contents of the /src_c:

change_parameter.c*     jw_main.c*              segmentation.c*
checkpoints.c*          jw_ptv.c*               sortgrid.c*
correspondences.c*      libtiff.lib*            testvrml.c*
demo.c*                 lsqadj.c*               tiff.h*
draw.c*                 mousefunction.c*        tiffio.h*
epi.c*                  multimed.c*             tiffvers.h*
globals.h*              orientation.c*          tools.c*
homemakefile*           peakfitting.c*          track.c*
homemakefile.mak*       pointpos.c*             trafo.c*
image_processing.c*     ptv.c*                  ttools.c*
imgcoord.c*             ptv.h*                  typedefs.h*
intersect.c*            ray_tracing.c*          unixmakefile*
jw_ImgFmtTIF.c*         rotation.c*             vrml.c*


The contents of the /src_tcl:

button.tcl*     display.tcl*    mainpar.tcl
calpar.tcl*     draw.tcl*       trackpar.tcl*


The source code is written in C in combination with Tcl/Tk.
The directory /src_c contains a makefile (homemakefile.mak)
which can be open with Microsoft Visual C++. During opening
this file, a new workspace will be generated.

Notice: The paths to the libs in the makefile have to be adjusted:

INC_DIR1 = C:\Tcl\include\

TCL_LIB = C:\Tcl\lib\tcl84.lib 
TK_LIB = C:\Tcl\lib\tk84.lib
TIFF_LIB = H:\tk84ptv\src_c\libtiff.lib

For compilation in the DOS-prompt first perform vcvars32.bat for 
initalization, after that nmake -f homemakefile.mak


IMPORTANT!
==========

Before running the software some paths have to be set.
In start.bat, may contain:

G:/tk84ptv/src_c/jw_prog G:/tk84ptv/ptv.tcl

has to be modified to the actual position on the PC.
First the path to the jw_prog.exe followed by the
path to the ptv.tcl script file.

In addition, change path in first line of ptv.tcl

first line:

set auto_path "G:/tk84ptv . $auto_path"

change to according path on PC!

The start.bat should be copied to the project data directory.
Start with double click!
