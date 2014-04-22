#----------------------------------------------------------------------#
#    Makefile for PTV                                    	       #
#----------------------------------------------------------------------#


INC_DIR1 = C:\Tcl\include\

TCL_LIB = C:\Tcl\lib\tcl84.lib 
TK_LIB = C:\Tcl\lib\tk84.lib
TIFF_LIB = H:\tk84ptv\src_c\libtiff.lib

OBJ = jw_main.obj \
	jw_ImgFmtTIF.obj \
	jw_ptv.obj \
	tools.obj \
	change_parameter.obj \
	segmentation.obj \
	image_processing.obj \
	draw.obj \
	peakfitting.obj \
	multimed.obj \
	trafo.obj \
	correspondences.obj \
	ray_tracing.obj \
	lsqadj.obj \
	epi.obj \
	imgcoord.obj \
	intersect.obj \
	rotation.obj \
	orientation.obj \
	pointpos.obj \
	sortgrid.obj \
	checkpoints.obj \
	mousefunction.obj \
	demo.obj \
	track.obj \
	ttools.obj \
	vrml.obj \
	ptv.obj


all:	$(OBJ)
	cl -o jw_prog $(OBJ) $(TCL_LIB) $(TK_LIB) $(TCL_LIB) $(TK_LIB) $(TIFF_LIB)

$(OBJ):  
	cl -c -I$(INC_DIR1) $*.c
