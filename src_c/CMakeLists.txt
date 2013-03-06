# This is the 3dptv/src_c  CMakeLists.txt file

message(STATUS "gathering ptv project files from 3dptv/src_c/:")

set(TARGET_NAME "3dptv")

# gather the header files

file ( GLOB ptv_hdrs 
globals.h
ptv.h
tcl.h
tiff.h
tiffio.h
tiffvers.h
tk.h
typedefs.h
)

# gather the source files

file ( GLOB ptv_srcs 
change_parameter.c
checkpoints.c
correspondences.c
demo.c
draw.c
epi.c
image_processing.c
imgcoord.c
intersect.c
jw_ImgFmtTIF.c
jw_main.c
jw_ptv.c
lsqadj.c
mousefunction.c
multimed.c
orientation.c
peakfitting.c
pointpos.c
ptv.c
ray_tracing.c
resource1.h
resource2.h
rotation.c
segmentation.c
sortgrid.c
tools.c
track.c
trafo.c
ttools.c
vrml.c
)

# organize the files into folders inside the project

source_group("Header Files" FILES ${ptv_hdrs})
source_group("Source Files" FILES ${ptv_srcs})
source_group("Tcl Files" FILES ${ptv_tcl})

add_executable(${TARGET_NAME} ${ptv_srcs} ${ptv_hdrs} ${ptv_tcl} )

if(WIN32)
  set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE C)  
  set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libc")
endif()

set(EXTRA_LIBS ${EXTRA_LIBS} ${TCL_LIBRARY} ${TK_LIBRARY} ${TIFF_LIBRARY})
target_link_libraries(${TARGET_NAME} ${EXTRA_LIBS} )