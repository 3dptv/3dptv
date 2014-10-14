//////////////////////////////////////////////////////////////////////////////
//
// this code was written by Beat Luthi at IfU, ETH Zürich, Okt 2007
// and updatet together with Marc Wolf to have flexible kernel length for the
// cubic polynomial fits along the trajectories, it used to be fixed at 21
//
// is represents an attempt to have ONE clean non-GUI version of the postPorcessing codes
// that float around in various Borland versions
//
// luethi@photrack.ch
//
// last update/change: August 2011 Marc Wolf
//
//////////////////////////////////////////////////////////////////////////////

1) (simple)
put the file "input.inp" in the C: folder
adjust input.inp to your needs
put the excecutable "post_process.exe" anywhere you like
run "post_process.exe" 

enjoy the output files xuap.xxx and trajPoint.xxx
the format for "xuap.xxx" are 14 columns:
link_past, link_future, 
x_raw, y_raw, z_raw,
x, y, z, u, v, w, ax, ay, az
control flag
the format for the "trajPOint.xxx" are 33 columns, described in the picture 
"list_of_output_columns.jpg"

2) (more involved)
if you want to change the location and/or name of your input file go the line 71 of
"post_process.cpp", rebuild, run and enjoy
