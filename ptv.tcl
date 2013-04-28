# --- Starting the Tcl GUI for the 3DPTV package ----
package require tooltip

lappend auto_path [file dirname [info script]]
lappend auto_path "."
# for debugging:
# puts $auto_path

wm title . "Measurement of Particles in flows"
wm iconname . "PTV"
wm geometry . +200+50

# Global Tcl variables 
global method methodmsg subt_ethz subt_poly cp mp tp sel clicki fy fx examine tbuf zoomfactor
set fx {}
set fy {}

set debug 1		;# generate extra information if set

# get examine value if specified
if {$argc  == "1"}  {set examine [lindex $argv 0]} else {set examine ""}

# Initialize the 3dPTV program to read the parameter files
init_proc_cmd

# -- Message strings at the top of the client window -- #
set subt_ethz "Using camera orientations for particle mapping (ETHZ method)"
set subt_poly "Using polynomials for particle mapping" 
if {$cp(method) == "ETHZ"} {set subtitle $subt_ethz} else {set subtitle $subt_poly}

label .msg1 -text "Multimedia Particle Positioning and Tracking" -font {Helvetica 12 bold}
label .msg2 -text $subtitle -font {Helvetica 10}
pack .msg1 -side top -pady {10 0}
pack .msg2 -side top -pady {0 10}

tooltip::tooltip .msg1 "Code developed by IPG and IfU ETH Zürich \n extended by TUE, last update: Dec. 2012"

# -- Create the main menu bar and add menu buttons -- #
# (the build functions are defined in the file ptv_menu.tcl)
frame .mbar -relief raised -bd 2
if { $cp(method) == "ETHZ" } {
	build_ETHZ_menubar .mbar 
} else {
	build_POLY_menubar .mbar
}
focus .mbar
pack .mbar -side top

# create image data structures for images
image create photo temp -width [ expr $cp(imx) ] -height [expr $cp(imy) ]

for {set i 1} {$i <= $mp(ncam)} {incr i} {
	image create photo image$i -width [ expr $cp(imx) ] -height [ expr $cp(imy) ]
	image create photo ori$i -width [ expr $cp(imx) ] -height [expr $cp(imy) ]
}

# messages during processing
listbox .text -height 5 
.text configure -font {Helvetica 13 bold} 
pack .text -pady 10 -expand yes -fill both



# -- function to switch between different mapping modes and calibration objects -- #
proc set_method {} {
    global cp method methodmsg subt_ethz subt_poly

	if {$cp(method) == "ETHZ" } { .msg2 configure -text $subt_ethz } \
	else 				        { .msg2 configure -text $subt_poly }
	
    foreach cw [winfo children .mbar] { destroy $cw }
    if {$cp(method) == "ETHZ"}	{ build_ETHZ_menubar .mbar } \
    else 						{ build_POLY_menubar .mbar }

    # Re-initialize 3dptv.exe by calling change_method_c() in jw_main.c
    change_method_cmd $cp(method) $cp(multi);
}
