
set auto_path "~willneff/prog/tk84ptv . $auto_path"
#set auto_path "G:/prog/tk84ptv . $auto_path"
 
wm title . "Measurement of Particles in flows"
wm iconname . "PTV"
wm geometry . +200+50

# Global tclvariables 
global cp mp tp sel clicki fy fx examine tbuf zoomfactor
set fx {}
set fy {}

# get examine value if specified

if { $argc  == "1" }  { set examine [lindex $argv 0] } else { set examine "" }

# #######################################
label .msg -wraplength 4i -justify left -text "Multimedia Particle Positioning and Tracking"  \
-font {Helvetica 12 bold}
pack .msg -side top -pady 8

# Initialize Programm to read parameterfiles

init_proc_cmd

# Mainmenu ###################################

frame .mbar -relief raised -bd 2

button .mbar.start -text "Start" -command "start_proc_cmd; bindingsstart "

menubutton .mbar.pre -text "Pretracking" -relief raised \
    -underline 0 -menu .mbar.pre.menu

button .mbar.3d -text "3D-Coordinates" -command "determination_cmd"

menubutton .mbar.track -text "Tracking" -relief raised \
    -underline 0 -menu .mbar.track.menu

menubutton .mbar.calib -text "Calibration" -relief raised \
    -underline 0 -menu .mbar.calib.menu

menubutton .mbar.seq -text "Sequence" -relief raised \
    -underline 0 -menu .mbar.seq.menu

menubutton .mbar.demo -text "Demos" -relief raised \
    -underline 0 -menu .mbar.demo.menu

menubutton .mbar.change -text "Change Parameters" -relief raised \
    -underline 7 -menu .mbar.change.menu

menubutton .mbar.options -text "Options" -relief raised \
    -underline 0 -menu .mbar.options.menu

button .mbar.quit -text "Quit" -command "quit_cmd;destroy ."

menu .mbar.pre.menu
.mbar.pre.menu add command -label "High Pass" -command pre_processing_cmd
.mbar.pre.menu add command -label "Image Coord" -command "bindings0;detection_proc_cmd"
.mbar.pre.menu add command -label "Correspondences" -command " bindings2;correspondences_cmd"

menu .mbar.seq.menu
.mbar.seq.menu add command -label "Sequence with display" -command "sequence_cmd 1"
.mbar.seq.menu add command -label "Sequence without display" -command "sequence_cmd 0"

menu .mbar.track.menu
.mbar.track.menu add command -label "Detected Particles" -command "marktrack_cmd"
.mbar.track.menu add command -label "Tracking with display" -command "trackcorr_cmd 1"
.mbar.track.menu add command -label "Tracking without display" -command "trackcorr_cmd 0"
.mbar.track.menu add command -label "Tracking backwards" -command "trackback_cmd 0"
.mbar.track.menu add command -label "Sequence/Tracking" -command "sequence_cmd 0;trackcorr_cmd 0"
.mbar.track.menu add command -label "Show Trajectories" -command "trajectories_cmd"
.mbar.track.menu add command -label "VRML Tracks" -command "VRMLtracks"
.mbar.track.menu add command -label "VRML Detection" -command "VRMLdetections"
.mbar.track.menu add command -label "VRML Detection + Tracks" -command "VRMLdettracks"
#.mbar.track.menu add command -label "PTV" -command "ptv_cmd"

menu .mbar.calib.menu
.mbar.calib.menu add command -label "Show Calib. Image" -command "set sel 1;calib_cmd;bindings1;"
.mbar.calib.menu add command -label "Detection" -command " set sel 2;bindings0;calib_cmd"
.mbar.calib.menu add command -label "Manual orientation" -command " set sel 3;calib_cmd"
.mbar.calib.menu add command -label "Orientation with file" -command " set sel 4;calib_cmd"
.mbar.calib.menu add command -label "Sortgrid" -command " set sel 5;calib_cmd"
.mbar.calib.menu add command -label "Orientation" -command " set sel 6;calib_cmd"
.mbar.calib.menu add command -label "Checkpoints" -command " set sel 7;calib_cmd"
.mbar.calib.menu add command -label "Ap figures" -command "set sel 8;calib_cmd"
  
menu .mbar.demo.menu
.mbar.demo.menu add command -label "Sequence 1" -command "flow_cmd 0"
.mbar.demo.menu add command -label "Sequence 2" -command "flow_cmd 1"
.mbar.demo.menu add command -label "Sequence 3" -command "flow_cmd 2"
.mbar.demo.menu add command -label "Sequence 4" -command "flow_cmd 3"

menu .mbar.change.menu
.mbar.change.menu add command -label "Change Main Parameters" \
    -command mainpar
.mbar.change.menu add command -label "Change Calibration Parameters" -command calpar
.mbar.change.menu add command -label "Change Tracking Parameters" -command trackpar


menu .mbar.options.menu
.mbar.options.menu add command -label "Refresh Images" -command "refresh"
.mbar.options.menu add command -label "Show original" -command "showori"
.mbar.options.menu add command -label "Clear Canvas" -command clearcam
.mbar.options.menu add command -label "Zoom  8x" -command "set zoomfactor 8;bindings4"
.mbar.options.menu add command -label "Zoom 16x" -command "set zoomfactor 16;bindings4"
.mbar.options.menu add command -label "Zoom 32x" -command "set zoomfactor 32;bindings4"
.mbar.options.menu add command -label "Zoom 64x" -command "set zoomfactor 64;bindings4"
.mbar.options.menu add command -label "Zoom 128x" -command "set zoomfactor 128;bindings4"
.mbar.options.menu add command -label "Zoom 256x" -command "set zoomfactor 256;bindings4"
.mbar.options.menu add command -label "Zoom corr  2x" -command "set zoomfactor 2;bindings3"
.mbar.options.menu add command -label "Zoom corr  4x" -command "set zoomfactor 4;bindings3"
.mbar.options.menu add command -label "Zoom corr  6x" -command "set zoomfactor 6;bindings3"
.mbar.options.menu add command -label "Zoom corr  8x" -command "set zoomfactor 8;bindings3"
.mbar.options.menu add command -label "Zoom corr 16x" -command "set zoomfactor 16;bindings3"
.mbar.options.menu add command -label "Zoom corr 32x" -command "set zoomfactor 32;bindings3"
.mbar.options.menu add command -label "Zoom corr 64x" -command "set zoomfactor 64;bindings3"
##################################################################
pack .mbar.start .mbar.pre \
        .mbar.3d .mbar.seq .mbar.track .mbar.calib .mbar.demo .mbar.change \
	.mbar.options .mbar.quit -side left -expand 1 -fill x

tk_menuBar .mbar .mbar.start .mbar.pre \
 .mbar.3d .mbar.seq .mbar.track .mbar.calib .mbar.demo .mbar.options
focus .mbar

pack .mbar -side top

# create image data structures for images

 image create photo temp -width [ expr $cp(imx) ] -height [expr $cp(imy) ]

for { set i 1} { $i <= $mp(ncam)} { incr i } {
    image create photo image$i -width [ expr $cp(imx) ] -height [ expr $cp(imy) ]
    image create photo ori$i -width [ expr $cp(imx) ] -height [expr $cp(imy) ]
}

# messages during processing
 
listbox .text -height 5 
.text configure -font {Helvetica 13 bold} 
pack .text -pady 10 -expand yes -fill both

