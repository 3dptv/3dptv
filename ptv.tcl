
package require tooltip

lappend auto_path [file dirname [info script]]
lappend auto_path "."
# for debugging:
# puts $auto_path

 
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


tooltip::tooltip .msg "code developed by IPG and IfU ETH Zürich \n last update: May 08"


# Initialize Programm to read parameterfiles

init_proc_cmd

# Mainmenu ###################################

frame .mbar -relief raised -bd 2

button .mbar.start -text "Start" -command "start_proc_cmd; bindingsstart "

tooltip::tooltip .mbar.start "Initiates all variables. Needs to be clicked first, except when doing only calibration."

menubutton .mbar.pre -text "Pretracking" -relief raised \
    -underline 0 -menu .mbar.pre.menu

tooltip::tooltip .mbar.pre "To get a first idea how and howmany 2d particles are detected and to see  
howmany triplets and quadruplets are found for the given calibration and settings."

button .mbar.3d -text "3D-Coordinates" -command "determination_cmd"

tooltip::tooltip .mbar.3d "Writes out the file rt_is.first containing the determined triplets and quadruplets, 
No., x, y, z, no. in cam1, no. in cam2, no. in cam3, no. in cam4."

menubutton .mbar.track -text "Tracking" -relief raised \
    -underline 0 -menu .mbar.track.menu
    
tooltip::tooltip .mbar.track "Tracks points of frame-sequence defined in Change Parameters->Change Main Parameters,
according to Change Tracking Parameters."

menubutton .mbar.calib -text "Calibration" -relief raised \
    -underline 0 -menu .mbar.calib.menu
    
tooltip::tooltip .mbar.calib "Opens the Calibration menu. 

A typical sequence is:
'Show Calib. Image' (see if the paths are correct.)
'Detection' (adjust 'Target recognition on plate' in 'Change Calibation Parameters')
'Manual orientation' (necessary only the first time, BEWARE: MOST TIME IS LOST HERE, WRONG CLICKING AND NOT NOTICING IT!)
'Orientation with file' (if 'Manual orientation' has been correctly performed)
'Sortgrid' (here the code tries to assign detected target points to points that are defined in the 'File on Coodrinates on Plate'. Search radius is defined in 'parameters/sortgrid.par')
'Orientation' (this is where Calibration happens, results are written to '<Calibration image>.ori', Orientation parameters are controlled from 'Change Calibration Parameters)

after you have detected  and linked a few triplets and quadruplets but you are not really happy yet with accuracy and number of points
you can replace target points with actual deteced particles,
i.e., you can re-calibrate for the situation during the actual experiment:
'Orientation with particle positions (Sequence/Tracking/Shaking)'
or just
'Orientation with particle positions (Shaking)'
this is controlled by 'Shaking parameters inside 'Change Calibration parameters'

'Restore previous Orientation' is an undo possibility for the last calibration."

menubutton .mbar.seq -text "Sequence" -relief raised \
    -underline 0 -menu .mbar.seq.menu
    
tooltip::tooltip .mbar.seq "Computes particle positions of frame-sequence defined in Change Parameters->Change Main Parameters."

menubutton .mbar.demo -text "Demos" -relief raised \
    -underline 0 -menu .mbar.demo.menu
tooltip::tooltip .mbar.demo "Don't worry about this, its only kept for nostalgic reasons."


menubutton .mbar.change -text "Change Parameters" -relief raised \
    -underline 7 -menu .mbar.change.menu
tooltip::tooltip .mbar.change "Important three sub-menus to adjust control parameters. More hints are given within sub-menus."


menubutton .mbar.options -text "Options" -relief raised \
    -underline 0 -menu .mbar.options.menu
tooltip::tooltip .mbar.options "There are not too many options, so don't raise your hopes too much."


button .mbar.quit -text "Quit" -command "quit_cmd;destroy ." 
tooltip::tooltip .mbar.quit "It is good practice to quit from time to time to be SURE that changed values are updated to everywhere in the code."


menu .mbar.pre.menu
.mbar.pre.menu add command -label "High Pass" -command pre_processing_cmd
.mbar.pre.menu add command -label "Image Coord" -command "bindings0;detection_proc_cmd"
.mbar.pre.menu add command -label "Correspondences" -command " bindings2;correspondences_cmd 0"

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
.mbar.calib.menu add command -label "Show initial guess" -command " set sel 9;calib_cmd"
.mbar.calib.menu add command -label "Sortgrid" -command " set sel 5;calib_cmd"
.mbar.calib.menu add command -label "Orientation" -command " set sel 6;calib_cmd"
.mbar.calib.menu add command -label "Orientation with particle positions (Sequence/Tracking/Shaking)" -command "sequence_cmd 2;trackcorr_cmd 2;  set sel 10;calib_cmd"
.mbar.calib.menu add command -label "Orientation with particle positions (Shaking)" -command " set sel 10;calib_cmd"
.mbar.calib.menu add command -label "Restore previous Orientation" -command "restore_cmd"
.mbar.calib.menu add command -label "Checkpoints" -command " set sel 7;calib_cmd"
.mbar.calib.menu add command -label "Ap figures" -command "set sel 8;calib_cmd"

tooltip::tooltip .mbar.calib.menu -index 0 "This is a menu tooltip"


  
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

