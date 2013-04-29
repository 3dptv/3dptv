
# -- Constructing the traditional menu bar (ETHZ method) -- #

proc build_ETHZ_menubar {mbar} {
    mnu_start       .mbar.start
    mnu_pretracking .mbar.pre
    mnu_3Dcoord     .mbar.3d
    mnu_sequence    .mbar.seq
    mnu_tracking    .mbar.track
    mnu_calibration .mbar.calib
    mnu_changeparms .mbar.change
    mnu_demo        .mbar.demo
    mnu_options     .mbar.options
    mnu_quit        .mbar.quit
    
    pack $mbar.start .mbar.pre \
        .mbar.3d .mbar.seq .mbar.track .mbar.calib .mbar.demo .mbar.change \
        .mbar.options .mbar.quit -side left -ipadx 2 -expand 1 -fill x
#    tk_menuBar .mbar .mbar.start .mbar.pre \
#        .mbar.3d .mbar.seq .mbar.track .mbar .mbar.calib .mbar.demo .mbar.options
}


# -- Constructing the menu bar if using polynomials -- #

proc build_POLY_menubar {mbar} {
	mnu_start        .mbar.start
    mnu_pretracking  .mbar.pre
    mnu_3Dcoord      .mbar.3d
    mnu_sequence     .mbar.seq
    mnu_tracking     .mbar.track
    mnu_poly_cal     .mbar.calib
    mnu_poly_change  .mbar.change
    mnu_demo		 .mbar.demo
    mnu_options      .mbar.options
    mnu_quit         .mbar.quit
    
    pack .mbar.start .mbar.pre .mbar.3d .mbar.seq .mbar.track .mbar.calib .mbar.change \
        .mbar.options .mbar.quit -side left -ipadx 2 -expand 1 -fill x
    
#    tk_menuBar .mbar .mbar.start .mbar.pre \
#        .mbar.3d .mbar.seq .mbar.calib .mbar.change .mbar.options
}

# -- menu items not depending of the mapping method -- #

proc mnu_start {w} {			;# .mbar.start
    button $w -text "Start" -command "start_proc_cmd; bindingsstart"
    tooltip::tooltip $w "Initiates all variables. Needs to be clicked first, except when doing only calibration."
}

proc mnu_pretracking {w} {		;# .mbar.pre
    menubutton $w -text "Pretracking" -relief raised -underline 0 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "High Pass"       -command pre_processing_cmd
    $w.menu add command -label "Image Coord"     -command "bindings0;detection_proc_cmd"
    $w.menu add command -label "Correspondences" -command "bindings2;correspondences_cmd 0"

    tooltip::tooltip $w "To get a first idea how and howmany 2d particles are detected and to see  
howmany triplets and quadruplets are found for the given calibration and settings."
}

proc mnu_3Dcoord {w} {			;# .mbar.3d
    button $w -text "3D-Coordinates" -command "determination_cmd"
    tooltip::tooltip $w "Writes out the file rt_is.first containing the determined triplets and quadruplets, 
No., x, y, z, no. in cam1, no. in cam2, no. in cam3, no. in cam4."
}

proc mnu_sequence {w} {			;# .mbar.seq
    menubutton $w -text "Sequence" -relief raised -underline 0 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "Sequence with display" -command "sequence_cmd 1"
    $w.menu add command -label "Sequence without display" -command "sequence_cmd 0"

    tooltip::tooltip $w "Computes particle positions of frame-sequence defined in Change Parameters->Change Main Parameters."
}

proc mnu_tracking {w} {			;# .mbar.track
    menubutton $w -text "Tracking" -relief raised -underline 0 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "Detected Particles" -command "marktrack_cmd"
    $w.menu add command -label "Tracking with display" -command "trackcorr_cmd 1"
    $w.menu add command -label "Tracking without display" -command "trackcorr_cmd 0"
    $w.menu add command -label "Tracking backwards" -command "trackback_cmd 0"
    $w.menu add command -label "Sequence/Tracking" -command "sequence_cmd 0;trackcorr_cmd 0"
    $w.menu add command -label "Show Trajectories" -command "trajectories_cmd"
    $w.menu add command -label "VRML Tracks" -command "VRMLtracks"
    $w.menu add command -label "VRML Detection" -command "VRMLdetections"
    $w.menu add command -label "VRML Detection + Tracks" -command "VRMLdettracks"
	$w.menu add command -label "Tracking (just) in 3d" -command "ptv_cmd"
    
    tooltip::tooltip $w "Tracks points of frame-sequence defined in Change Parameters->Change Main Parameters,
according to Change Tracking Parameters."
}

proc mnu_calibration {w} {		;# .mbar.calib (only visible if using the ETHZ method)
    menubutton $w -text "Calibration" -relief raised -underline 0 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "Show Calib. Image" -command "set sel 1;calib_cmd;bindings1;" 
    $w.menu add command -label "Detection" -command "set sel 2;bindings0;calib_cmd"
    $w.menu add command -label "Manual orientation" -command "set sel 3;calib_cmd"
    $w.menu add command -label "Orientation with file" -command "set sel 4;calib_cmd"
    $w.menu add command -label "Show initial guess" -command "set sel 9;calib_cmd"
    $w.menu add command -label "Sortgrid" -command "set sel 5;calib_cmd"
	$w.menu add command -label "Sortgrid = initial guess" -command " set sel 14;calib_cmd"
	$w.menu add command -label "Show number on detected points" -command " set sel 15;calib_cmd"
	$w.menu add command -label "Sortgrid with file" -command " set sel 16;calib_cmd"
    $w.menu add command -label "Orientation" -command "set sel 6;calib_cmd"
    $w.menu add command -label "Orientation with particle positions (Sequence/Tracking/Shaking)" -command "sequence_cmd 2;trackcorr_cmd 2;  set sel 10;calib_cmd"
    $w.menu add command -label "Orientation with particle positions (Shaking)" -command "set sel 10;calib_cmd"
	$w.menu add command -label "Orientation with particle positions (Shaking, discarding bad 3d points)" -command " set sel 20;calib_cmd"
	$w.menu add command -label "Orientation from dumbbell (Sequence/Correction/Shaking)" -command "sequence_cmd 3; set sel 12;calib_cmd"
    $w.menu add command -label "Restore previous Orientation" -command "restore_cmd"
    $w.menu add command -label "Checkpoints" -command "set sel 7;calib_cmd"
    $w.menu add command -label "Ap figures" -command "set sel 8;calib_cmd"
	$w.menu add command -label "map mm to pixel" -command " set sel 30;calib_cmd"

    tooltip::tooltip $w "Opens the Calibration menu. 

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
}

proc mnu_demo {w} {			;# .mbar.demo
    menubutton $w -text "Demos" -relief raised -underline 0 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "Sequence 1" -command "flow_cmd 0"
    $w.menu add command -label "Sequence 2" -command "flow_cmd 1"
    $w.menu add command -label "Sequence 3" -command "flow_cmd 2"
    $w.menu add command -label "Sequence 4" -command "flow_cmd 3"
    tooltip::tooltip $w "Don't worry about this, its only kept for nostalgic reasons."
}

proc mnu_changeparms {w} {	;# .mbar.change (only visible if using the ETHZ method)
    menubutton $w -text "Change Parameters" -relief raised -underline 7 -menu $w.menu
    menu $w.menu
    $w.menu add command -label "Change Main Parameters"        -command mainpar
    $w.menu add command -label "Change Calibration Parameters" -command calpar
    $w.menu add command -label "Change Tracking Parameters"    -command trackpar
    tooltip::tooltip $w "Important three sub-menus to adjust control parameters. More hints are given within sub-menus."
}

	proc add_zoomentries {mnu zoomfacname cmdstr} {
		$mnu add radio -label "zoom   25 %" -command $cmdstr -variable $zoomfacname -value -4
		$mnu add radio -label "zoom   33 %" -command $cmdstr -variable $zoomfacname -value -3
		$mnu add radio -label "zoom   50 %" -command $cmdstr -variable $zoomfacname -value -2
		$mnu add radio -label "zoom 100 %"  -command $cmdstr -variable $zoomfacname -value 1
		$mnu add radio -label "zoom 200 %"  -command $cmdstr -variable $zoomfacname -value 2
		$mnu add radio -label "zoom 300 %"  -command $cmdstr -variable $zoomfacname -value 3
		$mnu add radio -label "zoom 400 %"  -command $cmdstr -variable $zoomfacname -value 4
		$mnu add radio -label "zoom 600 %"  -command $cmdstr -variable $zoomfacname -value 6
		$mnu add radio -label "zoom 800 %"  -command $cmdstr -variable $zoomfacname -value 8
	}

	proc add_framesizes {mnu cmdstr} {
	 $mnu add checkbutton -label "size fixed" -onvalue 1 -offvalue 0 -variable zoompar(fixed)
	$mnu add sep
	$mnu add radio -label "frame   25 %" -command $cmdstr -variable framefactor -value -4
	$mnu add radio -label "frame   33 %" -command $cmdstr -variable framefactor -value -3
	$mnu add radio -label "frame   50 %" -command $cmdstr -variable framefactor -value -2
	$mnu add radio -label "frame 100 %"  -command $cmdstr -variable framefactor -value 1
	}

proc mnu_options {w} {		;# .mbar.options (only visible if using the ETHZ method)
    global cp method
    menubutton $w -text "Options" -relief raised -underline 7 -menu $w.menu
    menu $w.menu
    $w.menu add radiobutton -label "Original ETHZ method of calibration"    -value ETHZ -variable cp(method) -command "set_method"
    $w.menu add radiobutton -label "Using polynomials for particle mapping" -value POLY -variable cp(method) -command "set_method"
    $w.menu add separator
	if {$cp(method) != "ETHZ"} {
		$w.menu add radiobutton -label "3D-body calibration"      -value 0 -variable cp(multi) -command "set_method"
		$w.menu add radiobutton -label "Multi planes calibration" -value 1 -variable cp(multi) -command "set_method"
		$w.menu add separator
	}
	$w.menu add command -label "Refresh Images" -command "refresh"
	$w.menu add command -label "Show original" -command "showori"
	$w.menu add command -label "Clear Canvas" -command clearcam
	$w.menu add sep
	$w.menu add cascade -label "Zoom"        -menu $w.menu.zoommenu -underline 0
	$w.menu add cascade -label "Frame sizes" -menu $w.menu.framesize -underline 0

	proc add_zoomentries2 {mnu} {
		$mnu add radio -label "zoom   25 %" -command "test" -variable zoomfactor -value -4
		$mnu add radio -label "zoom   33 %" -command "test" -variable zoomfactor -value -3
	}
	
	# All views are zoomed immediately to the clicked value, so no bindings here. ad holten, 04-2013
	menu $w.menu.zoommenu -tearoff 0
	add_zoomentries $w.menu.zoommenu zoomfactor "mouse_cmd 10 \$zoomfactor \$zoompar\(fixed\)"

	menu $w.menu.framesize -tearoff 0
	add_framesizes $w.menu.framesize "new_framesizes \$framefactor"  
}

proc mnu_quit {w} {			;# mbar.quit
    button $w -text "Quit" -command "quit_cmd;destroy ." 
    tooltip::tooltip $w "It is good practice to quit from time to time to be SURE that changed values are updated to everywhere in the code."
}


# -----------  menu items if using polynomials ----------

proc change_plane {m step} {
	global calib_menu plane nplanes cp

	set plane [expr $plane + $step]
	set enbnext normal
	set enbprev normal
	
	if {$plane <= 1}        {set enbprev disabled; set plane 1}
	if {$plane >= $nplanes} {set enbnext disabled; set plane $nplanes}
	set cp(plane_id) $plane

	set stat [$m entrycget 6 -state]
	if {$stat != $enbnext} {$m entryconfigure 6 -state $enbnext}
	set stat [$m entrycget 7 -state]
	if {$stat != $enbprev} {$m entryconfigure 7 -state $enbprev}
	
	scan [$m entrycget 1 -label] "%\[^:\]" s
	$m entryconfigure 1 -label [format "%s: %d" $s $plane]
}

proc mnu_poly_cal {w} {
	global cp calib_menu plane nplanes
    menubutton $w -text "Calibration" -relief raised -underline 0 -menu $w.menu

	set m $w.menu
	menu $w.menu -tearoff 1
	set calib_menu $w.menu

	if {$cp(multi)} {
		set plane 1
		set nplanes 4
	
		$w.menu add command -label "Show Calib. Image of plane: $plane"\
															   -command "set sel 1; calpoly_cmd; bindings1;"
		$w.menu add command -label "    Detection"             -command "set sel 2; bindings0; calpoly_cmd"
		$w.menu add command -label "    Manual orientation"    -command "set sel 3; calpoly_cmd"
		$w.menu add command -label "    Orientation with file" -command "set sel 4; calpoly_cmd"
		$w.menu add command -label "    Sortgrid"              -command "set sel 5; calpoly_cmd"
		$w.menu add command -label "Go to next plane"          -command "change_plane $w.menu 1"
		$w.menu add command -label "Go to previous plane"      -command "change_plane $w.menu -1" -state disabled
	} else {
		$w.menu add command -label "Show Calib. Images"	   -command "set sel 1; calpoly_cmd; bindings1;"
		$w.menu add command -label "Detection"             -command "set sel 2; bindings0; calpoly_cmd"
		$w.menu add command -label "Manual orientation"    -command "set sel 3; calpoly_cmd"
		$w.menu add command -label "Orientation with file" -command "set sel 4; calpoly_cmd"
		$w.menu add command -label "Sortgrid"              -command "set sel 5; calpoly_cmd"
		$w.menu add command -label "Sortgrid with file"    -command "set sel 16; calpoly_cmd"
	}
	$w.menu add command -label "Show number on detected points" -command "set sel 15; calpoly_cmd"
	$w.menu add separator 

	$w.menu add command -label "Orientation"           -command "set sel 6; calpoly_cmd"
    $w.menu add command -label "Orientation with particle positions (Sequence/Tracking/Shaking)" \
                                                       -command "sequence_cmd 2; trackcorr_cmd 2; set sel 10; calpoly_cmd"
    $w.menu add command -label "Orientation with particle positions (Shaking)" \
                                                       -command "set sel 10; calpoly_cmd"
    # $w.menu add command -label "Restore previous Orientation" -command "restore_cmd"
    $w.menu add command -label "Checkpoints"           -command "set sel 7;calpoly_cmd"
	
    tooltip::tooltip $w "Calibration functions, using polynomials for mapping particle positions."
}

proc mnu_poly_multi_cal {w} {
	global calib_menu plane nplanes
    menubutton $w -text "Calibration" -relief raised -underline 0 -menu $w.menu

	set plane 1
	set nplanes 4
    set m $w.menu
    menu $w.menu -tearoff 1
	
	set calib_menu $w.menu
	
    $w.menu add command -label "Show Calib. Image of plane: $plane"\
												           -command "set sel 1; calpoly_cmd; bindings1;"
    $w.menu add command -label "    Detection"             -command "set sel 2; bindings0; calpoly_cmd"
    $w.menu add command -label "    Manual orientation"    -command "set sel 3; calpoly_cmd"
    $w.menu add command -label "    Orientation with file" -command "set sel 4; calpoly_cmd"
    $w.menu add command -label "    Sortgrid"              -command "set sel 5; calpoly_cmd"
    $w.menu add command -label "Go to next plane"          -command "change_plane $w.menu 1"
    $w.menu add command -label "Go to previous plane"      -command "change_plane $w.menu -1" -state disabled
    $w.menu add separator 
    $w.menu add command -label "Orientation"           -command "set sel 6; calpoly_cmd"
    $w.menu add command -label "Orientation with particle positions (Sequence/Tracking/Shaking)" \
                                                       -command "sequence_cmd 2;trackcorr_cmd 2;  set sel 10;calib_cmd"
    $w.menu add command -label "Orientation with particle positions (Shaking)" \
                                                       -command "set sel 10;calib_cmd"
    $w.menu add command -label "Restore previous Orientation" -command "restore_cmd"
    $w.menu add command -label "Checkpoints"           -command "set sel 7;calib_cmd"
    $w.menu add command -label "Ap figures"            -command "set sel 8;calib_cmd"
	
    tooltip::tooltip $w "Calibration functions, using polynomials for mapping particle positions."
}

proc mnu_poly_change {w} {
    global cp
    menubutton $w -text "Change Parameters" -relief raised -underline 0 -menu $w.menu

    menu $w.menu -tearoff 0
    $w.menu add command -label "Main parameters..."            -command mainpoly
    # $w.menu add separator
    $w.menu add command -label "Calibration (3D-body) ..."     -command poly_calpar
    $w.menu add command -label "Calibration (multi planes)..." -command multi_calpar
    # $w.menu add separator
    $w.menu add command -label "Tracking Parameters..."        -command trackpar
    $w.menu add separator
    $w.menu add command -label "Debugging..."                  -command examine3

    # Enable/Disable some of the menu items 
    if {$cp(multi) == 1} {
        $w.menu entryconfigure "Calibration (3D-body) ..." -state disabled 
        $w.menu entryconfigure "Calibration (multi planes)..." -state normal
    } else {
        $w.menu entryconfigure "Calibration (3D-body) ..." -state normal 
        $w.menu entryconfigure "Calibration (multi planes)..." -state disabled
    } 
    tooltip::tooltip $w "View/Edit all control parameters."
}

proc mnu_poly_cascade {w} {
    menubutton $w -text "Polynomials" -relief raised -underline 0 -menu $w.menu

    set m $w.menu
    menu $w.menu
    $w.menu add cascade -label "Change Parameters" -menu $w.menu.params -underline 0
    $w.menu add command -label "Calibration using a 3D-body" -command 3Dbody::create
    $w.menu add command -label "Calibration using Z-planes" -command Zcalib::create
    $w.menu add separator
    $w.menu add command -label "Run the shaking procedure once" -command "call_shaking"
    $w.menu add command -label "Start the 'Sequence, Tracking, Shaking' cycles" -command "fit_by_shaking 0"

    set m1 $w.menu.params
    menu $m1 -tearoff 0
    $m1 add command -label "Main parameters"           -command mainpoly
    $m1 add command -label "3D-Body"                   -command 3dbodypar
    $m1 add command -label "multidlg Planes"           -command multidlg::create
    $m1 add command -label "Tracking Parameters"       -command trackpar
    $m1 add command -label "Parameters shaking cycles" -command shakingcyclespar

    tooltip::tooltip $w "Using Polynomials for mapping particle positions."
}


