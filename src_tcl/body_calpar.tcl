#
# This dialog is used for changing the calibration parameters if polynomials are used
#

package require BWidget 1.9.1

proc poly_calpar {} {
    global cp mp pf

    set w .changemain
    catch {destroy $w}
    toplevel .changemain

    set bold12 {Helvetica 12 bold}
    set bold10 {Helvetica 10 bold}

	if {$cp(multi)} {
		wm title $w "Changing Calibration Parameters, Multi Planes"

		frame $w.title
		label $w.title.l -text "Calibration Parameters for plane: " -font $bold12

		# ------ Plane index -------------------------- #
		spinbox $w.title.plane -from 1 -to 10 -width 4 \
			-command {change_plane %s} -validate key -vcmd {string is integer %P}
		pack $w.title.l $w.title.plane -side left
		pack $w.title  -anchor center -pady 10
	} else	{
		wm title $w "Changing Calibration Parameters"
		label $w.title -text "Calibration Parameters" -font $bold12
		pack $w.title  -anchor center -pady 10
	}

    # -- Calibration images ------------------------ #
    set title "Calibration images"
    set tf1 [TitleFrame $w.tf1 -text $title -font $bold10]
    
    calibration_images [$tf1 getframe]
    pack $tf1 -fill x -padx 10 -pady 3

    # -- Target coordinates file  ------------------ #
    set title "Coordinates on plate"
    set tf2 [TitleFrame $w.tf2 -text $title -font $bold10]
    
    coordinates_on_plate [$tf2 getframe]
    pack $tf2 -fill x -padx 10 -pady 3

	# -- Point numbers for manual pre-orientation -- #
    set title "Point numbers for manual pre-orientation"
    set tf3 [TitleFrame $w.tf3 -text $title -font $bold10]
    
    pointnumbers [$tf3 getframe]
    pack $tf3 -fill x -padx 10 -pady 3

    # -- Target recognition on plate --------------- #
	if {$cp(multi)} {
		set title "Target recognition on plate, (same for all planes)"
	} else {
		set title "Target recognition on plate"
	}
    set tf4 [TitleFrame $w.tf4 -text $title -font $bold10]
    
    target_recognition [$tf4 getframe]
    pack $tf4 -fill x -padx 10 -pady 3
    
	# -- Degree of polynomials --------------------- #
    set title "Degree of polynomials"
    set tf5 [TitleFrame $w.tf5 -text $title -font $bold10]
    
    degree_of_polynomials [$tf5 getframe]
    pack $tf5 -fill x -padx 10 -pady 3

    # -------- Shaking parameters ------------------ #
    set title "Shaking parameters"
    set tf6 [TitleFrame $w.tf6 -text $title -font $bold10]

    shaking_parameters [$tf6 getframe]
    pack $tf6 -fill x -padx 10 -pady 3
    
	button $w.ok -text OK -command "update_parms; done_proc_cmd; destroy $w"
    pack $w.ok -side bottom -pady 5
}

proc update_parms {} {
	global cp
	
	# if multi planes, save values to the current plane parameters
	if {$cp(multi)} {
		set i $cp(plane_id)
		set cp(platecoord$i) $cp(platecoord) 
		for {set j 1} {$j<=4} {incr j} {
			set cp(cal$i$j) $cp(calp$j)
			for {set k 1} {$k<=4} {incr k} {
				set cp(p$i$j$k) $cp(p$j$k)
			}
		}
	}
}

proc change_plane {iplane} {
	global cp mp
	
	# save dialog values to the current plane parameters
	update_parms

	set i $iplane	;# show the values for the new plane
	set cp(plane_id) $i
	set cp(platecoord) $cp(platecoord$i)
	for {set j 1} {$j<=4} {incr j} {
		set cp(calp$j) $cp(cal$i$j)
		for {set k 1} {$k<=4} {incr k} {
			set cp(p$j$k) $cp(p$i$j$k) 
		}
	}
}


proc calibration_images { w } {
	global cp
	frame $w.col1
	frame $w.col2
	
	if {$cp(multi)} {
		set i $cp(plane_id)
		for {set j 1} {$j<=4} {incr j} {
			set cp(calp$j) $cp(cal$i$j)
		}
	}
	
    for {set i 1} {$i<=4} {incr i} {
        frame $w.$i
        set imgname $w.$i.imgname
        
        frame $imgname
        entry $imgname.e -width 30 -relief sunken -bd 1 -textvariable cp(calp$i)
        label $imgname.l -text "Camera $i:"
		pack  $imgname.l $imgname.e -side left
		if {$i==-1} {; # can be removed
			button $imgname.br_img -text "browse all" -command "browse_image_files $w"
			pack  $imgname.br_img -side left -padx 10
		}
		pack $imgname -side left -fill x
        pack $w.$i -pady 2 -fill x -in $w.col1
    }
	button $w.br_img -text "browse all" -command "browse_image_files $w"

	frame $w.framehor
	frame $w.framever
	label $w.framehor.l -text "hor. size in pixels: "
	entry $w.framehor.e -width 6 -relief sunken -bd 1 -textvariable cp(imx)
	label $w.framever.l -text "ver. size in pixels: "
	entry $w.framever.e -width 6 -relief sunken -bd 1 -textvariable cp(imy)
	
	pack $w.framehor.l $w.framehor.e $w.framever.l $w.framever.e -side left 
	pack $w.br_img -anchor nw -pady 12 -in $w.col2
	pack $w.framehor $w.framever -in $w.col2
	pack $w.col1 $w.col2 -side left -expand 1
}

proc browse_image_files {w} {
    global cp mcp
    set types {
        {{Image files} {.tif}     }
        {{Image files} {}     TIFF}
        {{All Files}   *          }
    }
	if {$cp(multi)} {
		set title "Select the calibration images for this plane" 
	} else {
		set title "Select the calibration images" 
	}
    set pathnames [tk_getOpenFile -multiple true -filetypes $types -title $title]

    set nfiles [llength $pathnames]
    if {$nfiles > 0} {
        # check nfiles
        if {$nfiles > 4} {set $nfiles 4}
        for {set i 1} {$i<=$nfiles} {incr i} {
            # remind the difference in indices
            set pathname [lindex $pathnames [expr $i -1]]
            set cp(calp$i) $pathname
        }
    }
}

proc coordinates_on_plate { w } {
	global cp
	
	if {$cp(multi)} {
		set i $cp(plane_id)
		set cp(platecoord) $cp(platecoord$i)
	}
	
	frame $w.targetname
    entry $w.targetname.e -width 30 -relief sunken -bd 1 -textvariable cp(platecoord)
    label $w.targetname.l -text "File with coordinates:"
	pack  $w.targetname.l $w.targetname.e -side left
	
	button $w.br_crd -text "browse" -command "browse_crd_file $w"
	pack $w.targetname $w.br_crd -side left -padx 10
}

proc browse_crd_file {w} {
    global cp
    set types {
        {{Coordinate file} {.txt}}
        {{All Files}       *     }
    }
	if {$cp(multi)} {
		set title "Select the coordinate file for this plane" 
	} else {
		set title "Select the coordinate file" 
	}	
	set pathname [tk_getOpenFile -initialdir 'd:\git\testing\fan\' -filetypes $types -title $title]
    set cp(platecoord) [file pathtype $pathname]
}


proc nog_niet {} {
    if {$nfiles > 0} {
        # check nfiles
        if {$nfiles > $mcp(maxlevels)} {set $nfiles $mcp(maxlevels)}
        for {set i 1} {$i<=$nfiles} {incr i} {
            # remind the difference in indices
            set pathname [lindex $pathnames [expr $i -1]]
            set cp(imgname$ic$i) [lindex [file split $pathname] end]
        }
        set mcp(imgdir) [relative_path [file dirname $pathname]]
        # check if we have to redraw the page
        if {$nfiles != $mcp(nfiles$ic)} {
            set mcp(nfiles$ic) $nfiles
            set thispage [$nb pages $ic]
            $nb delete $thispage
            set frame [$nb insert $ic PlanesCam$ic -text "Camera $ic"]
            multi_camerapage $nb $frame $ic $nfiles
            $nb raise $thispage
            $nb compute_size
        }
    }
}


proc target_recognition { w } {

    # grey value threshold
    frame $w.partgv
    label $w.partgv.l -text "Greyvalue threshold,"
    pack $w.partgv.l -side left
    for {set i 1} {$i<=4} {incr i} {
        label $w.partgv.l$i -text " $i. Img:"
        entry $w.partgv.e$i -width 5 -relief sunken -bd 1 -textvariable mp(partgv$i)
        pack $w.partgv.l$i $w.partgv.e$i -side left
    }
    pack $w.partgv

    # particle sizes, etc
    frame $w.col1
    frame $w.col2
    frame $w.col3
    frame $w.colls
	
    frame $w.minnp
    frame $w.minnpx
    frame $w.minnpy
    frame $w.maxnp
    frame $w.maxnpx
    frame $w.maxnpy
    frame $w.ppsumgv
    frame $w.toldisc
    frame $w.ppcross

    label $w.minnp.l  -text "min npix: "
    entry $w.minnp.e  -width 5 -relief sunken -bd 1 -textvariable mp(pminnpix)
    label $w.minnpx.l -text "min npix in x: "
    entry $w.minnpx.e -width 5 -relief sunken -bd 1 -textvariable mp(pminnpixx)
    label $w.minnpy.l -text "min npix in y: "
    entry $w.minnpy.e -width 5 -relief sunken -bd 1 -textvariable mp(pminnpixy)
      
    label $w.maxnp.l  -text "max npix: "
    entry $w.maxnp.e  -width 5 -relief sunken -bd 1 -textvariable mp(pmaxnpix)
    label $w.maxnpx.l -text "max npix in x: "
    entry $w.maxnpx.e -width 5 -relief sunken -bd 1 -textvariable mp(pmaxnpixx)
    label $w.maxnpy.l -text "max npix in y: "
    entry $w.maxnpy.e -width 5 -relief sunken -bd 1 -textvariable mp(pmaxnpixy)

    label $w.ppsumgv.l -text "Sum of grayvalue: "
    entry $w.ppsumgv.e -width 5 -relief sunken -bd 1 -textvariable mp(psumgv)
    label $w.toldisc.l -text "Tolerable discontinuity: "
    entry $w.toldisc.e -width 5 -relief sunken -bd 1 -textvariable mp(partdisc)
    label $w.ppcross.l -text "Size of crosses: "
    entry $w.ppcross.e -width 5 -relief sunken -bd 1 -textvariable mp(pcrossize)

	
	# -- Packing everything ---------------------#
    pack $w.minnp.e $w.minnp.l $w.minnpx.e $w.minnpx.l $w.minnpy.e $w.minnpy.l -side right
    pack $w.maxnp.e $w.maxnp.l $w.maxnpx.e $w.maxnpx.l $w.maxnpy.e $w.maxnpy.l -side right
    pack $w.ppsumgv.e $w.ppsumgv.l $w.toldisc.e $w.toldisc.l $w.ppcross.e $w.ppcross.l -side right

    pack $w.minnp  $w.maxnp  $w.ppsumgv -pady 1 -fill x -in $w.col1 
    pack $w.minnpx $w.maxnpx $w.toldisc -pady 1 -fill x -in $w.col2
    pack $w.minnpy $w.maxnpy $w.ppcross -pady 1 -fill x -in $w.col3 
    pack $w.col1 $w.col2 $w.col3 -side left -expand 1 -in $w.colls
	pack $w.colls
}

proc pointnumbers {w} {
    # -- Point numbers for manual pre-orientation -- #
	global cp
	
	if {$cp(multi)} {
		set i $cp(plane_id)
		for {set j 1} {$j<=4} {incr j} {
			for {set k 1} {$k<=4} {incr k} {
				set cp(p$j$k) $cp(p$i$j$k)
			}
		}
		set np 4
	} else {
		set np 6
	}
	
    for {set i 1} {$i<=4} {incr i} {
		set cam $w.cam$i
		frame $cam
        label $cam.name -text "Camera $i:"
		pack $cam.name -side left
        for {set j 1} {$j<=$np} {incr j} {
			frame $cam.pp$j
            label $cam.pp$j.l -text "  P$j:"
            entry $cam.pp$j.e -width 4 -relief sunken -bd 1 -textvariable cp(p$i$j)
			pack $cam.pp$j.l $cam.pp$j.e -side left
			pack $cam.pp$j -side left
        }
		pack $cam
    }
}

proc rest {w} {
    frame  $w.cp
    label  $w.cp.l   -text "Copy the first line to the other fields:"
    button $w.cp.but -text " copy " -command "copy_ori_lines $ic $nelem"
    
    # -- packing/placing the objects
	
	pack 
	   
	
    set yi 4
    for  {set i 1} {$i<=$nelem} {incr i} {
        place $w.imglab$i -x 20  -y $yi
        set xi 45 
        for {set j 1} {$j<=4} {incr j} {
            place $w.pl{$i}$j -x $xi -y $yi
            place $w.pe{$i}$j -x [expr $xi+28] -y $yi
            set xi [expr $xi+60]
        }
        set yi [expr $yi+22]
    }
    pack $w.cp.l $w.cp.but -side left
    place $w.cp -x 20 -y $yi
    $w configure -width 287 -height [expr $yi+22] 
        
    # -- Tooltips --
}

proc degree_of_polynomials {w} {
	frame $w.col1
	frame $w.col2
	frame $w.col3
	
	frame $w.order1
	frame $w.cross1
	frame $w.order2
	frame $w.cross2
	
    label $w.col1.l1  -text "3D-points to pixels      "
    label $w.col1.l2  -text "pixels to epilpolar lines"
	label $w.order1.l -text "order:"
	label $w.cross1.l -text "order cross terms:"
	label $w.order2.l -text "order:"
	label $w.cross2.l -text "order cross terms:"
	entry $w.order1.e -width 4 -relief sunken -bd 1 -textvariable cp(order1)
	entry $w.cross1.e -width 4 -relief sunken -bd 1 -textvariable cp(order2)
	entry $w.order2.e -width 4 -relief sunken -bd 1 -textvariable cp(order3)
	entry $w.cross2.e -width 4 -relief sunken -bd 1 -textvariable cp(order4)

	# -- Packing everything ---------------------#
	pack $w.col1.l1 $w.col1.l2
	pack $w.order1.e $w.order1.l  $w.cross1.e $w.cross1.l -side right
	pack $w.order2.e $w.order2.l  $w.cross2.e $w.cross2.l -side right

	pack $w.order1 $w.order2 -pady 1 -fill x -in $w.col2
	pack $w.cross1 $w.cross2 -pady 1 -fill x -in $w.col3
	pack $w.col1 $w.col2 $w.col3 -side left -fill x -padx 6
}

proc shaking_parameters { w } {
	frame $w.col1
	frame $w.col2

    frame $w.first
	frame $w.last
    frame $w.maxpoints
	frame $w.maxframes
	
    label $w.first.l -text "first frame:"
    entry $w.first.e -width 8 -relief sunken -bd 1 -textvariable cp(first_shake)
    label $w.last.l -text "last frame:"
    entry $w.last.e -width 8 -relief sunken -bd 1 -textvariable cp(last_shake)
	
	label $w.maxpoints.l -text "max # points used for calibration:"
	entry $w.maxpoints.e -width 5 -relief sunken -bd 2 -textvariable cp(maxPoints_shake)
	label $w.maxframes.l -text "max # frames used for calibration:"
	entry $w.maxframes.e -width 5 -relief sunken -bd 2 -textvariable cp(maxFrames_shake)

	# -- Packing everything ---------------------#
	pack $w.first.e $w.first.l  $w.maxpoints.e $w.maxpoints.l -side right
	pack $w.last.e $w.last.l    $w.maxframes.e $w.maxframes.l -side right

	pack $w.first $w.last -pady 1 -fill x -in $w.col1
	pack $w.maxpoints $w.maxframes -pady 1 -fill x -in $w.col2
	pack $w.col1 $w.col2 -side left -expand 1
}

