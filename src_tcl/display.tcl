##########################################
# Displaying procedures ##################
##########################################

proc camcanvas { i } {
	global cp clicki zoompar
	set clicki 0

	# -- Camera -- #
    set w .cam$i
    toplevel $w
    wm title $w "Images of Camera $i"
	# -- Canvas for Image i -- #
    canvas $w.pic -yscrollcommand "$w.y set" -xscrollcommand "$w.x set" -cursor dotbox 
    scrollbar $w.y -command "$w.pic yview" -orient vertical 
    scrollbar $w.x -command "$w.pic xview" -orient horizontal 
    grid $w.pic -sticky news
    grid $w.y 
    grid $w.x
    
	# -- Adjust canvassize -- #
    $w.pic configure -scrollregion [$w.pic bbox all] -height [expr $cp(imy)] -width [expr $cp(imx)]

	# -- pack Image i -- #
    pack $w.y -side right -fill y 
    pack $w.x -side bottom -fill x
    pack $w.pic -side right -fill both -expand yes -anchor nw
}

proc refresh { } {
    global mp
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
		.cam$i.pic create image 0 0 -anchor nw -image image$i
    }
    
}

proc clearcam { } {
    global mp
    clear_cmd		;# call the clearmarkers_c function first
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
	    .cam$i.pic delete all
	    .cam$i.pic config -bg grey
    }
}

proc keepori { i } {
    global cp
    ori$i copy image$i
}

proc showori { } {
    global mp
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
        .cam$i.pic create image 0 0 -anchor nw -image ori$i
    }
}

# --- functions added/changed to improve zooming of images --- #
#     ad holten, March 2013

# proc newimage { i } {
#     set w .cam$i.pic
#     set p image$i
# 
#     $p copy temp -compositingrule set
#     $w create image 0 0 -anchor nw -image $p
#     update idletasks
# }

proc newimage {i xpos ypos zfac lockframesize} {
    # newimage function with zoom/subsample functionality
    # zfac > 1  = zoom out
    # zfac < -1 = zoom in   (subsample)
    # xpos and ypos are the new centre position of the view in fractions of the image
    global zoom_f zoompar zoom

    set w .cam$i.pic
    set p image$i
    set zoom_f($i) $zfac
	set zoompar(fixed) $lockframesize
	set zoom(fixed$i)    $lockframesize
    # delete the old image
    $w delete all
    $w config -bg grey
    
    # -- create the zoomed image --
    image create photo $p
    if {$zfac >= -1 && $zfac <= 1} {	;# not zooming or bad zfac
        $p copy temp -compositingrule set
    } elseif {$zfac > 1} {              ;# zooming in
        $p copy temp -zoom $zfac -compositingrule set
    } elseif {$zfac < -1} {             ;# zooming out
        $p copy temp -subsample [expr -$zfac] -compositingrule set
    }
    $w create image 0 0 -anchor nw -image $p

    # -- get the current frame dimensions --
	scan [wm geometry .cam$i] "%dx%d+%d%d" frm_width frm_height frm_left frm_top
	if {$lockframesize == 0} {
		# -- adjust the frame size to the new image --
		# but keep the frame visible inside the display
		if {$frm_left < 0} {set frm_left 0}
		if {$frm_top  < 0} {set frm_top 0}
		set nw_width  [expr [image width  $p]+25]
		set nw_height [expr [image height $p]+25]
		wm geometry .cam$i ${nw_width}x${nw_height}+$frm_left+$frm_top
	}
    update idletasks

    # -- configure the scroll region and move the scroll region to the clicked position --
    $w configure -scrollregion [$w bbox all]
    set displ_width  [lindex [$w bbox all] 2]
	set displ_height [lindex [$w bbox all] 3]
    set frame_width  [winfo width .cam$i]		;# also scan [wm geometry .cam$i] ... will work
    set frame_height [winfo height .cam$i]
    if {$displ_width > $frame_width || $displ_height > $frame_height} {
        set dx [expr double($frame_width)/$displ_width]
        set dy [expr double($frame_height)/$displ_height]
        $w xview moveto [expr $xpos - $dx/2]
        $w yview moveto [expr $ypos - $dy/2]
    }
}

proc panimage {i xpos ypos} {
    # move the centre of the image to the position xpos, ypos
    set w .cam$i.pic

	set x1 [lindex [$w xview] 0]	; # left side position (fraction) of the view
	set x2 [lindex [$w xview] 1]	; # right side position (fraction) of the view
	set y1 [lindex [$w yview] 0]	; # bottom side position (fraction) of the view
	set y2 [lindex [$w yview] 1]	; # top side position (fraction) of the view
    $w xview moveto [expr $xpos - ($x2 - $x1)/2]
    $w yview moveto [expr $ypos - ($y2 - $y1)/2]
}

proc new_framesizes {fac} {
    global mp cp zoompar
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
		# -- get the current frame dimensions --
		scan [wm geometry .cam$i] "%dx%d+%d%d" frm_width frm_height frm_left frm_top
		if {$fac < 0} {
			set nw_width  [expr -$cp(imx)/$fac+25]
			set nw_height [expr -$cp(imy)/$fac+25]
		} else {
			set nw_width  [expr $cp(imx)*$fac +25]
			set nw_height [expr $cp(imy)*$fac +25]
		}
		wm geometry .cam$i ${nw_width}x${nw_height}+$frm_left+$frm_top
	}
	update idletasks
	set zoompar(fixed) 1
}

proc get_viewsize {i} {
    # get the current status of the view
	global zoom zoompar zoom_f;
    set w .cam$i.pic

	set x1 [lindex [$w xview] 0]	; # left side position (fraction) of the view
	set x2 [lindex [$w xview] 1]	; # right side position (fraction) of the view
	set y1 [lindex [$w yview] 0]	; # bottom side position (fraction) of the view
	set y2 [lindex [$w yview] 1]	; # top side position (fraction) of the view
	set imx [image width image$i];
	set imy [image height image$i];

	set zoom(xc$i)  [expr ($x2 + $x1)/2]; 
	set zoom(yc$i)  [expr ($y2 + $y1)/2];
	set zoom(vwx$i) [expr $imx*($x2 - $x1)];
	set zoom(vwy$i) [expr $imy*($y2 - $y1)];
	set zoom(fac$i) $zoom_f($i)
	set zoom(fixed$i) $zoompar(fixed)
}

# --- create the popup menus ---
proc zoom_to {fac} {
    global zoompar
    mouse_cmd 1 $zoompar(x) $zoompar(y) $zoompar(nr) $fac $zoompar(fixed)
}

proc create_popup_zoomed {} {			;# used in the function popup_zoommenu in button.tcl
    set m [menu .popupMenu -tearoff 0]

	$m add cascade -label "Zoom All..."   -menu .mbar.options.menu.zoommenu -underline 0
    $m add cascade -label "Frame size..." -menu .mbar.options.menu.framesize -underline 0
	$m add separator
	add_zoomentries $m zoompar(fac) "zoom_to \$zoompar\(fac\)"
}

