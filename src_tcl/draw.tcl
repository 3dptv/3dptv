##########################################
# Drawing procedures    ##################
##########################################


proc markparticle { x y size imgnr color } {
    
    set w .cam$imgnr.pic
	image_to_viewcrd xs ys $x $y $imgnr
    $w create line  [expr $xs -$size]  $ys [expr $xs +$size] $ys -width 2 -fill $color
    $w create line $xs [expr $ys -$size] $xs [expr $ys +$size] -width 2 -fill $color
}

proc drawline { x0 y0 x1 y1 b imgnr color } {
	
    set w .cam$imgnr.pic
	image_to_viewcrd xa ya $x0 $y0 $imgnr
	image_to_viewcrd xb yb $x1 $y1 $imgnr
    $w create line $xa $ya $xb $yb -width $b -fill $color
}

proc drawtext { x y pnr imgnr color } {
  
    set w .cam$imgnr.pic
	image_to_viewcrd vwx vwy $x $y $imgnr
    $w create text [expr $vwx + 10 ] [expr $vwy + 10 ] -text "$pnr" -font normal -fill $color   
}

proc mark { w x y imgnr } {				;# called from C-functions
	image_to_viewcrd vwx vwy $x $y $imgnr
	mark_vwcrd $w $vwx $vwy $imgnr
}

proc mark_vwcrd { w x y imgnr } {
    global k xs ys fx fy

    set k $w
    set xs [expr $x]
    set ys [expr $y]
    $w create line  [expr $xs -6]  $ys [expr $xs +6] $ys -width 1 -fill red 
    $w create line $xs [expr $ys -6] $xs [expr $ys +6] -width 1 -fill red
}

proc draw_oval { x y size imgnr color } {

    set w .cam$imgnr.pic
	image_to_viewcrd vx vy $x $y $imgnr
	$w create oval [expr {$vx-$size}] [expr {$vy-$size}] \
		[expr {$vx+$size}] [expr {$vy+$size}] -width 1 -outline $color
}

proc measure { p1 p2 p3 p4 imgnr } {
    global clicki xs ys fx fy k tbuf
    global px0 px1 px2 px3 py0 py1 py2 py3 

    for {set i 1} {$i < 5} {incr i 1} {
		.text delete 1
		.text insert 1 "measure points $p1 $p2 $p3 $p4 in image $imgnr"
		
		set point [set p$i]
		.text delete 3
		.text insert 3 "measure point $point >> with MS Left,   one step back >> with MS Right"

		tkwait variable clicki

		if { $clicki == 0 } {
			.text delete 2
			.text insert 2 "one step back"
			set i [expr $i - 2]
		} else {
			set whatimage .cam$imgnr.pic
			if {$k == $whatimage} {
				# xs and ys are in view coordinates
				view_to_imagecrd imx imy $xs $ys $imgnr
				set fx  [concat $fx $imx]
				set fy  [concat $fy $imy]
			} else { 
				.text delete 2
				.text insert 2 "not image $imgnr"
				set i [expr $i - 1]
			}
		}
    }
    set px0 [lindex $fx 0] 
    set py0 [lindex $fy 0]
    set px1 [lindex $fx 1] 
    set py1 [lindex $fy 1]
    set px2 [lindex $fx 2] 
    set py2 [lindex $fy 2]
    set px3 [lindex $fx 3] 
    set py3 [lindex $fy 3]

    .text delete 3
    .text insert 3 " "
    killfxfy
}

proc killfxfy { } {
    global fx fy
    set fx {}
    set fy {}
}

proc view_to_imagecrd {imagex imagey viewx viewy iimg} {
    global zoom_f
    upvar $imagex imx $imagey imy
    
    set zf $zoom_f($iimg)
    if {$zf >= -1 && $zf <= 1} {
        set imx $viewx
        set imy $viewy;
    } elseif {$zf > 1} {        ;# image zoomed
        set imx [expr double($viewx)/$zf]
        set imy [expr double($viewy)/$zf]
    } else {                    ;# image subsampled
        set zf -$zf;
        set imx [expr double($viewx)*$zf]
        set imy [expr double($viewy)*$zf]
    }
}

proc image_to_viewcrd {viewx viewy imagex imagey iimg} {
    global zoom_f
    upvar $viewx vwx $viewy vwy
    
    set zf $zoom_f($iimg)
    if {$zf >= -1 && $zf <= 1} {
        set vwx $imagex
        set vwy $imagey
    } elseif {$zf > 1} {        ;# image zoomed
        set vwx [expr double($imagex)*$zf]
        set vwy [expr double($imagey)*$zf]
    } else {                    ;# image subsampled
        set zf -$zf;
        set vwx [expr double($imagex)/$zf]
        set vwy [expr double($imagey)/$zf]
    }
}

