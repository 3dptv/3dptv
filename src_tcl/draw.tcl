##########################################
# Drawing procedures    ##################
##########################################


proc markparticle { x y size imgnr color } {

    global mp
    
    set w .cam$imgnr.pic
    set xs [expr $x ]
    set ys [expr $y ]
    $w create line  [expr $xs -$size]  $ys [expr $xs +$size] $ys -width 2 -fill $color
    $w create line $xs [expr $ys -$size] $xs [expr $ys +$size] -width 2 -fill $color
}

proc drawline { x0 y0 x1 y1 b imgnr color } {
    
    set w .cam$imgnr.pic
    set xa [expr $x0]
    set ya [expr $y0]
    set xb [expr $x1]
    set yb [expr $y1]
    
    $w create line $xa $ya $xb $yb -width $b -fill $color
    
}

proc drawtext { x y pnr imgnr color } {
    
    set w .cam$imgnr.pic
    set xt [expr $x ]
    set yt [expr $y ]
    
    $w create text [expr $xt + 10 ] [expr $yt + 10 ] -text "$pnr" -font normal -fill $color   
}


proc markier { w x y image} {

    global k xs ys fx fy

    set k $w
    set xs [expr $x]
    set ys [expr $y]
    $w create line  [expr $xs -6]  $ys [expr $xs +6] $ys -width 1 -fill red 
    $w create line $xs [expr $ys -6] $xs [expr $ys +6] -width 1 -fill red

}

proc measure { p1 p2 p3 p4 imgnr } {

    global clicki xs ys fx fy k tbuf
    global px0 px1 px2 px3 py0 py1 py2 py3 

    for {set i 1} {$i < 5} {incr i 1} {

	.text delete 1
	.text insert 1 "measure points $p1 $p2 $p3 $p4 in image $imgnr"
	
	set point [set p$i]
	.text delete 3
	.text insert 3 "measure point $point >> with MS Left"

	tkwait variable clicki
	set whatimage .cam$imgnr.pic

	if { $k == $whatimage } {
	    set fx  [concat $fx $xs]
	    set fy  [concat $fy $ys]
        } else { 
	    .text delete 2
	    .text insert 2 "not image $imgnr"
	    set i [expr $i - 1 ]
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

