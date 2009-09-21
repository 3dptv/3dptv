# Bindungs 
    
proc bindings0 { } {
# Bindings after imgcoord: left > measure case 5, middle > zoom normal case1, right > delete case 4

    global mp
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	
	bind .cam$n.pic <ButtonPress-1> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      set bnr [expr $nr +1];\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      set gv [ image$bnr get  $x $y];\
					      markier %W %x %y image$bnr; incr clicki 1;\
					      mouse_cmd $x $y $nr 5 [ lindex $gv 1] }
	
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1] ;\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      mouse_cmd $x $y $nr 1 }
	
	bind .cam$n.pic <ButtonPress-3> { set fn %W ; set nr [expr [ string index $fn 4 ] -1] ;\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      mouse_cmd $x $y $nr 4 }
	
	bind . <Enter> {focus %W} 
    }
    .text delete 0
    .text insert 0 "Mouse Buttons: left > measure, middle > zoom, right > delete "
    .text delete 1
    .text insert 1 " "
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " "
}



proc bindingsstart { } {

# Bindings at the beginning: left > measure case 5, middle > zoom normal case1

    global mp clicki
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	
	bind .cam$n.pic <ButtonPress-1> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      set bnr [expr $nr +1];\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      set gv [ image$bnr get  $x $y];\
					      markier %W %x %y image$bnr; incr clicki 1;\
					      mouse_cmd $x $y $nr 5 [ lindex $gv 1] }
	
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1] ;\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      mouse_cmd $x $y $nr 1 }	
	bind . <Enter> {focus %W} 
    }
    .text insert 0 "Mouse Buttons: left > measure, middle > zoom"
    .text insert 1 " "
    .text insert 2 " "
    .text insert 3 " "
}

proc bindings1 { } {
# Bindings at the beginning: left > measure case 5, middle > zoom normal case1

    global mp clicki
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	
	bind .cam$n.pic <ButtonPress-1> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      set bnr [expr $nr +1];\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      set gv [ image$bnr get  $x $y];\
					      markier %W %x %y image$bnr; incr clicki 1;\
					      mouse_cmd $x $y $nr 5 [ lindex $gv 1] }
	
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1] ;\
					      set x [expr int(%x + [%W canvasx 0 ])];\
					      set y [expr int( %y + [%W canvasy 0 ])];\
					      mouse_cmd $x $y $nr 1 }
	
	bind . <Enter> {focus %W} 

    }
    .text delete 0
    .text insert 0 "Mouse Buttons: left > measure, middle > zoom"
    .text delete 1
    .text insert 1 " "
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " "

}


proc bindings2 { } {
#Bindings after correspondences: right > show epiline case 3, middle > zoom corresp case 2
    
    global mp
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	set p image$n
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      mouse_cmd %x %y $nr 2 } 
	bind .cam$n.pic <ButtonPress-3>  { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					       mouse_cmd %x %y $nr 3 } 
    }

    .text delete 0
    .text insert 0 "Mouse Buttons:"
    .text delete 1
    .text insert 1 "left> measure, middle> zoom corresp areas, right> show epipolarline"
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " " 
}

proc bindings3 { } {
#Bindings for tracking middle > zoom case 4
    
    global mp zoomfactor
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	set p image$n
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      mouse_cmd %x %y $nr 2 $zoomfactor } 
    }

    .text delete 0
    .text insert 0 "Mouse Buttons:left> measure, middle> zoom with factor: $zoomfactor"
    .text delete 1
    .text insert 1 " "
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " " 
}

proc bindings4 { } {
#Bindings for tracking middle > zoom case 4
    
    global mp zoomfactor
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set n $i
	set p image$n
	bind .cam$n.pic <ButtonPress-2> { set fn %W ; set nr [expr [ string index $fn 4 ] -1];\
					      mouse_cmd %x %y $nr 6 $zoomfactor } 
    }

    .text delete 0
    .text insert 0 "Mouse Buttons:left> measure, middle> zoom with factor: $zoomfactor"
    .text delete 1
    .text insert 1 " "
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " " 
}
