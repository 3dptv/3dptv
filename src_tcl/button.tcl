# Bindungs 
#    

# -- helper function -- #
proc do_measure {w clickx clicky} {		;# Demands for using Tcl_GlobalEval(interp, val)
	global clicki						;# in mouse_proc_c() (mousefunctions.c)

	set nr  [string index $w 4]
	set x   [expr int($clickx + [$w canvasx 0])]
	set y   [expr int($clicky + [$w canvasy 0])]
	set gv  [image$nr get $x $y]
	mark_vwcrd $w $x $y image$nr
	view_to_imagecrd imx imy $x $y $nr
	incr clicki 1
	mouse_cmd 5 $imx $imy $nr [lindex $gv 1]
}

proc popup_zoommenu {clickx clicky w xpos ypos} {
	# first set the zoompar values before opening the popup menu.
    global zoom_f zoompar
	set nr [string index $w 4]
	set x  [expr int($clickx + [$w canvasx 0])]
	set y  [expr int($clicky + [$w canvasy 0])]
	view_to_imagecrd zoompar(x) zoompar(y) $x $y $nr
	set zoompar(w)  $w
	set zoompar(nr) $nr
	set zoompar(fac)  $zoom_f($nr)
	tk_popup .popupMenu $xpos $ypos
}
  
proc bindings0 { } {
# Bindings after imgcoord:  left   > measure,     case 5
#                           middle > zoom normal, case 1
#                           right  > delete,      case 4
    global mp
    if {! [winfo exists .popupMenu]} { create_popup_zoomed }
	
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
        bind .cam$i.pic <ButtonPress-1> {		;# measure
			do_measure %W %x %y
        } 
        bind .cam$i.pic <ButtonPress-2> {
			popup_zoommenu %x %y %W %X %Y
        } 
		bind .cam$i.pic <ButtonPress-3> {
			set nr [string index %W 4]
			set x  [expr int(%x + [%W canvasx 0])]
			set y  [expr int(%y + [%W canvasy 0])]
			view_to_imagecrd imx imy $x $y $nr
			mouse_cmd 4 $imx $imy $nr
		}
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
# Bindings at the beginning:  left   > measure,     case 5
#                             middle > zoom normal, case 1
    global mp clicki
    if {! [winfo exists .popupMenu]} { create_popup_zoomed }
	
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
        bind .cam$i.pic <ButtonPress-1> {		;# measure
			do_measure %W %x %y
        } 
		bind .cam$i.pic <ButtonPress-2> {
			popup_zoommenu %x %y %W %X %Y
		}	
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

proc bindings1 { } {		;# same as bindingsstart{} till now
# Bindings at the beginning: left   > measure      case 5	(start of calibration)
#                            middle > zoom normal  case1
    global mp clicki zoompar zoom_f
    if {! [winfo exists .popupMenu]} { create_popup_zoomed }
	
    for {set i 1} {$i <= $mp(ncam)} {incr i} {
        bind .cam$i.pic <ButtonPress-1> {		;# measure
			do_measure %W %x %y
        } 
		bind .cam$i.pic <ButtonPress-2> {
			popup_zoommenu %x %y %W %X %Y
        } 
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
#Bindings after correspondences: left   > center corresp areas      case 2
#                                middle > popup zoom menu           case 5
#                                right  > center and show epiline   case 3
    global mp clicki zoompar zoom_f
	if {! [winfo exists .popupMenu]} { create_popup_zoomed }

    for {set i 1} {$i <= $mp(ncam)} {incr i} {
		set p image$i
		bind .cam$i.pic <ButtonPress-1> {
			set nr  [string index %W 4]
			set x   [expr int(%x + [%W canvasx 0])]
			set y   [expr int(%y + [%W canvasy 0])]
			view_to_imagecrd imx imy $x $y $nr
			mouse_cmd 2 $imx $imy $nr $zoompar(fixed) $clicki
		}
		bind .cam$i.pic <ButtonPress-2> {
			popup_zoommenu %x %y %W %X %Y
        } 
		bind .cam$i.pic <ButtonPress-3> {
			set nr [string index %W 4]
			set x  [expr int(%x + [%W canvasx 0])]
			set y  [expr int(%y + [%W canvasy 0])]
			view_to_imagecrd imx imy $x $y $nr
			mouse_cmd 3 $imx $imy $nr
		}
    }
    .text delete 0
    .text insert 0 "Mouse Buttons: le.> centre corr. areas, mi.> zoom, ri.> show epipolar lines"
    .text delete 1  
	.text insert 1 " "
    .text delete 2
    .text insert 2 " "
    .text delete 3
    .text insert 3 " " 
}

# -- This binding function is not used anymore, ad holten, 04-2013 --
# proc bindings3 { } {
# #Bindings for tracking middle > zoom case 4
#     global mp zoomfactor
# 	
#     for {set i 1} {$i <= $mp(ncam)} {incr i} {
# 		set n $i
# 		set p image$n
# 		bind .cam$n.pic <ButtonPress-2> {
# 			set nr [string index %W 4]
# 			mouse_cmd 2 %x %y $nr $zoomfactor
# 		} 
# 	}
#     .text delete 0
#     .text insert 0 "Mouse Buttons:left> measure, middle> zoom with factor: $zoomfactor"
#     .text delete 1
#     .text insert 1 " "
#     .text delete 2
#     .text insert 2 " "
#     .text delete 3
#     .text insert 3 " " 
# }

# -- This binding function is not used anymore, ad holten, 04-2013 --
# proc bindings4 { } {
# #Bindings for tracking middle > zoom case 4
#     
#     global mp zoomfactor zoompar
#     for {set i 1} {$i <= $mp(ncam)} {incr i} {
# 		set n $i
# 		set p image$n
# 		bind .cam$n.pic <ButtonPress-2> {
# 			set nr [string index %W 4]
# 			mouse_cmd 6 %x %y $nr $zoomfactor
# 		} 
#     }
#     .text delete 0
#     .text insert 0 "Mouse Buttons:left> measure, middle> zoom with factor: $zoomfactor"
#     .text delete 1
#     .text insert 1 " "
#     .text delete 2
#     .text insert 2 " "
#     .text delete 3
#     .text insert 3 " " 
# }

# -- This binding function is not used anymore, ad holten, 04-2013 --
# proc bindings5 { } {
# # Bindings for manual detection: left > measure case 5, middle > zoom normal case1, right > delete measeured point case 7
# 
#     global mp
#     for {set i 1} {$i <= $mp(ncam)} {incr i} {
# 		set n $i
# 		set clicki 1
# 		bind .cam$n.pic <ButtonPress-1> {
# 			set nr  [string index %W 4]
# 			set x   [expr int(%x + [%W canvasx 0])]
# 			set y   [expr int(%y + [%W canvasy 0])]
# 			set gv  [image$nr get  $x $y]
# 			markier %W %x %y image$nr
# 			incr clicki 1
# 			mouse_cmd 5 $x $y $nr [lindex $gv 1]
# 		}
# 		bind .cam$n.pic <ButtonPress-2> {
# 			popup_zoommenu %x %y %W %X %Y
# #            set nr [expr [string index %W 4] -1]
# #            set x  [expr int(%x + [%W canvasx 0])]
# #            set y  [expr int(%y + [%W canvasy 0])]
# #            view_to_imagecrd zoompar(x) zoompar(y) $x $y [expr $nr +1]
# #            set zoompar(w)  %W
# #            set zoompar(nr) $nr
# #            tk_popup .popupMenu %X %Y
# 		}
# 		bind .cam$n.pic <ButtonPress-3> {
# 			set nr [string index %W 4]
# 			set x  [expr int(%x + [%W canvasx 0])]
# 			set y  [expr int(%y + [%W canvasy 0])]
# 			set clicki 0
# 			mouse_cmd 7 $x $y $nr
# 		}
# 		bind . <Enter> {focus %W} 
#     }
#     .text delete 0
#     .text insert 0 "Mouse Buttons: left > measure, middle > zoom, right > one step back"
#     .text delete 1
#     .text insert 1 " "
#     .text delete 2
#     .text insert 2 " "
#     .text delete 3
#     .text insert 3 " "
# }

