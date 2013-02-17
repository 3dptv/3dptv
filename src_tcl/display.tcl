##########################################
# Displaying procedures ##################
##########################################

proc camcanvas { i } {

global cp clicki
set clicki 0

# Camera  ################################
    set w .cam$i
    toplevel $w
    wm title $w "Images of Camera $i"
# Canvas for Image i #####################
    canvas $w.pic -yscrollcommand "$w.y set" -xscrollcommand "$w.x set" -cursor dotbox 
    scrollbar $w.y -command "$w.pic yview" -orient vertical 
    scrollbar $w.x -command "$w.pic xview" -orient horizontal 
    grid $w.pic -sticky news
    grid $w.y 
    grid $w.x
    
# Adjust canvassize ######################
    $w.pic configure -scrollregion [$w.pic bbox all] -height [ expr $cp(imy) ] -width [ expr $cp(imx) ]

# pack Image i ###########################
    pack $w.y -side right -fill y 
    pack $w.x -side bottom -fill x
    pack $w.pic -side right -fill both -expand yes -anchor nw

}


proc newimage { i } {

    set nr $i
    set w .cam$nr.pic
    set p image$nr

    $p copy temp -compositingrule set
    $w create image 0 0 -anchor nw -image $p
    update idletasks
}

proc refresh { } {
    
    global mp
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set nr $i
	.cam$nr.pic create image 0 0 -anchor nw -image image$nr
    }
    
}

proc clearcam { } {
    
    global mp
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set nr $i
	.cam$nr.pic delete all
	.cam$nr.pic config -bg grey
    }  
}

proc keepori { i } {
    
    global cp
    set nr $i

    ori$nr copy image$nr
}

proc showori { } {
    
    global mp
    for { set i 1} { $i <= $mp(ncam)} { incr i } {
	set nr $i
	.cam$nr.pic create image 0 0 -anchor nw -image ori$nr
    }
    
}
