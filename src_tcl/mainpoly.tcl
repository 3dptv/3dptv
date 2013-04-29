#
# Dialog for changing the main parameters in the case 
# of using polynomials for particle mapping
#

package require BWidget 1.9.1

proc mainpoly {} {

    set w .changemain
    catch {destroy $w}
    toplevel .changemain
    wm title $w "Main Parameters"

    global mp pf
    set bold12 {Helvetica 12 bold}
    set bold10 {Helvetica 10 bold}


    label $w.title -text "Main Parameters" -font $bold12 -width 20 -anchor center
    pack $w.title -pady 3

    # ------ Number of cameras -----------------------------------------
    frame $w.numcam
    set ncam $w.numcam.ncam
    set all4 $w.numcam.all4

    frame $ncam
    label $ncam.label -text "Number of cameras:   "
    entry $ncam.entry -width 2 -relief sunken -bd 2 -textvariable mp(ncam)
    pack $ncam.label $ncam.entry -side left

    frame $all4
    checkbutton $all4.cb -text "accept only points seen from all cameras?" -variable mp(allCam)
    pack $all4.cb -side left

    pack $ncam $all4 -padx 10 -side left
    pack $w.numcam -pady {2 4}

    # ----- Images to be opened at start, and calibration files ------------
    # set s1 "Starting Images to be opened on start, and calibration files"
    set s1 "Starting images and calibration files"
    set tf1 [TitleFrame $w.tf1 -text $s1 -font $bold10]
    
    images_and_calibration_filenames [$tf1 getframe]
    pack $tf1 -fill x -padx 10 -pady 3

    # -------- Parameters for particle recognition -----------------------------------
    set s2 "Parameters for particle recognition"
    set tf2 [TitleFrame $w.tf2 -text $s2 -font $bold10]
    
    particle_recognition [$tf2 getframe]
    pack $tf2 -fill x -padx 10 -pady 3
    
    # -------- Parameters for sequence processing ------------------------------------
    set s3 "Parameters for sequence processing"
    set tf3 [TitleFrame $w.tf3 -text $s3 -font $bold10]

    sequence_processing [$tf3 getframe]
    pack $tf3 -fill x -padx 10 -pady 3
    
    # -------- Illuminated layer data ------------------------------------------------
    set s4 "Illuminated volume data"
    set tf4 [TitleFrame $w.tf4 -text $s4 -font $bold10]

    illum_volume [$tf4 getframe]
    pack $tf4 -fill x -padx 10 -pady 3


    # -------- Criteria for correspondences-------------------------------------------
    set s5 "Criteria for correspondences"
    set tf5 [TitleFrame $w.tf5 -text $s5 -font $bold10]

    crit_correspondences [$tf5 getframe]
    pack $tf5 -fill x -padx 10 -pady 3
    
    button $w.ok -text OK -command "done_proc_cmd; destroy $w"
    pack $w.ok -side bottom -pady 5
}


proc images_and_calibration_filenames { w } {
    for {set i 1} {$i<=4} {incr i} {
        frame $w.$i
        set imgname $w.$i.imgname
        set calname $w.$i.calname
        
        frame $imgname
        entry $imgname.en -width 30 -relief sunken -bd 2 -textvariable mp(fimg$i)
        label $imgname.lb -text "Camera $i:  image:"

        frame $calname
        entry $calname.en -width 18 -relief sunken -bd 2 -textvariable mp(camcal$i)
        label $calname.lb -text " calib. file:"
                
        pack  $imgname.lb $imgname.en -side left
        pack  $calname.lb $calname.en -side left
        pack  $imgname $calname -side left
 
        pack  $w.$i -pady 2
    }
	checkbutton $w.highpass -text "Highpass-Filter" -variable mp(highpass)
	pack  $w.highpass -side left
}

proc particle_recognition { w } {

    # grey value threshold
    frame $w.partgv
    label $w.partgv.l -text "Greyvalue threshold,"
    pack $w.partgv.l -side left
    for {set i 1} {$i<=4} {incr i} {
        label $w.partgv.l$i -text " $i. Img:"
        entry $w.partgv.e$i -width 5 -relief sunken -bd 2 -textvariable mp(partgv$i)
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

	# -- use relative tolerable discontinuity -- #
	frame $w.rel_disc
	checkbutton $w.rel_disc.button -text "use relative tol. discontinuity (%)?" -variable mp(rel_disc)
	
	# -- read from target files -----------------#
	frame $w.target
	checkbutton $w.target.button -text "use existing _target files?" -variable mp(target)

	# -- Subtract Mask --------------------------#
	frame $w.mask
	label $w.mask.label
	checkbutton $w.mask.button -text "Subtract mask     " -variable mp(mask)
	$w.mask.label config -text "Basename for the first mask:"
	entry $w.mask.entry -width 25 -relief sunken -bd 2 -textvariable mp(maskname)
	
	# -- Packing everything ---------------------#
    pack $w.minnp.e $w.minnp.l $w.minnpx.e $w.minnpx.l $w.minnpy.e $w.minnpy.l -side right
    pack $w.maxnp.e $w.maxnp.l $w.maxnpx.e $w.maxnpx.l $w.maxnpy.e $w.maxnpy.l -side right
    pack $w.ppsumgv.e $w.ppsumgv.l $w.toldisc.e $w.toldisc.l $w.ppcross.e $w.ppcross.l -side right

    pack $w.minnp  $w.maxnp  $w.ppsumgv -pady 1 -fill x -in $w.col1 
    pack $w.minnpx $w.maxnpx $w.toldisc -pady 1 -fill x -in $w.col2
    pack $w.minnpy $w.maxnpy $w.ppcross -pady 1 -fill x -in $w.col3 
    pack $w.col1 $w.col2 $w.col3 -side left -expand 1 -in $w.colls
	pack $w.colls
	
	pack $w.rel_disc.button -side left
	pack $w.rel_disc

	pack $w.mask.button $w.mask.label $w.mask.entry -side left
	pack $w.mask -pady 4

	pack $w.target.button -side left
	pack $w.target -fill x
}

proc sequence_processing { w } {
    frame $w.seqp
    frame $w.seqp.f
    frame $w.seqp.l

    label $w.seqp.f.l -text "Sequence images,    First:  "
    entry $w.seqp.f.e -width 8 -relief sunken -bd 2 -textvariable mp(seqfirst)
    label $w.seqp.l.l -text "     Last: "
    entry $w.seqp.l.e -width 8 -relief sunken -bd 2 -textvariable mp(seqlast)
    
    pack  $w.seqp.f.l $w.seqp.f.e $w.seqp.l.l $w.seqp.l.e -side left -fill x
    pack  $w.seqp.f $w.seqp.l -in $w.seqp -side left  -fill x
    pack  $w.seqp -pady 1

    for {set i 1} {$i<=4} {incr i} {
        frame $w.basename$i
        label $w.basename$i.l -text "Basename for $i. sequence: "
        entry $w.basename$i.e -width 25 -relief sunken -bd 2 \
            -textvariable mp(basename$i)
    
        pack $w.basename$i.l $w.basename$i.e -side left
        pack $w.basename$i
    }
}

proc illum_volume { w } {
    frame  $w.illu
    set mn $w.illu.min
    set mx $w.illu.max
    frame $mn
    frame $mx

    frame $mn.xmin
    frame $mn.ymin
    frame $mn.zmin
    frame $mx.xmax
    frame $mx.ymax
    frame $mx.zmax

    # $w.illu.title config -text "Illuminated volume data" -font {Helvetica 13 bold}
    label $mn.xmin.lb -text "Xmin:  "
    label $mn.ymin.lb -text "   Ymin:  "
    label $mn.zmin.lb -text "   Zmin:  "
    label $mx.xmax.lb -text "Xmax: "
    label $mx.ymax.lb -text "   Ymax: "
    label $mx.zmax.lb -text "   Zmax: "

    entry $mn.xmin.en -width 6 -relief sunken -bd 2 -textvariable mp(xmin)
    entry $mn.ymin.en -width 6 -relief sunken -bd 2 -textvariable mp(ymin)
    entry $mn.zmin.en -width 6 -relief sunken -bd 2 -textvariable mp(zmin)
    entry $mx.xmax.en -width 6 -relief sunken -bd 2 -textvariable mp(xmax)
    entry $mx.ymax.en -width 6 -relief sunken -bd 2 -textvariable mp(ymax)
    entry $mx.zmax.en -width 6 -relief sunken -bd 2 -textvariable mp(zmax)

    pack $mn.xmin.lb $mn.xmin.en $mn.ymin.lb $mn.ymin.en $mn.zmin.lb $mn.zmin.en -side left
    pack $mx.xmax.lb $mx.xmax.en $mx.ymax.lb $mx.ymax.en $mx.zmax.lb $mx.zmax.en -side left

    pack $mn.xmin $mn.ymin $mn.zmin -side left
    pack $mx.xmax $mx.ymax $mx.zmax -side left
    pack $mn $mx -side top -pady 1

    pack $w.illu -pady 4
}

proc crit_correspondences { w } {
    frame $w.col1
    frame $w.col2
	
    frame $w.minnx
    frame $w.minny
    frame $w.minnp
    frame $w.mingv
    frame $w.minwcor
    frame $w.tolepi 
	
    label $w.minnx.l -text "min ratio nx: "
    entry $w.minnx.e -width 6 -relief sunken -bd 2 -textvariable mp(nx)
    label $w.minny.l -text "min ratio ny: "
    entry $w.minny.e -width 6 -relief sunken -bd 2 -textvariable mp(ny)

    label $w.minnp.l -text "min ratio npix: "
    entry $w.minnp.e -width 6 -relief sunken -bd 2 -textvariable mp(npix)
    label $w.mingv.l -text "min ratio grayvalues: "
    entry $w.mingv.e -width 6 -relief sunken -bd 2 -textvariable mp(sgv)

    label $w.minwcor.l -text "min for we. corr (depreciated): "
    entry $w.minwcor.e -width 6 -relief sunken -bd 2 -textvariable mp(mincorr)
    label $w.tolepi.l -text "Tol. to epipolar (frac. of particle size): "
    entry $w.tolepi.e -width 6 -relief sunken -bd 2 -textvariable mp(tolepi)

	# -- Packing everything --------------------- #
    pack $w.minnx.e $w.minnx.l $w.minny.e $w.minny.l -side right
    pack $w.minnp.e $w.minnp.l $w.mingv.e $w.mingv.l -side right
    pack $w.minwcor.e $w.minwcor.l $w.tolepi.e $w.tolepi.l -side right

    pack $w.minnx $w.minnp $w.minwcor -pady 1 -fill x -in $w.col1 
    pack $w.minny $w.mingv $w.tolepi  -pady 1 -fill x -in $w.col2
    pack $w.col1 $w.col2 -side left -padx 4
}

