################################################
# Changing of mainparameters

proc mainpar {} {

set w .changemain
catch {destroy $w}
toplevel .changemain
wm title $w "Changing Parameters"

global mp

label $w.title -text "Main Parameters" -font {Helvetica 12 bold} -width 20 -anchor center
pack $w.title -pady 3

# Imagenames ###############################################

frame $w.ncam
frame $w.fname
frame $w.fname.img
frame $w.fname.cal

frame $w.fname.img.1
frame $w.fname.cal.1
label $w.fname.img.1.l
label $w.fname.cal.1.l
label $w.ncam.label


$w.ncam.label config -text "Number of cameras:   "
entry $w.ncam.entry -width 2 -relief sunken -bd 2 -textvariable mp(ncam)
pack $w.ncam.label $w.ncam.entry -side left
pack $w.ncam -pady 2
$w.fname.img.1.l config -text "Name of 1. image:"
entry $w.fname.img.1.e -width 25 -relief sunken -bd 2 -textvariable mp(fimg1)
$w.fname.cal.1.l config -text " Calibration data for 1. camera: "
entry $w.fname.cal.1.e -width 10 -relief sunken -bd 2 -textvariable mp(camcal1)
pack $w.fname.img.1.l $w.fname.img.1.e $w.fname.cal.1.l $w.fname.cal.1.e  -side left
pack  $w.fname.img.1 $w.fname.cal.1

frame $w.fname.img.2
frame $w.fname.cal.2
label $w.fname.img.2.l
label $w.fname.cal.2.l
$w.fname.img.2.l config -text "Name of 2. image:"
entry $w.fname.img.2.e -width 25 -relief sunken -bd 2 -textvariable mp(fimg2)
$w.fname.cal.2.l config -text " Calibration data for 2. camera: "
entry $w.fname.cal.2.e -width 10 -relief sunken -bd 2 -textvariable mp(camcal2)
pack $w.fname.img.2.l $w.fname.img.2.e $w.fname.cal.2.l $w.fname.cal.2.e -side left
pack $w.fname.img.2 $w.fname.cal.2  

frame $w.fname.img.3
frame $w.fname.cal.3
label $w.fname.img.3.l
label $w.fname.cal.3.l
$w.fname.img.3.l config -text "Name of 3. image:"
entry $w.fname.img.3.e -width 25 -relief sunken -bd 2 -textvariable mp(fimg3)
$w.fname.cal.3.l config -text " Calibration data for 3. camera: "
entry $w.fname.cal.3.e -width 10 -relief sunken -bd 2 -textvariable mp(camcal3)
pack $w.fname.img.3.l $w.fname.img.3.e $w.fname.cal.3.l $w.fname.cal.3.e -side left
pack $w.fname.img.3 $w.fname.cal.3  

frame $w.fname.img.4
frame $w.fname.cal.4
label $w.fname.img.4.l
label $w.fname.cal.4.l
$w.fname.img.4.l config -text "Name of 4. image:"
entry $w.fname.img.4.e -width 25 -relief sunken -bd 2 -textvariable mp(fimg4)
$w.fname.cal.4.l config -text " Calibration data for 4. camera: "
entry $w.fname.cal.4.e -width 10 -relief sunken -bd 2 -textvariable mp(camcal4)
pack $w.fname.img.4.l $w.fname.img.4.e $w.fname.cal.4.l $w.fname.cal.4.e  -side left
pack $w.fname.img.4 $w.fname.cal.4  

pack $w.fname.img $w.fname.cal -side left
pack $w.fname -padx 10

# Imageheader Tiff, frame or field

frame $w.pic

checkbutton $w.pic.highpass -text "Highpass-Filter     " -variable mp(highpass)
checkbutton $w.pic.tiff -text "TIFF-Header     " -variable mp(tiff)
radiobutton $w.pic.fram -variable mp(type) -text "Frame   " -value 0
radiobutton $w.pic.fieldodd -variable mp(type) -text "Field odd   " -value 1
radiobutton $w.pic.fieldeven -variable mp(type) -text "Field even" -value 2

pack $w.pic.highpass $w.pic.tiff $w.pic.fram $w.pic.fieldodd $w.pic.fieldeven -side left
pack $w.pic -pady 3

# Refractive indices ###################################

frame $w.refrac

frame $w.refrac.air
frame $w.refrac.glass
frame $w.refrac.water
frame $w.refrac.thickness

label $w.refrac.t
label $w.refrac.air.l
label $w.refrac.glass.l
label $w.refrac.water.l
label $w.refrac.thickness.l

$w.refrac.t config -text "Refractive indices:" -font {Helvetica 13 bold} 
$w.refrac.air.l config -text "air: "
entry $w.refrac.air.e -width 6 -relief sunken -bd 2 -textvariable mp(air)
$w.refrac.glass.l config -text "    glass: "
entry $w.refrac.glass.e -width 6 -relief sunken -bd 2 -textvariable mp(glass)
$w.refrac.water.l config -text "    water: "
entry $w.refrac.water.e -width 6 -relief sunken -bd 2 -textvariable mp(water)
$w.refrac.thickness.l config -text "    thickness of glass (mm): "
entry $w.refrac.thickness.e -width 6 -relief sunken -bd 2 \
    -textvariable mp(thicknessglass)

pack $w.refrac.t -pady 2
pack $w.refrac.air.l $w.refrac.air.e $w.refrac.glass.l $w.refrac.glass.e \
    $w.refrac.water.l $w.refrac.water.e $w.refrac.thickness.l \
    $w.refrac.thickness.e -side left

pack $w.refrac.air $w.refrac.glass $w.refrac.water $w.refrac.thickness \
    -side left
pack $w.refrac -pady 2

# Parameters for particle recognition ##################


frame $w.particle
frame $w.partmax
frame $w.partmin
frame $w.part
frame $w.partgv
frame $w.pminnp
frame $w.pminnpx
frame $w.pminnpy
frame $w.pmaxnp
frame $w.pmaxnpx
frame $w.pmaxnpy
frame $w.m
frame $w.m.psumgv
frame $w.m.partdisc
frame $w.m.partcross

label $w.part.l
label $w.partgv.l
label $w.partgv.l1
label $w.partgv.l2
label $w.partgv.l3
label $w.partgv.l4
label $w.pminnp.l
label $w.pminnpx.l
label $w.pminnpy.l
label $w.pmaxnp.l
label $w.pmaxnpx.l
label $w.pmaxnpy.l
label $w.m.psumgv.l
label $w.m.partdisc.l
label $w.m.partcross.l


$w.part.l config -text "Parameters for particle recognition" -font {Helvetica 13 bold}

$w.partgv.l config -text "Greyvalue threshold,"
$w.partgv.l1 config -text " 1. Img:"
entry $w.partgv.e1 -width 5 -relief sunken -bd 2 -textvariable mp(partgv1)
$w.partgv.l2 config -text " 2. Img:"
entry $w.partgv.e2 -width 5 -relief sunken -bd 2 -textvariable mp(partgv2)
$w.partgv.l3 config -text " 3. Img:"
entry $w.partgv.e3 -width 5 -relief sunken -bd 2 -textvariable mp(partgv3)
$w.partgv.l4 config -text "  4. Img:"
entry $w.partgv.e4 -width 5 -relief sunken -bd 2 -textvariable mp(partgv4)

$w.pminnp.l config -text "               min npix: "
entry $w.pminnp.e -width 5 -relief sunken -bd 2 -textvariable mp(pminnpix)
$w.pminnpx.l config -text "                  min npix in x: "
entry $w.pminnpx.e -width 5 -relief sunken -bd 2 -textvariable mp(pminnpixx)
$w.pminnpy.l config -text "       min npix in y: "
entry $w.pminnpy.e -width 5 -relief sunken -bd 2 -textvariable mp(pminnpixy)
$w.pmaxnp.l config -text "              max npix: "
entry $w.pmaxnp.e -width 5 -relief sunken -bd 2 -textvariable mp(pmaxnpix)
$w.pmaxnpx.l config -text "                 max npix in x: "
entry $w.pmaxnpx.e -width 5 -relief sunken -bd 2 -textvariable mp(pmaxnpixx)
$w.pmaxnpy.l config -text "      max npix in y: "
entry $w.pmaxnpy.e -width 5 -relief sunken -bd 2 -textvariable mp(pmaxnpixy)



pack $w.part.l -pady 5
pack $w.part -in $w.particle
pack $w.partgv.l $w.partgv.l1 $w.partgv.e1 $w.partgv.l2 $w.partgv.e2 $w.partgv.l3 \
    $w.partgv.e3  $w.partgv.l4 $w.partgv.e4  -side left
pack $w.partgv  -pady 3 -side top -in $w.particle
pack $w.particle -pady 5 
pack $w.pminnp.l $w.pminnp.e $w.pminnpx.l $w.pminnpx.e \
    $w.pminnpy.l $w.pminnpy.e -side left
pack $w.pminnp $w.pminnpx $w.pminnpy -side left -in $w.partmin
pack $w.partmin
pack $w.pmaxnp.l $w.pmaxnp.e $w.pmaxnpx.l $w.pmaxnpx.e \
    $w.pmaxnpy.l $w.pmaxnpy.e -side left
pack $w.pmaxnp $w.pmaxnpx $w.pmaxnpy -side left -in $w.partmax
pack $w.partmax

$w.m.partdisc.l config -text "  Tolerable discontinuity: "
entry $w.m.partdisc.e -width 5 -relief sunken -bd 2 -textvariable mp(partdisc)
$w.m.psumgv.l config -text "Sum of greyvalue: "
entry $w.m.psumgv.e -width 5 -relief sunken -bd 2 -textvariable mp(psumgv)
$w.m.partcross.l config -text "  Size of crosses: "
entry $w.m.partcross.e -width 5 -relief sunken -bd 2 -textvariable mp(pcrossize)

pack  $w.m.psumgv.l $w.m.psumgv.e $w.m.partdisc.l $w.m.partdisc.e $w.m.partcross.l $w.m.partcross.e -side left
pack  $w.m.psumgv $w.m.partdisc $w.m.partcross  -side left -in $w.m
pack $w.m -pady 4


# Parameters for sequence processing #########################

frame $w.seqp
frame $w.seqp.f
frame $w.seqp.l

frame $w.basename1
frame $w.basename2
frame $w.basename3
frame $w.basename4

label $w.basename1.l
label $w.basename2.l
label $w.basename3.l
label $w.basename4.l

$w.basename1.l config -text "Basename for 1. sequence: "
entry $w.basename1.e -width 25 -relief sunken -bd 2 -textvariable mp(basename1)
$w.basename2.l config -text "Basename for 2. sequence: "
entry $w.basename2.e -width 25 -relief sunken -bd 2 -textvariable mp(basename2)
$w.basename3.l config -text "Basename for 3. sequence: "
entry $w.basename3.e -width 25 -relief sunken -bd 2 -textvariable mp(basename3)
$w.basename4.l config -text "Basename for 4. sequence: "
entry $w.basename4.e -width 25 -relief sunken -bd 2 -textvariable mp(basename4)
label $w.seqp.text
label $w.seqp.f.l
label $w.seqp.l.l
$w.seqp.text config -text "Parameters for sequence processing" -font {Helvetica 13 bold}
$w.seqp.f.l config -text "Sequence images,      First:  "
entry $w.seqp.f.e -width 8 -relief sunken -bd 2 -textvariable mp(seqfirst)
$w.seqp.l.l config -text "     Last: "
entry $w.seqp.l.e -width 8 -relief sunken -bd 2 -textvariable mp(seqlast)

pack $w.seqp.text -pady 4
pack  $w.seqp.f.l $w.seqp.f.e $w.seqp.l.l $w.seqp.l.e -side left -fill x
pack $w.seqp.f $w.seqp.l -in $w.seqp -side left  -fill x
pack $w.seqp -pady 1

pack $w.basename1.l $w.basename1.e -side left
pack $w.basename1
pack $w.basename2.l $w.basename2.e -side left
pack $w.basename2
pack $w.basename3.l $w.basename3.e -side left
pack $w.basename3
pack $w.basename4.l $w.basename4.e -side left
pack $w.basename4

# Illuminated layer data ############################

frame $w.illu
frame $w.illu.xmin
frame $w.illu.xmax
frame $w.illu.crit
frame $w.illu.crit.a
frame $w.illu.crit.b
frame $w.illu.crit.c
frame $w.illu.crit.d

label $w.illu.t
label $w.illu.xmin.lx
label $w.illu.xmin.lzmin
label $w.illu.xmin.lzmax
label $w.illu.xmax.lx
label $w.illu.xmax.lzmin
label $w.illu.xmax.lzmax
label $w.illu.crit.lt
label $w.illu.crit.lnx
label $w.illu.crit.lny
label $w.illu.crit.lnpix
label $w.illu.crit.lsumgv
label $w.illu.crit.lminwcorr
label $w.illu.crit.ltolepi

$w.illu.t config -text "Illuminated layer data" -font {Helvetica 13 bold}
$w.illu.xmin.lx config -text "Xmin:  "
entry $w.illu.xmin.ex -width 6 -relief sunken -bd 2 -textvariable mp(xmin)
$w.illu.xmin.lzmin config -text "   Zmin: "
entry $w.illu.xmin.ezmin -width 6 -relief sunken -bd 2 \
    -textvariable mp(xminzmin)
$w.illu.xmin.lzmax config -text "   Zmax: "
entry $w.illu.xmin.ezmax -width 6 -relief sunken -bd 2 \
    -textvariable mp(xminzmax)
$w.illu.xmax.lx config -text "Xmax: "
entry $w.illu.xmax.ex -width 6 -relief sunken -bd 2 -textvariable mp(xmax)
$w.illu.xmax.lzmin config -text "   Zmin: "
entry $w.illu.xmax.ezmin -width 6 -relief sunken -bd 2 \
    -textvariable mp(xmaxzmin)
$w.illu.xmax.lzmax config -text "   Zmax: "
entry $w.illu.xmax.ezmax -width 6 -relief sunken -bd 2 \
    -textvariable mp(xmaxzmax)
$w.illu.crit.lt config -text "Criteria for correspondences" -font {Helvetica 13 bold}
$w.illu.crit.lnx config -text "             min corr for ratio nx: "
entry $w.illu.crit.enx -width 6 -relief sunken -bd 2 -textvariable mp(nx)
$w.illu.crit.lny config -text "                           min corr for ratio ny: "
entry $w.illu.crit.eny -width 6 -relief sunken -bd 2 -textvariable mp(ny)
$w.illu.crit.lnpix config -text "          min corr for ratio npix: "
entry $w.illu.crit.enpix -width 6 -relief sunken -bd 2 -textvariable mp(npix)
$w.illu.crit.lsumgv config -text "                                         sum of gv: "
entry $w.illu.crit.esumgv -width 6 -relief sunken -bd 2 -textvariable mp(sgv)
$w.illu.crit.lminwcorr config -text "min for weighted correlation: "
entry $w.illu.crit.eminwcorr -width 6 -relief sunken -bd 2 \
    -textvariable mp(mincorr)
$w.illu.crit.ltolepi config -text "      Tolerance to epipolar band (mm): "
entry $w.illu.crit.etolepi -width 6 -relief sunken -bd 2 \
    -textvariable mp(tolepi)

pack $w.illu.t -pady 4

pack $w.illu.xmin.lx $w.illu.xmin.ex $w.illu.xmin.lzmin $w.illu.xmin.ezmin \
    $w.illu.xmin.lzmax $w.illu.xmin.ezmax -side left
pack $w.illu.xmin -side top -pady 2
pack $w.illu.xmax.lx $w.illu.xmax.ex $w.illu.xmax.lzmin $w.illu.xmax.ezmin \
    $w.illu.xmax.lzmax $w.illu.xmax.ezmax -side left
pack $w.illu.xmax -side top -pady 2

pack $w.illu.crit.lt -pady 2
pack $w.illu.crit.lnx $w.illu.crit.enx $w.illu.crit.lny $w.illu.crit.eny -side left -in $w.illu.crit.a
pack $w.illu.crit.lnpix $w.illu.crit.enpix $w.illu.crit.lsumgv $w.illu.crit.esumgv -side left -in $w.illu.crit.b


pack $w.illu.crit.lminwcorr $w.illu.crit.eminwcorr $w.illu.crit.ltolepi $w.illu.crit.etolepi \
 -side left  -in $w.illu.crit.c

pack $w.illu.crit.a  $w.illu.crit.b $w.illu.crit.c -side top -in $w.illu.crit
pack $w.illu.crit -pady 6

pack $w.illu -pady 6

button $w.ok -text OK -command "done_proc_cmd;destroy $w"
pack $w.ok -side bottom -pady 5
}

