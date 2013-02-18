# Changing of Calibration Parameter
########################################

proc calpar { } {
set w .changecalib
catch {destroy $w}
toplevel .changecalib
wm title $w "Changing Parameters"

global cp

label $w.title -text "Calibration Parameters"  -font {Helvetica 12 bold} -width 20 -anchor center

# Pictures for calibration and platepointcoordinatefile

frame $w.c1
frame $w.c2
frame $w.c3
frame $w.c4
frame $w.calp1
frame $w.cori1
frame $w.calp2
frame $w.cori2
frame $w.calp3
frame $w.cori3
frame $w.calp4
frame $w.cori4
frame $w.plate
label $w.plate.l
label $w.calp1.l
label $w.calp2.l
label $w.calp3.l
label $w.calp4.l
label $w.cori1.l
label $w.cori2.l
label $w.cori3.l
label $w.cori4.l

frame $w.pic


$w.plate.l config -text "                                                                        File of Coordinates on Plate: "
entry $w.plate.e -width 20 -relief sunken -bd 2 \
    -textvariable cp(platecoord)

$w.calp1.l config -text "Camera 1, Calibration image: "
entry $w.calp1.e -width 20 -relief sunken -bd 2 -textvariable cp(calp1)
$w.calp2.l config -text "Camera 2, Calibration image: "
entry $w.calp2.e -width 20 -relief sunken -bd 2 -textvariable cp(calp2)
$w.calp3.l config -text "Camera 3, Calibration image: "
entry $w.calp3.e -width 20 -relief sunken -bd 2 -textvariable cp(calp3)
$w.calp4.l config -text "Camera 4, Calibration image: "
entry $w.calp4.e -width 20 -relief sunken -bd 2 -textvariable cp(calp4)
$w.cori1.l config -text "Orientation data: "
entry $w.cori1.e -width 20 -relief sunken -bd 2 -textvariable cp(ori1)
$w.cori2.l config -text "Orientation data: "
entry $w.cori2.e -width 20 -relief sunken -bd 2 -textvariable cp(ori2)
$w.cori3.l config -text "Orientation data: "
entry $w.cori3.e -width 20 -relief sunken -bd 2 -textvariable cp(ori3)
$w.cori4.l config -text "Orientation data: "
entry $w.cori4.e -width 20 -relief sunken -bd 2 -textvariable cp(ori4)


# Image- and Pixelsize ****************

frame $w.isizeha
frame $w.isizehb
frame $w.isize
frame $w.isizeh
frame $w.isizev
frame $w.pixsize
frame $w.pixsizeh
frame $w.pixsizev
label $w.isize.l
label $w.isizeh.l
label $w.isizev.l
label $w.pixsize.l 
label $w.pixsizeh.l
label $w.pixsizev.l
$w.isize.l config -text "Imagesize,"
$w.isizeh.l config -text " horizontal: "
entry $w.isizeh.e -width 7 -relief sunken -bd 2 -textvariable  cp(imx)
$w.isizev.l config -text " vertical: "
entry $w.isizev.e -width 7 -relief sunken -bd 2 -textvariable cp(imy)
$w.pixsize.l config -text " Pixelsize, "
$w.pixsizeh.l config -text " horizontal: "
entry $w.pixsizeh.e -width 7 -relief sunken -bd 2  -textvariable cp(pix_x)
$w.pixsizev.l config -text " vertical: "
entry $w.pixsizev.e -width 7 -relief sunken -bd 2  -textvariable cp(pix_y)

# Target recognition

frame $w.targeta
frame $w.target
frame $w.gvthres
frame $w.minnpa
frame $w.minnp
frame $w.minnpx
frame $w.minnpy
frame $w.maxnpa
frame $w.maxnp
frame $w.maxnpx
frame $w.maxnpy
frame $w.m
frame $w.m.ppsumgv
frame $w.m.toldisc
frame $w.m.ppcross

label $w.target.l
label $w.gvthres.l
label $w.gvthres.l1
label $w.gvthres.l2
label $w.gvthres.l3
label $w.gvthres.l4
label $w.minnp.l
label $w.minnpx.l
label $w.minnpy.l
label $w.maxnp.l
label $w.maxnpx.l
label $w.maxnpy.l
label $w.m.toldisc.l
label $w.m.ppsumgv.l
label $w.m.ppcross.l

$w.target.l config -text "Target recognition on plate"  -font {Helvetica 12 bold}
$w.gvthres.l config -text "Grayvalue threshold, "
$w.gvthres.l1 config -text "  1.: "
$w.gvthres.l2 config -text "  2.: "
$w.gvthres.l3 config -text "  3.: "
$w.gvthres.l4 config -text "  4.: "
entry $w.gvthres.e1 -width 5 -relief sunken -bd 2 -textvariable cp(partgv1)
entry $w.gvthres.e2 -width 5 -relief sunken -bd 2 -textvariable cp(partgv2)
entry $w.gvthres.e3 -width 5 -relief sunken -bd 2 -textvariable cp(partgv3)
entry $w.gvthres.e4 -width 5 -relief sunken -bd 2 -textvariable cp(partgv4)

$w.minnp.l config -text "               min npix: "
entry $w.minnp.e -width 5 -relief sunken -bd 2 -textvariable cp(minnpix)
$w.minnpx.l config -text "                     min npix in x: "
entry $w.minnpx.e -width 5 -relief sunken -bd 2 -textvariable cp(minnpixx)
$w.minnpy.l config -text "          min npix in y: "
entry $w.minnpy.e -width 5 -relief sunken -bd 2 -textvariable cp(minnpixy)
  
$w.maxnp.l config -text "              max npix: "
entry $w.maxnp.e -width 5 -relief sunken -bd 2 -textvariable cp(maxnpix)
$w.maxnpx.l config -text "                    max npix in x: "
entry $w.maxnpx.e -width 5 -relief sunken -bd 2 -textvariable cp(maxnpixx)
$w.maxnpy.l config -text "         max npix in y: "
entry $w.maxnpy.e -width 5 -relief sunken -bd 2 -textvariable cp(maxnpixy)

$w.m.ppsumgv.l config -text "Sum of grayvalue: "
entry $w.m.ppsumgv.e -width 5 -relief sunken -bd 2 -textvariable cp(ppsumgv)
$w.m.toldisc.l config -text "     Tolerable discontinuity: "
entry $w.m.toldisc.e -width 5 -relief sunken -bd 2 -textvariable cp(toldisc)
$w.m.ppcross.l config -text "     Size of crosses: "
entry $w.m.ppcross.e -width 5 -relief sunken -bd 2 -textvariable cp(ppcrossize)


# Point number for manual pre-orientation
# Calibpicture 1

frame $w.img1a
frame $w.img1
frame $w.p11
frame $w.p12
frame $w.p13
frame $w.p14

label $w.mori
label $w.img1.l
label $w.p11.l
label $w.p12.l
label $w.p13.l
label $w.p14.l

$w.mori config -text "Point number for manual pre-orientation" -font {Helvetica 12 bold}

$w.img1.l config -text "Image 1:"
$w.p11.l config -text "  P1:"
entry $w.p11.e -width 4 -relief sunken -bd 2 -textvariable cp(p11)
$w.p12.l config -text "  P2:"
entry $w.p12.e -width 4 -relief sunken -bd 2 -textvariable cp(p12)
$w.p13.l config -text "  P3:"
entry $w.p13.e -width 4 -relief sunken -bd 2 -textvariable cp(p13)
$w.p14.l config -text "  P4:"
entry $w.p14.e -width 4 -relief sunken -bd 2 -textvariable cp(p14)

# Calibpicture 2

frame $w.img2a
frame $w.img2
frame $w.p21
frame $w.p22
frame $w.p23
frame $w.p24

label $w.img2.l
label $w.p21.l
label $w.p22.l
label $w.p23.l
label $w.p24.l

$w.img2.l config -text "Image 2:"
$w.p21.l config -text "  P1:"
entry $w.p21.e -width 4 -relief sunken -bd 2 -textvariable cp(p21)
$w.p22.l config -text "  P2:"
entry $w.p22.e -width 4 -relief sunken -bd 2 -textvariable cp(p22)
$w.p23.l config -text "  P3:"
entry $w.p23.e -width 4 -relief sunken -bd 2 -textvariable cp(p23)
$w.p24.l config -text "  P4:"
entry $w.p24.e -width 4 -relief sunken -bd 2 -textvariable cp(p24)


# Calibpicture 3

frame $w.img3a
frame $w.img3
frame $w.p31
frame $w.p32
frame $w.p33
frame $w.p34

label $w.img3.l
label $w.p31.l
label $w.p32.l
label $w.p33.l
label $w.p34.l

$w.img3.l config -text "Image 3:"
$w.p31.l config -text "  P1:"
entry $w.p31.e -width 4 -relief sunken -bd 2 -textvariable cp(p31)
$w.p32.l config -text "  P2:"
entry $w.p32.e -width 4 -relief sunken -bd 2 -textvariable cp(p32)
$w.p33.l config -text "  P3:"
entry $w.p33.e -width 4 -relief sunken -bd 2 -textvariable cp(p33)
$w.p34.l config -text "  P4:"
entry $w.p34.e -width 4 -relief sunken -bd 2 -textvariable cp(p34)


# Calibpicture 4

frame $w.img4a
frame $w.img4
frame $w.p41
frame $w.p42
frame $w.p43
frame $w.p44

label $w.img4.l
label $w.p41.l
label $w.p42.l
label $w.p43.l
label $w.p44.l

$w.img4.l config -text "Image 4:"
$w.p41.l config -text "  P1:"
entry $w.p41.e -width 4 -relief sunken -bd 2 -textvariable cp(p41)
$w.p42.l config -text "  P2:"
entry $w.p42.e -width 4 -relief sunken -bd 2 -textvariable cp(p42)
$w.p43.l config -text "  P3:"
entry $w.p43.e -width 4 -relief sunken -bd 2 -textvariable cp(p43)
$w.p44.l config -text "  P4:"
entry $w.p44.e -width 4 -relief sunken -bd 2 -textvariable cp(p44)

####################################################################
# Shaking parameters

frame $w.sh_img1a
frame $w.sh_img1
frame $w.sh_img2a
frame $w.sh_img2
frame $w.first_shake
frame $w.last_shake
frame $w.maxPoints_shake
frame $w.maxFrames_shake

label $w.sh_img1.l
label $w.sh_img2.l
label $w.mshake
label $w.first_shake.l
label $w.last_shake.l
label $w.maxPoints_shake.l
label $w.maxFrames_shake.l

$w.mshake config -text "Shaking parameters" -font {Helvetica 12 bold}

$w.first_shake.l config -text "  first frame:"
entry $w.first_shake.e -width 7 -relief sunken -bd 2 -textvariable cp(first_shake)
$w.last_shake.l config -text "  last frame:"
entry $w.last_shake.e -width 7 -relief sunken -bd 2 -textvariable cp(last_shake)
$w.maxPoints_shake.l config -text "  max # points used for orientation:"
entry $w.maxPoints_shake.e -width 4 -relief sunken -bd 2 -textvariable cp(maxPoints_shake)
$w.maxFrames_shake.l config -text "  max # frames used for orientation:"
entry $w.maxFrames_shake.e -width 4 -relief sunken -bd 2 -textvariable cp(maxFrames_shake)
######################################################################

# Orientation parameters

frame $w.oripar
frame $w.oripara
frame $w.pnrori
#################
frame $w.examineFlag
frame $w.combineFlag
#################
frame $w.oripp
frame $w.orilens
frame $w.oriaff
frame $w.oriinter

label $w.oripar.l
label $w.pnrori.l
#################
label $w.examineFlag.l
label $w.combineFlag.l
#################
label $w.lens
label $w.affin
label $w.interfaces

$w.oripar.l config -text "Orientation parameters"  -font {Helvetica 12 bold}
################
$w.examineFlag.l config -text "Calibrate with different z-positions?"
entry $w.examineFlag.e -width 2 -relief sunken -bd 2 -textvariable cp(examineFlag)
$w.combineFlag.l config -text "Combine preprocessed planes?"
entry $w.combineFlag.e -width 2 -relief sunken -bd 2 -textvariable cp(combineFlag)
################
$w.pnrori.l config -text "Point number for orientation"
entry $w.pnrori.e -width 2 -relief sunken -bd 2 -textvariable cp(pnrori)
checkbutton $w.pdist -text "Principle distance   " -variable cp(pdist)
checkbutton $w.ppointx -text "xp   " -variable cp(xp)
checkbutton $w.ppointy -text "yp   " -variable cp(yp)
$w.lens config -text "Lens distortion (Brown):  " 
checkbutton $w.lensk1 -text "k1  " -variable cp(k1)
checkbutton $w.lensk2 -text "k2  " -variable cp(k2)
checkbutton $w.lensk3 -text "k3  " -variable cp(k3)
checkbutton $w.lensp1 -text "p1  " -variable cp(p1)
checkbutton $w.lensp2 -text "p2  " -variable cp(p2)
$w.affin config -text "Affin transformation:   "
checkbutton $w.afftrafoscx -text "scx   " -variable cp(scx)
checkbutton $w.afftrafoshe -text "she   " -variable cp(she)
$w.interfaces config -text "Interfaces:   "
checkbutton $w.interfacesy -text "are variable   " -variable cp(interf)


# Plate of testcoord-, calib.picture- and pre-orientationfiles
pack $w.title -pady 8


pack $w.calp1.l $w.calp1.e $w.cori1.l $w.cori1.e -side left
pack $w.calp2.l $w.calp2.e $w.cori2.l $w.cori2.e -side left
pack $w.calp3.l $w.calp3.e $w.cori3.l $w.cori3.e -side left 
pack $w.calp4.l $w.calp4.e $w.cori4.l $w.cori4.e -side left

pack $w.calp1 $w.cori1 -side left -in $w.c1
pack $w.calp2 $w.cori2 -side left -in $w.c2
pack $w.calp3 $w.cori3 -side left -in $w.c3
pack $w.calp4 $w.cori4 -side left -in $w.c4

pack $w.c1 $w.c2 $w.c3 $w.c4 -side top -padx 8 -pady 2

pack $w.plate.l $w.plate.e -side left
pack $w.plate -pady 5

tooltip::tooltip $w.calp1.e "Path to the file for the calibration image of camera 1. 
Path and filename is also used for the result of orientation, i.e. the extension '.ori' is appended to this name.
Note that if this name + .ori is equal to 'Orientation data' to the left, then the initial guess will be automatically updated, but also overwritten."
tooltip::tooltip $w.calp2.e "Path to the file for the calibration image of camera 2. 
Path and filename is also used for the result of orientation, i.e. the extension '.ori' is appended to this name.
Note that if this name + .ori is equal to 'Orientation data' to the left, then the initial guess will be automatically updated, but also overwritten."
tooltip::tooltip $w.calp3.e "Path to the file for the calibration image of camera 3. 
Path and filename is also used for the result of orientation, i.e. the extension '.ori' is appended to this name.
Note that if this name + .ori is equal to 'Orientation data' to the left, then the initial guess will be automatically updated, but also overwritten."
tooltip::tooltip $w.calp4.e "Path to the file for the calibration image of camera 4. 
Path and filename is also used for the result of orientation, i.e. the extension '.ori' is appended to this name.
Note that if this name + .ori is equal to 'Orientation data' to the left, then the initial guess will be automatically updated, but also overwritten."

tooltip::tooltip $w.cori1.e "Path to the file with the initial guess for orientation of camera 1. The format of the file is:
x, y, z of camera position in (mm)
omega, phi, kappa of angles around x, y, z (rad). BEWARE for anglees 0,0,0 the camera looks into direction 0,0,-1 !!!
3 x 3 matrix entries don't matter, since omega, phi, kappa define this rotation matrix
dx, dy (microns) the chip is not perfectly centered around the optical axis
focal distance (mm), roughly you can guess it as chipsize / imagesize * dist from camera to image
vec_x, vec_y, vec_z (mm), orientation and distance of air-side of interface measured from origo defined by 'File of Coordinates on Plate''."
tooltip::tooltip $w.cori2.e "Path to the file with the initial guess for orientation of camera 2. The format of the file is:
x, y, z of camera position in (mm)
omega, phi, kappa of angles around x, y, z (rad). BEWARE for anglees 0,0,0 the camera looks into direction 0,0,-1 !!!
3 x 3 matrix entries don't matter, since omega, phi, kappa define this rotation matrix
dx, dy (microns) the chip is not perfectly centered around the optical axis
focal distance (mm), roughly you can guess it as chipsize / imagesize * dist from camera to image
vec_x, vec_y, vec_z (mm), orientation and distance of air-side of interface measured from origo defined by 'File of Coordinates on Plate''."
tooltip::tooltip $w.cori3.e "Path to the file with the initial guess for orientation of camera 3. The format of the file is:
x, y, z of camera position in (mm)
omega, phi, kappa of angles around x, y, z (rad). BEWARE for anglees 0,0,0 the camera looks into direction 0,0,-1 !!!
3 x 3 matrix entries don't matter, since omega, phi, kappa define this rotation matrix
dx, dy (microns) the chip is not perfectly centered around the optical axis
focal distance (mm), roughly you can guess it as chipsize / imagesize * dist from camera to image
vec_x, vec_y, vec_z (mm), orientation and distance of air-side of interface measured from origo defined by 'File of Coordinates on Plate''."
tooltip::tooltip $w.cori4.e "Path to the file with the initial guess for orientation of camera 4. The format of the file is:
x, y, z of camera position in (mm)
omega, phi, kappa of angles around x, y, z (rad). BEWARE for anglees 0,0,0 the camera looks into direction 0,0,-1 !!!
3 x 3 matrix entries don't matter, since omega, phi, kappa define this rotation matrix
dx, dy (microns) the chip is not perfectly centered around the optical axis
focal distance (mm), roughly you can guess it as chipsize / imagesize * dist from camera to image
vec_x, vec_y, vec_z (mm), orientation and distance of air-side of interface measured from origo defined by 'File of Coordinates on Plate''."

tooltip::tooltip $w.plate.e "Defines the path to the file where the coordinates of your calibration target are defined.
This is a 4 column ASCII file with point number, x,y,z positions in mm units. This file DEFINES the orientation of the coordinate system."


# Imageheader Tiff, frame or field for calibration
checkbutton $w.pic.tiff -text "TIFF-Header    " -variable cp(tiff)
radiobutton $w.pic.fram -variable cp(type) -text "Frame   " -value 0
radiobutton $w.pic.fieldodd -variable cp(type) -text "Field odd   " -value 1
radiobutton $w.pic.fieldeven -variable cp(type) -text "Field even" -value 2

tooltip::tooltip $w.pic.tiff "Should be checked."
tooltip::tooltip $w.pic.fram "Should be checked."
tooltip::tooltip $w.pic.fieldodd "Should NOT be checked."
tooltip::tooltip $w.pic.fieldeven "Should NOT be checked."


pack $w.pic.tiff $w.pic.fram $w.pic.fieldodd $w.pic.fieldeven -side left
pack $w.pic -pady 1

####################################################
pack $w.isize.l $w.isizeh.l  $w.isizeh.e $w.isizev.l $w.isizev.e   -side left
pack $w.pixsize.l $w.pixsizeh.l $w.pixsizeh.e $w.pixsizev.l  $w.pixsizev.e -side left

tooltip::tooltip $w.isizeh.e "number of pixels in horizontal direction."
tooltip::tooltip $w.isizev.e "number of pixels in vertical direction."
tooltip::tooltip $w.pixsizeh.e "The pixel size (mm) is provided by the camera manufacturer. You can also work it out by dividing chip width by horizontal pixel number."
tooltip::tooltip $w.pixsizev.e "The pixel size (mm) is provided by the camera manufacturer. You can also work it out by dividing chip height by vertical pixel number."


pack $w.isize $w.isizeh $w.isizev -in $w.isizeha -side left
pack $w.pixsize $w.pixsizeh $w.pixsizev -in $w.isizehb -side left

pack $w.isizeha $w.isizehb -pady 1 -padx 5

pack $w.target.l
pack $w.gvthres.l $w.gvthres.l1 $w.gvthres.e1 $w.gvthres.l2 $w.gvthres.e2 $w.gvthres.l3 $w.gvthres.e3 \
    $w.gvthres.l4 $w.gvthres.e4  -side left -pady 4
pack $w.minnp.l $w.minnp.e $w.minnpx.l $w.minnpx.e \
    $w.minnpy.l $w.minnpy.e -side left
pack $w.maxnp.l $w.maxnp.e $w.maxnpx.l $w.maxnpx.e \
    $w.maxnpy.l $w.maxnpy.e -side left



pack $w.target -side top -pady 6
pack $w.gvthres -side top -pady 2

tooltip::tooltip $w.gvthres "With these parameters you control the grayvalues your target points need to have in order to be recognized." 
tooltip::tooltip $w.minnp.e "How many pixels of a targetpoint neeed to be above grayvalue threshold in order to be recognized?"
tooltip::tooltip $w.minnpx.e "How many pixels in x direction of a targetpoint need to be above grayvalue threshold in order to be recognized?"
tooltip::tooltip $w.minnpy.e "How many pixels in y direction of a targetpoint need to be above grayvalue threshold in order to be recognized?"
tooltip::tooltip $w.maxnpx.e "Maximum number of pixels in x direction of a targetpoint above grayvalue threshold in order to still be recognized."
tooltip::tooltip $w.maxnpy.e "Maximum number of pixels in y direction of a targetpoint above grayvalue threshold in order to still be recognized."


pack $w.minnpa  $w.minnp $w.minnpx $w.minnpy -side top
pack $w.minnp $w.minnpx $w.minnpy $w.minnpy -in $w.minnpa -side left 
pack $w.minnpa -pady 2


pack $w.maxnpa  $w.maxnp $w.maxnpx $w.maxnpy -side top
pack $w.maxnp $w.maxnpx $w.maxnpy $w.maxnpy -in $w.maxnpa -side left 
pack $w.maxnpa -pady 2


pack  $w.m.ppsumgv.l $w.m.ppsumgv.e  $w.m.toldisc.l $w.m.toldisc.e  $w.m.ppcross.l $w.m.ppcross.e -side left
pack $w.m.ppsumgv $w.m.toldisc $w.m.ppcross -side left  -side left -in $w.m
pack $w.m -pady 3

tooltip::tooltip $w.m.ppsumgv.e "Minimum sum of all grayvalues of a traget point above threshold in order to be recognized."
tooltip::tooltip $w.m.toldisc.e "If - within a target point - grayvalues jump more than this value, then it is split into 2 target points."
tooltip::tooltip $w.m.ppcross.e "Influences only rendering."


# packing for point number for manual pre-orientation

pack $w.mori -pady 6
pack $w.img1.l $w.p11.l  $w.p11.e $w.p12.l $w.p12.e $w.p13.l \
	$w.p13.e $w.p14.l $w.p14.e  -side left
pack $w.img1 $w.p11 $w.p12 $w.p13 $w.p14  -in $w.img1a -side left
pack $w.img1a -pady 2
tooltip::tooltip $w.img1 "Points that need to be visible and manually clicked in the Calibration image 1."
tooltip::tooltip $w.p11 "First point that needs to be visible and manually clicked in Calibration image 1."
tooltip::tooltip $w.p12 "Second point that needs to be visible and manually clicked in Calibration image 1."
tooltip::tooltip $w.p13 "Third point that needs to be visible and manually clicked in Calibration image 1."
tooltip::tooltip $w.p14 "Fourth point that needs to be visible and manually clicked in Calibration image 1."

tooltip::tooltip $w.img2 "Points that need to be visible and manually clicked in the Calibration image 2."
tooltip::tooltip $w.p21 "First point that needs to be visible and manually clicked in Calibration image 2."
tooltip::tooltip $w.p22 "Second point that needs to be visible and manually clicked in Calibration image 2."
tooltip::tooltip $w.p23 "Third point that needs to be visible and manually clicked in Calibration image 2."
tooltip::tooltip $w.p24 "Fourth point that needs to be visible and manually clicked in Calibration image 2."

tooltip::tooltip $w.img3 "Points that need to be visible and manually clicked in the Calibration image 3."
tooltip::tooltip $w.p31 "First point that needs to be visible and manually clicked in Calibration image 3."
tooltip::tooltip $w.p32 "Second point that needs to be visible and manually clicked in Calibration image 3."
tooltip::tooltip $w.p33 "Third point that needs to be visible and manually clicked in Calibration image 3."
tooltip::tooltip $w.p34 "Fourth point that needs to be visible and manually clicked in Calibration image 3."

tooltip::tooltip $w.img4 "Points that need to be visible and manually clicked in the Calibration image 4."
tooltip::tooltip $w.p41 "First point that needs to be visible and manually clicked in Calibration image 4."
tooltip::tooltip $w.p42 "Second point that needs to be visible and manually clicked in Calibration image 4."
tooltip::tooltip $w.p43 "Third point that needs to be visible and manually clicked in Calibration image 4."
tooltip::tooltip $w.p44 "Fourth point that needs to be visible and manually clicked in Calibration image 4."


pack $w.img2.l $w.p21.l  $w.p21.e $w.p22.l $w.p22.e $w.p23.l \
	$w.p23.e $w.p24.l $w.p24.e  -side left
pack $w.img2 $w.p21 $w.p22 $w.p23 $w.p24 -in $w.img2a -side left
pack $w.img2a -pady 2

pack $w.img3.l $w.p31.l  $w.p31.e $w.p32.l $w.p32.e $w.p33.l \
	$w.p33.e $w.p34.l $w.p34.e  -side left
pack $w.img3 $w.p31 $w.p32 $w.p33 $w.p34 -in $w.img3a -side left
pack $w.img3a -pady 2

pack $w.img4.l $w.p41.l  $w.p41.e $w.p42.l $w.p42.e $w.p43.l \
	$w.p43.e $w.p44.l $w.p44.e -side left 
pack $w.img4 $w.p41 $w.p42 $w.p43 $w.p44 -in $w.img4a -side left
pack $w.img4a -pady 2

########################################################################
# packing for Shaking parameters

pack $w.mshake -pady 6
pack $w.sh_img1.l $w.first_shake.l  $w.first_shake.e $w.last_shake.l $w.last_shake.e -side left
pack $w.sh_img1 $w.first_shake $w.last_shake -in $w.sh_img1a -side left
pack $w.sh_img1a -pady 2

tooltip::tooltip $w.first_shake.e "Sequence for 'Orientation with particle position' begins at this frame."
tooltip::tooltip $w.last_shake.e "Sequence for 'Orientation with particle position' ends at this frame."

pack $w.mshake -pady 6
pack $w.sh_img2.l $w.maxPoints_shake.l  $w.maxPoints_shake.e $w.maxFrames_shake.l $w.maxFrames_shake.e -side left
pack $w.sh_img2 $w.maxPoints_shake $w.maxFrames_shake -in $w.sh_img2a -side left
pack $w.sh_img2a -pady 2

tooltip::tooltip $w.maxPoints_shake.e "If the number of linked(!) quadruplets of the first frame is not anyway larger than this number, 
also triplets will be included - possibly of more than one frame."
tooltip::tooltip $w.maxFrames_shake.e "If triplets are used to reach 'max # points used for orientation' the max number of frames is confined. 
The frames are spread equally between first and last frame, defined one line above."
########################################################################

# packing of checkbuttons for orientation parameters

pack $w.oripar.l -pady 8
##################
pack $w.examineFlag.l $w.examineFlag.e -side left
tooltip::tooltip $w.examineFlag.e "If you do a multiplane calibration then set to 1 otherwise to 0."
pack $w.combineFlag.l $w.combineFlag.e -side left
tooltip::tooltip $w.combineFlag.e "If you do a multiplane calibration AND if you have already treated all your individual planes (defined in the file parameters/multi_planes.par) 
and now you want to combine them, then set to 1 otherwise to 0."

##################
pack $w.pnrori.l $w.pnrori.e -side left
tooltip::tooltip $w.pnrori.e "Leave this number set to 0."
pack $w.oripar -pady 2
##################
pack $w.examineFlag -in $w.oripara
pack $w.combineFlag -in $w.oripara
##################
pack $w.pnrori -in $w.oripara
pack $w.oripara 
pack $w.pdist  $w.ppointx  $w.ppointy -in $w.oripp -side left
tooltip::tooltip $w.pdist "Principale distance is set free to be adjusted. Keep checked."
tooltip::tooltip $w.ppointx "chip offset dx is set free to be adjusted. Keep checked."
tooltip::tooltip $w.ppointy "chip offset dy  is set free to be adjusted. Keep checked."
pack $w.lens $w.lensk1 $w.lensk2 $w.lensk3 $w.lensp1 $w.lensp2 \
    -in $w.orilens -side left  -fill x
tooltip::tooltip $w.lens "It is recommended NOT TO USE all these (higher order) lens parameters." 
tooltip::tooltip $w.lensk1 "It is recommended NOT TO USE all these (higher order) lens parameters."
tooltip::tooltip $w.lensk2 "It is recommended NOT TO USE all these (higher order) lens parameters."
tooltip::tooltip $w.lensk3 "It is recommended NOT TO USE all these (higher order) lens parameters."
tooltip::tooltip $w.lensp1 "It is recommended NOT TO USE all these (higher order) lens parameters."
tooltip::tooltip $w.lensp2 "It is recommended NOT TO USE all these (higher order) lens parameters."
pack $w.affin $w.afftrafoscx $w.afftrafoshe -in $w.oriaff -side left
tooltip::tooltip $w.affin "It is recommended NOT TO USE all these (higher order) lens parameters." 
tooltip::tooltip $w.afftrafoscx "It is recommended NOT TO USE all these (higher order) lens parameters."
tooltip::tooltip $w.afftrafoshe "It is recommended NOT TO USE all these (higher order) lens parameters."
pack $w.interfaces  $w.interfacesy -in $w.oriinter -side left
tooltip::tooltip $w.interfaces "It is recommended to make an as good as possible guess for the orientation of the interfaces ('Orientation data') but then leave them as they are."
tooltip::tooltip $w.interfacesy "It is recommended to make an as good as possible guess for the orientation of the interfaces ('Orientation data') but then leave them as they are."
pack $w.oripp $w.orilens $w.oriaff $w.oriinter -side top


# OK-button ##################################
button $w.ok -text OK -command "done_proc_cmd ;destroy $w"
pack $w.ok -side bottom -pady 5
}