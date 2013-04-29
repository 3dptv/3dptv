# Changing of Tracking Parameters
########################################

proc trackpar { } {
set w .changetrack
catch {destroy $w}
toplevel .changetrack
wm title $w "Changing Parameters"

global tp

label $w.title -text "Tracking Parameters"  -font {Helvetica 12 bold} -width 50 -anchor center
pack $w.title -pady 10

#---------------------Test with 2 different Tracking posibilities-------------#


frame $w.option

radiobutton $w.option.2d -text 2D-Tracking -variable tp(choice) -value 2d
radiobutton $w.option.3d -text 3D-Tracking -variable tp(choice) -value 3d

pack $w.option.2d $w.option.3d -side left
pack $w.option



# Tracking Parameters

frame $w.dv
frame $w.dx
frame $w.dy
frame $w.dz
frame $w.da
frame $w.df
frame $w.dvxmin
frame $w.dvxmax
frame $w.dvymin
frame $w.dvymax
frame $w.dvzmin
frame $w.dvzmax
frame $w.dacc
frame $w.dangle
frame $w.dsumg
frame $w.dn
frame $w.dnx
frame $w.dny

label $w.feature -text "Differences of Particle Features" -width 20 -anchor center
label $w.dvxmin.label
label $w.dvxmax.label
label $w.dvymin.label
label $w.dvymax.label
label $w.dvzmin.label
label $w.dvzmax.label
label $w.dacc.label
label $w.dangle.label

label $w.dn.label
label $w.dnx.label
label $w.dny.label
label $w.dsumg.label


$w.dvxmin.label config -text "       dvxmin: "
entry $w.dvxmin.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvxmin)
$w.dvxmax.label config -text "  dvxmax: "
entry $w.dvxmax.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvxmax)
$w.dvymin.label config -text "       dvymin: "
entry $w.dvymin.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvymin)
$w.dvymax.label config -text "  dvymax: "
entry $w.dvymax.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvymax)
$w.dvzmin.label config -text "       dvzmin: "
entry $w.dvzmin.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvzmin)
$w.dvzmax.label config -text "  dvzmax: "
entry $w.dvzmax.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dvzmax)
$w.dangle.label config -text " angle \[gon\]: "
entry $w.dangle.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dangle)
$w.dacc.label config -text "      dacc: "
entry $w.dacc.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dacc)
$w.dsumg.label config -text "Sum of greyvalues: "
entry $w.dsumg.entry -width 6 -relief sunken -bd 2 \
    -textvariable tp(dsumg)
$w.dn.label config -text "Size in pixel:        "
entry $w.dn.entry -width 4 -relief sunken -bd 2 \
    -textvariable tp(dn)
$w.dnx.label config -text "Size in x:           "
entry $w.dnx.entry -width 3 -relief sunken -bd 2 \
    -textvariable tp(dnx)
$w.dny.label config -text "Size in y:           "
entry $w.dny.entry -width 3 -relief sunken -bd 2 \
    -textvariable tp(dny)

checkbutton $w.add -text "Add new particles position" -variable tp(add)

pack $w.title -pady 10

pack $w.dvxmin.label $w.dvxmin.entry -pady 2 -side left
pack $w.dvxmax.label $w.dvxmax.entry -pady 2 -side left
pack $w.dvymin.label $w.dvymin.entry -pady 2 -side left
pack $w.dvymax.label $w.dvymax.entry -pady 2 -side left 
pack $w.dvzmin.label $w.dvzmin.entry -pady 2 -side left
pack $w.dvzmax.label $w.dvzmax.entry -pady 2 -side left 
pack $w.dangle.label $w.dangle.entry -pady 6 -side left
pack $w.dacc.label $w.dacc.entry -pady 6 -side left  

pack $w.dvxmin $w.dvxmax -in $w.dx -side left
pack $w.dvymin $w.dvymax -in $w.dy -side left
pack $w.dvzmin $w.dvzmax -in $w.dz -side left
pack $w.dangle $w.dacc -in $w.da -side left
pack $w.dx $w.dy $w.dz $w.dz $w.da -in $w.dv -side top

pack $w.dsumg.label $w.dsumg.entry -pady 2 -side left  
pack $w.dn.label $w.dn.entry -pady 2 -side left  
pack $w.dnx.label $w.dnx.entry -pady 2 -side left  
pack $w.dny.label $w.dny.entry -pady 2 -side left
  
pack $w.dv -padx 10 -pady 5 -side top
pack $w.add -pady 5



#-------------------new tracking parameters octobre 08-------------------------#
#----textboxes-----#

frame $w.text
frame $w.text.a
frame $w.text.b
frame $w.text.a.maxnum
frame $w.text.a.deltat
frame $w.text.a.maxvel
frame $w.text.b.linktol
frame $w.text.b.jumptol

label $w.text.t
label $w.text.a.maxnum.label
label $w.text.a.deltat.label
label $w.text.a.maxvel.label
label $w.text.b.linktol.label
label $w.text.b.jumptol.label


$w.text.t config -text "New tracking algorithm" -font {Helvetica 13 bold}

$w.text.a.maxnum.label config -text "Max number per frame: " 
entry $w.text.a.maxnum.entry -width 5 -relief sunken -bd 2 -textvariable tp(maxnum)
$w.text.a.deltat.label config -text "  Delta t: "
entry $w.text.a.deltat.entry -width 5 -relief sunken -bd 2 -textvariable tp(deltat)
$w.text.a.maxvel.label config -text "  Max velocitiy: "
entry $w.text.a.maxvel.entry -width 5 -relief sunken -bd 2 -textvariable tp(vmax)
$w.text.b.linktol.label config -text "Linking tolerance: "
entry $w.text.b.linktol.entry -width 5 -relief sunken -bd 2 -textvariable tp(linktol)
$w.text.b.jumptol.label config -text "              Jumping tolerance: "
entry $w.text.b.jumptol.entry -width 5 -relief sunken -bd 2 -textvariable tp(jumptol)

pack $w.text.t -pady 10

pack $w.text.a.maxnum.label $w.text.a.maxnum.entry $w.text.a.deltat.label $w.text.a.deltat.entry $w.text.a.maxvel.label $w.text.a.maxvel.entry -side left
pack $w.text.a.maxnum $w.text.a.deltat $w.text.a.maxvel -side left -in $w.text.a
pack $w.text.a -pady 5

pack $w.text.b.linktol.label $w.text.b.linktol.entry $w.text.b.jumptol.label $w.text.b.jumptol.entry -side left
pack $w.text.b.linktol $w.text.b.jumptol -side left -in $w.text.b
pack $w.text.b -pady 5

pack $w.text.a $w.text.b -side top
pack $w.text -pady 10


#----checkboxes-----#

frame $w.check
frame $w.check.1
frame $w.check.2
frame $w.check.3
frame $w.check.4

checkbutton $w.check.1.gluing -text "  Gluing           " -variable tp(gluing)
checkbutton $w.check.1.track -text "  Track                       " -variable tp(track)
checkbutton $w.check.2.iterate -text " Iterate            " -variable tp(iterate)
checkbutton $w.check.2.remkin -text "  Remove kints           " -variable tp(remkin)
checkbutton $w.check.3.writeb -text "  Write binary  " -variable tp(wbin)
checkbutton $w.check.3.writet -text "  Write trajectories      " -variable tp(wtraj)
checkbutton $w.check.4.output -text "  Output          " -variable tp(output)
checkbutton $w.check.4.camconfig -text "  Same camera config" -variable tp(camconfig)


pack $w.check.1.gluing $w.check.1.track -side left
pack $w.check.1 -pady 5

pack $w.check.2.iterate $w.check.2.remkin -side left
pack $w.check.2 -pady 5

pack $w.check.3.writeb $w.check.3.writet -side left
pack $w.check.3 -pady 5

pack $w.check.4.output $w.check.4.camconfig -side left
pack $w.check.4 -pady 5

pack $w.check.1 $w.check.2 $w.check.3 $w.check.4 -side top
pack $w.check -pady 10


# OK-button ##################################
button $w.ok -text OK -command "done_proc_cmd ;destroy $w"
pack $w.add -pady 5
pack $w.ok -side bottom -pady 10
}
