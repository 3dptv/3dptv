# Changing of Tracking Parameters
########################################

proc trackpar { } {
set w .changetrack
catch {destroy $w}
toplevel .changetrack
wm title $w "Changing Parameters"

global tp

label $w.title -text "Tracking Parameters"  -font {Helvetica 12 bold} -width 20 -anchor center

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

label $w.feature -text "Differences of Particle Features" -width 30 -anchor center
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
  
#pack $w.dsumg $w.dn $w.dnx $w.dny -in $w.df -side top

#pack $w.dv $w.feature $w.df -pady 5 -side top
pack $w.dv -padx 10 -pady 5 -side top

# OK-button ##################################
button $w.ok -text OK -command "done_proc_cmd ;destroy $w"
pack $w.add -pady 5
pack $w.ok -side bottom -pady 10
}
