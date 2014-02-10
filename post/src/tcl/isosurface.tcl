#/*****************************************************************************
# *
# *  Elmer, A Finite Element Software for Multiphysical Problems
# *
# *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
# * 
# *  This program is free software; you can redistribute it and/or
# *  modify it under the terms of the GNU General Public License
# *  as published by the Free Software Foundation; either version 2
# *  of the License, or (at your option) any later version.
# * 
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program (in file fem/GPL-2); if not, write to the 
# *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
# *  Boston, MA 02110-1301, USA.
# *
# *****************************************************************************/

#*******************************************************************************
#*
#* Isosurfaces display parameter settings
#*
#*******************************************************************************
#*
#*                     Author:       Juha Ruokolainen
#*
#*                    Address: CSC - IT Center for Science Ltd.
#*                                Keilaranta 14, P.O. BOX 405
#*                                  02101 Espoo, Finland
#*                                  Tel. +358 0 457 2723
#*                                Telefax: +358 0 457 2302
#*                              EMail: Juha.Ruokolainen@csc.fi
#*
#*                       Date: 26 Sep 1995
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#*******************************************************************************
#
#

#
#
# 26 Apr 1996
#

set IsosurfaceStyle       0
set IsosurfaceLineStyle   0
set IsosurfaceQuality     1
set IsosurfaceRadius      1

set IsosurfaceContour     "none"
set IsosurfaceColor       "none"
set IsosurfaceNormal      "none"

set IsosurfaceContourMin   0.0
set IsosurfaceContourMax   1.0
set IsosurfaceContourSetMinMax 0

set IsosurfaceContours     1
set CurrentIsoContours     0
set IsosurfaceRecompute    1


set IsosurfaceColorMin  0.0
set IsosurfaceColorMax  1.0
set IsosurfaceColorSetMinMax 0


proc isosurface_update {} {

   global IsosurfaceColor IsosurfaceColorMin IsosurfaceColorMax

   UpdateVariable IsosurfaceColor

   .isosurface.cset.min delete 0 end
   .isosurface.cset.min insert end [format %-10.5g $IsosurfaceColorMin]

   .isosurface.cset.max delete 0 end
   .isosurface.cset.max insert end [format %-10.5g $IsosurfaceColorMax]
}


proc isosurface_set_value_array { lines ColorMin ColorMax } {
    global IsosurfaceValues

    do i 0 [@ $lines-1] {
       set t [@ ($i+1.0)/($lines+1.0)];
       set IsosurfaceValues($i) [@ (1-$t)*$ColorMin + $t*$ColorMax];
    }
}

proc isosurface_set_values { win lines ColorMin ColorMax } {
    global CurrentIsoContours IsosurfaceValues

    do i 0 [@ $CurrentIsoContours-1] {
        if { [winfo exists $win.value$i] } { destroy $win.value$i }
    }

    isosurface_set_value_array $lines $ColorMin $ColorMax

    do i 0 [@ $lines-1] {
        entry $win.value$i -textvariable IsosurfaceValues($i) -width 12
        pack $win.value$i
#        $win create window 0 [@ $i*20] -window $win.value$i
    }

     set CurrentIsoContours $lines
}

proc isosurface_edit { } {
    global IsosurfaceContours IsosurfaceLineStyle IsosurfaceQuality IsosurfaceRadius
    global IsosurfaceContour IsosurfaceNormal IsosurfaceColor CurrentIsoContours
    global IsosurfaceColorMin IsosurfaceColorMax IsosurfaceColorSetMinMax
    global IsosurfaceContourMin IsosurfaceContourMax IsosurfaceContourSetMinMax

    if { [winfo exists .isosurface] } {
        wm iconify .isosurface
        wm deiconify .isosurface
        return
    }

    toplevel .isosurface
    place_window .isosurface

    frame .isosurface.cont
    label .isosurface.cont.label -text "Number Of Isosurfaces: "
    entry .isosurface.cont.entry -width 5 -textvariable IsosurfaceContours -relief sunken

    bind .isosurface.cont.entry <Return> { isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax }

    frame .isosurface.cont.values
#     -yscrollcommand ".isosurface.cont.values.scroll set" -width 200 -height 200
#    scrollbar .isosurface.cont.values.scroll -command ".isosurface.cont.values yview"

    pack .isosurface.cont -side top
    pack .isosurface.cont.label -side left
    pack .isosurface.cont.entry -side left -fill x

    isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax
    pack .isosurface.cont.values -side top
#    pack .isosurface.cont.values.scroll -side left -expand 1 -fill both

#
# Generate ...
#
    frame .isosurface.set
    label .isosurface.set.min_lab -text "Min: "

    entry .isosurface.set.min -width 10 -textvariable IsosurfaceContourMin
    bind .isosurface.set.min <Return> { isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax }

    label .isosurface.set.max_lab -text "Max: "

    entry .isosurface.set.max -width 10 -textvariable IsosurfaceContourMax
    bind .isosurface.set.max <Return> { isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax }

    checkbutton .isosurface.set.keep -text "Keep" -variable IsosurfaceContourSetMinMax -command { \
         isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax }

    pack .isosurface.set.min_lab -side left
    pack .isosurface.set.min -side left
    pack .isosurface.set.max_lab -side left
    pack .isosurface.set.max -side left
#   pack .isosurface.set.gen -side left
    pack .isosurface.set.keep -side left
    pack .isosurface.set -side top

#
#
#
    frame .isosurface.style
    label .isosurface.style.label -text "Surface Style: "
    radiobutton .isosurface.style.line -value 0 -variable IsosurfaceStyle -text "Line"
    radiobutton .isosurface.style.surf -value 1 -variable IsosurfaceStyle -text "Surface"
    radiobutton .isosurface.style.both -value 2 -variable IsosurfaceStyle -text "Both"

    pack .isosurface.style -side top
    pack .isosurface.style.label -side left
    pack .isosurface.style.line -side left -fill x
    pack .isosurface.style.surf -side left  -fill x
    pack .isosurface.style.both -side left  -fill x

    frame .isosurface.line
    label .isosurface.line.label -text "Line Style: "
    radiobutton .isosurface.line.line -value 0 -variable IsosurfaceLineStyle -text "Line"
    radiobutton .isosurface.line.cyli -value 1 -variable IsosurfaceLineStyle -text "Solid"

    pack .isosurface.line -side top
    pack .isosurface.line.label -side left
    pack .isosurface.line.line -side left -fill x
    pack .isosurface.line.cyli -side left  -fill x

    frame .isosurface.qual
    label .isosurface.qual.label -text "Line Quality: "
    entry .isosurface.qual.entry -relief sunken -width 5 -textvariable IsosurfaceQuality

    pack .isosurface.qual -side top
    pack .isosurface.qual.label -side left
    pack .isosurface.qual.entry -side left -fill x

    frame .isosurface.radi
    label .isosurface.radi.label -text "Width Scale: "
    entry .isosurface.radi.entry -relief sunken -width 5 -textvariable IsosurfaceRadius

    pack .isosurface.radi -side top
    pack .isosurface.radi.label -side left
    pack .isosurface.radi.entry -side left -fill x

#
# contour
#
    frame .isosurface.contours
    label .isosurface.contours.label -text "Contour Variable: "
    button .isosurface.contours.but -textvariable IsosurfaceContour       \
              -command { set IsosurfaceContour [make_scalar_list]; \
                        UpdateVariable "IsosurfaceContour";        \
                        .isosurface.set.min configure -textvariable IsosurfaceContourMin; \
                        .isosurface.set.max configure -textvariable IsosurfaceContourMax; \
                        isosurface_set_values .isosurface.cont.values $IsosurfaceContours  $IsosurfaceContourMin $IsosurfaceContourMax }

    UpdateVariable "IsosurfaceContour"
    isosurface_set_values .isosurface.cont.values $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax

    pack .isosurface.contours -side top
    pack .isosurface.contours.label -side left
    pack .isosurface.contours.but -side left -fill x
#
# color
#
    frame .isosurface.vari
    label .isosurface.vari.label -text "Color Variable: "
    button .isosurface.vari.but -textvariable IsosurfaceColor     \
              -command { set IsosurfaceColor [make_scalar_list]; \
                        UpdateVariable "IsosurfaceColor";  } 

    pack .isosurface.vari -side top
    pack .isosurface.vari.label -side left
    pack .isosurface.vari.but -side left -fill x

#
# set color min max
#
    frame .isosurface.cset

    label .isosurface.cset.min_lab -text "Min: "
    entry .isosurface.cset.min -width 10 -textvariable IsosurfaceColorMin
    bind .isosurface.cset.min <Return> isosurface_update

    label .isosurface.cset.max_lab -text "Max: "
    entry .isosurface.cset.max -width 10 -textvariable IsosurfaceColorMax
    bind .isosurface.cset.max <Return> isosurface_update

    checkbutton .isosurface.cset.keep -text "Keep" -variable IsosurfaceColorSetMinMax \
          -command isosurface_update

    pack .isosurface.cset.min_lab -side left
    pack .isosurface.cset.min -side left
    pack .isosurface.cset.max_lab -side left
    pack .isosurface.cset.max -side left
    pack .isosurface.cset.keep -side left
    pack .isosurface.cset -side top
    
#
#
#
    frame .isosurface.norm
    label .isosurface.norm.label -text "Surface Normal Variable: "
    button .isosurface.norm.but -textvariable IsosurfaceNormal     \
              -command { set IsosurfaceNormal [make_vector_list]; }

    pack .isosurface.norm -side top
    pack .isosurface.norm.label -side left
    pack .isosurface.norm.but -side left -fill x

#
#
#
    frame .isosurface.buttons
    button .isosurface.buttons.apply -text "Apply" -command {
         if { $CurrentIsoContours != $IsosurfaceContours } { \
            isosurface_set_values .isosurface.cont.values $IsosurfaceContours \
            $IsosurfaceContourMin $IsosurfaceContourMax }; set IsosurfaceRecompute 1; UpdateObject; play;}
    button .isosurface.buttons.close -text "Close" -command "destroy .isosurface"

    pack .isosurface.buttons -side top
    pack .isosurface.buttons.apply -side left
    pack .isosurface.buttons.close -side left -fill x
}

isosurface_set_value_array $IsosurfaceContours $IsosurfaceContourMin $IsosurfaceContourMax

