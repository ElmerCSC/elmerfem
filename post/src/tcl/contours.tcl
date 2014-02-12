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
#* Contours display parameter settings
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
# 23 Apr 1996
#

set ContourLineStyle   0
set ContourQuality     1
set ContourRadius      1
set ContourColor       "none"
set ContourContour     "none"
set ContourColorMin    0.0
set ContourColorMax    1.0

set ContourLines         5
set CurrentLines         0

set ContourActive        0
set ContourColorMap(0,R) 0
set ContourColorMap(0,G) 0
set ContourColorMap(0,B) 0

proc contour_set_color { win args } {
    global ContourColorMap ContourActive

    set R [$win.red get]
    set G [$win.grn get]
    set B [$win.blu get]

    set R [@ int($R*2.55+0.5)]
    set G [@ int($G*2.55+0.5)]
    set B [@ int($B*2.55+0.5)]

    set value [format "#%02x%02x%02x" $R $G $B]
    .contour.cont.values.fr$ContourActive.valuecolor configure -back $value
}

proc contour_set_value_array { lines ColorMin ColorMax } {
    global ContourValues

    if {  $lines > 0 } {
       do i 0 [@ $lines-1] {
         set t [@ ($i+1.0)/($lines+1.0)]
         set ContourValues($i) [@ (1-$t)*$ColorMin + $t*$ColorMax]
       }
    }
}

proc contour_set_values { win lines ColorMin ColorMax } {
    global CurrentLines ContourValues colmap colmap_size ContourColorMap
    global ContourValues ContourColorMin ContourColorMax

    if {  $lines > 0 } {

       do i 0 [@ $CurrentLines-1] {
           if { [winfo exists $win.fr$i.value] } { destroy $win.fr$i.value }
           if { [winfo exists $win.fr$i.valuecolor] } { destroy $win.fr$i.valuecolor }
           if { [winfo exists $win.fr$i] } { destroy $win.fr$i }
       }

       contour_set_value_array $lines $ColorMin $ColorMax

       do i 0 [@ $lines-1] {
           set a [@ $ContourColorMax - $ContourColorMin]
           set b [@ $ContourValues($i)-$ContourColorMin]

           if { $a==0 } { set a 1.0 }
           set t [@ int(($colmap_size-1.0)*$b/$a+0.5)]

           frame $win.fr$i
           entry $win.fr$i.value -textvariable ContourValues($i) -width 12
           pack $win.fr$i.value -side left

           button $win.fr$i.valuecolor \
              -back $colmap([@ ($t<0)?0:($t>=$colmap_size)?$colmap_size-1:$t]) \
              -command "set ContourActive $i"
           pack $win.fr$i.valuecolor -side left

           pack $win.fr$i

           bind $win.fr$i.value <Return> {
                    set val [%W get];
                    set a [@ $ContourColorMax - $ContourColorMin];
                    set b [@ $val - $ContourColorMin];
                    set t [@  int(($colmap_size-1.0)*$b/$a+0.5)];
                    if { $t < 0 } { set t 0 }
                    if { $t >= $colmap_size } { set t [@ $colmap_size-1] }
                    %Wcolor configure -back $colmap($t);
                }
       }

       set CurrentLines $lines
    }
}

proc contour_edit { } {
    global ContourLines ContourLineStyle ContourQuality ContourRadius
    global ContourColor ContourContour ContourColorMin ContourColorMax
    global ContourColorSetMinMax

    if { [winfo exists .contour] } {
        wm iconify .contour
        wm deiconify .contour
        return
    }

    toplevel .contour
    place_window .contour

    frame .contour.cont
    label .contour.cont.label -text "Number Of Contours: "
    entry .contour.cont.entry -width 5 -textvariable ContourLines -relief sunken

    bind .contour.cont.entry <Return> { contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

    frame .contour.cont.values
# -yscrollcommand ".contour.cont.values.scroll set" -width 200 -height 200
#    scrollbar .contour.cont.values.scroll -command ".contour.cont.values yview"

    pack .contour.cont -side top
    pack .contour.cont.label -side left
    pack .contour.cont.entry -side left -fill x

    contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax
    pack .contour.cont.values -side top
#   pack .contour.cont.values.scroll -side left -expand 1 -fill both

#
# Generate ...
#
    frame .contour.set
    label .contour.set.min_lab -text "Min: "

    entry .contour.set.min -width 10 -textvariable ContourColorMin
    bind .contour.set.min <Return> { contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

    label .contour.set.max_lab -text "Max: "

    entry .contour.set.max -width 10 -textvariable ContourColorMax
    bind .contour.set.max <Return> { contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

#    button .contour.set.gen -text "Generate" -command { \
#         contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

    checkbutton .contour.set.keep -text "Keep" -variable ContourColorSetMinMax -command { \
         contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

    pack .contour.set.min_lab -side left
    pack .contour.set.min -side left
    pack .contour.set.max_lab -side left
    pack .contour.set.max -side left
#    pack .contour.set.gen -side left
    pack .contour.set.keep -side left
    pack .contour.set -side top
    
# color sliders
#
#    frame .contour.rgb
#    slider .contour.rgb.red -orient horizontal -command { contour_set_color .contour.rgb } \
#             -from 0 -to 100 -troughcolor red -digit 4 -resol 0.5
#    slider .contour.rgb.grn -orient horizontal -command { contour_set_color .contour.rgb } \
#             -from 0 -to 100 -troughcolor green -digit 4 -resol 0.5
#    slider .contour.rgb.blu -orient horizontal -command { contour_set_color .contour.rgb } \
#            -from 0 -to 100 -troughcolor blue -digit 4 -tick 25 -resol 0.5
#
#    pack .contour.rgb.red -side left -expand 1 -fill x
#    pack .contour.rgb.grn -side left -expand 1 -fill x
#    pack .contour.rgb.blu -side left -expand 1 -fill x
#
#    pack .contour.rgb.red -side top -fill x
#    pack .contour.rgb.grn -side top -fill x
#    pack .contour.rgb.blu -side top  -fill x
#    pack .contour.rgb -side top -expand 1 -fill both
#
#
#
    frame .contour.line
    label .contour.line.label -text "Line Style: "
    radiobutton .contour.line.line -value 0 -variable ContourLineStyle -text "Line"
    radiobutton .contour.line.cyli -value 1 -variable ContourLineStyle -text "Solid"

    pack .contour.line -side top
    pack .contour.line.label -side left
    pack .contour.line.line -side left -fill x
    pack .contour.line.cyli -side left  -fill x

    frame .contour.qual
    label .contour.qual.label -text "Line Quality: "
    entry .contour.qual.entry -relief sunken -width 5 -textvariable ContourQuality

    pack .contour.qual -side top
    pack .contour.qual.label -side left
    pack .contour.qual.entry -side left -fill x

    frame .contour.radi
    label .contour.radi.label -text "Width Scale: "
    entry .contour.radi.entry -relief sunken -width 5 -textvariable ContourRadius

    pack .contour.radi -side top
    pack .contour.radi.label -side left
    pack .contour.radi.entry -side left -fill x

#
#
#
#    frame .contour.iso
#    label .contour.iso.label -text "Contour Variable: "
#    button .contour.iso.but -textvariable ContourContour       \
#              -command { set ContourContour [make_scalar_list]; \
#                        UpdateVariable "ContourContour";        \
#                        contour_set_values .contour.cont.values $ContourLines }
#
#    UpdateVariable "ContourContour"
#    contour_set_values .contour.cont.values $ContourLines
#
#    pack .contour.iso -side top
#    pack .contour.iso.label -side left
#    pack .contour.iso.but -side left -fill x
#
#
#
    frame .contour.vari
    label .contour.vari.label -text "Color Variable: "
    button .contour.vari.but -textvariable ContourColor       \
              -command { set ContourColor [make_scalar_list]; \
                         UpdateVariable "ContourColor";       \
                         .contour.set.min configure -textvariable ContourColorMin; \
                         .contour.set.max configure -textvariable ContourColorMax; \
                         contour_set_values .contour.cont.values $ContourLines $ContourColorMin $ContourColorMax }

    UpdateVariable "ContourColor"

    pack .contour.vari -side top
    pack .contour.vari.label -side left
    pack .contour.vari.but -side left -fill x
#
#
#


    frame .contour.buttons
    button .contour.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .contour.buttons.close -text "Close" -command "destroy .contour"

    pack .contour.buttons -side top
    pack .contour.buttons.apply -side left
    pack .contour.buttons.close -side left -fill x
}

contour_set_value_array $ContourLines $ContourColorMin $ContourColorMax
