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
#* Vectors display parameter settings
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
# 23 Apr 1996
#

set VectorLineStyle         0
set VectorQuality           1
set VectorRadius            1
set VectorLength            "none"
set VectorColor             "none"
set VectorArrow             "none"
set VectorThreshold         "none"
set VectorFloor             "0.0"
set VectorCeiling           "1.0"
set VectorLengthScale       1.0

set VectorThresholdMin      0.0
set VectorThresholdMax      1.0

set VectorVariableNames(0)  "none"
set NumberOfVectorVariables 1

set scalar "none"
set vector "none"

set vceil  100.0
set vfloor   0.0

proc get_scalar_variable { w } {
    global scalar

    set scalar [$w get [$w curselection]]
}

proc get_vector_variable { w } {
    global vector

    set vector [$w get [$w curselection]]
}

proc vector_set_ceil { vceil } {
    global VectorCeiling VectorThresholdMin VectorThresholdMax

    set a [@ ($VectorThresholdMax-$VectorThresholdMin)*$vceil/100.0+$VectorThresholdMin]
    set VectorCeiling $a
}

proc vector_set_floor { vfloor } {
    global VectorFloor VectorThresholdMin VectorThresholdMax

    set a [@ ($VectorThresholdMax-$VectorThresholdMin)*$vfloor/100.0+$VectorThresholdMin]
    set VectorFloor $a
}

proc vector_set_vceil { } { 
    global VectorCeiling VectorThresholdMin VectorThresholdMax vceil

    set vceil [@ 100*($VectorCeiling-$VectorThresholdMin)/($VectorThresholdMax-$VectorThresholdMin)]
    set vceil [@ $vceil<0.0?0.0:$vceil]
    set vceil [@ $vceil>100.0?100.0:$vceil]
}

proc vector_set_vfloor { } { 
    global VectorFloor VectorThresholdMin VectorThresholdMax vfloor

    set vfloor [@ 100*($VectorFloor-$VectorThresholdMin)/($VectorThresholdMax-$VectorThresholdMin)]
    set vfloor [@ $vfloor<0.0?0.0:$vfloor]
    set vfloor [@ $vfloor>100.0?100.0:$vfloor]
}

#
# List box containing scalar variable names
#
proc make_scalar_list { } {
   global ScalarVariableNames NumberOfScalarVariables scalar savescalar

   if { ![info exists ScalarVariableNames] } { return }

   set savescalar $scalar

   toplevel .vlist
   place_window .vlist

   frame .vlist.vari -relief sunken -bg lightblue
   listbox .vlist.vari.list -yscroll ".vlist.vari.scroll set"
   scrollbar .vlist.vari.scroll -command ".vlist.vari.list yview"

   pack .vlist.vari.list -side left -fill y
   pack .vlist.vari.scroll -side left -expand 1 -fill both
   pack .vlist.vari -side top -expand 1 -fill both

   bind .vlist.vari.list <Double-1> { set scalar [get_scalar_variable %W] }


   frame .vlist.equ
   label .vlist.equ.label -text "math: "
   entry .vlist.equ.entry -width 20 -relief sunken -textvariable mathcmd
   pack .vlist.equ -side top
   pack .vlist.equ.label -side left
   pack .vlist.equ.entry -side left -fill x

   bind .vlist.equ.entry <Return> { math $mathcmd; set scalar "tryagain"; }

   frame .vlist.close
   button .vlist.close.ok  -text  "OK" -command { set scalar [get_scalar_variable .vlist.vari.list]  }
   button .vlist.close.cancel -text  "Cancel" -command { set scalar $savescalar; }

   pack .vlist.close -side top
   pack .vlist.close.ok -side right
   pack .vlist.close.cancel -side right

   set oldfocus [focus]
   set scalar "tryagain"
   while { $scalar == "tryagain" } {
       .vlist.vari.list delete 0 end

       .vlist.vari.list insert end none
       do i 0 [@ $NumberOfScalarVariables-1] {
           set val $ScalarVariableNames($i)
           .vlist.vari.list insert end $val
       }

       vwait scalar
   }
   focus $oldfocus
   destroy .vlist

   return $scalar
}

#
# List box containing vector variable names
#
proc make_vector_list { } {
   global VectorVariableNames NumberOfVectorVariables vector savevector

   if { ![info exists VectorVariableNames] } { return }

   set savevector $vector

   toplevel .vlist
   place_window .vlist

   frame .vlist.vari -relief sunken -bg lightblue
   listbox .vlist.vari.list -yscroll ".vlist.vari.scroll set"
   scrollbar .vlist.vari.scroll -command ".vlist.vari.list yview"

   pack .vlist.vari.list -side left -fill y
   pack .vlist.vari.scroll -side left -expand 1 -fill both
   pack .vlist.vari -side top -expand 1 -fill both

   bind .vlist.vari.list <Double-1> { get_vector_variable %W }

    frame .vlist.equ
    label .vlist.equ.label -text "math: "
    entry .vlist.equ.entry -width 20 -relief sunken -textvariable mathcmd
    pack .vlist.equ -side top
    pack .vlist.equ.label -side left
    pack .vlist.equ.entry -side left -fill x

    bind .vlist.equ.entry <Return> { math $mathcmd; set vector "tryagain" }

    frame .vlist.close
    button .vlist.close.ok  -text "OK" -command { set vector [get_vector_variable .vlist.vari.list]  }
    button .vlist.close.cancel -text  "Cancel" -command { set vector $savevector; }

    pack .vlist.close -side top
    pack .vlist.close.ok -side right
    pack .vlist.close.cancel -side right

    set oldfocus [focus]
    set vector "tryagain"
    while { $vector == "tryagain" } {
        .vlist.vari.list delete 0 end

       .vlist.vari.list insert end none
        do i 0 [@ $NumberOfVectorVariables-1] {
            set val $VectorVariableNames($i)
            .vlist.vari.list insert end $val
        }

        vwait vector
    }
    focus $oldfocus

    destroy .vlist

    return $vector
}

proc vector_edit { } {

    global VectorLines VectorLineStyle VectorQuality VectorRadius
    global VectorColor VectorArrow VectorLength VectorLengthScale
    global VectorThreshold VectorFloor VectorCeiling

    global vceil vfloor

    if { [winfo exists .vector] } {
        wm iconify .vector
        wm deiconify .vector
        return
    }

    toplevel .vector
    place_window .vector

#
# length style
#
    frame .vector.scale
    label .vector.scale.label -text "Vector Length Scale: "
    slider .vector.scale.slider -relief raised -bd 2 -orient horizontal \
            -from 0.0 -to 10.0 -resol 0.02 -digits 4 -variable VectorLengthScale

    pack .vector.scale -side top
    pack .vector.scale.label -side left
    pack .vector.scale.slider -side left -fill x

#
# line style
#
    frame .vector.line
    label .vector.line.label -text "Line Style: "
    radiobutton .vector.line.line -value 0 -variable VectorLineStyle -text "Line"
    radiobutton .vector.line.cyli -value 1 -variable VectorLineStyle -text "Solid"

    pack .vector.line -side top
    pack .vector.line.label -side left
    pack .vector.line.line -side left -fill x
    pack .vector.line.cyli -side left  -fill x

#
# cyl qual
#
    frame .vector.qual
    label .vector.qual.label -text "Line Quality: "
    entry .vector.qual.entry -relief sunken -width 5 -textvariable VectorQuality

    pack .vector.qual -side top
    pack .vector.qual.label -side left
    pack .vector.qual.entry -side left -fill x

#
# cyl rad
#
    frame .vector.radi
    label .vector.radi.label -text "Width Scale: "
    entry .vector.radi.entry -relief sunken -width 5 -textvariable VectorRadius

    pack .vector.radi -side top
    pack .vector.radi.label -side left
    pack .vector.radi.entry -side left -fill x

#
# threshold
#
    frame .vector.thres
    label .vector.thres.label -text "Threshold Variable: "
    button .vector.thres.but -textvariable VectorThreshold \
        -command { set VectorThreshold [make_scalar_list]; \
                       UpdateVariable "VectorThreshold";   \
             vector_set_floor $vfloor; vector_set_ceil $vceil  }

    UpdateVariable "VectorThreshold"
    vector_set_floor $vfloor;
    vector_set_ceil $vceil

    pack .vector.thres -side top
    pack .vector.thres.label -side left -fill x
    pack .vector.thres.but -side left -fill x

    frame .vector.floor
    label .vector.floor.floorlab -text "Min: "
    entry .vector.floor.floor -relief sunken -width 12 -textvariable VectorFloor
    bind .vector.floor.floor <Return> { vector_set_vfloor }

    slider .vector.floor.slider -relief raised -bd 2 -orient horizontal \
        -from 0.0 -to 100.0 -resol 0.5 -digits 4 \
             -variable vfloor -command { vector_set_floor }

    frame .vector.ceil
    label .vector.ceil.ceillab -text "Max: "
    entry .vector.ceil.ceil -relief sunken -width 12 -textvariable VectorCeiling
    bind .vector.ceil.ceil  <Return> { vector_set_vceil }

    slider .vector.ceil.slider -relief raised -bd 2 -orient horizontal \
         -from 0.0 -to 100.0 -resol 0.5 -digits 4 \
              -variable vceil -command { vector_set_ceil  }

    pack .vector.floor -side top
    pack .vector.floor.floorlab -side left -fill x
    pack .vector.floor.floor -side left -fill x
    pack .vector.floor.slider -side left

    pack .vector.ceil -side top
    pack .vector.ceil.ceillab -side left -fill x
    pack .vector.ceil.ceil -side left -fill x
    pack .vector.ceil.slider -side left

#
# color
#
    frame .vector.color
    label .vector.color.label -text "Color Variable: "
    button .vector.color.but -textvariable VectorColor \
           -command { set VectorColor [make_scalar_list]; UpdateVariable "VectorColor" }
    UpdateVariable "VectorColor"

    pack .vector.color -side top
    pack .vector.color.label -side left
    pack .vector.color.but -side left -fill x

#
# length
#
    frame .vector.length
    label .vector.length.label -text "Length Variable: "
    button .vector.length.but -textvariable VectorLength \
              -command { set VectorLength [make_scalar_list]; UpdateVariable "VectorLength" }
    UpdateVariable "VectorLength"

    pack .vector.length -side top
    pack .vector.length.label -side left
    pack .vector.length.but -side left -fill x

#
# arrow
#
    frame .vector.arrow
    label .vector.arrow.label -text "Arrow Variable: "
    button .vector.arrow.but -textvariable VectorArrow \
               -command { set VectorArrow [make_vector_list]; }

    pack .vector.arrow -side top
    pack .vector.arrow.label -side left
    pack .vector.arrow.but -side left -fill x

#
# close & apply
#
    frame .vector.buttons
    button .vector.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .vector.buttons.close -text "Close" -command "destroy .vector"

    pack .vector.buttons -side top
    pack .vector.buttons.apply -side left
    pack .vector.buttons.close -side left -fill x
}
