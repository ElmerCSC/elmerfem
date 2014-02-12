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
#* Sphere display parameter settings
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
#
# 22 Apr 1996
#

set SphereStyle         0
set SphereLineStyle     0
set SphereQuality       1
set SphereRadiusScale   1
set SphereColor         "none"
set SphereRadius        "none"
set SphereThreshold     "none"
set SphereFloor         0.0
set SphereCeiling       0.0

set sceil  100.0
set sfloor 0.0

set SphereThresholdMin 0.0
set SphereThresholdMax 1.0

proc sphere_set_ceil { sceil } {
    global SphereCeiling SphereThresholdMin SphereThresholdMax

    set a [@ ($SphereThresholdMax-$SphereThresholdMin)*$sceil/100.0+$SphereThresholdMin]
    set SphereCeiling $a
}

proc sphere_set_floor { sfloor } {
    global SphereFloor SphereThresholdMin SphereThresholdMax

    set a [@ ($SphereThresholdMax-$SphereThresholdMin)*$sfloor/100.0+$SphereThresholdMin]
    set SphereFloor $a
}

proc sphere_set_sceil { } { 
    global SphereCeiling SphereThresholdMin SphereThresholdMax sceil

    set sceil [@ 100*($SphereCeiling-$SphereThresholdMin)/($SphereThresholdMax-$SphereThresholdMin)]
    set sceil [@ $sceil<0.0?0.0:$sceil]
    set sceil [@ $sceil>100.0?100.0:$sceil]
}

proc sphere_set_sfloor { } { 
    global SphereFloor SphereThresholdMin SphereThresholdMax sfloor

    set sfloor [@ 100*($SphereFloor-$SphereThresholdMin)/($SphereThresholdMax-$SphereThresholdMin)]
    set sfloor [@ $sfloor<0.0?0.0:$sfloor]
    set sfloor [@ $sfloor>100.0?100.0:$sfloor]
}


proc sphere_edit { } {
    global SphereStyle SphereLineStyle SphereQuality SphereRadius SphereColor SphereRadiusScale
    global SphereThreshold SphereCeiling SphereFloor

    global sceil sfloor

    if { [winfo exists .sphere] } {
        wm iconify .sphere
        wm deiconify .sphere
        return
    }

    toplevel .sphere
    place_window .sphere

#    frame .sphere.style
#    label .sphere.style.label -text "Sphere Style: "
#    radiobutton .sphere.style.line -value 0 -variable SphereStyle -text "Line"
#    radiobutton .sphere.style.surf -value 1 -variable SphereStyle -text "Surface"
#    radiobutton .sphere.style.both -value 2 -variable SphereStyle -text "Both"
#
#    pack .sphere.style -side top
#    pack .sphere.style.label -side left
#    pack .sphere.style.line -side left -fill x
#    pack .sphere.style.surf -side left  -fill x
#    pack .sphere.style.both -side left  -fill x
#
#    frame .sphere.line
#    label .sphere.line.label -text "Line Style: "
#    radiobutton .sphere.line.line -value 0 -variable SphereLineStyle -text "Line"
#    radiobutton .sphere.line.cyli -value 1 -variable SphereLineStyle -text "Solid"
#
#    pack .sphere.line -side top
#    pack .sphere.line.label -side left
#    pack .sphere.line.line -side left -fill x
#    pack .sphere.line.cyli -side left  -fill x
#
#    frame .sphere.edge
#    label .sphere.edge.label -text "Edge Style: "
#    radiobutton .sphere.edge.all  -value 0 -variable SphereEdgeStyle -text "All"
#    radiobutton .sphere.edge.free -value 1 -variable SphereEdgeStyle -text "Free"
#
#    pack .sphere.edge -side top
#    pack .sphere.edge.label -side left
#    pack .sphere.edge.all -side left -fill x
#    pack .sphere.edge.free -side left  -fill x
#
    frame .sphere.qual
    label .sphere.qual.label -text "Sphere Quality: "
    entry .sphere.qual.entry -relief sunken -width 5 -textvariable SphereQuality

    pack .sphere.qual -side top
    pack .sphere.qual.label -side left
    pack .sphere.qual.entry -side left -fill x

    frame .sphere.radi
    label .sphere.radi.label -text "Radius Scale: "
    entry .sphere.radi.entry -relief sunken -width 5 -textvariable SphereRadiusScale

    pack .sphere.radi -side top
    pack .sphere.radi.label -side left
    pack .sphere.radi.entry -side left -fill x

#
# sphere thresholding
#
    frame .sphere.thres
    label .sphere.thres.label -text "Threshold Variable: "
    button .sphere.thres.but -textvariable SphereThreshold \
        -command { set SphereThreshold [make_scalar_list]; \
                       UpdateVariable "SphereThreshold";   \
             sphere_set_floor $sfloor; sphere_set_ceil $sceil  }

    UpdateVariable "SphereThreshold"
    sphere_set_floor $sfloor;
    sphere_set_ceil $sceil

    pack .sphere.thres -side top
    pack .sphere.thres.label -side left -fill x
    pack .sphere.thres.but -side left -fill x

    frame .sphere.floor
    label .sphere.floor.floorlab -text "Min: "
    entry .sphere.floor.floor -relief sunken -width 12 -textvariable SphereFloor
    bind .sphere.floor.floor <Return> { sphere_set_sfloor }

    slider .sphere.floor.slider -relief raised -bd 2 -orient horizontal \
        -from 0.0 -to 100.0 -resol 0.5 -digits 4 \
             -variable sfloor -command { sphere_set_floor }

    frame .sphere.ceil
    label .sphere.ceil.ceillab -text "Max: "
    entry .sphere.ceil.ceil -relief sunken -width 12 -textvariable SphereCeiling
    bind .sphere.ceil.ceil  <Return> { sphere_set_sceil }

    slider .sphere.ceil.slider -relief raised -bd 2 -orient horizontal \
         -from 0.0 -to 100.0 -resol 0.5 -digits 4 \
              -variable sceil -command { sphere_set_ceil  }

    pack .sphere.floor -side top
    pack .sphere.floor.floorlab -side left -fill x
    pack .sphere.floor.floor -side left -fill x
    pack .sphere.floor.slider -side left

    pack .sphere.ceil -side top
    pack .sphere.ceil.ceillab -side left -fill x
    pack .sphere.ceil.ceil -side left -fill x
    pack .sphere.ceil.slider -side left



#
# sphere color variable
#
    frame .sphere.colvari
    label .sphere.colvari.label -text "Color Variable: "
    button .sphere.colvari.but -textvariable SphereColor -command { set SphereColor [make_scalar_list] }

    pack .sphere.colvari -side top
    pack .sphere.colvari.label -side left
    pack .sphere.colvari.but -side left -fill x

#
# sphere radius variable
#
    frame .sphere.radvari
    label .sphere.radvari.label -text "Radius Variable: "
    button .sphere.radvari.but -textvariable SphereRadius -command { set SphereRadius [make_scalar_list] }

    pack .sphere.radvari -side top
    pack .sphere.radvari.label -side left
    pack .sphere.radvari.but -side left -fill x

#
# buttons
#
    frame .sphere.buttons
    button .sphere.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .sphere.buttons.close -text "Close" -command "destroy .sphere"

    pack .sphere.buttons -side top
    pack .sphere.buttons.apply -side left
    pack .sphere.buttons.close -side left -fill x
}
