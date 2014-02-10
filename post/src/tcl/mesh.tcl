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
#* Mesh display parameter settings
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

set MeshStyle     0
set MeshLineStyle 0
set MeshEdgeStyle 1
set MeshQuality   1
set MeshRadius    1
set MeshColor     "none"
set MeshColorMin  0.0
set MeshColorMax  1.0
set MeshColorSetMinMax 0

set ScalarVariableNames(0) none
set NumberOfScalarVariables 1

proc mesh_update {} {

   global MeshColor MeshColorMin MeshColorMax

   UpdateVariable MeshColor

   .mesh.set.min delete 0 end
   .mesh.set.min insert end [format %-10.5g $MeshColorMin]

   .mesh.set.max delete 0 end
   .mesh.set.max insert end [format %-10.5g $MeshColorMax]
}

proc mesh_edit { } {

    global MeshStyle MeshLineStyle MeshEdgeStyle MeshQuality MeshRadius
    global MeshColor MeshColorMin MeshColorMax MeshColorSetMinMax

    if { [winfo exists .mesh] } {
        wm iconify .mesh
        wm deiconify .mesh
        return
    }

    toplevel .mesh
    wm title .mesh "Color Mesh Edit"

    place_window .mesh

    frame .mesh.style
    label .mesh.style.label -text "Mesh Style: "
    radiobutton .mesh.style.line -value 0 -variable MeshStyle -text "Line"
    radiobutton .mesh.style.surf -value 1 -variable MeshStyle -text "Surface"
    radiobutton .mesh.style.both -value 2 -variable MeshStyle -text "Both"

    pack .mesh.style -side top
    pack .mesh.style.label -side left
    pack .mesh.style.line -side left -fill x
    pack .mesh.style.surf -side left  -fill x
    pack .mesh.style.both -side left  -fill x


    frame .mesh.line
    label .mesh.line.label -text "Line Style: "
    radiobutton .mesh.line.line -value 0 -variable MeshLineStyle -text "Line"
    radiobutton .mesh.line.cyli -value 1 -variable MeshLineStyle -text "Solid"

    pack .mesh.line -side top
    pack .mesh.line.label -side left
    pack .mesh.line.line -side left -fill x
    pack .mesh.line.cyli -side left  -fill x

    frame .mesh.edge
    label .mesh.edge.label -text "Edge Style: "
    radiobutton .mesh.edge.all  -value 0 -variable MeshEdgeStyle -text "All"
    radiobutton .mesh.edge.free -value 1 -variable MeshEdgeStyle -text "Free"

    pack .mesh.edge -side top
    pack .mesh.edge.label -side left
    pack .mesh.edge.all -side left -fill x
    pack .mesh.edge.free -side left  -fill x

    frame .mesh.qual
    label .mesh.qual.label -text "Line Quality: "
    entry .mesh.qual.entry -relief sunken -width 5 -textvariable MeshQuality

    pack .mesh.qual -side top
    pack .mesh.qual.label -side left
    pack .mesh.qual.entry -side left -fill x

    frame .mesh.radi
    label .mesh.radi.label -text "Width Scale: "
    entry .mesh.radi.entry -relief sunken -width 5 -textvariable MeshRadius

    pack .mesh.radi -side top
    pack .mesh.radi.label -side left
    pack .mesh.radi.entry -side left -fill x
#
# mesh color
#
    frame .mesh.vari
    label .mesh.vari.label -text "Color Variable: "
    button .mesh.vari.but -textvariable MeshColor -command { set MeshColor [make_scalar_list]; mesh_update; }

    pack .mesh.vari -side top
    pack .mesh.vari.label -side left
    pack .mesh.vari.but -side left -fill x
#
# Generate ...
#
    frame .mesh.set

    label .mesh.set.min_lab -text "Min: "
    entry .mesh.set.min -width 10 -textvariable MeshColorMin
    bind .mesh.set.min <Return> mesh_update

    label .mesh.set.max_lab -text "Max: "
    entry .mesh.set.max -width 10 -textvariable MeshColorMax
    bind .mesh.set.max <Return> mesh_update

    checkbutton .mesh.set.keep -text "Keep" -variable MeshColorSetMinMax \
          -command mesh_update

    pack .mesh.set.min_lab -side left
    pack .mesh.set.min -side left
    pack .mesh.set.max_lab -side left
    pack .mesh.set.max -side left
    pack .mesh.set.keep -side left
    pack .mesh.set -side top
    
#
# buttons
#
    frame .mesh.buttons
    button .mesh.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .mesh.buttons.close -text "Close" -command "destroy .mesh"

    pack .mesh.buttons -side top
    pack .mesh.buttons.apply -side left
    pack .mesh.buttons.close -side left -fill x
}
