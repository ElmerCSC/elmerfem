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

#*******************************************************************************
#*
#*    JPG save utility routines
#*
#*******************************************************************************
#*
#*                     Author:       Mikko Lyly
#*
#*                    Address: CSC - IT Center for Science Ltd.
#*                                Keilaranta 14, P.O. BOX 405
#*                                  02101 Espoo, Finland
#*                                  Tel. +358 0 457 2723
#*                                Telefax: +358 0 457 2302
#*                              EMail: Juha.Ruokolainen@csc.fi
#*
#*                       Date: 04 Oct 2007
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#******************************************************************************
set savejpg_quality 80
set JPGFileName "elmerpost.jpg"

proc savejpg.Control { } {
    global savejpg_control savejpg_quality JPGFileName

    set savejpg_control .savejpg_control
    
    if { [winfo exists $savejpg_control] } {
	destroy $savejpg_control.title
	destroy $savejpg_control.quality
	destroy $savejpg_control.file
	destroy $savejpg_control.save_button
	destroy $savejpg_control.buttons
    } else {
	toplevel $savejpg_control
	place_window $savejpg_control
    }

    wm title $savejpg_control "Savejpg control"
    #
    #   Quality control:
    #
    frame $savejpg_control.quality

    label $savejpg_control.quality.label1 -width 8 -text "Quality:   "
    entry $savejpg_control.quality.value -width 5 \
	-textvariable savejpg_quality
    label $savejpg_control.quality.label2 -text "(1=low ... 100=best)"
    pack $savejpg_control.quality.label1 $savejpg_control.quality.value \
	$savejpg_control.quality.label2 -side left 
    
    pack $savejpg_control.quality -expand 1 -fill both -side top
    #
    # File name:
    #
    frame $savejpg_control.file
    
    label $savejpg_control.file.label -width 8 -text "File name:"
    entry $savejpg_control.file.name -width 30 -textvariable JPGFileName
    button $savejpg_control.file.button -text "Browse.." \
	-command { set JPGFileName [tk_getSaveFile -parent .savejpg_control \
					-title "Save Picture To File"]; }
    pack $savejpg_control.file.label $savejpg_control.file.name \
	$savejpg_control.file.button -side left -expand 1

    pack $savejpg_control.file -expand 1 -fill both -side top
    #
    #   Buttons:
    #
    frame $savejpg_control.buttons

    button $savejpg_control.buttons.close -text "Close" \
	-command { destroy $savejpg_control }

    frame $savejpg_control.save_button
    button $savejpg_control.save_button.save \
	-text "Save" -command { savejpg $JPGFileName $savejpg_quality }
    pack $savejpg_control.save_button.save \
	-side left -expand 1 -fill x

    pack $savejpg_control.buttons.close -side right

    pack $savejpg_control.save_button -side top -expand 1 -fill x
    pack $savejpg_control.buttons -side top -fill x -expand 1
}
