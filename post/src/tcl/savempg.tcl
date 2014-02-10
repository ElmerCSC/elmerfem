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
#*    MPG save utility routines
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
#*******************************************************************************
set savempg_start_stop_state "Start"
set savempg_append_state "-"
set savempg_bitrate 1000000
set savempg_method 0
set savempg_bframes 2
set savempg_gopsize 20
set MPGFileName "elmerpost.es"
#
# Helper procs:
#
proc savempg_start_stop {  } {
    global savempg_start_stop_state savempg_append_state \
	savempg_bitrate MPGFileName savempg_method savempg_bframes \
	savempg_gopsize

    if { $savempg_start_stop_state == "Start" } {

	if { $savempg_method == 0 } {
	    savempg codec mpg1
	}
	if { $savempg_method == 1 } {
	    savempg codec mpg2
	}
	if { $savempg_method == 2 } {
	    savempg codec mpg4
	}
	if { $savempg_method == 3 } {
	    savempg codec yuv4
	}

	savempg bitrate $savempg_bitrate

	savempg gopsize $savempg_gopsize

	savempg bframes $savempg_bframes

	savempg start $MPGFileName

	set savempg_start_stop_state "Stop"
        set savempg_append_state "Append"
    } else {
	savempg stop
	set savempg_start_stop_state "Start"
        set savempg_append_state "-"
    }
}

proc savempg_append { } {
    global savempg_start_stop_state

    if { $savempg_start_stop_state == "Stop" } {
	savempg append
    }
}

proc savempg_close { } {
    global savempg_start_stop_state savempg_append_state

    if { $savempg_start_stop_state == "Stop" } {
	savempg stop
	set savempg_start_stop_state "Start"
	set savempg_append_state "-"
    }
}
#
# Main control:
#
proc savempg.Control { } {
    global savempg_control savempg_start_stop_state savempg_append_state \
	savempg_bitrate MPGFileName savempg_method savempg_bframes \
	savempg_gopsize

    set savempg_control .savempg_control
    
    if { [winfo exists $savempg_control] } {
	destroy $savempg_control.title
	destroy $savempg_control.bitrate 
	destroy $savempg_control.file
	destroy $savempg_control.start_stop_button
	destroy $savempg_control.append_button
	destroy $savempg_control.buttons
	destroy $savempg_control.method_buttons
	destroy $savempg_control.gop
	destroy $savempg_control.bf
    } else {
	toplevel $savempg_control
	place_window $savempg_control
    }

    wm title $savempg_control "Savempg control"
    #
    # Codec:
    #
    frame $savempg_control.method_buttons
    label $savempg_control.method_buttons.label -width 8 -text "Codec:"
    radiobutton $savempg_control.method_buttons.mpg1 -value 0 \
	-variable savempg_method -text "MPG1"
    radiobutton $savempg_control.method_buttons.mpg2 -value 1 \
	-variable savempg_method -text "MPG2"
    radiobutton $savempg_control.method_buttons.mpg4 -value 2 \
	-variable savempg_method -text "MPG4"
    radiobutton $savempg_control.method_buttons.yuv4 -value 3 \
	-variable savempg_method -text "YUV4"
    pack $savempg_control.method_buttons.label \
	$savempg_control.method_buttons.mpg1 \
	$savempg_control.method_buttons.mpg2 \
	$savempg_control.method_buttons.mpg4 \
	$savempg_control.method_buttons.yuv4 \
	-side left -fill x

    pack $savempg_control.method_buttons -side top -expand 1 -fill x
    #
    #   Bitrate control:
    #
    frame $savempg_control.bitrate

    label $savempg_control.bitrate.label1 -width 8 -text "Bitrate:"
    entry $savempg_control.bitrate.value -width 10 \
	-textvariable savempg_bitrate
    label $savempg_control.bitrate.label2 \
	-text "bps (25 fps)"
    pack $savempg_control.bitrate.label1 $savempg_control.bitrate.value \
	$savempg_control.bitrate.label2 -side left 
    
    pack $savempg_control.bitrate -expand 1 -fill both -side top
    #
    #   GOP-size:
    #
    frame $savempg_control.gop

    label $savempg_control.gop.label1 -width 8 -text "GOP-size:"
    entry $savempg_control.gop.value1 -width 10 \
	-textvariable savempg_gopsize
    label $savempg_control.gop.label2 -text "frames max."

    pack $savempg_control.gop.label1 $savempg_control.gop.value1 \
	$savempg_control.gop.label2 -side left 
    
    pack $savempg_control.gop -fill both -side top
    #
    #   B-frames:
    #
    frame $savempg_control.bf

    label $savempg_control.bf.label1 -width 8 -text "B-size:"
    entry $savempg_control.bf.value1 -width 10 \
	-textvariable savempg_bframes
    label $savempg_control.bf.label2 -text "frames max."

    pack $savempg_control.bf.label1 $savempg_control.bf.value1 \
	$savempg_control.bf.label2 -side left 
    
    pack $savempg_control.bf -fill both -side top
    #
    # File name:
    #
    frame $savempg_control.file
    
    label $savempg_control.file.label -width 8 -text "File name:"
    entry $savempg_control.file.name -width 32 -textvariable MPGFileName
    button $savempg_control.file.button -text "Browse.." \
	-command { set MPGFileName [tk_getSaveFile -parent .savempg_control \
					-title "Save Video Clip To File"]; }
    pack $savempg_control.file.label $savempg_control.file.name \
	$savempg_control.file.button -side left -expand 1

    pack $savempg_control.file -expand 1 -fill both -side top
    #
    #   Buttons:
    #
    frame $savempg_control.buttons
    #
    # Close button:
    #
    button $savempg_control.buttons.close -text "Close" \
	-command { savempg_close; destroy $savempg_control }
    #
    # Start/stop:
    #
    frame $savempg_control.start_stop_button
    button $savempg_control.start_stop_button.start_stop \
	-textvariable savempg_start_stop_state  \
	-command { savempg_start_stop }
    pack $savempg_control.start_stop_button.start_stop \
	-side left -expand 1 -fill x
    #
    # Append:
    #
    frame $savempg_control.append_button
    button $savempg_control.append_button.append \
	-textvariable savempg_append_state \
	-command savempg_append
    pack $savempg_control.append_button.append -side left -expand 1 -fill x

    #
    # Pack:
    #
    pack $savempg_control.buttons.close -side right


    pack $savempg_control.start_stop_button -side top -expand 1 -fill x
    pack $savempg_control.append_button -side top -expand 1 -fill x
    pack $savempg_control.buttons -side top -fill x -expand 1
}
