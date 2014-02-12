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
#
#/*******************************************************************************
# *
# *     Model coordinate clip plane editor.
# *
# *******************************************************************************
# *
# *                     Author:       Juha Ruokolainen
# *
# *                    Address: CSC - IT Center for Science Ltd.
# *                                Keilaranta 14, P.O. BOX 405
# *                                  02101 Espoo, Finland
# *                                  Tel. +358 0 457 2723
# *                                Telefax: +358 0 457 2302
# *                              EMail: Juha.Ruokolainen@csc.fi
# *
# *                       Date: 27 Sep 1995
# *
# *                Modified by:
# *
# *       Date of modification:
# *
# ******************************************************************************/


set ClipLowX  -1.25
set ClipHighX  1.25
set ClipLowY  -1.25
set ClipHighY  1.25
set ClipLowZ  -1.25
set ClipHighZ  1.25

proc clip_Edit {} {

  global ClipLowX ClipHighX ClipLowY ClipHighY ClipLowZ ClipHighZ

  if { [winfo exists .clip_edit] } {
    wm iconify .clip_edit
	wm deiconify .clip_edit
  } else {
    toplevel .clip_edit
    place_window .clip_edit

	frame  .clip_edit.xlow -relief sunken -bd 2
	label  .clip_edit.xlow.lab -text "Low X Plane: "
	entry  .clip_edit.xlow.ent -width 7 -textvariable ClipLowX
	slider .clip_edit.xlow.sld -from -1.25 -to 1.25 -variable ClipLowX -resol 0.005 -orient horizontal
    pack .clip_edit.xlow.lab .clip_edit.xlow.ent .clip_edit.xlow.sld -side left
	pack .clip_edit.xlow -side top

	frame  .clip_edit.xhigh -relief sunken -bd 2
	label  .clip_edit.xhigh.lab -text "High X Plane: "
	entry  .clip_edit.xhigh.ent -width 7 -textvariable ClipHighX
	slider .clip_edit.xhigh.sld -from -1.25 -to 1.25 -variable ClipHighX -resol 0.005 -orient horizontal
    pack .clip_edit.xhigh.lab .clip_edit.xhigh.ent .clip_edit.xhigh.sld -side left
	pack .clip_edit.xhigh -side top


	frame  .clip_edit.ylow -relief sunken -bd 2
	label  .clip_edit.ylow.lab -text "Low Y Plane: "
	entry  .clip_edit.ylow.ent -width 7 -textvariable ClipLowY
	slider .clip_edit.ylow.sld -from -1.25 -to 1.25 -variable ClipLowY -resol 0.005 -orient horizontal
    pack .clip_edit.ylow.lab .clip_edit.ylow.ent .clip_edit.ylow.sld -side left
	pack .clip_edit.ylow -side top

	frame  .clip_edit.yhigh -relief sunken -bd 2
	label  .clip_edit.yhigh.lab -text "High Y Plane: "
	entry  .clip_edit.yhigh.ent -width 7 -textvariable ClipHighY
	slider .clip_edit.yhigh.sld -from -1.25 -to 1.25 -variable ClipHighY -resol 0.005 -orient horizontal
    pack .clip_edit.yhigh.lab .clip_edit.yhigh.ent .clip_edit.yhigh.sld -side left
	pack .clip_edit.yhigh -side top


	frame  .clip_edit.zlow -relief sunken -bd 2
	label  .clip_edit.zlow.lab -text "Low Z Plane: "
	entry  .clip_edit.zlow.ent -width 7 -textvariable ClipLowZ
	slider .clip_edit.zlow.sld -from -1.25 -to 1.25 -variable ClipLowZ -resol 0.005 -orient horizontal
    pack .clip_edit.zlow.lab .clip_edit.zlow.ent .clip_edit.zlow.sld -side left
	pack .clip_edit.zlow -side top

	frame  .clip_edit.zhigh -relief sunken -bd 2
	label  .clip_edit.zhigh.lab -text "High Z Plane: "
	entry  .clip_edit.zhigh.ent -width 7 -textvariable ClipHighZ
	slider .clip_edit.zhigh.sld -from -1.25 -to 1.25 -variable ClipHighZ -resol 0.005 -orient horizontal
    pack .clip_edit.zhigh.lab .clip_edit.zhigh.ent .clip_edit.zhigh.sld -side left
	pack .clip_edit.zhigh -side top

	frame .clip_edit.buttons
	button .clip_edit.buttons.apply -text Apply -command clip_apply
	button .clip_edit.buttons.ok -text OK -command "clip_apply; destroy .clip_edit"
	button .clip_edit.buttons.close -text Close -command "destroy .clip_edit"
	pack .clip_edit.buttons.apply .clip_edit.buttons.ok .clip_edit.buttons.close -side left
    pack .clip_edit.buttons -side top


    bind .clip_edit <Return> clip_apply


  }
}


proc clip_apply { args } {
    global ClipLowX ClipLowY ClipLowZ ClipHighX ClipHighY ClipHighZ

	clip 0  1  0  0 [@ 0-$ClipLowX]
	clip 1 -1  0  0 $ClipHighX
    clip 2  0  1  0 [@ 0-$ClipLowY]
	clip 3  0 -1  0 $ClipHighY
	clip 4  0  0 -1 [@ 0-$ClipLowZ]
    clip 5  0  0 1 $ClipHighZ
	display
}
