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
#* Text display utility widget.
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
# User level command for setting text parametrs
#

set txt_x 0
set txt_y 0
set txt_z -1
set txt_str ""
set txt_fontname ""

#
# Create a text editor toplevel window for a user to use.
# There can be several editors running at any time provided 
# a different name for editor is given at each invokation.
#
# 12 Sep 1995
#
proc txt_Edit { name {UpdateProc text} } {
   global txt_x txt_y txt_z txt_str txt_fontname

   set txt_edit .txt_edit$name

   if { [winfo exists $txt_edit] } { destroy $txt_edit }

   toplevel $txt_edit
   place_window $txt_edit

   wm minsize $txt_edit 500 100
   wm title $txt_edit $name

   label $txt_edit.title -height 1;
   pack $txt_edit.title -side top -expand 1 -fill both
#
# Font name
#
   frame $txt_edit.font
   entry $txt_edit.font.str -relief sunken -width 40 -textvariable txt_fontname
   pack $txt_edit.font
   pack $txt_edit.font.str
#  
#
# List box containing font names
#
   frame $txt_edit.lbox -relief sunken -bg lightblue
   listbox $txt_edit.lbox.list -yscroll "$txt_edit.lbox.scroll set"
   scrollbar $txt_edit.lbox.scroll -command "$txt_edit.lbox.list yview"

   pack $txt_edit.lbox.scroll -side left -fill y
   pack $txt_edit.lbox.list -side left -expand 1 -fill both
   pack $txt_edit.lbox -side top -expand 1 -fill both

   bind $txt_edit.lbox.list <Double-1> "txt_SetListName %W; destroy $txt_edit"

#
# button box at window bottom (Apply,OK,Cancel)
#
   frame $txt_edit.bbox  -bg lightblue

   button $txt_edit.bbox.ok     -text "OK"     -command "destroy $txt_edit;"
   button $txt_edit.bbox.cancel -text "Cancel" -command "destroy $txt_edit"

   pack $txt_edit.bbox.ok  -side left -expand 1 -fill x
   pack $txt_edit.bbox.cancel  -side left -expand 1 -fill x

   pack $txt_edit.bbox -side top -expand 1 -fill both
#
#  add font names to listbox
#
   txt_AddFontNames $txt_edit.lbox.list

   set oldfocus [focus]

   grab  $txt_edit
   focus $txt_edit

   tkwait window $txt_edit

   if { [winfo exists $txt_edit] } {
       bind $txt_edit <Destroy> {};
       destroy $txt_edit
   }
   focus $oldfocus

   return $txt_fontname
}

#
# Get fontr name from listbox
#
# 12 Sep 1995
#
proc txt_SetListName {w} {
    global txt_fontname

    set txt_fontname [$w get [$w curselection]]
}

#
# Add font names to a listbox given as argument
#
# 13 Sep 1995
#
proc txt_AddFontNames {w} {
    global NumberOfFonts FontNames

    do i 0 [@ $NumberOfFonts-1] { $w insert end $FontNames($i) }
}
