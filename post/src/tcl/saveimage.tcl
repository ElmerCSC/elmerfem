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
#*     Model file read utility routines
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

set OUTPS      0
set PSFileName ""

proc saveimage { File outps }   {
  global GlobalOptions PSFileName

  if { $outps == 1 } {
      puts "ps"
      set PSFileName $File
      set GlobalOptions(OutputPS) 1
      display
      set GlobalOptions(OutputPS) 0
      display
  } else {
      if { $outps == 0 } {
	  puts "ppm"
	  screensave $File
      } else {
	  puts "jpg"
	  savejpg $File
      }
  }
}


proc SaveImage {} {
    global PSFileName OutputPS OUTPS

    set w .screensave

    if { [winfo exists $w] } {
      wm iconify $w
      wm deiconify $w
      return
    } else {
      toplevel $w
      place_window $w
      wm title $w "Save Screen"
    }

#
# File name
#
    label $w.sp1 -text "Save as:"
    radiobutton $w.rb2 -value 0 -variable "OUTPS" -text "PPM Image"
    checkbutton $w.rb3  -variable "GlobalOptions(FitToPagePS)" -text "Fit PS to page"
    radiobutton $w.rb1 -value 1 -variable "OUTPS" -text "Postscript"
    radiobutton $w.rb4 -value 2 -variable "OUTPS" -text "JPG Image"
    pack $w.sp1 $w.rb1 $w.rb3 $w.rb2 $w.rb4

#    pack $w.rb3 -side top

    label $w.sp3 -text "\nSelect file:\n" -font "Helvetica-Bold 12"
    pack $w.sp3 -side top

    frame $w.file -relief sunken
    label $w.file.lab -text "File Name: "
    entry $w.file.name -width 40 -textvariable PSFileName

    button $w.file.but -text "Browse..." -command \
       { set PSFileName [tk_getSaveFile -parent .screensave -title "Save Image To File"]; }

    pack $w.file.lab $w.file.name $w.file.but -side left
    pack $w.file -side top

    label $w.sp4 -text \n
    pack $w.sp4 -side top

    frame $w.but
    button  $w.but.exit -text "Close" -command "destroy $w"
    button  $w.but.save -text "Save" -command { saveimage $PSFileName $OUTPS }

    pack $w.but $w.but.exit $w.but.save -side right
}
