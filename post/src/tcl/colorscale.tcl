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
#* Colorscale display parameter settings
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

set ColorScaleStyle         1
set ColorScaleEntries       6
set ColorScaleDecimals      4
set ColorScaleX          -0.8
set ColorScaleY          -0.8
set ColorScaleThickness   0.1
set ColorScaleLength      1.5
set ColorScaleFont       "-adobe-helvetica-bold-r-normal--17-120-100-100-p-88-iso8859-1"
set ColorScaleFontColor  [@ 256*256*255+256*255+255]
set ColorScaleColor      "none"
set ColorScaleColorMin   0.0
set ColorScaleColorMax   1.0

proc colscale_update {} {

   global ColorScaleColor ColorScaleColorMin ColorScaleColorMax

   UpdateVariable ColorScaleColor

   .colscale.set.min delete 0 end
   .colscale.set.min insert end [format %-10.5g $ColorScaleColorMin]

   .colscale.set.max delete 0 end
   .colscale.set.max insert end [format %-10.5g $ColorScaleColorMax]
}

proc colscale_edit { } {
    global ColorScaleStyle ColorScaleEntries ColorScaleDecimals ColorScaleX ColorScaleY
    global ColorScaleThickness ColorScaleLength
    global ColorScaleFont ColorScaleFontColor ColorScaleColor

    global ColorScaleColorMin ColorScaleColorMax ColorScaleColorSetMinMax

    if { [winfo exists .colscale] } {
        wm iconify .colscale
        wm deiconify .colscale
        return
    }

    toplevel .colscale
    wm title .colscale "Color Scale Edit"
    place_window .colscale

    frame .colscale.style
    label .colscale.style.label -text "Colorscale Style: "
    radiobutton .colscale.style.vert -value 0 -variable ColorScaleStyle -text "Vertical"
    radiobutton .colscale.style.hori -value 1 -variable ColorScaleStyle -text "Horizontal"

    pack .colscale.style -side top
    pack .colscale.style.label -side left
    pack .colscale.style.vert -side left -fill x
    pack .colscale.style.hori -side left  -fill x


    frame .colscale.entr
    label .colscale.entr.label -text "Labels: "
    entry .colscale.entr.entry -relief sunken -width 3 -textvariable ColorScaleEntries

    label .colscale.entr.ldecim -text "Decimals: "
    entry .colscale.entr.edecim -relief sunken -width 3 -textvariable ColorScaleDecimals

    pack .colscale.entr -side top
    pack .colscale.entr.label -side left
    pack .colscale.entr.entry -side left -fill x

    pack .colscale.entr.ldecim -side left
    pack .colscale.entr.edecim -side left -fill x

#
# colscale color
#
    frame .colscale.vari
    label .colscale.vari.label -text "Color Variable: "
    button .colscale.vari.but -textvariable ColorScaleColor -command { set ColorScaleColor [make_scalar_list]; colscale_update }

    pack .colscale.vari -side top
    pack .colscale.vari.label -side left
    pack .colscale.vari.but -side left -fill x

#
# Generate ...
#
    frame .colscale.set

    label .colscale.set.min_lab -text "Min: "
    entry .colscale.set.min -width 10 -textvariable ColorScaleColorMin
    bind .colscale.set.min <Return> colscale_update

    label .colscale.set.max_lab -text "Max: "
    entry .colscale.set.max -width 10 -textvariable ColorScaleColorMax
    bind .colscale.set.max <Return> colscale_update

    checkbutton .colscale.set.keep -text "Keep" -variable ColorScaleColorSetMinMax \
          -command colscale_update

    pack .colscale.set.min_lab -side left
    pack .colscale.set.min -side left
    pack .colscale.set.max_lab -side left
    pack .colscale.set.max -side left
    pack .colscale.set.keep -side left
    pack .colscale.set -side top
#
#
#
    frame .colscale.x
    label .colscale.x.label -text "X Position: "
    entry .colscale.x.entry -relief sunken -width 5 -textvariable ColorScaleX

    frame .colscale.y
    label .colscale.y.label -text "Y Position: "
    entry .colscale.y.entry -relief sunken -width 5 -textvariable ColorScaleY

    frame .colscale.l
    label .colscale.l.label -text "Length: "
    entry .colscale.l.entry -relief sunken -width 5 -textvariable ColorScaleLength

    frame .colscale.t
    label .colscale.t.label -text "Thickness: "
    entry .colscale.t.entry -relief sunken -width 5 -textvariable ColorScaleThickness

    pack .colscale.x -side top
    pack .colscale.x.label -side left
    pack .colscale.x.entry -side left -fill x

    pack .colscale.y -side top
    pack .colscale.y.label -side left
    pack .colscale.y.entry -side left -fill x

    pack .colscale.l -side top
    pack .colscale.l.label -side left
    pack .colscale.l.entry -side left -fill x

    pack .colscale.t -side top
    pack .colscale.t.label -side left
    pack .colscale.t.entry -side left -fill x
#
#
#
      frame .colscale.font
      button .colscale.font.font  -text "Label Font"  -command { \
          setfont [txt_Edit "FontSelect"]; }
      button .colscale.font.color -text "Label Color" -command { \
         set value [tk_chooseColor -title "Choose Font Color" -parent .colscale.font.color -initialcolor white]; \
          scan $value "#%02x%02x%02x" R G B; set ColorScaleFontColor [@ 256*256*$R+256*$G+$B]; }
 
      pack .colscale.font -side top
      pack .colscale.font.font  -side left
      pack .colscale.font.color -side left -fill x

#
# buttons
#
    frame .colscale.buttons
    button .colscale.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .colscale.buttons.close -text "Close" -command "destroy .colscale"

    pack .colscale.buttons -side top
    pack .colscale.buttons.apply -side left
    pack .colscale.buttons.close -side left -fill x
}
