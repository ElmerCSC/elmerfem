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
#  Dialog box widget .
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

proc dl_dialog {title text bitmap default args} {
    global button

    toplevel .dl_dialog -class dialog

    place_window .dl_dialog
    wm title .dl_dialog $title
    wm iconname .dl_dialog Dialog

    frame .dl_dialog.top -relief raised -bd 1
    pack .dl_dialog.top -side top -fill both

    frame .dl_dialog.bot -relief raised -bd 1
    pack .dl_dialog.bot -side bottom -fill both

    message .dl_dialog.top.msg -width 3i -text $text -font -Adobe-Times-Medium-R-Normal-*-140-*

    if { $bitmap != "" } {
        pack .dl_dialog.top.msg -side right -expand 1 -fill both -padx 5m -pady 5m

        label .dl_dialog.top.bitmap -bitmap $bitmap
        pack .dl_dialog.top.bitmap -side left -padx 5m -pady 5m
    } else {
        pack .dl_dialog.top.msg -expand 1 -fill both -padx 5m -pady 5m
    } 

    set i 0
    foreach but $args {
        if { $i == $default } {
            button .dl_dialog.bot.button$i -bg gray -text $but -command "set button $i"
            pack .dl_dialog.bot.button$i -side left -expand 1 -padx 5m -pady 5m -ipadx 2m -ipady 1m
        } else {
            button .dl_dialog.bot.button$i -text $but -command "set button $i"
            pack .dl_dialog.bot.button$i  -side left -expand 1 -padx 5m -pady 5m -ipadx 2m -ipady 1m
        }
        incr i
    }

    if { $default >= 0 } {
        bind .dl_dialog <Return> ".dl_dialog.bot.button$default flash; set button $default"
        bind .dl_dialog <Destroy> "set button $default"
    }

    set oldfocus [focus]
    grab .dl_dialog
    focus .dl_dialog

    tkwait variable button
    if { [winfo exists .dl_dialog] } { bind .dl_dialog <Destroy> {}; destroy .dl_dialog }
    focus $oldfocus

    return $button
}

proc dl_dialog_int_set { window inc min max } {
   global dl_dialog_ret

   set val [@ [$window get]+$inc];

   if { $val < $min } { set val $min };
   if { $val > $max } { set val $max };

   set dl_dialog_ret $val
}

#
# Get an integer value from user
# 
# 18 Sep 1995
#
proc dl_dialog_int { title text {default 0} {min 0} {max 512} } {
    global dl_dialog_ret dl_dialog_dummy

    set win .dl_dialog_int_win

    toplevel $win -class dialog
    place_window $win

    wm minsize $win 300 200
    wm title $win "Give an Integer"

    frame $win.top -relief raised -bd 2
    frame $win.mid -relief raised -bd 2
    frame $win.bot -relief raised -bd 2

    message $win.top.msg -text $text -font -Adobe-Times-Medium-R-Normal-*-180-* -width 300

    entry $win.mid.entry -relief sunken -textvariable dl_dialog_ret -width 6

    button $win.mid.down_slow  -text "<"  -command "dl_dialog_int_set $win.mid.entry -1 $min $max"
    button $win.mid.down_fast  -text "<<" -command "dl_dialog_int_set $win.mid.entry -5 $min $max"

    button $win.mid.up_slow    -text ">"  -command "dl_dialog_int_set $win.mid.entry +1 $min $max"
    button $win.mid.up_fast    -text ">>" -command "dl_dialog_int_set $win.mid.entry +5 $min $max"

    pack $win.top.msg

    pack $win.mid.down_slow -side left
    pack $win.mid.down_fast -side left

    pack $win.mid.entry -side left

    pack $win.mid.up_fast -side left
    pack $win.mid.up_slow -side left

    button $win.bot.cancel -text CANCEL -command "set dl_dialog_ret {}; destroy $win"
    button $win.bot.ok -text OK -command "destroy $win"

    pack $win.bot.cancel -side left -expand 1 -padx 5m -pady 5m -ipadx 2m -ipady 1m
    pack $win.bot.ok -side left -expand 1 -padx 5m -pady 5m -ipadx 2m -ipady 1m

    pack $win.bot -side bottom
    pack $win.mid -side bottom
    pack $win.top -side bottom -fill both -expand 1
     
    set dl_dialog_ret $default

    bind $win <Return> "destroy $win"

    set oldfocus [focus]

    grab  $win
    focus $win

    tkwait window $win

    if { [winfo exists $win] } {
        bind $win <Destroy> {};
        destroy $win
    }
    focus $oldfocus

    return $dl_dialog_ret
}
