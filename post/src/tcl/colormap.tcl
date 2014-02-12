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
# Colormap Editor Utility Widget
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

set ci_marks 0
set ci_mark_n 0
set ci_applyto "-default"

#
# User level command to give colormap from a file
#
proc colormap { file_name {opt "-default"} } {
    set Usage "Usage: colormap file -vectors -mesh -contour -isosurface -sphere -particles"
    global colmap colmap_size

#   if { [llength $args] != 1 } { return -code error "colormap: Argument number mismatch.\n\n$Usage" }
#    set file_name $args

    set ci_edit .ci_editColormap

    if { [winfo exists $ci_edit] } {
        ci_ReadColorMap $ci_edit $file_name
        ci_SendColorMap $ci_edit

        return
    }

    ci_ReadColorMapFile $file_name color knot maxindex
    if { $maxindex <= 0 } return

    do i 0 [@ $maxindex-1] { set colmap($i) #000000 }
 
    if { [info exists color] } {
        set sid [array startsearch color]
        while { [array anymore color $sid] != 0 } {
 
            set k [array nextelement color $sid] 
            set colmap($k) $color($k)

        }
        array donesearch color $sid
    }

    if { [info exists knot] } {
        set box_list ""
 
        set sid [array startsearch knot]
        while { [array anymore knot $sid] != 0 } {
            set n [array nextelement knot $sid] 
 
            if { $box_list == ""  } {
                set box_list $n
            } else {
                set box_list "$box_list $n"
            }
        }
        array donesearch knot $sid

        set box_list [split $box_list]
        set box_list [lsort -integer $box_list]

        set n0 [lindex $box_list 0]
        set color0 $knot($n0)

        do i 1 [@ [llength $box_list]-1] {
             set n1 [lindex $box_list $i]

             if { $n0 == $n1 } continue;
             set color1 $knot($n1)

             set k 0
             set t [split [GetInterpolate $n0 $n1 $color0 $color1] " "]
             do j $n0 [@ $n1-1] { set colmap($j) [lindex $t $k]; incr k }

             set n0 $n1
             set color0 $color1
         }
         set colmap($n1) $color1

    }

    set rgblist ""
    do i 0 [@ $maxindex-1] { set rgblist $rgblist$colmap($i) }

    set colmap_size $maxindex
    GetColorMap $opt $maxindex $rgblist
    UpdateObject;
}

#
# Create a colormap editor, optionally given filename
#
# 17 Sep 1995
proc ci_ColorMap { name {file_name ""} {nboxes 0} } {
    global ELMER_POST_HOME

    set ci_edit .ci_edit$name

    if { $file_name == "" } {
        if { $nboxes > 0 } {
            ci_CreateColorMap $ci_edit $nboxes 1
        } else {
            ci_ReadColorMap $ci_edit $ELMER_POST_HOME/lib/colormaps/default.cm
        }
    } else {
        ci_ReadColorMap $ci_edit $file_name
    }
}

#
# Create a new clormap
#
# 18 Sep 1995
proc ci_NewColorMap { ci_edit } {
    global global_limits

    set n [dl_dialog_int {NumberOfEntries} {Give number of colormap entries?} 32 0 $global_limits(MAX_COLORMAP_ENTRIES)]

    if { $n == "" } return
 
    if { $n <= 0 || $n >= $global_limits(MAX_COLORMAP_ENTRIES) } {
        dl_dialog {Error} "Invalid number of colormap entries: \[$n\]" error 0 "OK"
        return
    }

    ci_CreateColorMap $ci_edit $n 1
}

#
# Create a colormap editor (or modify current)
#
# 17 Sep 1995
#
proc ci_CreateColorMap { ci_edit {nboxes 32} {initialize 1} } {

    global ci_apply ci_applyto ci_list ci_marks ci_size ci_isize ci_nbox_x ci_nbox_y me_rgb_names ELMER_POST_HOME

    set ci_nbox_x 32
    while { [@ $nboxes % $ci_nbox_x] != 0 } { set ci_nbox_x [@ $ci_nbox_x-1] }

    set ci_nbox_y [@ $nboxes/$ci_nbox_x]

    set maxbox $ci_nbox_x
    if { $ci_nbox_y > $ci_nbox_x } { set maxbox $ci_nbox_y }

    set ci_size [@ int(20.0*(1.0+(32-$maxbox)/32.0))];
    set ci_isize [@ 1.0/$ci_size]

    set ci_len_x  [@ $ci_size*$ci_nbox_x]
    set ci_len_y  [@ $ci_size*$ci_nbox_y]

    if { [winfo exists $ci_edit] } {
        destroy $ci_edit.bbox
        destroy $ci_edit.c
        destroy $ci_edit.red
        destroy $ci_edit.green
        destroy $ci_edit.blue
        destroy $ci_edit.lbox
        destroy $ci_edit.menubar
    } else { toplevel $ci_edit; place_window $ci_edit }

    wm minsize $ci_edit $ci_len_x $ci_len_y

    frame $ci_edit.menubar
    menubutton $ci_edit.menubar.file -menu $ci_edit.menubar.file.menu -text "File" -underline 0
    menu $ci_edit.menubar.file.menu

    $ci_edit.menubar.file.menu add command 
    $ci_edit.menubar.file.menu entryconfigure last -label "  New  " -command "ci_NewColorMap $ci_edit"

    $ci_edit.menubar.file.menu add command 
    $ci_edit.menubar.file.menu entryconfigure last -label "  Read  " -command "ci_ReadColorMap $ci_edit"

    $ci_edit.menubar.file.menu add command
    $ci_edit.menubar.file.menu entryconfigure last -label "  Save  " -command "ci_WriteColorMap $ci_edit"

    $ci_edit.menubar.file.menu add command
    $ci_edit.menubar.file.menu entryconfigure last -label "  Quit  " -command "destroy $ci_edit"

    menubutton $ci_edit.menubar.edit -menu $ci_edit.menubar.edit.menu -text "Edit" -underline 0
    menu $ci_edit.menubar.edit.menu

    $ci_edit.menubar.edit.menu add command
    $ci_edit.menubar.edit.menu entryconfigure last -label "  Delete current knot (Shift-Middle Mouse)  " -command "ci_DeleteCurrent $ci_edit"

    $ci_edit.menubar.edit.menu add command
    $ci_edit.menubar.edit.menu entryconfigure last -label "  Add knot (Middle Mouse)  " -command "ci_CreateMark 0 0 $ci_edit"

    $ci_edit.menubar.edit.menu add command
    $ci_edit.menubar.edit.menu entryconfigure last -label "  Reset  " -command "ci_ColorMap Colormap"

######################################################################################################

    menubutton $ci_edit.menubar.apply -menu $ci_edit.menubar.apply.menu -text "Apply-To" -underline 0
    menu $ci_edit.menubar.apply.menu

    set ci_apply 0
    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Default" -command { set ci_applyto "-default" } \
                                -value 0 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Mesh" -command { set ci_applyto "-mesh" } \
                                -value 1 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Contours" -command { set ci_applyto "-contour" } \
                                -value 2 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Vectors" -command { set ci_applyto "-vectors" } \
                                -value 3 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Isosurfaces" -command { set ci_applyto "-isosurface" } \
                                -value 4 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Spheres" -command { set ci_applyto "-sphere" } \
                                -value 5 -variable ci_apply

    $ci_edit.menubar.apply.menu add radiobutton
    $ci_edit.menubar.apply.menu entryconfigure last -label "Particles" -command { set ci_applyto "-particles" } \
                                -value 6 -variable ci_apply


######################################################################################################

    menubutton $ci_edit.menubar.help -menu $ci_edit.menubar.help.menu -text "Help" -underline 0
    menu $ci_edit.menubar.help.menu

    $ci_edit.menubar.help.menu add command
    $ci_edit.menubar.help.menu entryconfigure last -label "  I'm lost  " -command "epHelp edit.html#Colormap"

    $ci_edit.menubar.help.menu add command
    $ci_edit.menubar.help.menu entryconfigure last -label "  About This Program  " -command "puts {Juha is responsible...}"

    pack $ci_edit.menubar.file -side left
    pack $ci_edit.menubar.edit -side left
    pack $ci_edit.menubar.apply -side left
    pack $ci_edit.menubar.help -side right

    pack $ci_edit.menubar -side top -fill x

    canvas $ci_edit.c -width $ci_len_x -height $ci_len_y
    pack  $ci_edit.c -side top

    set k 0
    do i 0 [@ $ci_nbox_y-1] {
       do j 0 [@ $ci_nbox_x-1] {
           set x0 [@ $ci_size*$j]
           set x1 [@ $ci_size*($j+1)]

           set y0 [@ $ci_size*$i]
           set y1 [@ $ci_size*($i+1)]

           $ci_edit.c create rect $x0 $y0 $x1 $y1 -fill #000000 -outline white
           incr k
       }
    }

    frame $ci_edit.red
    frame $ci_edit.green
    frame $ci_edit.blue

    slider $ci_edit.red.scl -relief raised -bd 5 -troughcolor red -orient horizontal  -from 0 -to 100 -resol 0.5 -digit 4
    slider $ci_edit.green.scl -relief raised -bd 5 -troughcolor green -orient horizontal -from 0 -to 100 -resol 0.5 -digit 4
    slider $ci_edit.blue.scl -relief raised -bd 5 -troughcolor blue -orient horizontal -from 0 -to 100 -resol 0.5 -digit 4 -tick 25

    pack $ci_edit.red.scl -side left -expand 1 -fill x
    pack $ci_edit.green.scl -side left -expand 1 -fill x
    pack $ci_edit.blue.scl -side left -expand 1 -fill x

    pack $ci_edit.red -side top -fill x
    pack $ci_edit.green -side top -fill x
    pack $ci_edit.blue -side top  -fill x
    frame $ci_edit.bbox

    button $ci_edit.bbox.apply -text "Apply" -command "ci_SendColorMap $ci_edit;"
    button $ci_edit.bbox.cancel -text "Cancel" -command "destroy $ci_edit"
    button $ci_edit.bbox.ok -text "OK" -command "ci_SendColorMap $ci_edit; destroy $ci_edit"

    pack $ci_edit.bbox.apply  -side left -expand 1 -fill x
    pack $ci_edit.bbox.ok  -side left -expand 1 -fill x
    pack $ci_edit.bbox.cancel  -side left -expand 1 -fill x
    pack $ci_edit.bbox -side top -expand 1 -fill both

    if { $initialize != 0 } {
        set ci_marks 0
        if { [info exists ci_list] } { unset ci_list }

        $ci_edit.c itemconfigure 1 -fill #000000
        $ci_edit.c itemconfigure [@ $k] -fill #ffffff

        ci_CreateMark 0 0 $ci_edit
        set x [@ $ci_size*($ci_nbox_x-1)+$ci_size/2]
        set y [@ $ci_size*($ci_nbox_y-1)+$ci_size/2]
        ci_CreateMark $x $y $ci_edit

        ci_InterpolateColor $ci_edit
    }

    $ci_edit.red.scl configure -command "ci_ShowColor $ci_edit"
    $ci_edit.green.scl configure -command "ci_ShowColor $ci_edit"
    $ci_edit.blue.scl configure -command "ci_ShowColor $ci_edit"

    frame $ci_edit.lbox -relief sunken -bg lightblue
    listbox $ci_edit.lbox.list -yscroll "$ci_edit.lbox.scroll set"
    scrollbar $ci_edit.lbox.scroll -command "$ci_edit.lbox.list yview"

    pack $ci_edit.lbox.scroll -side left -fill y
    pack $ci_edit.lbox.list -side left -expand 1 -fill both
    pack $ci_edit.lbox -side top -expand 1 -fill both

    bind $ci_edit.lbox.list <Double-1> {ci_SetListColor %W}

    if { ![info exists me_rgb_list] } me_GetColorNames
    me_AddColorNames $ci_edit.lbox.list

    bind $ci_edit.c <Button-2> { ci_CreateMark %x %y %W }
    bind $ci_edit.c <Shift-Button-2> { ci_DeleteMark %x %y %W }
}

proc ci_ReadColorMapFile { file_name color_array knot_array ncolors } {
   global global_limits

   upvar $color_array color
   upvar $knot_array  knot
   upvar $ncolors maxindex

   set file [open $file_name RDONLY]

   set maxindex 0
   while { ![eof $file] } {
       set str [string trim [gets $file]]

       if { [scan $str "%s %d %f %f %f" code index R G B] <= 0 && [eof $file] } break;
       set code [string tolower [string index $code 0]]

       if { $code == "c" || $code == "k" } {

           if { $index >= $global_limits(MAX_COLORMAP_ENTRIES) || $index < 0 } {
               return -code error "Colormap index \[$index\] out of range, aborting"
           }

           set R [@ int(2.55*$R+0.5)]
           set G [@ int(2.55*$G+0.5)]
           set B [@ int(2.55*$B+0.5)]

           set value [format "#%02x%02x%02x" $R $G $B]

           if { $code == "c" } {
               set color($index) $value
           } elseif { $code == "k" } {
               set knot($index) $value
           }

           if { $index > $maxindex } { set maxindex $index }
       } 
   }
   close $file

   incr maxindex
}

#
# Read a colormap definitiom from a given file
#
# 17 Sep 1995
#
proc ci_ReadColorMap { w {file_name ""} } {
   global ci_nbox_x ci_nbox_y ci_marks ci_size ci_list ELMER_POST_HOME

   if { $file_name == "" } { set file_name [fs_FileSelect $ELMER_POST_HOME/lib/colormaps *.cm read] }
   if { $file_name == "" } return

   ci_ReadColorMapFile $file_name color knot maxindex

   ci_CreateColorMap $w $maxindex 0

   if { [info exists color] } {
       set sid [array startsearch color]
       while { [array anymore color $sid] != 0 } {

           set k [array nextelement color $sid] 
           $w.c itemconfigure [@ $k+1]  -fill $color($k)

       }
       array donesearch color $sid
   }

   set ci_marks 0
   if { [info exists ci_list] } { unset ci_list }

   if { [info exists knot] } {
       set sid [array startsearch knot]
       while { [array anymore knot $sid] != 0 } {

           set k [array nextelement knot $sid] 

           set x [@ $ci_size*($k % $ci_nbox_x)+$ci_size/2]
           set y [@ $ci_size*(int($k/(1.0*$ci_nbox_x)))+$ci_size/2]

           $w.c itemconfigure [@ $k+1]  -fill $knot($k)
           ci_CreateMark $x $y $w

       }
       array donesearch knot $sid
   }

   ci_InterpolateColor $w  
}

#
# Write a colormap definitiom to a given file
#
# 17 Sep 1995
#
proc ci_WriteColorMap { w {file_name ""} } {
    global ci_nbox_x ci_nbox_y ci_list ci_isize env

    if { $file_name == "" } { set file_name [fs_FileSelect [pwd] *.cm write] }
    if { $file_name == "" } return

    if { [file exists $file_name] } {
        set ret [dl_dialog {File Exists} "File \[[file tail $file_name]\] already exists, overwrite ?" {warning} 1 {CANCEL} {OK}]
        if { $ret == 0 } return

        set outfile [open $file_name {WRONLY TRUNC}]
    } else {
        set outfile [open $file_name {WRONLY CREAT}]
    }

    set maxbox 0
    set minbox [@ $ci_nbox_x*$ci_nbox_y]

    puts $outfile  "!\n!\n! ElmerPost colormap definition file."
    puts $outfile  "!\n!\n! Created: [$ date] by $env(USER)\n!\n!"
    puts $outfile  "! code\tindex\tred(%)\tgreen(%)\tblue(%)"

    set sid [array startsearch ci_list]
    while { [array anymore ci_list $sid] != 0 } {
       set k [array nextelement ci_list $sid] 
       set n $ci_list($k)

       set l [$w.c coords $n]

       set x [lindex $l 0]
       set y [lindex $l 1]

       set box_j [@ int($x*$ci_isize)]
       set box_i [@ int($y*$ci_isize)]

       set box_n [@ $box_i*$ci_nbox_x+$box_j]

       scan [$w.c itemcget [@ $box_n+1] -fill] "#%02x%02x%02x" R G B

       set R [@ $R/2.55]
       set G [@ $G/2.55]
       set B [@ $B/2.55]

       puts $outfile "  knot\t$box_n\t$R\t$G\t\t$B"

       if { $box_n < $minbox } { set minbox $box_n }
       if { $box_n > $maxbox } { set maxbox $box_n }
    }
    array donesearch ci_list $sid

    do i 0 [@ $minbox-1] {
        scan [$w.c itemcget [@ $i+1] -fill] "#%02x%02x%02x" R G B
        set R [@ $R/2.55]
        set G [@ $G/2.55]
        set B [@ $B/2.55]
        puts $outfile " color $i $R $G $B"
    }

    do i [@ $maxbox+1] [@ $ci_nbox_x*$ci_nbox_y-1] {
        scan [$w.c itemcget [@ $i+1] -fill] "#%02x%02x%02x" R G B
        set R [@ $R/2.55]
        set G [@ $G/2.55]
        set B [@ $B/2.55]
        puts $outfile "  color\t$i\t$R\t$G\t$B"
    }

    close $outfile
}

#
# Send colormap definition to drawing
#
# 17 Sep 1995
#
proc ci_SendColorMap { w } {
    global ci_applyto ci_nbox_x ci_nbox_y colmap colmap_size

    set rgblist ""
    set n [@ $ci_nbox_x*$ci_nbox_y]
    do i 1 $n {
         set t [$w.c itemcget $i -fill]
         set colmap([@ $i-1]) $t
         set rgblist $rgblist$t
    }

    set colmap_size $n
    GetColorMap $ci_applyto $n $rgblist
}

#
# Set slider colors according to the current box color
#
# 17 Sep 1995
#
proc ci_SetColor { w box_n args } {
   global ci_nbox_x ci_nbox_y ci_isize ci_current

   set value [$w.c itemcget $box_n -fill]

   scan $value "#%02x%02x%02x" R G B
   
   set R [@ $R/2.55]
   set G [@ $G/2.55]
   set B [@ $B/2.55]

   $w.red.scl set   $R
   $w.green.scl set $G
   $w.blue.scl set  $B
}

#
# Set box color according to the current slider values
#
# 17 Sep 1995
#
proc ci_ShowColor {w args} {
   global ci_nbox_x ci_nbox_y ci_isize ci_current

   set R [$w.red.scl   get]
   set G [$w.green.scl get]
   set B [$w.blue.scl  get]

   set R [@ int(2.55*$R+0.5)]
   set G [@ int(2.55*$G+0.5)]
   set B [@ int(2.55*$B+0.5)]

   set value [format "#%02x%02x%02x" $R $G $B]

   set i $ci_current($w.c,obj)
   set l [$w.c coords $i]

   set x [lindex $l 0]
   set y [lindex $l 1]

   set box_j [@ int($x*$ci_isize)]
   set box_i [@ int($y*$ci_isize)]

   set box_n [@ $box_i*$ci_nbox_x+$box_j]

   incr box_n
   $w.c itemconfigure $box_n -fill $value

   ci_InterpolateColor $w

}

#
# Get color name from listbox and show the color with a button
#
# 12 Sep 1995
#
proc ci_SetListColor { w } {
    global me_rgb_list

    set top [winfo toplevel $w]

    set value [$w get [$w curselection]]

    $top.red.scl   set $me_rgb_list($value,R)
    $top.green.scl set $me_rgb_list($value,G)
    $top.blue.scl  set $me_rgb_list($value,B)

    ci_ShowColor $top
}


#
# interpolate colors given a set of knots and color values at the knots
#
# 17 Sep 1995
#
proc ci_InterpolateColor { w } {
   global ci_marks ci_list ci_isize ci_size ci_nbox_x ci_nbox_y

   set box_list ""

   set sid [array startsearch ci_list]
   while { [array anymore ci_list $sid] != 0 } {
       set n [array nextelement ci_list $sid] 
       set l [$w.c coords $ci_list($n)]

       set x [lindex $l 0]
       set y [lindex $l 1]

       set box_j [@ int($x*$ci_isize)]
       set box_i [@ int($y*$ci_isize)]

       set box_n [@ $box_i*$ci_nbox_x+$box_j]

       if { $box_list == ""  } {
           set box_list $box_n
       } else {
           set box_list "$box_list $box_n"
       }
   }
   array donesearch ci_list $sid

   set box_list [split $box_list]
   set box_list [lsort -integer $box_list]

   set i 0

   set n0 [lindex $box_list $i]
   set color0 [$w.c itemcget [@ $n0+1] -fill]

   do i 1 [@ [llength $box_list]-1] {
        set n1 [lindex $box_list $i]

        if { $n0 == $n1 } continue;
        set color1 [$w.c itemcget [@ $n1+1] -fill]

        set rgblist [split [GetInterpolate $n0 $n1 $color0 $color1] " "]
 
        set k 0
        do j [@ $n0+1] $n1 { 
             $w.c itemconfigure $j -fill [lindex $rgblist $k]
             incr k
         }

        set n0 $n1
        set color0 $color1
   }
}

#
# Create a now knot given x,y coordinates over canvas
#
# 17 Sep 1995
#
proc ci_CreateMark { x y w } { 
    global ci_list ci_marks ci_mark_n ci_size ci_isize ci_nbox_x ci_nbox_y

    set w [winfo toplevel $w]

    set x [@ $ci_size*int($x*$ci_isize) + $ci_size/2 ]
    set y [@ $ci_size*int($y*$ci_isize) + $ci_size/2 ]

    set ci_list($ci_mark_n) [$w.c create arc [@ $x-5] [@ $y-5] [@ $x+5] [@ $y+5] -fill red -start 0 -extent 359.99 -tag ci_mark_$ci_mark_n]

    $w.c bind $ci_list($ci_mark_n) <B1-Motion> { ci_Drag %x %y %W }
    $w.c bind $ci_list($ci_mark_n) <Button-1> { ci_GetClosest %x %y %W } 

    ci_GetClosest $x $y $w.c
 
    incr ci_marks
    incr ci_mark_n
}

#
# Delete current selected knot
#
# 18 Sep 1995
#
proc ci_DeleteCurrent { w } { 
   global ci_list ci_marks ci_size ci_isize ci_nbox_x ci_nbox_y ci_current

   if { $ci_marks <= 0 } return

   set w [winfo toplevel $w]

   set str [$w.c itemcget $ci_current($w.c,obj) -tag]

   scan $str "ci_mark_%d" ctag
   unset ci_list($ctag)
 
   $w.c delete $ci_current($w.c,obj)
   set ci_marks [@ $ci_marks-1]

   if { $ci_marks <= 0 } return

   set sid [array startsearch ci_list]

   set n [array nextelement ci_list $sid] 
   set l [$w.c coords $ci_list($n)]

   set x [lindex $l 0]
   set y [lindex $l 1]
 
   set x [@ $ci_size*int($x*$ci_isize) + $ci_size/2 ]
   set y [@ $ci_size*int($y*$ci_isize) + $ci_size/2 ]

   ci_GetClosest $x $y $w

   array donesearch ci_list $sid

   ci_InterpolateColor $w
}

proc ci_DeleteMark { x y w } { 
    global ci_list ci_marks ci_size ci_isize ci_nbox_x ci_nbox_y

    if { $ci_marks <= 0 } return

    set w [winfo toplevel $w]

    set x [@ $ci_size*int($x*$ci_isize) + $ci_size/2 ]
    set y [@ $ci_size*int($y*$ci_isize) + $ci_size/2 ]

    ci_GetClosest $x $y $w.c
 
    set a [$w.c find closest $x $y]
    if { $a <= [@ $ci_nbox_x*$ci_nbox_y] } return;

    set str [$w.c itemcget $a -tag]
    scan $str "ci_mark_%d" ctag
    unset ci_list($ctag)
 
    $w.c delete $a
 
    set ci_marks [@ $ci_marks-1]

    if { $ci_marks <= 0 } return

    set sid [array startsearch ci_list]

    set n [array nextelement ci_list $sid] 
    set l [$w.c coords $ci_list($n)]

    set x [lindex $l 0]
    set y [lindex $l 1]

    set x [@ $ci_size*int($x*$ci_isize) + $ci_size/2 ]
    set y [@ $ci_size*int($y*$ci_isize) + $ci_size/2 ]

    ci_GetClosest $x $y $w

    array donesearch ci_list $sid

    ci_InterpolateColor $w
}

#
# Get closest knot given x,y coordinates 
#
# 17 Sep 1995
#
proc ci_GetClosest { x y w } { 
    global ci_current ci_size ci_isize ci_nbox_x ci_nbox_y

    set w [winfo toplevel $w]

    set a [$w.c find closest $x $y]
    if { $a <= [@ $ci_nbox_x*$ci_nbox_y] } return;

    if { [info exists ci_current($w.c,obj)] } {
        $w.c itemconfigure $ci_current($w.c,obj) -fill red
    }

    set ci_current($w.c,obj) [$w.c find closest $x $y]
 
    $w.c itemconfigure $ci_current($w.c,obj) -fill green

    set box_i [@ int($y*$ci_isize)]
    set box_j [@ int($x*$ci_isize)]

    set x [@ $ci_size*$box_j + $ci_size/2]
    set y [@ $ci_size*$box_i + $ci_size/2]

    set ci_current($w.c,x) $x
    set ci_current($w.c,y) $y

    set box_n [@ $box_i*$ci_nbox_x+$box_j]

    incr box_n

    ci_SetColor $w $box_n

    set ci_current($w.c,color) [$w.c itemcget $box_n -fill]

puts "closest $box_n $ci_current($w.c,color)"
}

#
# Move a knot interavtively using leftmost mouse button
#
# 17 Sep 1995
#
proc ci_Drag { x y w } { 
    global ci_current ci_size ci_isize ci_nbox_x ci_nbox_y

    set w [winfo toplevel $w]

    if { $x > [@ $ci_size*$ci_nbox_x] } { set x [@ $ci_size*($ci_nbox_x-1)] }
    if { $y > [@ $ci_size*$ci_nbox_y] } { set y [@ $ci_size*($ci_nbox_y-1)] }

    if { $x < 0 }  { set x 0 }
    if { $y < 0 }  { set y 0 }

    set oldbox_j [@ int($ci_current($w.c,x)*$ci_isize)]
    set oldbox_i [@ int($ci_current($w.c,y)*$ci_isize)]

    set newbox_j [@ int($x*$ci_isize)]
    set newbox_i [@ int($y*$ci_isize)]

    set x [@ $ci_size*$newbox_j + $ci_size/2]
    set y [@ $ci_size*$newbox_i + $ci_size/2]

    set dx [@ $x - $ci_current($w.c,x)]
    set dy [@ $y - $ci_current($w.c,y)]      

    $w.c move $ci_current($w.c,obj) $dx $dy 

    set ci_current($w.c,x) $x
    set ci_current($w.c,y) $y

    set box_n [@ $newbox_i*$ci_nbox_x+$newbox_j]

    incr box_n
    $w.c itemconfigure $box_n -fill $ci_current($w.c,color)

    ci_InterpolateColor $w
}
