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
#* Material editor utility widget.
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
# User level command for setting material parametrs
#
proc material { args } {
     global me_color me_A me_S

     set Usage "\
Usage: material -diffuse r,g,b -specular r,g,b -shininess exp -opacity percent\n\n\
-(r,g,b) values are given in percents of the full intensity\n -exp can be range from  0 to 128,\n -opactity is in percents.\n"

     set opt(1,name)  -diffuse;     set opt(1,args)  3
     set opt(2,name)  -specular;    set opt(2,args)  3
     set opt(3,name)  -shininess;   set opt(3,args)  1
     set opt(4,name)  -opacity;     set opt(4,args)  1
     set opt(5,name)  -noupdate;    set opt(5,args)  0
     set opt(6,name)  -vectors;     set opt(6,args)  0
     set opt(7,name)  -contour;     set opt(7,args)  0
     set opt(8,name)  -mesh;        set opt(8,args)  0
     set opt(9,name)  -sphere;      set opt(9,args)  0
     set opt(10,name) -particles;   set opt(10,args) 0
     set opt(11,name) -isosurface;  set opt(11,args) 0
     set opt(12,name) -default;     set opt(12,args) 0

     check_args material $Usage 12 0 opt opt_val 0 0 arg_val $args

     if { ![info exists opt_val] } { return $Usage }

     set me_edit .me_editMaterial
     set update 1

     set cmd "material"

     set sid [array startsearch opt_val]

     set defopt "-default"

     while { [array anymore opt_val $sid] != 0 } {

         set n [array nextelement opt_val $sid]
         set value $opt_val($n)
         set n [string range $n 1 [string length $n]]

         switch $n {
            diffuse {
                 if { [info exists me_color($me_edit,diffuse,R)] } {
                     set R $me_color($me_edit,diffuse,R)
                     set G $me_color($me_edit,diffuse,G)
                     set B $me_color($me_edit,diffuse,B)
                 }

                 scan $value "%f,%f,%f" R G B

                 set me_color($me_edit,diffuse,R) $R
                 set me_color($me_edit,diffuse,G) $G
                 set me_color($me_edit,diffuse,B) $B

                 set cmd "$cmd -diffuse $R,$G,$B"
            }
            specular {
                 if { [info exists me_color($me_edit,specular,R)] } {
                     set R $me_color($me_edit,specular,R)
                     set G $me_color($me_edit,specular,G)
                     set B $me_color($me_edit,specular,B)
                 }
     
                 scan $value "%f,%f,%f" R G B
  
                 set me_color($me_edit,specular,R) $R
                 set me_color($me_edit,specular,G) $G
                 set me_color($me_edit,specular,B) $B

                 set cmd "$cmd -specular $R,$G,$B"
            }
            shininess {
                if { [info exists me_S($me_edit)] }  { set S $me_S($me_edit) }
                scan $value "%d" S
                set me_S($me_edit) $S

                set cmd "$cmd -shininess $S"
            } 
            opacity {
                if { [info exists me_A($me_edit)] }  { set A $me_A($me_edit) }
                scan $value "%d" A
                set me_A($me_edit) $A

                set cmd "$cmd -opacity $A"
            }
            noupdate {
                set update 0
            }
            vectors    { set defopt "-vectors" }
            mesh       { set defopt "-mesh"    }
            contour    { set defopt "-contour" }
            sphere     { set defopt "-sphere"  }
            particles  { set defopt "-particles" }
            isosurface { set defopt "-isosurface" }
         }
     }

     array donesearch opt_val $sid

     if { $cmd == "" } return

     if { $update && [winfo exists $me_edit] } { me_UpdateAll $me_edit }

     UpdateColor $defopt color diffuse $me_color($me_edit,diffuse,R) $me_color($me_edit,diffuse,G) $me_color($me_edit,diffuse,B)     \
                 color specular $me_color($me_edit,specular,R) $me_color($me_edit,specular,G) $me_color($me_edit,specular,B) \
                 opacity $me_A($me_edit) shininess $me_S($me_edit)
     UpdateObject;

     if { !$update } { history add $cmd }
}

#
# User level command to set background color
#
proc background { args } {
    global me_color

    set Usage \
"Usage: background red green blue\n\nRed,Green, and Blue are in percents of the full intensity.\n"

    set opt(1,name) -noupdate;  set opt(1,args) 0

    check_args background $Usage 1 1 opt opt_val 1 3 arg_val $args

    set me_edit .me_editBackground

    set update 1

    if { [info exists opt_val(-noupdate)] } { set update 0 }

    if { [info exists me_color($me_edit,diffuse,R)] } {
        set R $me_color($me_edit,diffuse,R)
        set G $me_color($me_edit,diffuse,G)
        set B $me_color($me_edit,diffuse,B)
    }

    scan $args "%f %f %f" R G B

    set me_color($me_edit,diffuse,R) $R
    set me_color($me_edit,diffuse,G) $G
    set me_color($me_edit,diffuse,B) $B

    if { $update && [winfo exists $me_edit] } { me_UpdateColor $me_edit }

    UpdateBackColor color diffuse $me_color($me_edit,diffuse,R) $me_color($me_edit,diffuse,G) $me_color($me_edit,diffuse,B)

    UpdateEdgeColor color diffuse [@ 100-$R] [@ 100-$G] [@ 100-$B]

    if { !$update } { history add "background $R $G $B" }
}

#
# Create a material editor toplevel window for a user to use.
# There can be several editors running at any time provided 
# a different name for editor is given at each invokation.
#
# 12 Sep 1995
#
proc me_Edit { name just_colors {UpdateProc material} } {
   global me_color me_A me_S me_edit_var me_rgb_list me_R me_G me_B me_apply_to me_defopt

   set me_edit .me_edit$name

   if { [winfo exists $me_edit] } { destroy $me_edit }

   if { ![info exists me_color($me_edit,diffuse,R)] } {
       set me_color($me_edit,diffuse,R) 100
       set me_color($me_edit,diffuse,G) 100
       set me_color($me_edit,diffuse,B) 100
   }

   if { ![info exists me_color($me_edit,specular,R)] } {
       set me_color($me_edit,specular,R) 0
       set me_color($me_edit,specular,G) 0
       set me_color($me_edit,specular,B) 0
   }

   if { ![info exists me_S($me_edit)] } { set me_S($me_edit)  20 }

   if { ![info exists me_A($me_edit)] } { set me_A($me_edit) 100 }

   if { ![info exists me_edit_var($me_edit)] } { set me_edit_var($me_edit) "diffuse" }

   set me_R($me_edit) $me_color($me_edit,$me_edit_var($me_edit),R)
   set me_G($me_edit) $me_color($me_edit,$me_edit_var($me_edit),G)
   set me_B($me_edit) $me_color($me_edit,$me_edit_var($me_edit),B)

   toplevel $me_edit
   place_window $me_edit
   wm minsize $me_edit 100 100
   wm title $me_edit $name

   label $me_edit.title -height 1;#-text "Elmer Post - Material" -height 3 -font -Adobe-Times-Medium-R-Normal-*-180-*
   pack $me_edit.title -side top -expand 1 -fill both

    frame $me_edit.menubar

    menubutton $me_edit.menubar.file -menu $me_edit.menubar.file.menu -text "Apply-To" -underline 0
    menu $me_edit.menubar.file.menu

    set me_apply_to 0
    set me_defopt   "-default"
    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Default" -value 0 -variable me_apply_to \
               -command { set me_defopt "-default" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Mesh" -value 1 -variable me_apply_to \
               -command { set me_defopt "-mesh" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Vectors" -value 2 -variable me_apply_to \
               -command { set me_defopt "-vectors" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Contour" -value 3 -variable me_apply_to \
               -command { set me_defopt "-contour" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Isosurface" -value 4 -variable me_apply_to \
               -command { set me_defopt "-isosurface" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Spheres" -value 5 -variable me_apply_to \
               -command { set me_defopt "-sphere" }

    $me_edit.menubar.file.menu add radiobutton
    $me_edit.menubar.file.menu entryconfigure last -label "Particles" -value 6 -variable me_apply_to \
               -command { set me_defopt "-particles" }

    pack $me_edit.menubar.file -side left
    pack $me_edit.menubar -side top -expand 1 -fill both

   frame $me_edit.red   -bd 10 -relief ridge -bg lightblue
   frame $me_edit.green -bd 10 -relief ridge -bg lightblue
   frame $me_edit.blue  -bd 10 -relief ridge -bg lightblue

   if { !$just_colors } {
#
# Radiobuttons selecting color (diffuse or specular) to be edited
#
       frame $me_edit.rbox -bd 5 -relief ridge -bg lightblue
       radiobutton $me_edit.rbox.diffuse  -variable me_edit_var($me_edit) -value "diffuse" \
                 -text "Ambient & Diffuse" -command "me_ChangeSpec $me_edit"

       radiobutton $me_edit.rbox.specular -variable me_edit_var($me_edit) -value "specular" \
                     -text "Specular" -command "me_ChangeSpec $me_edit"



       pack $me_edit.rbox.diffuse  -side left -expand 1 -fill x
       pack $me_edit.rbox.specular -side left -expand 1 -fill x
       pack $me_edit.rbox -side top -expand 1 -fill both

#
# slider for Shininess setting
#
       frame $me_edit.sbox -bd 5 -relief ridge -bg lightblue
       slider $me_edit.sbox.scl -orient horizontal -command "me_SetShine $me_edit" -from 0 -to 128 -digit 4 -resol 0.5 -tick 32
       label $me_edit.sbox.lab -text "Shininess" -relief raised -bd 5

       pack $me_edit.sbox.scl -side bottom -expand 1 -fill x
       pack $me_edit.sbox.lab -side bottom -expand 1 -fill x
       pack $me_edit.sbox -expand 1 -fill both

#
# slider for Opacity setting
#
       frame $me_edit.obox -bd 5 -relief ridge -bg lightblue
       slider $me_edit.obox.scl -orient horizontal -command "me_SetOpac $me_edit" -from 0 -to 100 -digit 4 -resol 0.5
       label $me_edit.obox.lab -text "Opacity (%)" -relief raised -bd 5

       pack $me_edit.obox.scl -side bottom -expand 1 -fill x
       pack $me_edit.obox.lab -side bottom -expand 1 -fill x
       pack $me_edit.obox -expand 1 -fill both
   }

#
# button, background showing current color, pressing it same as Apply-button (downwards)
#
   button $me_edit.value -bd 10 -relief raised -command \
           "$UpdateProc \$me_defopt -\$me_edit_var($me_edit) \$me_R($me_edit),\$me_G($me_edit),\$me_B($me_edit) \
                  -opacity \$me_A($me_edit) -shininess \$me_S($me_edit) -noupdate"


#
# rgb sliders
#
   slider $me_edit.red.scl -orient horizontal -command "me_ShowColor $me_edit" -from 0 -to 100 -troughcolor red -digit 4 -resol 0.5
   slider $me_edit.green.scl -orient horizontal -command "me_ShowColor $me_edit" -from 0 -to 100 -troughcolor green -digit 4 -resol 0.5
   slider $me_edit.blue.scl -orient horizontal -command "me_ShowColor $me_edit" -from 0 -to 100 -troughcolor blue -digit 4 -tick 25 -resol 0.5

   pack $me_edit.red.scl -side left -expand 1 -fill x
   pack $me_edit.green.scl -side left -expand 1 -fill x
   pack $me_edit.blue.scl -side left -expand 1 -fill x

   pack $me_edit.red -side top -fill x
   pack $me_edit.green -side top -fill x
   pack $me_edit.blue -side top  -fill x
   pack $me_edit.value -side top -expand 1 -fill both

#
# List box containing color names
#
   frame $me_edit.lbox -relief sunken -bg lightblue
   listbox $me_edit.lbox.list -yscroll "$me_edit.lbox.scroll set"
   scrollbar $me_edit.lbox.scroll -command "$me_edit.lbox.list yview"

   pack $me_edit.lbox.scroll -side left -fill y
   pack $me_edit.lbox.list -side left -expand 1 -fill both
   pack $me_edit.lbox -side top -expand 1 -fill both

   bind $me_edit.lbox.list <Double-1> {me_SetListColor %W}

#
# button box at window bottom (Apply,OK,Cancel)
#
   frame $me_edit.bbox  -bg lightblue

   button $me_edit.bbox.apply -text "Apply" -command \
           "$UpdateProc \$me_defopt -\$me_edit_var($me_edit) \$me_R($me_edit),\$me_G($me_edit),\$me_B($me_edit) \
                  -opacity \$me_A($me_edit) -shininess \$me_S($me_edit) -noupdate"

   button $me_edit.bbox.cancel -text "Cancel" -command "destroy $me_edit"

   button $me_edit.bbox.ok -text "OK" -command \
           "$UpdateProc \$me_defopt -\$me_edit_var($me_edit) \$me_R($me_edit),\$me_G($me_edit),\$me_B($me_edit) \
              -opacity \$me_A($me_edit) -shininess \$me_S($me_edit) -noupdate; destroy $me_edit"

   pack $me_edit.bbox.apply  -side left -expand 1 -fill x
   pack $me_edit.bbox.ok  -side left -expand 1 -fill x
   pack $me_edit.bbox.cancel  -side left -expand 1 -fill x

   pack $me_edit.bbox -side top -expand 1 -fill both

#
#  Get database of color names -> RGB-values, and add them to a listbox
#
   if { ![info exists me_rgb_list] } me_GetColorNames
   me_AddColorNames $me_edit.lbox.list

#
# Update the slider etc. values
#
   if { $name == "Material"   } { me_UpdateAll   $me_edit }
   if { $name == "Background" } { me_UpdateColor $me_edit }
}

#
# Update the color settings given current state
#
proc me_UpdateColor { top } {
   global me_edit_var me_R me_G me_B me_color

   if { $me_edit_var($top) == "diffuse" } {
       set me_R($top) $me_color($top,diffuse,R)
       set me_G($top) $me_color($top,diffuse,G)
       set me_B($top) $me_color($top,diffuse,B)

       $top.red.scl set $me_R($top)
       $top.green.scl set $me_G($top)
       $top.blue.scl set $me_B($top)
   } else {
       set me_R($top) $me_color($top,specular,R)
       set me_G($top) $me_color($top,specular,G)
       set me_B($top) $me_color($top,specular,B)

       $top.red.scl set $me_R($top)
       $top.green.scl set $me_G($top)
       $top.blue.scl set $me_B($top)
   }
}


#
# Update the interaction window, given current state
#
proc me_UpdateAll { top } {
    global me_S me_A

    me_UpdateColor $top

    $top.sbox.scl set $me_S($top)
    $top.obox.scl set $me_A($top)
}

#
# User pressed one of the radiobuttons changing the diffuse/specular setting
#
proc me_ChangeSpec { top } {
    global me_edit_var me_R me_G me_B me_color

    set me_R($top) $me_color($top,$me_edit_var($top),R)
    set me_G($top) $me_color($top,$me_edit_var($top),G)
    set me_B($top) $me_color($top,$me_edit_var($top),B)

    $top.red.scl   set $me_R($top)
    $top.green.scl set $me_G($top)
    $top.blue.scl  set $me_B($top)
}

#
# Get color name from listbox and show the color with a button
#
# 12 Sep 1995
#
proc me_SetListColor {w} {
    global me_edit_var me_rgb_list me_R me_G me_B me_color

    set top [winfo toplevel $w]

    set value [$w get [$w curselection]]

    set me_R($top) $me_rgb_list($value,R)
    set me_G($top) $me_rgb_list($value,G)
    set me_B($top) $me_rgb_list($value,B)

    $top.red.scl   set $me_R($top)
    $top.green.scl set $me_G($top)
    $top.blue.scl  set $me_B($top)

    set R $me_R($top)
    set G $me_G($top)
    set B $me_B($top)

    set me_color($top,$me_edit_var($top),R) $R
    set me_color($top,$me_edit_var($top),G) $G
    set me_color($top,$me_edit_var($top),B) $B

    set R [@ int(2.55*$R+0.5)]
    set G [@ int(2.55*$G+0.5)]
    set B [@ int(2.55*$B+0.5)]

    set value  [format "#%02x%02x%02x" $R $G $B]
    set value1 [format "#%02x%02x%02x" [@ 255-$R] [@ 255-$G] [@ 255-$B]]

    $top.value configure -bg $value -activeback $value -fg $value1
}

proc me_Nil { args } {}

#
# Get rgb-definition from sliders and set button background accordingly
#
# 12 Sep 1995
#
proc me_ShowColor {top args} {
    global me_R me_G me_B me_edit_var me_color

    set me_R($top) [$top.red.scl   get]
    set me_G($top) [$top.green.scl get]
    set me_B($top) [$top.blue.scl  get]

    set R $me_R($top)
    set G $me_G($top)
    set B $me_B($top)

    set me_color($top,$me_edit_var($top),R) $R
    set me_color($top,$me_edit_var($top),G) $G
    set me_color($top,$me_edit_var($top),B) $B

    set R [@ int(2.55*$R+0.5)]
    set G [@ int(2.55*$G+0.5)]
    set B [@ int(2.55*$B+0.5)]

    set value  [format "#%02x%02x%02x" $R $G $B]
    set value1 [format "#%02x%02x%02x" [@ 255-$R] [@ 255-$G] [@ 255-$B]]

    $top.value configure -bg $value -activeback $value -fg $value1
}

#
# Set material state for opacity according to slider value 
#
# 19 Sep 1995
#
proc me_SetOpac {top value} {
    global me_A

    set me_A($top) $value
}

#
# Set material state for shininess according to slider value 
#
# 19 Sep 1995
#
proc me_SetShine {top value} {
    global me_S

    set me_S($top) $value
}

#
# Add color names to a listbox given as argument
#
# 13 Sep 1995
#
proc me_AddColorNames {w} {
    global me_rgb_count me_rgb_names

    do i 0 [@ $me_rgb_count-1] { $w insert end $me_rgb_names($i) }
}

#
# Read color names and their rgb-definitions from a file rgb.txt
#
# 13 Sep 1995
#
proc me_GetColorNames {} {
   global me_rgb_list me_rgb_names me_rgb_count ELMER_POST_HOME

   set me_rgb_file [open "$ELMER_POST_HOME/lib/rgb.txt" RDONLY]

   set me_rgb_count 0
   while {![eof $me_rgb_file]} {
      set me_str [gets $me_rgb_file]
      set me_str [string trim $me_str]
      if { [string index $me_str 0] == "!" } continue

      set R 0
      set G 0
      set B 0

      set me_name1 ""
      scan $me_str "%d %d %d %s %s" R G B me_name me_name1
      set me_name "$me_name $me_name1"
      set me_name [string trim $me_name]

      set R [@ $R/2.55]
      set G [@ $G/2.55]
      set B [@ $B/2.55]

      set me_rgb_list($me_name,R) $R
      set me_rgb_list($me_name,G) $G
      set me_rgb_list($me_name,B) $B
      set me_rgb_names($me_rgb_count) $me_name

      incr me_rgb_count
   }

   close $me_rgb_file
}
