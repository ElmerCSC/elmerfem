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
# Camera settings
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
#*                       Date: 11 Oct 1995
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#*******************************************************************************
#

proc cam_CreateEdit { { current "" } } {
   global camera camera_names camera_n camera_current

   set camera_edit .camera_edit

   if { ![winfo exists $camera_edit] } {
      toplevel $camera_edit
      place_window $camera_edit
   } else {
      destroy $camera_edit.left
      destroy $camera_edit.right
	  wm iconify $camera_edit
	  wm deiconify $camera_edit
   }


   if { $current == "" } { set current $camera_names(0) }
   set camera_current [c_CurrentCamera $current]

   bind $camera_edit <Return> "c_SetCamera"

   frame $camera_edit.left  -bd 10 -relief ridge
   frame $camera_edit.right  -bd 10 -relief ridge

   frame $camera_edit.left.from  -bd 10 -relief ridge
   entry $camera_edit.left.from.ent_x -relief sunken -width 5 -textvariable camera(from_x)
   label $camera_edit.left.from.lab_x -text "From x"

   entry $camera_edit.left.from.ent_y -relief sunken -width 5 -textvariable camera(from_y)
   label $camera_edit.left.from.lab_y -text "From y" 

   entry $camera_edit.left.from.ent_z -relief sunken -width 5 -textvariable camera(from_z)
   label $camera_edit.left.from.lab_z -text "From z"

   frame $camera_edit.left.to -bd 10 -relief ridge
   entry $camera_edit.left.to.ent_x -relief sunken -width 5 -textvariable camera(to_x)
   label $camera_edit.left.to.lab_x -text "To x" 

   entry $camera_edit.left.to.ent_y -relief sunken -width 5 -textvariable camera(to_y)
   label $camera_edit.left.to.lab_y -text "To y" 

   entry $camera_edit.left.to.ent_z -relief sunken -width 5 -textvariable camera(to_z)
   label $camera_edit.left.to.lab_z -text "To z" 

   frame $camera_edit.left.up -bd 10 -relief ridge
   entry $camera_edit.left.up.ent_x -relief sunken -width 5 -textvariable camera(up_x)
   label $camera_edit.left.up.lab_x -text "Up x" 

   entry $camera_edit.left.up.ent_y -relief sunken -width 5 -textvariable camera(up_y)
   label $camera_edit.left.up.lab_y -text "Up y" 

   entry $camera_edit.left.up.ent_z -relief sunken -width 5 -textvariable camera(up_z)
   label $camera_edit.left.up.lab_z -text "Up z" 

   set camera(name) $current

   frame $camera_edit.right.status -bd 10 -relief ridge
   label $camera_edit.right.status.label -text "Current: $current, " 
   label $camera_edit.right.status.status -text "Status: " 
   radiobutton $camera_edit.right.status.on -variable camera(display) -text "On" -value "on"
   radiobutton $camera_edit.right.status.off -variable camera(display) -text "Off" -value "off"

   frame $camera_edit.right.projection -bd 10 -relief ridge

   label $camera_edit.right.projection.label -text "Projection: " 

   radiobutton $camera_edit.right.projection.perspective -variable camera(projection) \
        -text "Perspective" -value "perspective" -command { c_SetCamera; cam_CreateEdit $camera_names($camera_current) }

   radiobutton $camera_edit.right.projection.ortho -variable camera(projection) \
          -text "Ortho" -value "ortho" -command { c_SetCamera; cam_CreateEdit $camera_names($camera_current) }

   if { $camera(projection) == "perspective" } {

       frame $camera_edit.right.field -bd 10 -relief ridge
       entry $camera_edit.right.field.ent -relief sunken -width 5 -textvariable camera(field_angle)
       label $camera_edit.right.field.lab -text "Field Angle" 

       frame $camera_edit.right.near -bd 10 -relief ridge
       entry $camera_edit.right.near.ent -relief sunken -width 5 -textvariable camera(near)
       label $camera_edit.right.near.lab -text "Near Clip" 
   
       frame $camera_edit.right.far -bd 10 -relief ridge
       entry $camera_edit.right.far.ent -relief sunken -width 5 -textvariable camera(far)
       label $camera_edit.right.far.lab -text "Far Clip" 

   } else {

       frame $camera_edit.right.left -bd 10 -relief ridge
       entry $camera_edit.right.left.ent -relief sunken -width 5 -textvariable camera(left)
       label $camera_edit.right.left.lab -text "Left Clip" 

       frame $camera_edit.right.right -bd 10 -relief ridge
       entry $camera_edit.right.right.ent -relief sunken -width 5 -textvariable camera(right)
       label $camera_edit.right.right.lab -text "Right Clip" 

       frame $camera_edit.right.bottom -bd 10 -relief ridge
       entry $camera_edit.right.bottom.ent -relief sunken -width 5 -textvariable camera(bottom)
       label $camera_edit.right.bottom.lab -text "Bottom Clip" 

       frame $camera_edit.right.top -bd 10 -relief ridge
       entry $camera_edit.right.top.ent -relief sunken -width 5 -textvariable camera(top)
       label $camera_edit.right.top.lab -text "Top Clip" 

       frame $camera_edit.right.near -bd 10 -relief ridge
       entry $camera_edit.right.near.ent -relief sunken -width 5 -textvariable camera(near)
       label $camera_edit.right.near.lab -text "Near Clip" 

       frame $camera_edit.right.far -bd 10 -relief ridge
       entry $camera_edit.right.far.ent -relief sunken -width 5 -textvariable camera(far)
       label $camera_edit.right.far.lab -text "Far Clip" 

   }
   
   frame $camera_edit.right.menubar -relief raised

   menubutton $camera_edit.right.menubar.file -menu $camera_edit.right.menubar.file.menu -text File -underline 0
   menu $camera_edit.right.menubar.file.menu

   $camera_edit.right.menubar.file.menu add command
   $camera_edit.right.menubar.file.menu entryconfigure last -label "Open...Ctrl+O" \
       -command { cam_LoadCamera; cam_CreateEdit }

   $camera_edit.right.menubar.file.menu add command
   $camera_edit.right.menubar.file.menu entryconfigure last -label "New..." -command cam_CreateCamera

   $camera_edit.right.menubar.file.menu add command
   $camera_edit.right.menubar.file.menu entryconfigure last -label "Save..." -command cam_SaveCameras

   $camera_edit.right.menubar.file.menu add command
   $camera_edit.right.menubar.file.menu entryconfigure last -label "Quit..." -command "destroy $camera_edit"

   menubutton $camera_edit.right.menubar.edit -menu $camera_edit.right.menubar.edit.menu -text "Edit" -underline 0
   menu $camera_edit.right.menubar.edit.menu
 
   do i 0 [@ $camera_n-1] {
       $camera_edit.right.menubar.edit.menu add radiobutton
       $camera_edit.right.menubar.edit.menu entryconfigure last -variable camera_current \
                -command "cam_CreateEdit $camera_names($i)" -label $camera_names($i) -value $i
   }

   frame $camera_edit.right.viewport_x -bd 10 -relief ridge
   entry $camera_edit.right.viewport_x.low -relief sunken -width 5 -textvariable camera(view_low_x);
   label $camera_edit.right.viewport_x.low_lab  -text "Viewport Low X" 

   entry $camera_edit.right.viewport_x.high -relief sunken -width 5 -textvariable camera(view_high_x);
   label $camera_edit.right.viewport_x.high_lab -text "Viewport High X" 

   frame $camera_edit.right.viewport_y -bd 10 -relief ridge
   entry $camera_edit.right.viewport_y.low -relief sunken -width 5 -textvariable camera(view_low_y)
   label $camera_edit.right.viewport_y.low_lab  -text "Viewport Low Y" 

   entry $camera_edit.right.viewport_y.high -relief sunken -width 5 -textvariable camera(view_high_y)
   label $camera_edit.right.viewport_y.high_lab -text "Viewport High Y" 

   frame $camera_edit.right.buttons -bd 10 -relief ridge

   button $camera_edit.right.buttons.apply -text "Apply" -relief raised \
          -command "c_SetCamera"

   button $camera_edit.right.buttons.ok -text "OK" -relief raised \
          -command "c_SetCamera; destroy $camera_edit"

   button $camera_edit.right.buttons.cancel -text "Cancel" -relief raised \
          -command "destroy $camera_edit" 

   pack $camera_edit.left.from.lab_x -side top -fill x
   pack $camera_edit.left.from.ent_x -side top

   pack $camera_edit.left.from.lab_y -side top -fill x
   pack $camera_edit.left.from.ent_y -side top

   pack $camera_edit.left.from.lab_z -side top -fill x
   pack $camera_edit.left.from.ent_z -side top

   pack $camera_edit.left.to.lab_x  -side top   -fill x
   pack $camera_edit.left.to.ent_x  -side top

   pack $camera_edit.left.to.lab_y  -side top   -fill x
   pack $camera_edit.left.to.ent_y  -side top

   pack $camera_edit.left.to.lab_z  -side top   -fill x
   pack $camera_edit.left.to.ent_z  -side top

   pack $camera_edit.left.up.lab_x  -side top   -fill x
   pack $camera_edit.left.up.ent_x  -side top

   pack $camera_edit.left.up.lab_y  -side top   -fill x
   pack $camera_edit.left.up.ent_y  -side top

   pack $camera_edit.left.up.lab_z  -side top   -fill x
   pack $camera_edit.left.up.ent_z  -side top

   pack $camera_edit.left.from -side top -fill x
   pack $camera_edit.left.to   -side top -fill x
   pack $camera_edit.left.up   -side top -fill x

   pack $camera_edit.right.projection.label -side left
   pack $camera_edit.right.projection.perspective -side left
   pack $camera_edit.right.projection.ortho       -side left

   if { $camera(projection) == "perspective" } {

      pack $camera_edit.right.field.ent  -side left
      pack $camera_edit.right.field.lab  -side left

      pack $camera_edit.right.near.ent   -side left
      pack $camera_edit.right.near.lab   -side left

      pack $camera_edit.right.far.ent    -side left
      pack $camera_edit.right.far.lab    -side left

   } else {

      pack $camera_edit.right.left.ent   -side left
      pack $camera_edit.right.left.lab   -side left

      pack $camera_edit.right.right.ent  -side left
      pack $camera_edit.right.right.lab  -side left

      pack $camera_edit.right.bottom.ent -side left
      pack $camera_edit.right.bottom.lab -side left

      pack $camera_edit.right.top.ent    -side left
      pack $camera_edit.right.top.lab    -side left

      pack $camera_edit.right.near.ent   -side left
      pack $camera_edit.right.near.lab   -side left

      pack $camera_edit.right.far.ent    -side left
      pack $camera_edit.right.far.lab    -side left

   }

   pack $camera_edit.right.menubar.file  -side left -fill x
   pack $camera_edit.right.menubar.edit  -side left -fill x

   pack $camera_edit.right.status.label -side left -fill x
   pack $camera_edit.right.status.status -side left -fill x
   pack $camera_edit.right.status.on -side left -fill x
   pack $camera_edit.right.status.off -side left -fill x

   pack $camera_edit.right.viewport_x.low -side left
   pack $camera_edit.right.viewport_x.low_lab -side left

   pack $camera_edit.right.viewport_x.high -side left
   pack $camera_edit.right.viewport_x.high_lab -side left

   pack $camera_edit.right.viewport_y.low -side left
   pack $camera_edit.right.viewport_y.low_lab -side left

   pack $camera_edit.right.viewport_y.high -side left
   pack $camera_edit.right.viewport_y.high_lab -side left

   pack $camera_edit.right.menubar -side top -fill x
   pack $camera_edit.right.status -side top -fill x
   pack $camera_edit.right.projection -side top -fill x

   if { $camera(projection) == "perspective" } {

       pack $camera_edit.right.field -side top -fill x
       pack $camera_edit.right.near  -side top -fill x
       pack $camera_edit.right.far   -side top -fill x

   } else {

       pack $camera_edit.right.left   -side top -fill x
       pack $camera_edit.right.right  -side top -fill x
       pack $camera_edit.right.bottom -side top -fill x
       pack $camera_edit.right.top    -side top -fill x
       pack $camera_edit.right.near   -side top -fill x
       pack $camera_edit.right.far    -side top -fill x

   }

   pack $camera_edit.right.viewport_x -side top -fill x
   pack $camera_edit.right.viewport_y -side top -fill x

   pack $camera_edit.right.buttons.cancel -side right -fill x
   pack $camera_edit.right.buttons.ok     -side right -fill x
   pack $camera_edit.right.buttons.apply  -side right -fill x

   pack $camera_edit.right.buttons -side top -fill both

   pack $camera_edit.left  -side left -fill y
   pack $camera_edit.right -side left -fill y
}

proc cam_CreateCamera { name } {

   global camera
  
   set current $name

   set camera(name) $current

   set camera(from_x) 0.0
   set camera(from_y) 0.0
   set camera(from_z) 5.0

   set camera(to_x)   0.0
   set camera(to_y)   0.0
   set camera(to_z)   0.0

   set camera(up_x)   0.0
   set camera(up_y)   1.0
   set camera(up_z)   0.0

   set camera(str_t)  0.03
   set camera(str_r)  5.00

   set camera(projection) perspective

   set camera(field_angle)       30.0
   set camera(near)   0.1
   set camera(far)   20.0

   if { $camera(projection) == "ortho" } {
	  set camera(left)        -1.0
	  set camera(right)        1.0
	  set camera(bottom)      -1.0
	  set camera(top)          1.0
      set camera(near)      -100.0
      set camera(far)        100.0
   }

   set camera(view_low_x)         0.0
   set camera(view_low_y)         0.0
   set camera(view_high_x)        1.0
   set camera(view_high_y)        1.0


   cam_CreateEdit "$current"
}

proc cam_LoadCamera {} {
    global ELMER_POST_HOME log

    #
    # set name [fs_FileSelect $ELMER_POST_HOME/lib/cameras *.cam]
    #

    #   Type names              Extension(s)    Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
        {"Camera files"            {.cam}          TEXT}
        {"All files"            *}
    }

    set name [tk_getOpenFile -filetypes $types -initialdir $ELMER_POST_HOME/lib/cameras \
	    -parent .camera_edit -defaultextension .cam]

    if { $name != "" } { c_LoadCamera $name }
}

proc cam_SaveCameras {} {
   global camera camera_n camera_names camera_current

   set types {
       {"Camera files"            {.cam}          TEXT}
       {"All files"            *}
   }

   set name [tk_getSaveFile -filetypes $types -parent .camera_edit -defaultextension .cam]
   set file [open $name w]

   do i 0 [@ $camera_n-1] {
	  echo "Saving camera: $camera_names($i)"

      c_CurrentCamera $camera_names($i)

      puts $file "camera $camera_names($i)"
      puts $file "projection $camera(projection)"
      puts $file "from $camera(from_x) $camera(from_y) $camera(from_z)"
      puts $file "up $camera(up_x) $camera(up_y) $camera(up_z)"
      puts $file "viewport $camera(view_low_x) $camera(view_low_y) $camera(view_high_x) $camera(view_high_y)"
      puts $file "field angle $camera(field_angle)"
      puts $file "clip $camera(near) $camera(far)"
      puts $file " "
   }

   c_CurrentCamera $camera_names($camera_current)

   close $file

   echo "Created camera file $name"
}
