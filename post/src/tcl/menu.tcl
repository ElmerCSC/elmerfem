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
#* Command line interface & main menu definitions for ElmerPost
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

if { [winfo exists .menubar] } { puts uuh; destroy .menubar }
if { [winfo exists .tbox] } { destroy .tbox }
if { [winfo exists .lbox] } { destroy .lbox }
if { [winfo exists .cbox] } { destroy .cbox }

#wm minsize . 10 3

frame .menubar

menubutton .menubar.file -menu .menubar.file.menu -text "File" -underline 0
menu .menubar.file.menu

set ModelFileName ""
.menubar.file.menu add command  -underline 0
.menubar.file.menu entryconfigure last -label "Open ... Ctrl+O" -command { ReadFile }
bind . <Control-o> ReadFile

.menubar.file.menu add command  -underline 0
.menubar.file.menu entryconfigure last -label "Save Image" -command { SaveImage }

#.menubar.file.menu add command
#.menubar.file.menu entryconfigure last -label "  Write  " -command {fs_FileSelect [pwd] *}

.menubar.file.menu add command  -underline 0
.menubar.file.menu entryconfigure last -label "Load Sicopolis ... " -command { LoadSicopolis }

.menubar.file.menu add command -underline 0
.menubar.file.menu entryconfigure last -label "Quit ... Ctrl+Q  " -command "Exit"
bind . <Control-q> Exit

pack .menubar.file -side left

#
# EDIT MENU 
#

menubutton .menubar.edit -menu .menubar.edit.menu -text "Edit"
menu .menubar.edit.menu

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "ColorMap ..." -command "ci_ColorMap Colormap"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Grouping ..." -command "group_edit"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Material ..." -command "me_Edit Material 0 material"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Background ..." -command "me_Edit Background 1 background"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Timestep Control ..." -command "time_Edit"

.menubar.edit.menu add command
.menubar.edit.menu entryconfigure last -label "Math Module Window ..." -command "math_win"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Stop Processing ..." -command "StopProcessing"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Camera Settings ..." -command "cam_CreateEdit"

.menubar.edit.menu add command -underline 0
.menubar.edit.menu entryconfigure last -label "Object Clip Planes..." -command clip_Edit

pack .menubar.edit -side left

#
#  DISPLAY MENU
#
menubutton .menubar.disp -menu .menubar.disp.menu -text "Display" -underline 0
menu .menubar.disp.menu

.menubar.disp.menu add checkbutton -underline 0
.menubar.disp.menu entryconfigure last -label "Mesh Lines ..." -variable DisplayStyle(MeshLines) -command "UpdateObject"

.menubar.disp.menu add checkbutton -underline 0
.menubar.disp.menu entryconfigure last -label "Color Mesh ... Ctrl+M" -variable DisplayStyle(ColorMesh)         \
         -command { UpdateObject; if { $DisplayStyle(ColorMesh) } mesh_edit }
bind . <Control-m> { UpdateObject; if { $DisplayStyle(ColorMesh) } { set DisplayStyle(ColorMesh) 0 } \
      else { set DisplayStyle(ColorMesh) 1; mesh_edit } }

.menubar.disp.menu add checkbutton -underline 8
.menubar.disp.menu entryconfigure last -label "Contour Lines ..." -variable DisplayStyle(Contours)       \
         -command { UpdateObject; if { $DisplayStyle(Contours) } contour_edit }

.menubar.disp.menu add checkbutton -underline 8
.menubar.disp.menu entryconfigure last -label "Contour Surfaces ..." -variable DisplayStyle(Isosurfaces) \
         -command { UpdateObject; if { $DisplayStyle(Isosurfaces) } isosurface_edit }

.menubar.disp.menu add checkbutton -underline 0
.menubar.disp.menu entryconfigure last -label "Vectors ..." -variable DisplayStyle(Vectors)     \
         -command { UpdateObject; if { $DisplayStyle(Vectors) } vector_edit }

.menubar.disp.menu add checkbutton -underline 0
.menubar.disp.menu entryconfigure last -label "Particles ..." -variable DisplayStyle(Particles) \
         -command { UpdateObject; if { $DisplayStyle(Particles) } particle_edit }

.menubar.disp.menu add checkbutton
.menubar.disp.menu entryconfigure last -label "Spheres ..." -variable DisplayStyle(Spheres)      \
         -command { UpdateObject; if { $DisplayStyle(Spheres) } sphere_edit }

.menubar.disp.menu add checkbutton
.menubar.disp.menu entryconfigure last -label "ColorScale ..." -variable DisplayStyle(ColorScale)      \
         -command { UpdateObject; if { $DisplayStyle(ColorScale) } colscale_edit }
pack .menubar.disp -side left

.menubar.disp.menu add checkbutton
.menubar.disp.menu entryconfigure last -label "Mode Display ..." -variable DisplayStyle(ModeDisplay)      \
         -command { UpdateObject; if { $DisplayStyle(ModeDisplay) } ModeDisplay }
pack .menubar.disp -side left


set DisplayStyle(MeshLines)   0
set DisplayStyle(ColorMesh)   0
set DisplayStyle(Contours)    0
set DisplayStyle(Vectors)     0
set DisplayStyle(Isosurface)  0
set DisplayStyle(Spheres)     0
set DisplayStyle(ColorScale)  0
set DisplayStyle(ModeDisplay)  0

#menubutton .menubar.options -menu .menubar.options.menu -text "Options" -underline 0
#menu .menubar.options.menu
#.menubar.options.menu add checkbutton
#.menubar.options.menu entryconfigure last -label "Construct Volume Element Sides..." \
#          -variable GlobalOptions(VolumeSides)
#pack .menubar.options -side left

menubutton .menubar.help -menu .menubar.help.menu -text "Help" -underline 0
menu .menubar.help.menu 

.menubar.help.menu add command -underline 0
.menubar.help.menu entryconfigure last -label "Help..." -command { help }

#.menubar.help.menu add command -underline 0
#.menubar.help.menu entryconfigure last -label "About Elmerpost..." \
#   -command { message "\nElmerpost Version 0.9 (12/2/97)\n\nAuthor: Juha Ruokolainen\n\nCopyright: CSC - IT Center for Science Ltd. (CSC)\n                  P.O.Box. 405, 02320 Espoo, Finland\n" "About ElmerPost" }
#.menubar.help.menu add command -underline 0
#.menubar.help.menu entryconfigure last -label "About the Elmer project..." -command {  help http://www.csc.fi/CFD/elmer/ }
pack .menubar.help -side left
# -side right -fill x

frame .buts

#image create bitmap QuitImage -file \
#    /mnt/mds/proj2/proj/Nov97/elmer/jpr/elmerpost_src/tcl7.6/tk4.2/bitmaps/questhead.bmp
#button .buts.quit -image QuitImage -command Exit
#pack .buts.quit -side left -fill x

image create photo ReadImage -file $ELMER_POST_HOME/lib/images/file.gif
button .buts.read -image ReadImage -command ReadFile
pack .buts.read -side left -fill x

image create photo ColorMeshImage -file $ELMER_POST_HOME/lib/images/cmesh.gif
checkbutton .buts.cmesh -image ColorMeshImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(ColorMesh) -command { UpdateObject; if { $DisplayStyle(ColorMesh) } mesh_edit }
pack .buts.cmesh -side left -fill x

image create photo ContourImage -file $ELMER_POST_HOME/lib/images/cont.gif
checkbutton .buts.cont -image ContourImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(Contours) -command { UpdateObject; if { $DisplayStyle(Contours) } contour_edit }
pack .buts.cont -side left -fill x

image create photo IsoContourImage -file $ELMER_POST_HOME/lib/images/isocont.gif
checkbutton .buts.isocont -image IsoContourImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(Isosurfaces) -command { UpdateObject; if { $DisplayStyle(Isosurfaces) } isosurface_edit }
pack .buts.isocont -side left -fill x

image create photo VectorImage -file $ELMER_POST_HOME/lib/images/vect.gif
checkbutton .buts.vect -image VectorImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(Vectors) -command { UpdateObject; if { $DisplayStyle(Vectors) } vector_edit }
pack .buts.vect -side left -fill x

image create photo ParticleImage -file $ELMER_POST_HOME/lib/images/part.gif
checkbutton .buts.part -image ParticleImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(Particles) -command { UpdateObject; if { $DisplayStyle(Particles) } particle_edit }
pack .buts.part -side left -fill x

image create photo ColscaleImage -file $ELMER_POST_HOME/lib/images/colscale.gif
checkbutton .buts.colscale -image ColscaleImage -indicatoron false -border 5 -selectcolor blue \
 -variable DisplayStyle(ColorScale) -command { UpdateObject; if { $DisplayStyle(ColorScale) } colscale_edit }
pack .buts.colscale -side left -fill x

#set PlayStatus "Play"
#button .buts.play -textvariable PlayStatus -command "UpdateObject; display"
image create photo PlayImage -file $ELMER_POST_HOME/lib/images/play.gif
button .buts.play -image PlayImage -back green -command "UpdateObject; display"
pack .buts.play -side left -fill x

image create photo RMX -file $ELMER_POST_HOME/lib/images/rmx.gif
image create photo RPX -file $ELMER_POST_HOME/lib/images/rpx.gif
image create photo RMY -file $ELMER_POST_HOME/lib/images/rmy.gif
image create photo RPY -file $ELMER_POST_HOME/lib/images/rpy.gif

frame .rots
button .rots.rmx -image RMX -command "rotate -r -x -5"
button .rots.rpx -image RPX -command "rotate -r -x  5"

button .rots.rmy -image RMY -command "rotate -r -y -5"
button .rots.rpy -image RPY -command "rotate -r -y  5"

button .rots.rmz -text "-Z" -command "rotate -r -z -5"
button .rots.rpz -text "+Z" -command "rotate -r -z  5"

frame .trans

image create photo Left      -file $ELMER_POST_HOME/lib/images/left.gif
image create photo Right     -file $ELMER_POST_HOME/lib/images/right.gif
image create photo Up        -file $ELMER_POST_HOME/lib/images/up.gif
image create photo Down      -file $ELMER_POST_HOME/lib/images/down.gif
image create photo ScaleUp   -file $ELMER_POST_HOME/lib/images/scaleup.gif
image create photo ScaleDown -file $ELMER_POST_HOME/lib/images/scaledown.gif

button .trans.tmx -image Left -command "translate -r -x -0.1"
button .trans.tpx -image Right -command "translate -r -x  0.1"

button .trans.tmy -image Down   -command "translate -r -y -0.1"
button .trans.tpy -image Up -command "translate -r -y  0.1"

button .trans.tmz -text "-Z" -command "translate -r -z -0.1"
button .trans.tpz -text "+Z" -command "translate -r -z  0.1"

button .rots.scd -image ScaleDown -command "scale -r -0.1"
button .rots.scu -image ScaleUp   -command "scale -r  0.1"
button .rots.res -text "RESET" -command "rotate 0 0 0; scale 1; translate 0 0 0"

menubutton .trans.rpri -menu .trans.rpri.menu -text "Rot PRI" -relief raised
menu .trans.rpri.menu
set RotationPriority 6
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "xyz" -value 0 -variable RotationPriority -command { rotpriority xyz }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "xzy" -value 1 -variable RotationPriority -command { rotpriority xzy }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "yxz" -value 2 -variable RotationPriority -command { rotpriority yxz }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "yzx" -value 3 -variable RotationPriority -command { rotpriority yzx }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "zxy" -value 4 -variable RotationPriority -command { rotpriority zxy }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "zyx" -value 5 -variable RotationPriority -command { rotpriority zyx }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "local" -value 6 -variable RotationPriority -command { rotpriority local }
.trans.rpri.menu add radiobutton
.trans.rpri.menu entryconfigure last -label "global" -value 7 -variable RotationPriority -command { rotpriority global }

menubutton .trans.tpri -menu .trans.tpri.menu -text "Trn PRI" -relief raised
menu .trans.tpri.menu
set TransformationPriority 0
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "trs" -value 0 -variable TransformationPriority -command { trnpriority trs }
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "tsr" -value 1 -variable TransformationPriority -command { trnpriority tsr }
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "rts" -value 2 -variable TransformationPriority -command { trnpriority rts }
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "rst" -value 3 -variable TransformationPriority -command { trnpriority rst }
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "str" -value 4 -variable TransformationPriority -command { trnpriority str }
.trans.tpri.menu add radiobutton
.trans.tpri.menu entryconfigure last -label "srt" -value 5 -variable TransformationPriority -command { trnpriority srt }



pack .rots.rmx .rots.rpx .rots.rmy .rots.rpy .rots.rmz .rots.rpz -side left
pack .trans.tmx .trans.tpx .trans.tmy .trans.tpy .trans.tmz .trans.tpz -side left
pack .rots.scd .rots.scu .rots.res -side left
pack .trans.rpri -side left
pack .trans.tpri -side left

frame .options
checkbutton .options.keep -indicatoron false -border 5 -selectcolor green -variable KeepScale -text "Freeze Scaling"
checkbutton .options.norm -indicatoron false -border 5 -selectcolor green -variable NormalUpdate -text "Update Normals"
pack .options.keep .options.norm -side left

frame .lbox -bd 5 -bg lightblue
frame .tbox -bd 5 -bg lightblue
frame .cbox
label .cbox.lab -text "Elmer-Post: "
entry .cbox.cmd -relief sunken

text .tbox.text -yscroll ".tbox.scroll set" -wrap none
scrollbar .tbox.scroll -orient vertical -command ".tbox.text yview"
pack .tbox.scroll -side left -fill y
pack .tbox.text -side left -fill both -expand 1 

listbox .lbox.list -yscroll ".lbox.scroll set"
scrollbar .lbox.scroll -orient vertical -command ".lbox.list yview"
pack .lbox.scroll -side left -fill y
pack .lbox.list -side left -fill both -expand 1 

set main_out .tbox.text
set main_log .lbox.list
set main_cmd .cbox.cmd

bind $main_log <Button-1> {
    $main_cmd delete 0 end
    set main_cur [$main_log curselection]
    if { $main_cur != "" } {
        $main_cmd insert 0 [$main_log get $main_cur]
    }
}

bind $main_log <Double-Button-1> {
    $main_cmd delete 0 end

    set main_cur [$main_log curselection]
    if { $main_cur != "" } {
        $main_cmd insert 0 [$main_log get $main_cur]
        CommandReturn $main_cmd $main_log $main_out
    }
}

proc CommandReturn {w l t} {
    if { $l != "" } {
        $l insert end [$w get]
        history add [$w get]

        set n [$l size]

        set k [@ $n-5]
        if { $k < 0 } { set k 0}
        $l yview $k

        $l select clear 0 end
        $l select set end end

        if { $t != "" } {
            $t insert end [uplevel [$w get]]\n

            set n [split [$t index end] .]
            set row [lindex $n 0]

            set h [$t cget -height]
            set row [@ $row-$h]

            $t yview $row
        }
    } else {
        uplevel [$w get]]\n
    }

    $w delete 0 end
}

proc CommandUpKey {w l} {
   if { $l != "" } {
       set n [$l curselection]
       if { $n == "" } {
           $l select clear 0 end
           $l select set end end
           set n [$l curselection]
       }

       set n [@ $n-1]
       if { $n < 0 } { set n 0 }

        $l select clear 0 end
        $l select set $n $n

        set k [@ $n-5]
        if { $k < 0 } { set k 0 }
        $l yview $k
   }

   $w delete 0 end
   if { $l != "" } { $w insert 0 [$l get $n] }
}

proc CommandDownKey {w l} {
   if { $l != "" } {
       set n [$l curselection]
       if { $n == "" } {
           $l select clear 0 end
           $l select set end end
           $l yview  end
           set n [$l curselection]
       }
       set n [@ $n+1]
       if { $n >= [$l size] } { set $n [@ [$l size]-1] }

       $l select clear 0 end
       $l select set $n $n

       set k [@ $n-5]
       if { $k < 0 } { set k 0 }
       $l yview $k
   }

   $w delete 0 end
   if { $l != "" } { $w insert 0 [$l get $n] }
}

proc CommandLeftKey {w} {
   set t [$w index insert]
   set t [@ $t-1]

   if { $t < 0 } { set t 0 }

   $w icursor $t
}

proc CommandRightKey {w} {
   set t [$w index insert]
   set t [@ $t+1]

   if {$t >= [$w index end]} {set t [$w index end]}

   $w icursor $t
}

proc CommandToEnd {w} { $w icursor end }

proc CommandToBegin {w} { $w icursor 0 }

proc CommandDelete {w} { $w delete [$w index insert] }

proc CommandDeleteEnd {w} { $w delete [$w index insert] end }

bind $main_cmd <Return>    {CommandReturn %W $main_log $main_out}

bind $main_cmd   <Up>        {CommandUpKey %W $main_log}
bind $main_cmd   <Down>      {CommandDownKey %W $main_log}
#bind Entry <Left>      {CommandLeftKey %W}
#bind Entry <Right>     {CommandRightKey %W}

#bind Entry <Control-b> {CommandLeftKey %W}
#bind Entry <Control-f> {CommandRightKey %W}
bind $main_cmd <Control-p>   {CommandUpKey %W $main_log}
bind $main_cmd <Control-n>   {CommandDownKey %W $main_log}

#bind Entry <Control-a> {CommandToBegin %W}
#bind Entry <Control-e> {CommandToEnd %W}
#bind Entry <Control-d> {CommandDelete %W}
#bind Entry <Control-k> {CommandDeleteEnd %W}

proc TextUpKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set row [@ $row-1]

   $w yview -pickplace $row
   $w mark set insert $row.$col
}

proc TextDownKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set row [@ $row+1]

   $w yview -pickplace $row
   $w mark set insert $row.$col
}

proc TextLeftKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set col [@ $col-1]
   $w mark set insert $row.$col
}

proc TextRightKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set col [@ $col+1]
   $w mark set insert $row.$col
}

proc TextHomeKey {w} {
    $w yview -pickplace 0
    $w mark set insert 0.0
}

proc TextEndKey {w} {
    $w yview -pickplace end
    $w mark set insert end
}

proc TextPageDownKey {w} {
    set h [$w cget -height]
    set h [@ $h-1]

    set t [split [$w index insert] .]

    set row [lindex $t 0]
    set col [lindex $t 1]

    set row [@ $row+$h]

    $w yview -pickplace $row
    $w mark set insert $row.$col
}

proc TextPageUpKey {w} {
    set h [$w cget -height]
    set h [@ $h-1]

    set t [split [$w index insert] .]

    set row [lindex $t 0]
    set col [lindex $t 1]

    set row [@ $row-$h]

    $w yview -pickplace $row
    $w mark set insert $row.$col
}

proc TextToBegin {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col 0

   $w mark set insert $row.$col
}

proc TextToEnd {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col 0
   $w mark set here $row.$col

   set r [@ $row+1]
   $w mark set insert $r.$col

   set str [$w get here insert]

   set col [@ [string length $str]-1]

   $w mark set insert $row.$col
   $w mark unset here
}

proc TextDelete {w} { $w delete insert }

proc Nil {} {}

#bind Text <Up>        {TextUpKey %W}
#bind Text <Down>      {TextDownKey %W}
#bind Text <Left>      {TextLeftKey %W}
#bind Text <Right>     {TextRightKey %W}
#bind Text <Home>      {TextHomeKey %W}
#bind Text <End>       {TextEndKey %W}

#bind Text <Control-b> {TextLeftKey %W}
#bind Text <Control-f> {TextRightKey %W}
#bind Text <Control-p> {TextUpKey %W}
#bind Text <Control-n> {TextDownKey %W}

bind Text <Escape> {Nil}
#bind Text <Escape>"\<" {TextEndKey %W}
#bind Text <Escape>"\>" {TextHomeKey %W}
bind Text <Escape>v     {TextPageUpKey %W}
bind Text <Control-v>   {TextPageDownKey %W}

#bind Text <Control-a> {TextToBegin %W}
#bind Text <Control-e> {TextToEnd %W}
#bind Text <Control-d> {TextDelete %W}

$main_out configure -setgrid 1 -width 55 -height 5 -background white
$main_log configure -setgrid 1 -width 55 -height 5
$main_cmd configure -width 55

pack .menubar -side top  -fill x

pack .buts  -side top
pack .rots  -side top
pack .options -side top
pack .trans -side top

pack .lbox -side top -fill both -expand 1
pack .tbox -fill both -expand 1

pack .cbox.lab -side left -expand 1
pack .cbox.cmd -side right -expand 1 -fill both
pack .cbox  -side bottom

#
#tk_menuBar .menubar .menubar.file .menubar.edit
#bind all <Control-c>  { puts here; return }
#

$main_out insert end "\n\nWelcome to Elmer Post Processing\n"
$main_out insert end "Have Fun! \n\n\n"

focus $main_cmd

