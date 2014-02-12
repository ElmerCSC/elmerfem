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

set ModelInfo(StartTime) 1
set ModelInfo(EndTime) 1
set ModelInfo(IncTime) 1

proc ConvertFidap { name } {
    global reply ModelFileName ELMER_POST_HOME FidapReadDone

    if { $ModelFileName != "" } {
       set reply yes;

       if { [file exists $ModelFileName] } {
         set reply [tk_messageBox -icon warning -type yesno \
             -message "File \"$ModelFileName\" already exists.\nDo you want to overwrite it?"]
       }

      if { $reply == "yes" } {
          math "$ exec $ELMER_POST_HOME/bin/fdneut2ep $name $ModelFileName; tcl(\"set FidapReadDone 1\");"
      }
    }
}

proc ReadGetFidapFile { name } {
    global ModelFileName TempName FidapReadDone
 
    set TempName $name

    set w .fidap

    toplevel $w
    place_window $w

    label $w.text1 -text "This seems to be a FIDAP NEUTRAL file!" 
    label $w.text2 -text "Wanna make conversion ?"

    pack $w.text1 -side top
    pack $w.text2 -side top

    frame $w.f
    label $w.f.filelab -text "Save as: "
    entry $w.f.filenam -width 40 -textvariable ModelFileName
    button $w.f.filebut -text "Browse..." -command { set ModelFileName [tk_getSaveFile -parent .fidap] }

    pack $w.f.filelab $w.f.filenam $w.f.filebut -side left
    pack $w.f -side top

    frame  $w.buttons
    button $w.buttons.ok    -text "OK" -command "ConvertFidap $TempName"
    button $w.buttons.close -text "Close" -command { set FidapReadDone 1 }

    pack $w.buttons.ok $w.buttons.close -side left
    pack $w.buttons -side top

#   listbox $w.lbox
#   pack $w.lbox -side top

    bind $w <Return> "ConvertFidap $TempName"

    tkwait variable FidapReadDone
    destroy $w
}

proc ReadFileHeader { FileName } {
   global ModelInfo ModelFileStatus ModelFileName ELMER_POST_HOME

   if { $FileName == "" } {
       set ModelFileName [tk_getOpenFile -parent .read -title "Open File"];
       set FileName $ModelFileName

       if { $FileName == "" } return 
    }

    if { [file  extension $FileName] == ".FDNEUT" } {

       set ModelFileStatus "Converting FIDAP Neutral File..."

       ReadGetFidapFile $FileName
       set FileName $ModelFileName

       if { $FileName == "" } return 
    }

   set file [open $FileName "r"] 
   set header [gets $file]

   scan $header "%d %d %d %d" ModelInfo(Nodes) \
     ModelInfo(Elements) ModelInfo(DOFS) ModelInfo(Timesteps) 

   set str [gets $file]

   set ModelInfo(Description) ""
   if { [string first "#title " $str] >= 0 } {
      set ModelInfo(Description) $str
   }

#
#  for { set i 0 } { $i<$ModelInfo(Nodes) } { incr i } {
#     gets $file
#  }
#
#  set ModelInfo(GroupNames) ""
#  for { set i 0 } { $i<$ModelInfo(Elements) } { incr i } {
#     set str [gets $file]
#     if { [string index $str 0] == "#" } {
#       set name ""
#       scan $str "#group %s" name
#       if { $name != "" } {
#         set ModelInfo(GroupNames) $ModelInfo(GroupNames)$name
#         set ModelInfo(GroupNames) $ModelInfo(GroupNames)\n
#       }
#       set i [@ $i-1]
#     }
#  }
#

   close $file

   set ModelInfo(DOFNames) ""
   set v [string first "vector:" $header]
   set s [string first "scalar:" $header]

   for { set i 0 } { $v>=0 || $s>=0 } { incr i } {
     if { $v >= 0 && ($v < $s || $s < 0) } {
       set header [string range $header [@ $v+7] [@ [string length $header]-1]]
       scan $header %s name
       set name "Vector:\t$name"
     } else {
       set header [string range $header [@ $s+7] [@ [string length $header]-1]]
       scan $header %s name
       set name "Scalar:\t$name"
     }
     set ModelInfo(DOFNames) $ModelInfo(DOFNames)$name
     set ModelInfo(DOFNames) $ModelInfo(DOFNames)\n

     set v [string first "vector:" $header]
     set s [string first "scalar:" $header]
   }

   if { $ModelInfo(StartTime) < 1 } { set ModelInfo(StartTime) 1 }

   if { $ModelInfo(StartTime) > $ModelInfo(Timesteps) } {
      set ModelInfo(StartTime) $ModelInfo(Timesteps)
   }

   if { $ModelInfo(EndTime) < $ModelInfo(StartTime) } {
      set ModelInfo(EndTime) $ModelInfo(StartTime)
   }

   if { $ModelInfo(EndTime) > $ModelInfo(Timesteps) } {
      set ModeInfo(EndTime) $ModelInfo(Timesteps)
   }

   if { $ModelInfo(EndTime) > $ModelInfo(Timesteps) } {
      set ModelInfo(EndTime) $ModelInfo(Timesteps)
   }

   if { $ModelInfo(IncTime) < 1 } { set ModelInfo(IncTime) 1 }

   set ModelFileStatus "Header Read"
}

proc readfile { File {StartTime 1} {EndTime 1} {IncTime 1}  }   {
  global ModelFileName ModelFileStatus ModelInfo

  set ModelInfo(StartTime) $StartTime
  set ModelInfo(EndTime)   $EndTime
  set ModelInfo(IncTime)   $IncTime

  ReadFileHeader $File

  set ModelFileStatus "Reading"
  update

  cReadFile $File $StartTime $EndTime $IncTime
  catch colscale_update;
  catch mesh_update;
  catch isosurface_update;

  set ModelFileName $File
  set ModelFileStatus "Read Done"
  update

  if { [winfo exists .time_edit] } { time_Edit }

  math who
}

proc DoReadFile { FileName } {

    global ModelFileStatus ModelInfo ModelFileName

    if { $FileName == "" } {
       set ModelFileName [tk_getOpenFile -parent .read -title "Open File"];
       set FileName $ModelFileName

       if { $FileName == "" } return 
    }

    ReadFileHeader $FileName

    readfile $FileName $ModelInfo(StartTime) $ModelInfo(EndTime) \
                       $ModelInfo(IncTime)
}

proc ReadFile {} {
    global main_out main_log GlobalOptions ModelFileName ModelFileStatus
    global ModelInfo 

    set w .read

    if { [winfo exists $w] } {
      wm iconify $w
      wm deiconify $w
      return
    } else {
      toplevel $w
      place_window $w
      wm title $w "Read Model File"
    }

    frame $w.top



    set ModelFileStatus "Not Done"
    label $w.top.lab -text "Status: " -font "Helvetica-Bold 12"
#   -font -Adobe-Times-Medium-R-Normal-*-180-*
    label $w.top.stat -textvariable ModelFileStatus -relief sunken -font "Helvetica-Bold 12"
#    -Adobe-Times-Medium-R-Normal-*-180-*

    pack $w.top.stat $w.top.lab -side right
    pack $w.top -side top -expand 1 -fill x

    label $w.sp0 -text "\nOptions:\n" -font "Helvetica-Bold 12"
#   -font -Adobe-Times-Medium-R-Normal-*-180-*
    pack $w.sp0 -side top

    frame $w.sides -relief sunken

    checkbutton $w.sides.surf -text "Generate Surface Element Sides" \
             -variable GlobalOptions(SurfaceSides)

    checkbutton $w.sides.vol -text "Generate Volume Element Sides" \
             -variable GlobalOptions(VolumeSides)

    checkbutton $w.sides.vol_edge -text "Generate Volume Element Edges" \
             -variable GlobalOptions(VolumeEdges)

    pack $w.sides.surf -side top -expand 1 -fill x
    pack $w.sides.vol  -side top -expand 1 -fill x
    pack $w.sides.vol_edge  -side top -expand 1 -fill x
    pack $w.sides -side top 

    label $w.sp1 -text "\nFile Information:\n" -font "Helvetica-Bold 12"
# -font -Adobe-Times-Medium-R-Normal-*-180-*
    pack $w.sp1 -side top 

    label $w.desc -textvariable ModelInfo(Description)
    pack $w.desc -side top 

    frame $w.nodeinfo -relief sunken
    label $w.nodeinfo.nodl  -text "Nodes:\t\t"
    entry $w.nodeinfo.nodes -textvariable ModelInfo(Nodes) -state disabled -relief flat -back white

    pack $w.nodeinfo.nodl $w.nodeinfo.nodes -side left
    pack $w.nodeinfo -side top

    frame $w.eleminfo -relief sunken
    label $w.eleminfo.elel  -text "Elements:\t"
    entry $w.eleminfo.elems -textvariable ModelInfo(Elements) -state disabled -relief flat -back white

    pack $w.eleminfo.elel $w.eleminfo.elems -side left
    pack $w.eleminfo -side top

    frame $w.timesinfo -relief sunken
    label $w.timesinfo.dofl  -text "Timestps:\t"
    entry $w.timesinfo.times -textvariable ModelInfo(Timesteps) -state disabled -relief flat -back white

    pack $w.timesinfo.dofl $w.timesinfo.times -side left
    pack $w.timesinfo -side top 

    frame $w.dofsinfo -relief sunken
    label $w.dofsinfo.dofl  -text "DOFS:\t\t"
    entry $w.dofsinfo.dofs -textvariable ModelInfo(DOFS) -state disabled -relief flat -back white
    label $w.dofnames -textvariable ModelInfo(DOFNames)

    pack $w.dofsinfo.dofl $w.dofsinfo.dofs -side left
    pack $w.dofsinfo $w.dofnames -side top 

    label $w.sp2 -text "\nSelect timesteps:\n" -font "Helvetica-Bold 12"
#    -font -Adobe-Times-Medium-R-Normal-*-180-*
    pack $w.sp2 -side top 
#
# Timestep
#
    frame $w.time -relief raised
    label $w.time.start_lab -text "First: "
    entry $w.time.start_ent -width 5 -textvariable ModelInfo(StartTime)

    label $w.time.end_lab -text "Last: "
    entry $w.time.end_ent -width 5 -textvariable ModelInfo(EndTime)

    label $w.time.inc_lab -text "Increment "
    entry $w.time.inc_ent -width 5 -textvariable ModelInfo(IncTime)

    button $w.time.but -text "All" -command \
       { set ModelInfo(StartTime) 1; set ModelInfo(EndTime) $ModelInfo(Timesteps); \
                               set ModelInfo(IncTime) 1; }

    pack $w.time.start_lab $w.time.start_ent -side left
    pack $w.time.end_lab $w.time.end_ent -side left
    pack $w.time.inc_lab $w.time.inc_ent -side left
    pack $w.time.but -side left
    pack $w.time -side top

#
# File name
#
    label $w.sp3 -text "\nSelect file:\n" -font "Helvetica-Bold 12"
# -font -Adobe-Times-Medium-R-Normal-*-180-*
    pack $w.sp3 -side top

    frame $w.file -relief sunken
    label $w.file.lab -text "Model file: "
    entry $w.file.name -width 40 -textvariable ModelFileName

    button $w.file.but -text "Browse..." -command \
       { set ModelFileName [tk_getOpenFile -parent .read -title "Read Model File"]; ReadFileHeader $ModelFileName  }

    pack $w.file.lab $w.file.name $w.file.but -side left
    pack $w.file -side top

    label $w.sp4 -text \n
    pack $w.sp4 -side top

    frame $w.buttons
    button $w.buttons.head -text "Read header" -command \
             { ReadFileHeader $ModelFileName }

    button $w.buttons.read -text "Read file" -command \
            { DoReadFile $ModelFileName }

    button $w.buttons.ok -text "OK" -command \
     { DoReadFile $ModelFileName; destroy .read; }

    button $w.buttons.close -text "Close" -command "destroy $w"

    pack $w.buttons.head $w.buttons.read $w.buttons.ok \
               $w.buttons.close -side left

    pack $w.buttons -side top

    bind $w <Return> { DoReadFile $ModelFileName }
}
