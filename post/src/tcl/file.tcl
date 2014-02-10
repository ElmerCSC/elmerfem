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
#* File Selector Utility Widget
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

set fs_file_pattern ""
set fs_file_dir ""
set fs_file_ok 0
set fs_file_read read

#
# main rutine of the widget, give initial directory and 
# file matching pattern as arguments
#
# 11 Sep 95
#
proc fs_FileSelect {directory pattern {readwrite read} } {
    global fs_file_pattern fs_file_dir fs_file_select fs_file_ok fs_file_read

    if { $directory != "" } { set fs_file_dir $directory }
    if { $pattern != "" }   { set fs_file_pattern $pattern }

    set fs_file_read $readwrite 

    if {[winfo exists .fs_file]} { destroy .fs_file }

    toplevel .fs_file
    place_window .fs_file
    wm minsize .fs_file 50 10
    wm title .fs_file "File Selection"

    frame .fs_file.pat_box  -bg lightblue; # frame for pattern
    frame .fs_file.dir_box  -bg lightblue; # frame for directory
    frame .fs_file.but_box  -bg lightblue; # frame for cancel,ok buttons
    frame .fs_file.name_box -bg lightblue; # frame for file name
    frame .fs_file.list_box -bg lightblue; # frame for directory list

    label .fs_file.name_box.lab -text "File name: "
    entry .fs_file.name_box.name -relief sunken

    label .fs_file.pat_box.lab -text "Pattern   : "
    entry .fs_file.pat_box.pattern -relief sunken

    label .fs_file.dir_box.lab -text  "Directory: "
    entry .fs_file.dir_box.directory -relief sunken

    button .fs_file.but_box.fs_file_ok -text "OK" -command {fs_NameGetSelect .fs_file.name_box.name} -bd 5
    button .fs_file.but_box.fs_file_cancel -text "CANCEL" -command {destroy .fs_file} -bd 5

    listbox .fs_file.list_box.list -yscroll ".fs_file.list_box.scroll set"
    scrollbar .fs_file.list_box.scroll  -command ".fs_file.list_box.list yview"

    pack .fs_file.list_box.scroll -side left -fill y
    pack .fs_file.list_box.list -side left -expand 1 -fill both

    pack .fs_file.name_box.lab -side left
    pack .fs_file.name_box.name -side right -expand 1 -fill x

    pack .fs_file.dir_box.lab -side left
    pack .fs_file.dir_box.directory -side right -expand 1 -fill x

    pack .fs_file.pat_box.lab -side left
    pack .fs_file.pat_box.pattern -side right -expand 1 -fill x

    pack .fs_file.but_box.fs_file_ok -side right
    pack .fs_file.but_box.fs_file_cancel -side right

    pack .fs_file.pat_box  -side top -fill x
    pack .fs_file.dir_box  -side top -fill x
    pack .fs_file.name_box -side top -fill x
    pack .fs_file.list_box -side top -fill x
    pack .fs_file.but_box  -side top -fill x

    bind .fs_file.list_box.list <Double-Button-1> {fs_ShowAndExecuteSelect %W}
    bind .fs_file.list_box.list <Button-1> {fs_ShowSelect %W %y}

    bind .fs_file.name_box.name <Return>  {fs_NameGetSelect %W}
    bind .fs_file.dir_box.directory <Return>  {fs_DirectoryGetSelect %W}

    bind .fs_file.pat_box.pattern <Return>  {fs_PatternGetSelect %W}

    bind .fs_file.name_box.name <Up> {CommandUpKey %W .fs_file.list_box.list}
    bind .fs_file.name_box.name <Down> {CommandDownKey %W .fs_file.list_box.list}
    bind .fs_file.name_box.name <Control-p> {CommandUpKey %W .fs_file.list_box.list}
    bind .fs_file.name_box.name <Control-n> {CommandDownKey %W .fs_file.list_box.list}

    .fs_file.dir_box.directory insert end $fs_file_dir
    .fs_file.pat_box.pattern insert end $fs_file_pattern

    .fs_file.list_box.list configure -setgrid 1

    fs_GetFileList $directory $pattern

    bind .fs_file <Destroy> {set fs_file_ok 0}

    set oldfocus [focus] 
    set fs_file_ok 0
    focus .fs_file.name_box.name

    tkwait variable fs_file_ok

    set retval ""
    if { $fs_file_ok != 0 } { set retval $fs_file_dir/$fs_file_select }

    if {[winfo exists .fs_file]} {destroy .fs_file}

    focus $oldfocus

    return $retval
}

#
# get list of files matching given pattern and directory
#
# 11 Sep 95 
#
proc fs_GetFileList {directory pattern} {
    .fs_file.list_box.list delete 0 end

    set a  [lsort [glob -nocomplain {.*} $directory/$pattern]]

    do i 0 [llength $a] {
         set b [lindex $a $i]
         if { $b != "" } {
             if { [file type $b] == "directory" } {
                 .fs_file.list_box.list insert end [list "D" [file tail $b]]
             } else {
                 .fs_file.list_box.list insert end [file tail $b]
             }
         }
    }

    .fs_file.list_box.list select clear 0 end
    .fs_file.list_box.list select set 0
}

#
# get the selection from listbox and update view and/or give singal
# for completition. called by signal.
#
# 11 Sep 95 
#
proc fs_ShowAndExecuteSelect {w} {
    global fs_file_select fs_file_pattern fs_file_dir  fs_file_ok fs_file_read

    set fs_file_select [$w get [$w curselect]]

    set a [split $fs_file_select " "]
    if { [llength $a] > 1 } {set fs_file_select [lindex $a 1]}

    if { $fs_file_read == "read" && ![file exists $fs_file_dir/$fs_file_select] } {
        dl_dialog {FileError} "File doesn't exist: $fs_file_select" {} 0 {OK}
        return
    }

    if { [file type $fs_file_dir/$fs_file_select] == "directory" } {

        set fs_file_dir $fs_file_dir/$fs_file_select
        set temp [pwd]
        cd $fs_file_dir
        set fs_file_dir [pwd]
        cd $temp

        .fs_file.dir_box.directory delete 0 end
        .fs_file.dir_box.directory insert end $fs_file_dir

        $w delete 0 end

        fs_GetFileList $fs_file_dir $fs_file_pattern

        .fs_file.name_box.name delete 0 end
    } else { set fs_file_ok 1 }
}

#
# get directory name from entry. called by  signal.
#
# 11 Sep 95 
#
proc fs_DirectoryGetSelect {w} {
    global fs_file_select fs_file_pattern fs_file_dir fs_file_ok
    global ELMER_POST_HOME env
    
    eval "set fs_file_dir [$w get]"

    if { ![file exists $fs_file_dir] } {
        dl_dialog {FileError} "Directory doesn't exist: $fs_file_dir" {} 0 {OK}
        return
    }

    fs_GetFileList $fs_file_dir $fs_file_pattern
}

#
# Get pattern string from entry. Called by  signal.
#
# 11 Sep 95 
#
proc fs_PatternGetSelect {w} {
    global fs_file_select fs_file_pattern fs_file_dir fs_file_ok
    
    set fs_file_pattern [$w get]

    fs_GetFileList $fs_file_dir $fs_file_pattern
}

#
# Get name string from entry. Called by  signal.
#
# 11 Sep 95 
#
proc fs_NameGetSelect {w} {
    global fs_file_select fs_file_dir fs_file_pattern fs_file_ok fs_file_read

    set fs_file_select [$w get]

    set a [split $fs_file_select " "]
    if { [llength $a] > 1 } {set fs_file_select [lindex $a 1]}

    if {$fs_file_select =="" || $fs_file_read == "read" && ![file exists $fs_file_dir/$fs_file_select]} {
        dl_dialog {FileError} "File doesn't exist: $fs_file_select" {} 0 {OK}
        return
    }


#    if { [catch [file type $fs_file_dir/$fs_file_select] == "directory"] } {
#        .fs_file.dir_box.directory insert end "/$fs_file_select"
#
#        set fs_file_dir $fs_file_dir/$fs_file_select
#
#        $w delete 0 end
#
#        fs_GetFileList $fs_file_dir $fs_file_pattern
#
#        .fs_file.name_box.name delete 0 end
#    } else { set fs_file_ok 1 }
     set fs_file_ok 1
}

#
# Get name string from listbox. Called by  signal.
#
# 11 Sep 95 
#
proc fs_ShowSelect {w y} {
    global fs_file_select fs_file_ok
 
    $w select clear 0 end
    $w select set  [$w nearest $y]

    set fs_file_select [$w get [$w curselect]]

    .fs_file.name_box.name delete 0 end
    .fs_file.name_box.name insert end $fs_file_select
}
