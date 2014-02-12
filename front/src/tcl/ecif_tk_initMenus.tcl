#/*****************************************************************************/
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

#***********************************************************************
#Program:   ELMER Front 
#Module:    ecif_tk_initMenus.tcl
#Language:  Tcl
#Date:      20.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Menus.
#
#************************************************************************


#NOTE:
# Menu widget names are stored in the global variable gMW.
# Individual button LABELS within menus are stored in gMT.
# Properties of the buttons can be set using these two variables like:
# $gMW(File) entryconfigure $gMT(File,modfile) -state disabled
# if model file import button in file-menu should be disabled (greyed)
#
 
# NOTE: This is a huge proc
################################
### Start of INIT_MENUS proc ###
################################
proc INIT_MENUS {} {

global globMenu Info Status Model ModelFlags gMW gMT gBW gBT UserSetting


# Set screen parameters
Screen::setParameters

# Init user settings data
# =======================
# NOTE This is done only once per session!
# NOTE: This is needed here because using UserSetting variables
# as command button variables will create these array fields and
# Panel::initFields does not any set the initial values because
# fields exists already (we cannot use force_initial option because
# default value may come from a user setting file !!!
#
Panel::initFields UserSetting $UserSetting(initialFields)

#======================
# MAIN WINDOW structure
#======================
set mw $Info(mainWindow)

set menuFrame [frame $mw.menuFrame]
set rowButtonsFrame [frame $mw.rowButtonsFrame -relief groove -bd 2]
set centerFrame [frame $mw.centerFrame]
set messageFrame [frame $centerFrame.messageFrame]
set colButtonsFrame [frame $centerFrame.colButtonsFrame -relief groove -bd 2]
set statusFrame [frame $mw.statusFrame]


#================================
# CONSTRUCTING Menu-bar and menus 
#================================
#
#NOTE:
# Menu widget names are stored in the global variable gMW.
# Individual button LABELS within menus are stored in gMT.
# Properties of the buttons can be set using these two variables like:
# $gMW(File) entryconfigure $gMT(File,modelfile) -state disabled
# if model file button in file-menu should be disabled (greyed)
#
#---Menu bar is created
set menuBar [append $menuFrame .menu]
menu $menuBar -relief groove -bd 2 -tearoff 0 ;#NOTE: -relief & -bd does not work!

# FILE-menu
#
set gMW(File) [MenuBuild::createMenuBarItem $menuBar file 0]
MenuBuild::addMenuBarItem $menuBar $gMW(File) "File" 0
MenuBuild::FileMenu


# EDIT-menu
#
set gMW(Edit) [MenuBuild::createMenuBarItem $menuBar edit 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Edit) "Edit" 0
MenuBuild::EditMenu


# DISPLAY-menu
#
if { $Info(ELMER_FRONT_THETIS_SUPPORT) } {
  set gMW(Display) [MenuBuild::createMenuBarItem $menuBar display 1]
} else {
  set gMW(Display) [MenuBuild::createMenuBarItem $menuBar display 0]
}
MenuBuild::addMenuBarItem $menuBar $gMW(Display) "Display" 0
MenuBuild::DisplayMenu


# PROBLEM-menu
#
set gMW(Problem) [MenuBuild::createMenuBarItem $menuBar problem 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Problem) "Problem" 0
MenuBuild::ProblemMenu


# MODEL-menu
#
set gMW(Model) [MenuBuild::createMenuBarItem $menuBar model 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Model) "Model" 0
MenuBuild::ModelMenu


# MESH-menu
#
set gMW(Mesh) [MenuBuild::createMenuBarItem $menuBar mesh 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Mesh) "Mesh" 0
MenuBuild::MeshMenu


# SOLVER-menu
#
set gMW(Solver) [MenuBuild::createMenuBarItem $menuBar solver 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Solver) "Solver" 0
MenuBuild::SolverMenu


# RUN-menu
#
set gMW(Run) [MenuBuild::createMenuBarItem $menuBar run 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Run) "Run" 0
MenuBuild::RunMenu


# WINDOW-menu
#
set gMW(Window) [MenuBuild::createMenuBarItem $menuBar window 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Window) "Window" 0
MenuBuild::WindowMenu


# HELP-menu
#
set gMW(Help) [MenuBuild::createMenuBarItem $menuBar help 0]
MenuBuild::addMenuBarItem $menuBar $gMW(Help) "Help" 0
MenuBuild::HelpMenu


#---Packing menu widgets
$Info(mainWindow) configure -menu $menuBar 

#===========================
# CONSTRUCTING Buttons Areas
#===========================

# ===========
# Row buttons
# ===========

set rbf $rowButtonsFrame
#-Rows (two rows currently)
set rbf1 [frame $rbf.r1]
set rbf2 [frame $rbf.r2]

#-Separator frames
frame $rbf1.sep1 -width 4
frame $rbf1.sep2 -width 6
frame $rbf2.sep1 -width 4
frame $rbf2.sep2 -width 6

set ipath $Info(imagePath)
set hg1 20

#---Move buttons
set globMenu(moveButtons) { moveLeft moveRight moveUp moveDown scaleUp scaleDown}
set fname move
set m [frame $rbf1.$fname]

set gBW(moveLeft)   [ button $m.moveLeft ]
set gBW(moveRight)  [ button $m.moveRight ]
set gBW(moveUp)     [ button $m.moveUp ]
set gBW(moveDown)   [ button $m.moveDown ]
set gBW(scaleUp)    [ button $m.scaleLeft ]
set gBW(scaleDown)  [ button $m.scaleDown ]
set gBT(moveLeft) [ image create photo moveLeft -file [file join $ipath left.gif] ]
set gBT(moveRight) [ image create photo moveRight -file [file join $ipath right.gif] ]
set gBT(moveUp) [ image create photo moveUp -file [file join $ipath up.gif] ]
set gBT(moveDown) [ image create photo moveDown -file [file join $ipath down.gif] ]
set gBT(scaleUp) [ image create photo scaleUp -file [file join $ipath scaleup.gif] ]
set gBT(scaleDown) [ image create photo scaleDown -file [file join $ipath scaledown.gif] ]

$gBW(moveLeft) configure -height $hg1 -image $gBT(moveLeft) \
      -command "MenuExec::rendererMove -1 X"
$gBW(moveRight) configure -height $hg1 -image $gBT(moveRight) \
      -command "MenuExec::rendererMove 1 X"
$gBW(moveUp) configure -height $hg1 -image $gBT(moveUp) \
      -command "MenuExec::rendererMove 1 Y"
$gBW(moveDown) configure -height $hg1 -image $gBT(moveDown) \
      -command "MenuExec::rendererMove -1 Y"
$gBW(scaleUp) configure -height $hg1 -image $gBT(scaleUp) \
      -command "MenuExec::rendererScale 1"
$gBW(scaleDown) configure -height $hg1 -image $gBT(scaleDown) \
      -command "MenuExec::rendererScale -1"

foreach btn $globMenu(moveButtons) { 
  pack $gBW($btn) -side left -padx 1
}


#---Rotate buttons
set globMenu(rotateButtons)   {rotateLeft rotateRight axisLabel axisX axisY axisZ }
set globMenu(rotateButtons2D) {rotateLeft rotateRight axisLabel axisZ }
set fname rotate
set m [frame $rbf2.$fname]

set gBW(rotateLeft)  [ button $m.rotateLeft ]
set gBW(rotateRight) [ button $m.rotateRight ]
set gBW(axisLabel)  [ label $m.axisLabel]
set gBW(axisX)       [ checkbutton $m.axisX -selectcolor $Info(select_color) ]
set gBW(axisY)       [ checkbutton $m.axisY -selectcolor $Info(select_color) ]
set gBW(axisZ)       [ checkbutton $m.axisZ -selectcolor $Info(select_color) ]
set gBT(rotateLeft)  [ image create photo rotateLeft -file [file join $ipath rotateLeft.gif] ]
set gBT(rotateRight) [ image create photo rotateRight -file [file join $ipath rotateRight.gif] ]
set gBT(axisX) "X"
set gBT(axisY) "Y"
set gBT(axisZ) "Z"

$gBW(rotateLeft) configure -height $hg1 -image $gBT(rotateLeft) \
      -command "MenuExec::rendererRotate -1"

$gBW(rotateRight) configure -height $hg1 -image $gBT(rotateRight) \
      -command "MenuExec::rendererRotate 1"

$gBW(axisLabel) configure -height 1 -text " Rotate:" 

$gBW(axisX) configure -height 1 -text $gBT(axisX) \
      -variable Info(currentRotateAxis) \
      -command "MenuExec::rendererSetRotatePriorities" \
      -onvalue "X" -offvalue ""

$gBW(axisY) configure -height 1 -text $gBT(axisY) \
      -variable Info(currentRotateAxis) \
      -command "MenuExec::rendererSetRotatePriorities" \
      -onvalue "Y" -offvalue ""

$gBW(axisZ) configure -height 1 -text $gBT(axisZ) \
      -variable Info(currentRotateAxis) \
      -command "MenuExec::rendererSetRotatePriorities" \
      -onvalue "Z" -offvalue ""

# Pack buttons
foreach btn $globMenu(rotateButtons) {
  pack $gBW($btn) -side left -padx 1
}


#---Draw target buttons 
#
set globMenu(drawTargetButtons) {draw_target_bodies draw_target_surfaces draw_target_edges }
image create photo bodies     -file [file join $ipath bodies.gif]
image create photo boundaries -file [file join $ipath boundaries.gif]
image create photo mesh       -file [file join $ipath mesh.gif]

# Texts for the buttons
set gBT(draw_target_bodies) "Draw\nbodies"
set gBT(draw_target_surfaces) "Draw\nsurfaces"
set gBT(draw_target_edges) "Draw\nedges"

# Frame for the buttons
set fname draw_target
set m [frame $rbf1.$fname]

# Button sizes
 if { $Info(platform) == "windows" } {
  set heights { 1 1 1 }
  set wid 6
} else {
  set heights { 2 2 2 }
  set wid 8
}

# Create and pack buttons
foreach btn $globMenu(drawTargetButtons) \
        hig $heights {
  #-Control variable values
  set ModelFlags($btn,oldValue)   0
  #-Create button
  set gBW($btn) [checkbutton $m.$btn -selectcolor $Info(select_color)]
  set group [string toupper $fname]
  set field [string toupper $btn]
  $gBW($btn) configure -height $hig -width $wid -text $gBT($btn) -indicatoron 0 \
        -command "Util::setFlagValue $group $field" \
        -variable ModelFlags($field) -font $Info(cmdBtnFont)
  #-Pack button
  pack $gBW($btn) -side left -padx 1
}


#---Draw source buttons
#
set globMenu(drawSourceButtons) { draw_source_cad draw_source_mesh }
image create photo mesh -file [file join $ipath mesh.gif]

# Texts for the buttons
set gBT(draw_source_cad) "Cad\ngeometry"
set gBT(draw_source_mesh) "Mesh\ngeometry"

# Button sizes
if {$Info(platform) == "windows" } {
  set heights { 1 1 }
  set wid 6
} else {
  set heights { 2 2 }
  set wid 8
}

set fname draw_source
set m [frame $rbf1.$fname]

foreach btn $globMenu(drawSourceButtons) \
        hig $heights {
  #-Control variable values
  set ModelFlags($btn,oldValue)   0
  #-Create button
  set gBW($btn) [checkbutton $m.$btn -selectcolor $Info(select_color)]
  set group [string toupper $fname]
  set field [string toupper $btn]
  $gBW($btn) configure -height $hig -width $wid -text $gBT($btn) -indicatoron 0 \
        -command "Util::setFlagValue $group $field" \
        -variable ModelFlags($field) -font $Info(cmdBtnFont)
  #-Pack button
  pack $gBW($btn) -side left -padx 1
}

#---User settings control buttons

set globMenu(settingsButtons) {auto_load_mesh}

# Texts for the buttons
set gBT(auto_load_mesh) "Auto load\nmesh"

# Frame for the buttons
set fname settings
set m [frame $rbf2.$fname]

# Button sizes
if { $Info(platform) == "windows" } {
  set heights { 1 }
  set wid 8
} else {
  set heights { 2 }
  set wid 8
}


# Create and pack buttons
foreach btn $globMenu(settingsButtons) \
        hig $heights {

  #-Create button
  set gBW($btn) [checkbutton $m.$btn -selectcolor $Info(select_color)]
  set group [string toupper $fname]
  set field [string toupper $btn]
  $gBW($btn) configure -height $hig -width $wid -text $gBT($btn) -indicatoron 0 \
                       -variable UserSetting($field) -font $Info(cmdBtnFont)

  #-Pack button
  pack $gBW($btn) -side left -padx 1
}


#---Display buttons
#
set globMenu(displayButtons) { selectBodies display reset }

# Texts and commands for the buttons
#
set texts {"Select\nbodies" "Display" "Reset" }
set commands {"Util::execProc BodyDisplay::openPanel" MenuExec::rendererDisplay MenuExec::rendererReset}

# Button sizes
if { $Info(platform) == "windows" } {
  set heights { 1 1 1 }
  set wid 6
} else {
  set heights { 2 2 2 }
  set wid 6
}

set fname display
set m [frame $rbf2.$fname]

foreach btn $globMenu(displayButtons) \
        txt $texts                    \
        cmd $commands                 \
        hig $heights {
  #-Create button
  set gBW($btn) [button $m.$btn]
  $gBW($btn) configure -height $hig -width $wid -text $txt -command $cmd  -font $Info(cmdBtnFont)

  #-Pack button
  pack $gBW($btn) -side left -padx 1 -anchor e
}

# Pack row button frames
pack $rbf1 -side top -anchor w -fill x -expand 0 -pady 2
pack $rbf2 -side top -anchor e -fill x -expand 0 -pady 2

pack $rbf1.move        \
     $rbf1.sep1        \
     $rbf1.draw_target \
     $rbf1.sep2        \
     $rbf1.draw_source -side left -fill x -expand 0
pack $rbf2.rotate      \
     $rbf2.sep1        \
     $rbf2.settings    \
     $rbf2.sep2        \
     $rbf2.display -side left -fill x -expand 0 -anchor e


# ==============
# Column buttons
# ==============

# NOTE One column currently
#
set cbf $colButtonsFrame

#-Separator frames
frame $cbf.sep1 -height 2
frame $cbf.sep2 -height 3
frame $cbf.sep3 -height 4

# These groups are used for widget state manipulation
#
set globMenu(commandButtons) { saveModel mesh loadMesh
                               solve results break info_}

set globMenu(browseButtons) { info_ clearMessageArea process }

# This group is used for widget creating
#
set globMenu(columnButtons) { saveModel mesh loadMesh
                              solve results break
                              clearMessageArea process info_
                            }

set texts    { "Save\nmodel" "Mesh" "Load\nmesh"
               "Solve" "Results" "Break"
               "Clear" "Process" "Info" }

set commands { "MenuExec::saveModelFile" "Util::execProc MeshDefine::openPanel" "MenuExec::loadMesh"
               "MenuExec::prepare_elmer_exec Solver"  "MenuExec::show_results"
               "MenuExec::doBreak"
               "Message::clearMessageArea" "Util::execProc ProcessTable::openPanel" "MenuExec::info_" }

# Column button sizes
if { $Info(platform) == "windows" } {
  set heights { 1 1 1 1 1 1 1 1 1}
  set wid 6
} else {
  set heights {2 1 2 1 1 1 1 1 1}
  set wid 6
}

set padys    {1 1 1 1 1 1 1 1 1}
set sep_wids {5 0 0 0 5 0 0 5 0}
set fname command
set m [frame $cbf.$fname]

foreach btn $globMenu(columnButtons) \
        txt $texts                   \
        cmd $commands                \
        hig $heights                 \
        pady $padys                  \
        swid $sep_wids {
  #-Create button
  set gBW($btn)  [ button $m.$btn ]
  $gBW($btn) configure -height $hig -width $wid -text $txt -command $cmd  -font $Info(cmdBtnFont)

  #-Separator frame
  frame $m.sep$btn -height $swid -width $wid

  #-Pack button
  pack $gBW($btn) -side top -pady $pady
  pack $m.sep$btn -side top 
}

# Pack command buttons frame
pack $cbf.command -side top -expand 1 -fill y


#=========================
# CONSTRUCTING Status Area
#=========================
set m $statusFrame

#Subframes
#
set frames {1 2 3}
set padys {0 4 0}

#Subframe fields
#
set id_sets { {1 2}
              {1 2 3}
              {1 2}
            }
set field_sets { {equations timesteps}
                 {materials initialConditions bodyForces}
                 {outerBoundaryConditions innerBoundaryConditions meshes }
               }
set text_sets  { {"Equations:" "Timesteps:"}
                 {"Materials:" "Initial cond:" "Body forces:"}
                 {"Outer bc:" "Inner bc:" "Meshes:" }
               }
set wid_sets   { {38 8}
                 {8 8 8}
                 {8 8 8}
               }

set padx_sets  { {0 0}
                 {0 21 0}
                 {0 21 0}
               }

set lwid_sets { {11 11}
                {11 13 11}
                {11 13 11}
              } 

set bg $Info(nonActiveBg)

#-Create sub-frames and their widgets and pack them
#
foreach frm $frames pady $padys \
        ids $id_sets fields $field_sets \
        texts $text_sets wids $wid_sets \
        padxs $padx_sets lwids $lwid_sets {

  # Form status area sub-frame
  set sub_frame [frame $m.$frm]
  pack $sub_frame -side top -fill x -anchor w -pady $pady

  # Form fields in the sub-frame
  foreach id $ids fld $fields txt $texts wid $wids padx $padxs lwid $lwids {
    set entry_var Status($fld)
    set frame_fld [frame $sub_frame.$id]
    set frame_lbl [frame $frame_fld.fl]
    set frame_ent [frame $frame_fld.fe]
    set fld_lbl [label $frame_lbl.l -text $txt -width $lwid]
    set fld_ent [entry $frame_ent.e -textvariable $entry_var \
                                    -width $wid -bg $bg -state disabled]
    set Status($fld,entryWidget) $fld_ent
    pack $frame_fld -side left -anchor w -padx $padx
    pack $frame_lbl $frame_ent -side left -anchor w
    pack $fld_lbl $fld_ent -anchor w
  }
}


#==========================
# CONSTRUCTING Message Area
#==========================
set m $messageFrame
set lb_frame [frame $m.lb_frame]

set bg  $Info(nonActiveBgLight)

#set lb [listbox $lb_frame.lb -selectmode browse -height 12 -width 45 
set lb [listbox $lb_frame.lb  -height 12 -width 45 -font $Info(msgAreaFont) \
        -xscrollcommand [list $lb_frame.sb_x set] \
        -yscrollcommand [list $lb_frame.sb_y set] \
        -bg $bg ]
set sb_x [scrollbar $lb_frame.sb_x -orient horizontal \
                                   -command [list $lb_frame.lb xview] ]
set sb_y [scrollbar $lb_frame.sb_y -orient vertical \
                                   -command [list $lb_frame.lb yview] ]

set Info(messageWindow,listBox) $lb

pack $lb_frame -expand 1 -fill both
pack $sb_x -side bottom -fill x -expand 0
pack $lb -side left -fill both -expand 1 
pack $sb_y -side left -fill y -expand 0 

  
#===============
# INITIALIZATION
#===============

#---Disable menubuttons which are either initially or completely not working.
MenuBuild::configureMenuButtons File Save 0
MenuBuild::configureMenuButtons File Browsers 0
MenuBuild::configureMenuButtons File LoadMesh 0
MenuBuild::configureMenuButtons File Mesh 0
MenuBuild::configureMenuButtons Edit NearlyAll 0
MenuBuild::configureMenuButtons Edit Editors 0
MenuBuild::configureMenuButtons Edit MatcGmtr 0
MenuBuild::configureMenuButtons Display All 0
#MenuBuild::configureMenuButtons Problem All 0
MenuBuild::configureMenuButtons Problem Model 0
#MenuBuild::configureMenuButtons Model All 0
MenuBuild::configureMenuButtons Model Model 0
MenuBuild::configureMenuButtons Mesh All 0
#MenuBuild::configureMenuButtons Solver All 0
MenuBuild::configureMenuButtons Solver Model 0
MenuBuild::configureMenuButtons Run NearlyAll 0
MenuBuild::configureButtons "all" 0
MenuBuild::configureButtons "break" 0

# Priority setting n Run menu available only in Win32!
if { $Info(platform) == "windows" } {
  MenuBuild::configureMenuButtons Run bgPriorityLevel 0
}


#===========================
# PACKING main window frames
#===========================
#
pack $menuFrame -side top -expand 0
pack $rowButtonsFrame -side top -anchor e -expand 0
pack $statusFrame -side bottom  -anchor e -padx 0  -pady 10 -expand 0
pack $centerFrame -side bottom  -pady 10 -fill both -expand 1
pack $colButtonsFrame -side right -anchor ne -pady 10 -expand 0 
pack $messageFrame -side right -padx 50  -pady 10 -anchor c -fill both -expand 1

#================
# WELCOME message
#================
Message::showMessage "Welcome to using $Info(FRONT_NAME)!" blue 1

}
#End of INIT_MENUS proc 
##############################
##############################
##############################




#############################
#### MenuBuild functions ####
#############################

# Generic
# =======

proc MenuBuild::applyUserLevel {} {
  global Info

  MenuBuild::FileMenu
  MenuBuild::EditMenu
}


proc MenuBuild::configurePanelState {menu panel state} {
  global gMW gMT
  $gMW($menu) entryconfigure $gMT($menu,$panel) -state $state
}


# Add one menu item to the menu bar
#
proc MenuBuild::addMenuBarItem {frame item menu_text acc_index} {
  global Info

  $frame add cascade -label $menu_text -menu $item -underline $acc_index
}


# Add a command button to a menu item
#
proc MenuBuild::addMenuCommand {menu name cmd acc_index {acc_cmd ""} } {
  global gMT gMW 

  set widget $gMW($menu)
  set state  $gMW($menu,$name,state)

  set label  $gMT($menu,$name)

  #-State given
  if { $state != "" } {
    $widget add command -label $label -command $cmd -state $state -underline $acc_index -accelerator $acc_cmd

  #-State not known currently
  } else {
    $widget add command -label $label -command $cmd -underline $acc_index -accelerator $acc_cmd
  }

  set gMW($menu,$name,exists) 1
}


# Add a cascade submenu to a menu item
#
proc MenuBuild::addMenuCascade {menu name acc_index {acc_cmd ""} } {
  global gMT gMW 

  set widget $gMW($menu)
  set cmenu  $gMW($menu,$name)
  set state  $gMW($menu,$name,state)

  set label  $gMT($menu,$name)
  
  #-State given
  if { $state != "" } {
    $widget add cascade -label $label -menu $cmenu -state $state -underline $acc_index -accelerator $acc_cmd

  #-State not known currently
  } else {
    $widget add cascade -label $label -menu $cmenu -underline $acc_index -accelerator $acc_cmd
  }

  set gMW($menu,$name,exists) 1
}


# Set menu command infos
#
proc MenuBuild::initMenuInfo {menu names} {
  global gMT gMW 

  foreach name $names {

    # Init state variable if not defined
    if { ![info exists gMW($menu,$name,state) ] } {
      set gMW($menu,$name,state) ""
    }

    # Mark not yet added
    set gMW($menu,$name,exists) 0
  }
}


#Construct one menu bar item
#
proc MenuBuild::createMenuBarItem {frame widget_name tear_off} {
  global Info
  
  set m $frame.$widget_name
  menu $m -tearoff $tear_off 
  
  return $m
}


# Delete one menu bar item
#
proc MenuBuild::deleteMenuBarItem {frame menu_text} {
  
  $frame delete $menu_text
}


# Insert one menu bar item after the item having menu_text
#
proc MenuBuild::insertMenuBarItem {frame item menu_text acc_index after_menu_text} {
  global Info

  $frame insert $after_menu_text cascade -label $menu_text -menu $item -underline $acc_index
}


#-Body colors list menu
#
proc MenuBuild::createBodyListMenu { menu_widget {command_proc ""} } {
  global Info Model ObjectTable

  set m $menu_widget
  $m delete 0 end
  set cmd $command_proc
  set bd_index ""

  set counter 0

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    set label $ObjectTable($id,nm)
    set color $ObjectTable($id,clr)
    
    if { $command_proc != "" } {
      $m add command -label $label -background $color -command "$command_proc $counter"
    } else {
      $m add command -label $label -background $color
    }
  }
}


#-Window list menu
#
proc MenuBuild::createWindowListMenu {menu_widget} {
  global Info

  set m $menu_widget
  $m delete 0 end

  # If no open windows
  if { $Info(windowList) == "" } {
    
    # Destroy possible mepty window list window    
    if { [info exists Info(windowListTearoff)] &&
         [winfo exists $Info(windowListTearoff)]
       } {
      destroy $Info(windowListTearoff)
    }

    return
  }

  # Update window list menu
  foreach row $Info(windowList) {

    set arr [lindex $row 0]
    set wid [lindex $row 1]
    set wnm [string trimlef [lindex $row 2] "*"]

    $m add command -label $wnm -command "Util::checkPanelWindow $arr $wid \"$wnm\"  \"\" 1"
  }
}


# Configure one option for one menu button
#
# NOTE: Menu button is identified by the menu text!!!
#
proc MenuBuild::configureMenuButtonOption { menu button {option state} {value normal} } {
  global gMW gMT

  set gMW($menu,$button,$option) $value
  $gMW($menu) entryconfigure $gMT($menu,$button) -$option $value
}


# Configure menu buttons state
#
proc MenuBuild::configureMenuButtons {menu buttons state} {
  global globMenu gMW gMT ModelFlags

  # This flag tells if we are done
  # already in the switch statement
  # ie. we have called some special
  # menu handling proc
  set done 0

  # Default is plain button name (ie. no group name)
  set btns $buttons

  # Check if buttons is for some menu group
  switch $menu {

    File { 
      switch $buttons {
        All { set btns { modelfile_save modelfile_saveAs copy_parameters} }
        Save { set btns { modelfile_save modelfile_saveAs copy_parameters} }
        Browsers { set btns {sif srf emf egf}; set menu "File,browsers" }
        LoadMesh { MenuBuild::configureFileLoadMeshMenus $state; set done 1 }
        Mesh { set btns { msavers save_elmer_meshfile save_elmerpost_meshfile save_thetis_meshfile } }
      }
    } 

    Edit {
     switch $buttons {
        NearlyAll { set btns { bodyProperties boundaries removeCadGeometry removeInactiveParameterData
                               editSolverInputFile editSolverReloadFile } }
        All { set btns { bodyProperties boundaries removeCadGeometry removeInactiveParameterData 
                         editSolverInputFile editSolverReloadFile userSettings } }
        Editors { set btns {sif srf egf}; set menu "Edit,editors" }
        MatcGmtr { set btns { updateCadGeometry dropMatcDefinitions } }
      }
    } 

    Display { 
      switch $buttons {
        All { set btns { draw } }
      }
    } 

    Problem { 
      switch $buttons {
        All { set btns { modelProperties coordinates timesteps datafiles constants equations equationOrder} }
        Model { set btns { modelProperties coordinates timesteps datafiles} }
      }
    } 

    Model { 
      switch $buttons {
        All { set btns { modelInfo bodyInfo bodyList bodyParam boundaryParam initial force mater bndr } }
        Model { set btns { modelInfo bodyInfo bodyList } }
      }
    } 

    Mesh { 
      switch $buttons {
        All { set btns { defineMesh selectMesh } }
      }
    } 

    Solver { 
      switch $buttons {
        All { set btns { parameters calculators processor } }
        Model { set btns { processors } }
      }
    } 

    Run { 
      switch $buttons {
        NearlyAll { set btns { {Mesh GEOMETRY_TYPE_CAD} Viewfactors GebhardtFactors
                                Solver resultFile processTable bgPriorityLevel} }
        All { set btns { {Mesh GEOMETRY_TYPE_CAD} Viewfactors GebhardtFactors
                          Solver resultFile processTable bgPriorityLevel Post } }
      }
    } 

  }

  if {$done} {
    return
  }

  # Handle all buttons given in the btns argument
  # check also is modelFlag forces the state into
  # disabled
  #
  foreach btn $btns {

    set bt [lindex $btn 0]
    set flag [lindex $btn 1]

    if {$state == 0} {
      set state_value disabled
    } else {
      set state_value normal
    }

    # If we have some flag dependency
    # which possibly force state to disabled
    #
    if {$flag != "" } {
      if { $ModelFlags($flag) == 0 } {
        set state_value disabled
      }
    }

    set gMW($menu,$bt,state) $state_value

    # Configure button widget (if it exists)
    #
    if { [info exists gMW($menu,$bt,exists)] &&
         $gMW($menu,$bt,exists)
       } {
      MenuBuild::configureMenuButtonOption $menu $bt state $state_value
    }
  }
}


# Menus (and button) for reading mesh geometry from the DB 
# Special handling is needed here!!!
# 
proc MenuBuild::configureFileLoadMeshMenus {mode} {
  global gMW gMT

  set load_state    disabled
  set unload_state  disabled

  # NOTE: You have to put the widget into the "normal" state
  # before you can configure it!
  # mode == 1 means "load": load --> active, unload --> inactive

  # Command button
  if {$mode == 1} {
    set load_state    normal
    set unload_state  disabled
    MenuBuild::configureButtonOption loadMesh state normal
    MenuBuild::configureButtonOption loadMesh text "Load\nmesh"
    MenuBuild::configureButtonOption loadMesh command "MenuExec::loadMesh"

  # mode == -1 means "unload":  load --> inactive, unload --> active
  } elseif {$mode == -1} {
    set load_state    disabled
    set unload_state  normal
    MenuBuild::configureButtonOption loadMesh state normal
    MenuBuild::configureButtonOption loadMesh text "Unload\nmesh"
    MenuBuild::configureButtonOption loadMesh command "MenuExec::unloadMesh"

  # mode == 0 means both inactive
  } elseif {$mode == 0 } {
    set load_state    disabled
    set unload_state  disabled
    MenuBuild::configureButtonOption loadMesh state disabled
  }

  # Menubuttons
  MenuBuild::configureMenuButtonOption File loadMesh state $load_state
  MenuBuild::configureMenuButtonOption File unloadMesh state $unload_state
}


# Configure buttons state
#
proc MenuBuild::configureButtons {button_group state} {
  global globMenu gBW gBT ModelFlags

  set groupNames { moveButtons rotateButtons
                   drawSourceButtons drawTargetButtons       
                   displayButtons commandButtons}

  # From all button group
  set allButtons ""
  foreach grp $groupNames {
    foreach btn $globMenu($grp) {
      lappend allButtons $btn
    }
  }
  set globMenu(allButtons) $allButtons

  switch $button_group {
    all                   { set buttons $globMenu(allButtons) }
    allButtons            { set buttons $globMenu(allButtons) }
    moveButtons           { set buttons $globMenu(moveButtons) }
    rotateButtons         { set buttons $globMenu(rotateButtons) }
    rotateButtons2D       { set buttons $globMenu(rotateButtons2D) }
    drawTargetButtons     { set buttons $globMenu(drawTargetButtons) }
    drawSourceButtons     { set buttons $globMenu(drawSourceButtons) }
    displayButtons        { set buttons $globMenu(displayButtons) }
    commandButtons        { set buttons $globMenu(commandButtons) }
    draw_target_bodies    { set buttons {draw_target_bodies} }
    draw_target_surfaces  { set buttons {draw_target_surfaces} }
    draw_target_edges     { set buttons {draw_target_edges} }
    draw_source_cad       { set buttons {draw_source_cad} }
    draw_source_mesh      { set buttons {draw_source_mesh} }
    saveModel             { set buttons {saveModel} }
    mesh                  { set buttons {mesh} }
    loadMesh              { set buttons {loadMesh} }
    settings              { set buttons {settings} }
    solve                 { set buttons {solve} }
    results               { set buttons {results} }
    info                  { set buttons {info} }
    default               { set buttons $button_group }
  }

  set state_value normal

  if {$state == 0} {
    set state_value disabled
  }

  # Loop all buttons in the group
  #
  foreach btn $buttons {
    
    # Check possible flag dependency
    #
    set bt [lindex $btn 0]
    set flag [lindex $btn 1]

    # State option is applicabale only
    # for button widgets
    set wclass [winfo class $gBW($bt)]
    if { ![regexp -nocase button $wclass] } {
      continue
    }

    if {$state == 0} {
      set state_value disabled
    } else {
      set state_value normal
    }

    # If we have some flag dependency
    # which possibly force state to disabled
    #
    if {$flag != "" } {
      if { $ModelFlags($flag) == 0 } {
        set state_value disabled
      }
    }

    $gBW($bt) configure -state $state_value
  }
}


# Configue one option for one button
# NOTE: Option must be applicable for the
# NOTE: Option WITHOUT minus (-) sign!
# button, this is not checked!!!
#
proc MenuBuild::configureButtonOption {button {option state} {value normal} } {
  global globMenu gBW gBT

  $gBW($button) configure -$option $value
}


# Toggle check-type menu-command label (add On/Off)
#
# NOTE: On/Off text is added after the 'base-text'
# NOTE: This is used mainly because check-mark is not well visible in Tcl-menus!!!
#
proc MenuBuild::setCheckMenuLabel { menu base_label indicator_var} {
  global gMW
  upvar #0 $indicator_var iv

  if { $iv } {
    set lbl "$base_label  (On)"
  } else {
    set lbl "$base_label  (Off)"
  }

  $gMW($menu) entryconfigure "$base_label*" -label $lbl
}


########################
#                      #
# MENU IMPLEMENTATIONS #
#                      #
########################

#===========#
# File menu #
#===========#

proc MenuBuild::FileMenu {} {
  global gMT gMW Info

  set dots $Info(menuDots)

  #--Cascade menu for file browser
  set gMW(File,browsers) $gMW(File).browsers
  set m $gMW(File,browsers)

  set gMT(File,browsers,sif) "Sif file"
  set gMT(File,browsers,srf) "Solver reload file"
  set gMT(File,browsers,emf) "Model file"
  set gMT(File,browsers,egf) "Egf file"
  set gMT(File,browsers,def) "Definitions"
  set gMT(File,browsers,any) "File$dots"

  if { ![winfo exists $m] } {

    set gMW(File,browsers,sif,exists) 1
    set gMW(File,browsers,srf,exists) 1
    set gMW(File,browsers,emf,exists) 1
    set gMW(File,browsers,egf,exists) 1
    set gMW(File,browsers,def,exists) 1
    set gMW(File,browsers,any,exists) 1

    menu $m  -tearoff 1
    $m add command -label $gMT(File,browsers,sif) -command "MenuExec::browseSifFile"
    $m add command -label $gMT(File,browsers,srf) -command "MenuExec::browseSolverReloadFile"
    $m add command -label $gMT(File,browsers,emf) -command "MenuExec::browseEmfFile"
    $m add command -label $gMT(File,browsers,egf) -command "MenuExec::browseEgfFile"
    $m add command -label $gMT(File,browsers,def) -command "MenuExec::browseCurrentDefinitions"
    $m add command -label $gMT(File,browsers,any) -command "MenuExec::browseFile"
  }

  #--Cascade menu for mesh file savers
  set gMW(File,msavers) $gMW(File).msavers
  set m $gMW(File,msavers)

  if { ![winfo exists $m] } {
    menu $m  -tearoff 0
    set as_external 1
    $m add command -label "Save mesh in Elmer format$dots"    -command "MenuExec::saveElmerMeshFile $as_external"
    $m add command -label "Save mesh in ElmerPost format$dots"  -command "MenuExec::saveElmerPostMeshFile"
    if { $Info(ELMER_FRONT_THETIS_SUPPORT) } {
      $m add command -label "Save mesh in Thetis format$dots" -command "MenuExec::saveThetisMeshFile"
    }
  }

  # All menu command names
  #defs_browse modelfile_browse solver_input_browse
  #save_elmer_meshfile save_elmerpost_meshfile save_thetis_meshfile
  set names { cadfile meshfile modelfile copy_parameters
              loadMesh unloadMesh
              modelfile_save modelfile_saveAs
              browsers
              deffile_load deffile_save
              msavers
              exit }

  set gMT(File,cadfile) "Open Cad file$dots"
  set gMT(File,meshfile) "Open mesh file$dots"
  set gMT(File,modelfile) "Open model file$dots"
  set gMT(File,copy_parameters) "Copy parameters$dots"
  set gMT(File,loadMesh) "Load mesh"
  set gMT(File,unloadMesh) "Unload mesh"
  set gMT(File,modelfile_save) "Save model file"
  set gMT(File,modelfile_saveAs) "Save model file As$dots"
  set gMT(File,browsers) "Browse$dots"
  set gMT(File,deffile_load) "Load definition file$dots"
  set gMT(File,deffile_save) "Save definition file As$dots"
  set gMT(File,msavers) "Save mesh$dots"
  set gMT(File,exit) "Exit"

  MenuBuild::initMenuInfo File $names
  $gMW(File) delete 0 end
  set m $gMW(File)

  MenuBuild::addMenuCommand File cadfile "MenuExec::openCadFile" 5 
  $m add separator
  MenuBuild::addMenuCommand File meshfile "MenuExec::openMeshFile" 8
  $m add separator
  MenuBuild::addMenuCommand File modelfile "MenuExec::openModelFile" 5
  $m add separator
  MenuBuild::addMenuCommand File copy_parameters "MenuExec::copyParameters"  5
  $m add separator
  MenuBuild::addMenuCommand File loadMesh "MenuExec::loadMesh" 0 Ctrl+L
  MenuBuild::addMenuCommand File unloadMesh "MenuExec::unloadMesh" 0 
  $m add separator
  MenuBuild::addMenuCommand File modelfile_save "MenuExec::saveModelFile" 0 Ctrl+S
  MenuBuild::addMenuCommand File modelfile_saveAs "MenuExec::saveModelFileAs"  16

  if { $Info(userLevel) >= $Info(advancedUser) } {
    $m add separator
    MenuBuild::addMenuCommand File deffile_load "MenuExec::loadDefinitionFile" 1
    MenuBuild::addMenuCommand File deffile_save "MenuExec::saveDefinitionFileAs" 5
  }

  $m add separator
  MenuBuild::addMenuCascade File browsers  0
  $m add separator
  MenuBuild::addMenuCascade File msavers  0

  if { 0 && $Info(ELMER_FRONT_THETIS_SUPPORT) } {
    $m add separator
    MenuBuild::addMenuCommand File save_thetis_meshfile "MenuExec::saveThetisMeshFile" 13
  }

  $m add separator
  MenuBuild::addMenuCommand File exit "MenuExec::cifExit" 1 Alt+X

  # File menu accelerators:
  set mw $Info(mainWindow)
  set m $gMW(File)
  bind $mw <Control-s> "$m invoke \"$gMT(File,modelfile_save)\""
  bind $mw <Control-S> "$m invoke \"$gMT(File,modelfile_save)\""
  bind $mw <Control-l> "$m invoke \"$gMT(File,loadMesh)\""
  bind $mw <Control-L> "$m invoke \"$gMT(File,loadMesh)\""
  bind $mw <Alt-x> "$m invoke \"$gMT(File,exit)\""
  bind $mw <Alt-X> "$m invoke \"$gMT(File,exit)\""
}


#===========#
# Edit menu #
#===========#

proc MenuBuild::EditMenu {} {
  global gMT gMW Info
  global BodyProperty Boundary UserSetting

  set dots $Info(menuDots)

  #--Cascade menu for file editors
  set gMW(Edit,editors) $gMW(Edit).editors
  set m $gMW(Edit,editors)

  set gMT(Edit,editors,sif) "Sif file"
  set gMT(Edit,editors,srf) "Solver reload file"
  set gMT(Edit,editors,egf) "Egf file"
  #set gMT(Edit,editors,def) "Definitions"
  set gMT(Edit,editors,any) "File$dots"

  if { ![winfo exists $m] } {

    set gMW(Edit,editors,sif,exists) 1
    set gMW(Edit,editors,srf,exists) 1
    set gMW(Edit,editors,egf,exists) 1
    #set gMW(Edit,editors,def,exists) 1
    set gMW(Edit,editors,any,exists) 1

    menu $m  -tearoff 1
    $m add command -label $gMT(Edit,editors,sif) -command "MenuExec::browseSifFile 1"
    $m add command -label $gMT(Edit,editors,srf) -command "MenuExec::browseSolverReloadFile 1"
    $m add command -label $gMT(Edit,editors,egf) -command "MenuExec::browseEgfFile 1"
    #$m add command -label $gMT(Edit,editors,def) -command "MenuExec::browseCurrentDefinitions 1"
    $m add command -label $gMT(Edit,editors,any) -command "MenuExec::browseFile 1"
  }

  #--A cascade for Matc input file names
  set gMW(Edit,setMatcInputFileNames) $gMW(Edit).setMatcInputFileNames
  set m $gMW(Edit,setMatcInputFileNames)

  if { ![winfo exists $m] } {
    menu $m  -tearoff 0
    $m add command -label "For model file$dots" -command "MenuExec::setMatcInputFile emf" 
    $m add command -label "For sif-file$dots" -command "MenuExec::setMatcInputFile sif" 
  }

  #--A cascade for user level setting
  set gMW(Edit,userLevel) $gMW(Edit).userLevel
  set m $gMW(Edit,userLevel)

  if { ![winfo exists $m] } {
    menu $m  -tearoff 1
    $m add checkbutton -label "Novice" -variable Info(userLevel,new) -onvalue $Info(noviceUser) \
                       -command MenuExec::checkUserLevel
    $m add checkbutton -label "Advanced" -variable Info(userLevel,new) -onvalue $Info(advancedUser) \
                       -command MenuExec::checkUserLevel
    $m add checkbutton -label "Power user" -variable Info(userLevel,new) -onvalue $Info(powerUser)  \
                       -command MenuExec::checkUserLevel
  }
    
  # All menu command names
  set names { bodyProperties boundaries
              removeCadGeometry removeInactiveParameterData
              setMeshInputUnit
              editors
              matcDefinitions setMatcInputFileNames 
              updateCadGeometry
              workingDirectory
              userSettings userLevel clearDefinitions}

  #--Menu texts
  set gMT(Edit,bodyProperties) $BodyProperty(winTitle)$dots
  set gMT(Edit,boundaries) $Boundary(winTitle)$dots
  set gMT(Edit,removeCadGeometry) "Remove Cad geometry$dots"
  set gMT(Edit,removeInactiveParameterData) "Remove inactive paramter data$dots"
  set gMT(Edit,setMeshInputUnit) "Set Mesh Input Unit$dots"
  set gMT(Edit,editors) "Edit$dots"
  set gMT(Edit,matcDefinitions) "MATC definitions$dots"
  set gMT(Edit,setMatcInputFileNames) "MATC input file names$dots"
  set gMT(Edit,updateCadGeometry) "Update Cad geometry$dots"
  set gMT(Edit,workingDirectory) "Working directory$dots"
  set gMT(Edit,userLevel) "User level$dots"
  set gMT(Edit,userSettings) $UserSetting(winTitle)$dots
  set gMT(Edit,clearDefinitions) "Clear definitions$dots"

  #--Pick current states and destroy possible old menu
  MenuBuild::initMenuInfo Edit $names
  $gMW(Edit) delete 0 end

  set gMT(Edit,dropMatcDefinitions) "Drop Matc definitions from model*"
  set gMW(Edit,dropMatcDefinitions,exists) 1

  set m $gMW(Edit)
  MenuBuild::addMenuCommand Edit bodyProperties "Util::execProc BodyProperty::openPanel" 0
  MenuBuild::addMenuCommand Edit boundaries "Util::execProc Boundary::openPanel" 3
  $m add separator
  if { $Info(userLevel) >= $Info(advancedUser) } {
    MenuBuild::addMenuCommand Edit removeCadGeometry "MenuExec::removeCadGeometry" 0
  }
  MenuBuild::addMenuCommand Edit removeInactiveParameterData "MenuExec::removeInactiveParameterData" 7
  $m add separator

  MenuBuild::addMenuCommand Edit setMeshInputUnit "MenuExec::setMeshInputUnit" 15
  $m add separator

  MenuBuild::addMenuCascade Edit editors 0
  $m add separator
  MenuBuild::addMenuCommand Edit matcDefinitions "Util::execProc MatcDefinitions::openPanel" 5
  $m add separator
  MenuBuild::addMenuCascade Edit setMatcInputFileNames 4
  $m add separator
  MenuBuild::addMenuCommand Edit updateCadGeometry "MenuExec::verifyUpdateCadGeometry" 0
  $m add separator
  $m add checkbutton -label "Drop Matc definitions from model  (Off)" \
                     -command "MenuExec::verifyDropMatcDefs" \
                     -variable Info(dropModelMatcDefinitions) -onvalue 1
  $m add separator
  MenuBuild::addMenuCommand Edit workingDirectory "MenuExec::setWorkingDirectory" 0
  $m add separator
  MenuBuild::addMenuCascade Edit userLevel 0        
  $m add separator
  MenuBuild::addMenuCommand Edit userSettings "Util::execProc UserSetting::openPanel" 2
  $m add separator
  MenuBuild::addMenuCommand Edit clearDefinitions "MenuExec::clearCurrentDefinitions" 0
}


#==============#
# Display menu #
#==============#

proc MenuBuild::DisplayMenu {} {
  global gMT gMW Info
  global BodyDisplay BoundaryDisplay VertexDisplay LabelDisplay

  set dots $Info(menuDots)

  set names { draw transform reset labels bodies boundaries vertices }

  set gMT(Display,draw) "Display model"
  set gMT(Display,transform) "Transform$dots"
  set gMT(Display,reset)  "Reset"
  set gMT(Display,labels) $LabelDisplay(winTitle)$dots
  set gMT(Display,bodies) $BodyDisplay(winTitle)$dots
  set gMT(Display,boundaries) $BoundaryDisplay(winTitle)$dots
  set gMT(Display,vertices) $VertexDisplay(winTitle)$dots
  #

  MenuBuild::initMenuInfo Display $names
  $gMW(Display) delete 0 end

  set m $gMW(Display)

  MenuBuild::addMenuCommand Display draw "MenuExec::rendererDisplay" 0
  $m add separator

  set gMW(Display,transform) $gMW(Display).transform
  set mt $gMW(Display,transform)

  if { ![winfo exists $mt] } {
    menu $mt  -tearoff 1
    $mt add command -label "Rot X+" -command "MenuExec::rendererRotate 1 X"
    $mt add command -label "Rot X-" -command "MenuExec::rendererRotate -1 X"
    $mt add command -label "Rot Y+" -command "MenuExec::rendererRotate 1 Y"
    $mt add command -label "Rot Y-" -command "MenuExec::rendererRotate -1 Y"
    $mt add command -label "Rot Z+" -command "MenuExec::rendererRotate 1 Z"
    $mt add command -label "Rot Z-" -command "MenuExec::rendererRotate -1 Z"
    $mt add separator
    $mt add command -label "Mov X+" -command "MenuExec::rendererMove 1 X"
    $mt add command -label "Mov X-" -command "MenuExec::rendererMove -1 X"
    $mt add command -label "Mov Y+" -command "MenuExec::rendererMove 1 Y"
    $mt add command -label "Mov Y-" -command "MenuExec::rendererMove -1 Y"
    $mt add separator
    $mt add command -label "Zoom+" -command "MenuExec::rendererScale 1"
    $mt add command -label "Zoom-" -command "MenuExec::rendererScale -1"
  }

  MenuBuild::addMenuCascade Display transform 0        

  $m add separator
  MenuBuild::addMenuCommand Display reset "MenuExec::rendererReset" 0 Ctrl+R
  $m add separator
  MenuBuild::addMenuCommand Display labels "Util::execProc LabelDisplay::openPanel" 0 Ctrl+L
  $m add separator
  MenuBuild::addMenuCommand Display bodies "Util::execProc BodyDisplay::openPanel" 8 Ctrl+B
  $m add separator
  MenuBuild::addMenuCommand Display boundaries "Util::execProc BoundaryDisplay::openPanel" 11 Ctrl+X
  $m add separator
  MenuBuild::addMenuCommand Display vertices "Util::execProc VertexDisplay::openPanel" 11 Ctrl+X

  # Display menu accelerators:
  set mw $Info(mainWindow)
  set m $gMW(Display)
  bind $mw <Control-r> "$m invoke \"$gMT(Display,reset)\""
  bind $mw <Control-R> "$m invoke \"$gMT(Display,reset)\""
  bind $mw <Control-l> "$m invoke \"$gMT(Display,labels)\""
  bind $mw <Control-L> "$m invoke \"$gMT(Display,labels)\""
  bind $mw <Control-b> "$m invoke \"$gMT(Display,bodies)\""
  bind $mw <Control-B> "$m invoke \"$gMT(Display,bodies)\""
  bind $mw <Control-h> "$m invoke \"$gMT(Display,boundaries)\""
  bind $mw <Control-H> "$m invoke \"$gMT(Display,boundaries)\""
  bind $mw <Control-x> "$m invoke \"$gMT(Display,vertices)\""
  bind $mw <Control-X> "$m invoke \"$gMT(Display,vertices)\""
}


#==============#
# Problem menu #
#==============#

proc MenuBuild::ProblemMenu {} {
  global gMT gMW Info
  global ModelParameter SimulationParameter
  global Constant Coordinate Datafile Equation Calculator
  global ModelProperty SifIncludeFile SolverOrder Timestep

  set dots $Info(menuDots)

  #--Command names
  set names { modelProperties datafiles modelParam simuParam coordinates timesteps constants
              equations calculators equationOrder}

  #--Menu
  set gMT(Problem,modelProperties) $ModelProperty(winTitle)$dots
  set gMT(Problem,datafiles) $Datafile(winTitle)$dots
  set gMT(Problem,modelParam) $ModelParameter(winTitle)$dots
  set gMT(Problem,simuParam) $SimulationParameter(winTitle)$dots
  set gMT(Problem,coordinates) $Coordinate(winTitle)$dots
  set gMT(Problem,timesteps) $Timestep(winTitle)$dots
  set gMT(Problem,constants) $Constant(winTitle)$dots
  set gMT(Problem,equations) $Equation(winTitle)$dots
  set gMT(Problem,calculators) $Calculator(winTitle)$dots
  set gMT(Problem,equationOrder) $SolverOrder(winTitle)$dots

  MenuBuild::initMenuInfo Problem $names
  $gMW(Problem) delete 0 end

  set m $gMW(Problem)
  MenuBuild::addMenuCommand Problem modelProperties "Util::execProc ModelProperty::openPanel" 0
  MenuBuild::addMenuCommand Problem datafiles "Util::execProc Datafile::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem modelParam "Util::execProc ModelParameter::openPanel" 6
  MenuBuild::addMenuCommand Problem simuParam "Util::execProc SimulationParameter::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem coordinates "Util::execProc Coordinate::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem timesteps "Util::execProc Timestep::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem constants "Util::execProc Constant::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem equations "StdPanelCreate::openPanel Equation" 0 Ctrl+E
  $m add separator
  MenuBuild::addMenuCommand Problem calculators "Util::execProc Calculator::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Problem equationOrder "Util::execProc SolverOrder::openPanel" 17

  # Problem menu accelerators:
  set mw $Info(mainWindow)
  set m $gMW(Problem)
  bind $mw <Control-e> "$m invoke \"$gMT(Problem,equations)\""
  #bind $mw <Control-E> "$m invoke \"$gMT(problem,equations)\""
}


#============#
# Model menu #
#============#

proc MenuBuild::ModelMenu {} {
  global gMT gMW Info
  global InputFileInfo ModelInfo BodyInfo 
  global BodyParameter BoundaryParameter
  global BodyForce InitialCondition Material BoundaryCondition

  set dots $Info(menuDots)

  #--A cascade menu for body-list
  set gMW(Model,bodyList) $gMW(Model).bodyList
  set m $gMW(Model,bodyList)

  if { ![winfo exists $m] } {
    menu $m
  }

  #--Command names
  set names { inputFileInfo modelInfo bodyInfo bodyList bodyParam boundaryParam initial force mater bndr }
  
  #--Menu
  set gMT(Model,inputFileInfo)  $InputFileInfo(winTitle)$dots
  set gMT(Model,modelInfo)  $ModelInfo(winTitle)$dots
  set gMT(Model,bodyInfo) $BodyInfo(winTitle)$dots
  set gMT(Model,bodyList) "Body list"
  set gMT(Model,bodyParam) $BodyParameter(winTitle)$dots
  set gMT(Model,boundaryParam) $BoundaryParameter(winTitle)$dots
  set gMT(Model,initial) $InitialCondition(winTitle)$dots
  set gMT(Model,force) $BodyForce(winTitle)$dots
  set gMT(Model,mater) $Material(winTitle)$dots
  set gMT(Model,bndr) $BoundaryCondition(winTitle)$dots

  MenuBuild::initMenuInfo Model $names
  $gMW(Model) delete 0 end

  set m $gMW(Model)
  MenuBuild::addMenuCommand Model inputFileInfo "InputFileInfo::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Model modelInfo "ModelInfo::openPanel" 0
  MenuBuild::addMenuCommand Model bodyInfo "BodyInfo::openPanel" 6
  $m add separator
  MenuBuild::addMenuCascade Model bodyList 5
  $m add separator
  MenuBuild::addMenuCommand Model bodyParam "StdPanelCreate::openPanel BodyParameter" 0
  MenuBuild::addMenuCommand Model boundaryParam "StdPanelCreate::openPanel BoundaryParameter" 3
  $m add separator
  MenuBuild::addMenuCommand Model initial "StdPanelCreate::openPanel InitialCondition" 0
  MenuBuild::addMenuCommand Model force "StdPanelCreate::openPanel BodyForce" 5
  MenuBuild::addMenuCommand Model mater "StdPanelCreate::openPanel Material" 0
  $m add separator
  MenuBuild::addMenuCommand Model bndr "StdPanelCreate::openPanel BoundaryCondition" 9
}


#===========#
# Mesh menu #
#===========#

proc MenuBuild::MeshMenu {} {
  global gMT gMW Info
  global MeshDefine

  set dots $Info(menuDots)

  #--Command names
  set names {defineMesh}
  
  #--Menu texts
  set gMT(Mesh,defineMesh) $MeshDefine(winTitle)$dots

  MenuBuild::initMenuInfo Mesh $names
  $gMW(Mesh) delete 0 end

  set m $gMW(Mesh)
  
  MenuBuild::addMenuCommand Mesh defineMesh "Util::execProc MeshDefine::openPanel" 0
}


#=============#
# Solver menu #
#=============#

proc MenuBuild::SolverMenu {} {
  global gMT gMW Info
  global Solver SolverControl Processor

  set dots $Info(menuDots)

  #--Command names
  set names { parameters control processor }

  #--Menu texts
  set gMT(Solver,parameters) $Solver(winTitle)$dots
  set gMT(Solver,control) $SolverControl(winTitle)$dots
  set gMT(Solver,processor) $Processor(winTitle)$dots

  MenuBuild::initMenuInfo Solver $names
  $gMW(Solver) delete 0 end

  set m $gMW(Solver)
  MenuBuild::addMenuCommand Solver parameters "Util::execProc Solver::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Solver control "Util::execProc SolverControl::openPanel" 7

  #
  # Processor panel not in use! MVe 02.09.01
  #
  #$m add separator
  #MenuBuild::addMenuCommand Solver processor "Util::execProc Processor::openPanel" 0
}


#==========#
# Run menu #
#==========#

proc MenuBuild::RunMenu {} {
  global gMT gMW Info
  global PostFileSelect ProcessTable

  set dots $Info(menuDots)

  #--A cascade for bg process priority level
  set gMW(Run,bgPriorityLevel) $gMW(Run).bgPriorityLevel
  set m $gMW(Run,bgPriorityLevel) 

  if { ![winfo exists $m] } {
    menu $m  -tearoff 1

    $m add checkbutton -label "Favor GUI" -variable Info(bgPriorityLevel) \
                       -onvalue $Info(ElmerRun,lowPriority)
    $m add checkbutton -label "Favor process" -variable Info(bgPriorityLevel) \
                       -onvalue $Info(ElmerRun,normalPriority)
  }

  #--Command names
  set names { Mesh Viewfactors GebhardtFactors Solver Results resultFile 
              processTable bgPriorityLevel }

  #--Menu texts
  set gMT(Run,Mesh)            "Generate mesh"
  set gMT(Run,Viewfactors)     "Calculate View factors"
  set gMT(Run,GebhardtFactors) "Calculate Gebhardt factors"
  set gMT(Run,Solver)          "Solver"
  set gMT(Run,Results)         "Postprocessor"
  set gMT(Run,resultFile)      $PostFileSelect(winTitle)$dots
  set gMT(Run,processTable)    $ProcessTable(winTitle)$dots
  set gMT(Run,bgPriorityLevel) "Process priority level$dots"

  MenuBuild::initMenuInfo Run $names
  $gMW(Run) delete 0 end

  set m $gMW(Run)
  MenuBuild::addMenuCommand Run Mesh "MenuExec::prepare_elmer_exec Mesh2D"  9 Ctrl+F1
  $m add separator
  MenuBuild::addMenuCommand Run Viewfactors "MenuExec::prepare_elmer_exec Viewfactors" 10 Ctrl+F2
  MenuBuild::addMenuCommand Run GebhardtFactors "MenuExec::prepare_elmer_exec GebhardtFactors" 10 Ctrl+F3
  $m add separator
  MenuBuild::addMenuCommand Run Solver "MenuExec::prepare_elmer_exec Solver" 0 Ctrl+F4
  $m add separator
  MenuBuild::addMenuCommand Run Results "MenuExec::prepare_elmer_exec Results" 0 Ctrl+F5
  MenuBuild::addMenuCommand Run resultFile "PostFileSelect::openPanel" 7 Ctrl+O
  $m add separator
  $m add checkbutton -label "Display run messages  (Off)" \
                     -command "MenuBuild::setCheckMenuLabel Run \"Display run messages\" Info(displayExecMessages)"  \
                     -variable Info(displayExecMessages) -onvalue 1
  $m add separator
  MenuBuild::addMenuCommand Run processTable "Util::execProc ProcessTable::openPanel" 8 Ctrl+T
  MenuBuild::addMenuCascade Run bgPriorityLevel 17 

  # run menu accelerators:
  set mw $Info(mainWindow)
  set m $gMW(Run)
  bind $mw <Control-F1> "$m invoke \"$gMT(Run,Mesh)\" "
  bind $mw <Control-F2> "$m invoke \"$gMT(Run,Viewfactors)\" "
  bind $mw <Control-F3> "$m invoke \"$gMT(Run,GebhardtFactors)\" "
  bind $mw <Control-F4> "$m invoke \"$gMT(Run,Solver)\" "
  bind $mw <Control-F5> "$m invoke \"$gMT(Run,Results)\" "
  bind $mw <Control-o>  "$m invoke \"$gMT(Run,resultFile)\" "
  bind $mw <Control-O>  "$m invoke \"$gMT(Run,resultFile)\" "
  bind $mw <Control-t>  "$m invoke \"$gMT(Run,processTable)\" "
  bind $mw <Control-T>  "$m invoke \"$gMT(Run,processTable)\" "
}


#=============#
# Window menu #
#=============#

proc MenuBuild::WindowMenu {} {
  global gMT gMW Info

  set dots $Info(menuDots)

  #--A cascade menu for window-list
  set gMW(Window,windowList) $gMW(Window).windowList
  set m $gMW(Window,windowList)
  
  if { ![winfo exists $m] } {
    menu $m
  }

  #--Command names
  set names { closeAll closeUnModified windowList }

  #--Menu texts
  set gMT(Window,closeAll) "Close all"
  set gMT(Window,closeUnModified) "Close unmodified"
  set gMT(Window,windowList) "Window list$dots"

  MenuBuild::initMenuInfo Window $names
  $gMW(Window) delete 0 end

  set m $gMW(Window)
  MenuBuild::addMenuCommand Window closeAll "MenuExec::closeAllWindows"  0
  MenuBuild::addMenuCommand Window closeUnModified "MenuExec::closeUnmodifiedWindows" 6
  $m add separator
  MenuBuild::addMenuCascade Window windowList 0

  set widget $gMW(Window,windowList)
  $widget configure -tearoffcommand "Util::setTearOffTitle windowListTearoff"

}


#===========#
# Help menu #
#===========#

proc MenuBuild::HelpMenu {} {
  global gMT gMW Info
  global SystemInfo

  set dots $Info(menuDots)

  #--Command names
  set names { index contents systemInfo graphicsInfo about
              sourceList tclCommand updateSource}

  #--Menu texts
  set gMT(Help,index) "Index"
  set gMT(Help,contents) "Contents"
  set gMT(Help,graphicsInfo) "Graphics info$dots"
  set gMT(Help,systemInfo) $SystemInfo(winTitle)$dots
  set gMT(Help,about) "About$dots"
  set gMT(Help,updateSource) "Update Tcl"
  set gMT(Help,sourceList) "Source list$dots"
  set gMT(Help,tclCommand) "Tcl command$dots"

  MenuBuild::initMenuInfo Help $names
  $gMW(Help) delete 0 end

  set m $gMW(Help)
  MenuBuild::addMenuCommand Help index "MenuExec::helpIndex" 0
  MenuBuild::addMenuCommand Help contents "MenuExec::helpContents" 0
  $m add separator
  MenuBuild::addMenuCommand Help graphicsInfo "RendererInfo::openPanel" 0
  MenuBuild::addMenuCommand Help systemInfo "SystemInfo::openPanel" 0
  $m add separator
  MenuBuild::addMenuCommand Help about "About::openPanel" 0

  #==================================
  # NOTE THIS IS FOR DEBUGGING TCL/TK
  #==================================
  if { $Info(ELMER_FRONT_DEBUG_TCL) == 1 } {
    $m add separator
    MenuBuild::addMenuCommand Help updateSource "updateSource"  -1 Ctrl+W
    MenuBuild::addMenuCommand Help sourceList "sourceListPanel" -1 Ctrl+M
    MenuBuild::addMenuCommand Help tclCommand "tclCommandPanel" -1 Ctrl+N

    # debug accelerators:
    set mw $Info(mainWindow)
    set m $gMW(Help)
    bind $mw <Control-m> "$m invoke \"$gMT(Help,sourceList)\""
    bind $mw <Control-M> "$m invoke \"$gMT(Help,sourceList)\""
    bind $mw <Control-n> "$m invoke \"$gMT(Help,tclCommand)\" "
    bind $mw <Control-N> "$m invoke \"$gMT(Help,tclCommand)\" "
    bind $mw <Control-w> "$m invoke \"$gMT(Help,updateSource)\" "
    bind $mw <Control-W> "$m invoke \"$gMT(Help,updateSource)\" "
  }
}


# ecif_tk_initMenus.tcl
# ********************
