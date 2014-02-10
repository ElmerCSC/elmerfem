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
#Module:    ecif_tk_meshDefinePanel.tcl
#Language:  Tcl
#Date:      31.05.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining mesh properties (name and mesh-H currently)
#
#************************************************************************
 
#------Mesh define proc------

proc MeshDefine::openPanel { {adding_new_mesh 0} } {
  ## This procedure displays the model level definition panel
  ## Global variables
  global Info Common MeshDefine Model ModelProperty ModelFlags

  set w $MeshDefine(winName)
  set wgeom $MeshDefine(winGeometry)

  set id [winfo atom $w]
  set MeshDefine(winId) $id

  set Info(thisWindow) $w

  if { !$ModelFlags(GEOMETRY_TYPE_CAD) } {
    set MeshDefine(winTitle)  "Define mesh  (No remeshing -only mesh geometry available!)"
  } else {
    set MeshDefine(winTitle) "Define mesh"
  }

  if { 1 == [Util::checkPanelWindow MeshDefine $id $MeshDefine(winTitle) $wgeom] } {
    return
  }  

  if { $ModelFlags(GEOMETRY_TYPE_CAD) } {
    set MeshDefine(acceptNewMesh) 1

  } else {
    set MeshDefine(acceptNewMesh) $adding_new_mesh
  }

  set MeshDefine(dataChanged) 0 
  set MeshDefine(dataModified) 0

  set MeshDefine(meshListChanged) 0
  set MeshDefine(addedMeshIndices) ""

  set MeshDefine(edited,GridH) 0
  set MeshDefine(edited,GridParameter) 0

  set MeshDefine(edit_mode) ""
  set MeshDefine(command_mode) ""
   
  toplevel $w
  set this $w
  focus $w 

  wm title $w $MeshDefine(winTitle)
  wm geometry $w $wgeom 

  # Copy Model-array lists to MeshDefine
  #
  MeshDefine::setMeshDefineLists

  MeshDefine::updateFields $Model(currentMeshIndex)

  set MeshDefine(MESH_INDEX) $Model(currentMeshIndex)

  Panel::initFields MeshDefine
  Panel::backupFields MeshDefine


  #-----WIDGET CREATION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  set bg $Info(nonActiveBg)

  #-Data is organized into three basic frames.
  #-Mesh H related data
  frame $w.f1 ;#---Mesh names lbox and new, rename, delete buttons frame
  frame $w.f1.lF  ;#--Mesh names lbox frame
  frame $w.f1.bF  ;#--Mesh names buttons outer frame
  frame $w.f1.bF1 ;#--Mesh names buttons frame1
  frame $w.f1.bF2 ;#--Mesh names buttons frame2

  frame $w.f2 ;#---Mesh data frame
  frame $w.f2.nF  ;#--Mesh name
  frame $w.f2.iF  ;#--Mesh info
  frame $w.f2.iF.fnn   ;#-Nof nodes
  frame $w.f2.iF.fne   ;#-Nof elements

  frame $w.f2.hF  ;#--Mesh-h related
  frame $w.f2.hF.fl    ;#-Mesh-h label frame
  frame $w.f2.hF.fH    ;#-Mesh-h entry frame
  frame $w.f2.hF.cB    ;#-Default, Current buttons

  frame $w.f2.fF  ;#--Mesh-factor related
  frame $w.f2.fF.fl    ;#-Mesh-f label frame
  frame $w.f2.fF.fF    ;#-Mesh-f entry frame
  frame $w.f2.fF.cB    ;#-Default, Current buttons

  frame $w.f2.fBG  ;#---Background mesh file related
  frame $w.f2.fBG.f1  ;#--Bg mesh file name realated
  frame $w.f2.fBG.f1.fl  ;#-Bg mesh label frame
  frame $w.f2.fBG.f1.fF  ;#-Bg mesh file entry frame
  frame $w.f2.fBG.f2  ;#--Bg mesh control buttons related
  frame $w.f2.fBG.f2.cB  ;#-Bg mesh control buttons frame

  frame $w.f2.cB  ;#--Mesh structure, density and generate buttons frame

  $w.f2 configure        -bd 2 -relief groove
  $w.f2.hF.cB configure  -bd 2 -relief groove
  $w.f2.fF.cB configure  -bd 2 -relief groove
  $w.f2.fBG configure    -bd 2 -relief groove
  $w.f2.cB configure     -bd 2 -relief groove

  #-Buttons (Ok etc.)
  frame $w.fB ;#--Buttons

  # If no cad-geometry, no remeshing!
  if { !$ModelFlags(GEOMETRY_TYPE_CAD) } {
    set def_state disabled
    set def_bg  $Info(nonActiveBg)

  # we have cad-geometry
  } else {
    set def_state normal
    set def_bg  white
  }

  # Mesh names frame
  # ================

  #--Create and pack names listbox frame
  #
  set bg $Info(nonActiveBg)
  
  frame $w.f1.lbF
  set m $w.f1.lbF

  set lb [listbox $m.lb  -height 8 -width 40 \
          -xscrollcommand [list $m.sb_x set]  \
          -yscrollcommand [list $m.sb_y set]  \
          -bg $bg ]
  set sb_x [scrollbar $m.sb_x -orient horizontal \
                              -command [list $m.lb xview] ]
  set sb_y [scrollbar $m.sb_y -orient vertical \
                              -command [list $m.lb yview] ]

  set MeshDefine(meshNamesLB) $lb

  MeshDefine::updateMeshNamesLB

  bind $MeshDefine(meshNamesLB) <ButtonRelease-1> "+MeshDefine::meshSelected"

  pack $sb_x -side bottom -fill x -expand 0
  pack $lb -side left -fill both -expand 1 
  pack $sb_y -side left -fill y -expand 0 

  #--Create and pack names buttons frame
  #
  set bw 6

  #-First column of buttons
  set m $w.f1.bF1
  button $m.new -text "New" -command "MeshDefine::meshNew" \
                -state $def_state -width $bw

  button $m.rename -text "Rename" -command "MeshDefine::meshRename" \
                   -state $def_state -width $bw

  button $m.delete -text "Delete" -command "MeshDefine::meshDelete" \
                   -state $def_state -width $bw

  pack $m.new $m.rename $m.delete -side top -padx $fpx1 -pady $fpy3 

  set MeshDefine(newButton)    $m.new
  set MeshDefine(renameButton) $m.rename
  set MeshDefine(deleteButton) $m.delete

  #-Second column of buttons
  set m $w.f1.bF2
  button $m.add -text "Add..." -command "MeshDefine::meshAdd" \
                   -state $def_state -width $bw

  button $m.remove -text "Remove" -command "MeshDefine::meshRemove" \
                   -state $def_state -width $bw

  button $m.display -text "Display" -command "MeshDefine::meshDisplay" \
                   -state $def_state -width $bw

  button $m.use -text "Use" -command "MeshDefine::meshUse" \
                   -state $def_state -width $bw

  pack $m.add  $m.remove $m.display $m.use -side top -padx $fpx1 -pady $fpy3 

  set MeshDefine(addButton) $m.add
  set MeshDefine(removeButton) $m.remove
  set MeshDefine(displayButton) $m.display
  set MeshDefine(useButton) $m.use


  # Pack mesh names subframes
  pack $w.f1.bF1 $w.f1.bF2 -in $w.f1.bF -side left -expand 0 -fill y -padx $fpx2 -pady $fpy3
  pack $w.f1.lbF -side left -expand 1 -fill both -padx $fpx1 -pady $fpy3
  pack $w.f1.bF -side left -expand 0 -fill y -padx $fpx1 -pady $fpy3


  # Create and pack mesh data frame
  # ===============================

  #--Create and pack mesh-name entry
  #
  set m $w.f2.nF

  #-Name label
  frame $m.lf
  label $m.lf.l -text "Mesh name:"
  pack $m.lf.l -side top -anchor w

  #-Name entry
  frame $m.ef

  #NOTE: This entry updates directly corresponding ModelProperty field!
  # 
  set wdg [entry $m.ef.e  -textvariable MeshDefine(MESH_NAME) -width 32 \
                 -font $Info(entryFont) -state disabled  -background $bg]

  set MeshDefine(nameEntry) $wdg
  set MeshDefine(allWidgets,MESH_NAME) $wdg

  bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 MeshDefine $wdg {%A %K}"
  Widget::setEntryBindings non_standard_single MeshDefine MESH_NAME $wdg
  bind  $wdg <KeyPress-Return> "MeshDefine::meshNameModified"

  #set action1  "$MeshDefine(deleteButton) configure -state disabled"
  #set action2  "$MeshDefine(renameButton) configure -state disabled"
  #bind $wdg <KeyRelease> "+Widget::bindEditingEvent \"$action1\" {%A %K}"
  #bind $wdg <KeyRelease> "+Widget::bindEditingEvent \"$action2\" {%A %K}"


  pack $m.ef.e -side left -anchor w -padx $fpx1 -pady $fpy3 -expand 0
  pack $m.lf $m.ef -side top -anchor w -padx $fpx1 -pady $fpy3 -expand 0

  #--Create and pack info frame (nof nodes/elements)
  #
  set m $w.f2.iF.fnn
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Current\nnof nodes: "
  frame $m.ef
  entry $m.ef.e  -textvariable MeshDefine(nofNodes) -width 12 -justify center \
                 -font $Info(entryFont) -state disabled -background $bg

  pack $m.lf $m.ef -side left
  pack $m.lf.l $m.ef.e -side left

  set m $w.f2.iF.fne
  frame $m.lf
  label $m.lf.l  -anchor nw -text "Current\nnof elements: "
  frame $m.ef
  entry $m.ef.e  -textvariable MeshDefine(nofElements) -width 12 -justify center \
                 -font $Info(entryFont) -state disabled -background $bg

  pack $m.lf $m.ef -side left
  pack $m.lf.l $m.ef.e -side left

  pack $w.f2.iF.fnn $w.f2.iF.fne -side left -padx $fpx1 -pady $fpy2


  #--Create and pack Mesh H entry
  #
  set m $w.f2.hF.fl
  set fld MESH_H

  #-Label
  label $m.l -text "Model mesh H \[m\]: "  -width 15
  pack $m.l -side left -anchor w
  set MeshDefine(allWidgets,label,$fld) $m.l

  set m $w.f2.hF.fH
  
  set wdg [entry $m.e  -textvariable MeshDefine($fld) -width 15 \
                 -font $Info(entryFont) -state $def_state -bg $def_bg]

  set MeshDefine(allWidgets,$fld) $wdg
  bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 MeshDefine $wdg {%A %K}"
  Widget::setEntryBindings non_standard_single MeshDefine $fld $wdg
  pack $m.e -side left -anchor w -padx $fpx1 -pady $fpy3 -expand 1

  #--Create and pack Mesh-H Buttons
  #
  set m $w.f2.hF.cB

  button $m.current -text Current -command "MeshDefine::setCurrentH" \
                       -state $def_state 

  button $m.default -text Default -command "MeshDefine::setDefaultH" \
                       -state $def_state  

  set MeshDefine(currentHButton) $m.current
  set MeshDefine(defaultHButton) $m.default

  pack $m.current $m.default -side left -padx $fpx1 -pady $fpy2

  # Pack mesh-h subframes
  pack $w.f2.hF.fl $w.f2.hF.fH $w.f2.hF.cB -side left -expand 0 -padx $fpx1 -pady $fpy3


  #--Create and pack Mesh F entry
  #
  set m $w.f2.fF.fl
  set fld MESH_F

  #-Label
  label $m.l -text "Mesh scaling factor:" -width 15
  pack $m.l -side left -anchor w
  set MeshDefine(allWidgets,label,$fld) $m.l

  set m $w.f2.fF.fF
  
  set wdg [entry $m.e  -textvariable MeshDefine($fld) -width 15 \
                 -font $Info(entryFont) -state $def_state -bg $def_bg]

  set MeshDefine(allWidgets,$fld) $wdg
  bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 MeshDefine $wdg {%A %K}"
  Widget::setEntryBindings non_standard_single MeshDefine $fld $wdg
  pack $m.e -side left -anchor w -padx $fpx1 -pady $fpy3 -expand 1

  #--Create and pack Mesh-F Buttons
  #
  set m $w.f2.fF.cB

  button $m.current -text Current -command "MeshDefine::setCurrentF" \
                       -state $def_state 

  button $m.default -text Default -command "MeshDefine::setDefaultF" \
                       -state $def_state  

  set MeshDefine(currentFButton) $m.current
  set MeshDefine(defaultFButton) $m.default

  pack $m.current $m.default -side left -padx $fpx1 -pady $fpy2

  # Pack mesh-f subframes
  pack $w.f2.fF.fl $w.f2.fF.fF $w.f2.fF.cB -side left -expand 0 -padx $fpx1 -pady $fpy3

  #--Create and pack background mesh file fields
  set fld MESH_BG_MESH_FILE  

  #
  #-Label
  set m $w.f2.fBG.f1.fl
  label $m.l -text "Bg mesh file"
  pack $m.l -side left -anchor w
  set MeshDefine(allWidgets,label,$fld) $m.l

  #-Bg file entry
  set m $w.f2.fBG.f1.fF
  
  frame $m.fFbb ;# Browse button frame

  set wdg [StdPanelCreate::createBrowsableFileField MeshDefine $fld $m $m.fFbb "" normal] 
  set MeshDefine(allWidgets,$fld) $wdg
  $wdg configure -width 32 -state $def_state

  pack $wdg -side left -anchor w


  #-Use file and Use as control file checkboxes
  set m $w.f2.fBG.f2.cB

  set fld MESH_BG_MESH_ACTIVE

  set wdg [checkbutton $m.bgActive -text "Use bg mesh file" -anchor w \
                 -variable MeshDefine($fld)            \
                 -command MeshDefine($fld,editProc)    \
                 -offvalue 0 -onvalue 1 -state $def_state                               \
                 -indicatoron 1  -selectcolor $Info(select_color) ]

  Widget::setCheckBoxBindings non_standard MeshDefine $fld $wdg
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 MeshDefine $fld $wdg"


  set fld MESH_BG_MESH_CONTROL

  set wdg [checkbutton $m.bgControl -text "Use as control file" -anchor w \
                 -variable MeshDefine($fld)      \
                 -offvalue 0 -onvalue 1 -state $def_state                          \
                 -indicatoron 1  -selectcolor $Info(select_color) ]

  Widget::setCheckBoxBindings non_standard MeshDefine $fld $wdg
  bind $wdg <ButtonPress-3> "Widget::genericButton-3 MeshDefine $fld $wdg"
  
  pack $m.bgActive $m.bgControl -side left -padx $fpx3

  # Pack subframes
  pack $w.f2.fBG.f1 $w.f2.fBG.f2 -side top
  pack $w.f2.fBG.f1.fl $w.f2.fBG.f1.fF -side left -expand 0 -padx $fpx1 -pady $fpy3
  pack $w.f2.fBG.f2.cB -side left


  #--Create and pack mesh structure buttons
  #
  set m $w.f2.cB

  button $m.structure -text "Mesh structure" \
                      -command "MeshDefine::callPanel GridParameter" \
                      -width 12 -state $def_state -justify center

  button $m.density -text "Mesh density" \
                    -command "MeshDefine::callPanel GridH" \
                    -width 15 -state $def_state 

  button $m.generate -text "Generate mesh" -command "MeshDefine::meshGenerate 1" \
                         -width 12 -state $def_state 

  set MeshDefine(structureButton) $m.structure
  set MeshDefine(densityButton)   $m.density
  set MeshDefine(generateButton)  $m.generate

  pack $m.structure  $m.density -side left -padx $fpx1 -pady $fpy2 
  pack $m.generate -side left -padx $fpx3 -pady $fpy2 

  # Pack subframes
  pack $w.f2.nF $w.f2.iF $w.f2.hF $w.f2.fF $w.f2.fBG -side top -anchor w -expand 0 -padx $fpx1 -pady $fpy1
  pack $w.f2.cB -side top -anchor c -expand 0 -padx $fpx1 -pady $fpy3


  # Create and pack buttons
  # =======================

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.fB.ok -text OK -command "MeshDefine::panelOk $this"]
  set cn_btn [button $w.fB.cancel -text Cancel -command "MeshDefine::panelCancel $this" \
                                  -state $ca]
  set ap_btn [button $w.fB.apply -text Apply -command MeshDefine::panelApply \
                                 -state $ap]
  
  focus $ok_btn
  set MeshDefine(applyButton)  $ap_btn
  set MeshDefine(cancelButton) $cn_btn
    
  pack $ok_btn $cn_btn $ap_btn -side left -expand 0 -padx $fpx2 -pady $fpy2


  # Pack all frames
  # ===============
  pack $w.f1 -side top -expand 1 -fill both -padx $fpx1 -pady $fpy3
  pack $w.f2 -side top -expand 0 -fill both -padx $fpx1 -pady $fpy3
  pack $w.fB -side top -expand 0 -fill y -padx $fpx1 -pady $fpy3


  # New mesh
  if { $MeshDefine(MESH_NAME) == "" } {
    MeshDefine::meshNew

  # Display info for the current mesh
  } else {
    MeshDefine::meshSelected
  }

  MeshDefine::setAllButtonsState normal

  if { $adding_new_mesh } {
    MeshDefine::meshNew
  }

  # Set field label bindings for right-button help
  Widget::setLabelBindings MeshDefine
}
# End MeshDefine::openPanel


#===========================================================#
#                   Button procs                            #
#===========================================================#

# Activate name-entry for entering a new mesh name
#
proc MeshDefine::meshNew {} {
  global Info MeshDefine ModelProperty

  set MeshDefine(edit_mode) "new"

  set MeshDefine(MESH_INDEX,old) $MeshDefine(MESH_INDEX)

  #set MeshDefine(MESH_INDEX) [llength $MeshDefine(meshNames)]
 
  Widget::configureS $MeshDefine(newButton) normal

  $MeshDefine(nameEntry) configure -state normal -bg white

  set MeshDefine(MESH_NAME) ""
  #set MeshDefine(MESH_INDEX) $Info(NO_INDEX)
  set MeshDefine(MESH_INDEX) [llength $MeshDefine(meshNames)]

  set MeshDefine(nofNodes) ""
  set MeshDefine(nofElements) ""

  MeshDefine::setDefaultH
  MeshDefine::setDefaultF

  focus $MeshDefine(nameEntry)
  $MeshDefine(nameEntry) icursor 0

  MeshDefine::setMeshEditButtonsState disabled
  MeshDefine::setCommandButtonsState normal
}


# Rename an existing mesh (ie. directory)
#
proc MeshDefine::meshRename {} {
  global Info MeshDefine

  if { $MeshDefine(edit_mode) == "rename" } {
    MeshDefine::meshNameModified
    set MeshDefine(edit_mode) ""
    return
  }

  set mn $MeshDefine(MESH_NAME)

  set msg [list "***WARNING: Mesh directory name will be changed\n" \
                "and model data will be updated accordingly!\n\n" \
                "Edit mesh name and press Enter to apply the changes\n\n" \
                $Info(continueOk) ]

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return 
  }


  set MeshDefine(edit_mode) "rename"

  MeshDefine::setAllButtonsState disabled
  Widget::configureS $MeshDefine(renameButton) normal

  $MeshDefine(nameEntry) configure -state normal -bg white
  $MeshDefine(nameEntry) icursor end

  focus $MeshDefine(nameEntry)
}


# Delete an existing mesh (directory and mesh files)
#
proc MeshDefine::meshDelete {} {
  global Info MeshDefine

  set mn $MeshDefine(MESH_NAME)

  if { [MeshDefine::checkSolverMeshUsage $mn] } {

    set msg [list "***WARNING: Mesh is used by solver(s)!\n\n" \
                  "Are you sure to remove the mesh!:\n$mn\n\n" \
                  $Info(continueOk) ]
  } else {
    set msg [list "Are you sure to delete the mesh!:\n$MeshDefine(MESH_NAME)\n\n" \
                   $Info(continueOk) ]
  }

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return 
  }

  set MeshDefine(edit_mode) "delete"

  MeshDefine::applyMeshDelete $mn
}


# Add new meshe(s) into model's mesh list (from MESHDIR)
#
proc MeshDefine::meshAdd {} {
  global Info MeshDefine ModelProperty

  if { $ModelProperty(MODEL_DIRECTORY,absolute) == "" } {
    return
  } 

  set mdir [file join $ModelProperty(MODEL_DIRECTORY,absolute) $Info(meshDirectoryName)] 

  # These meshes are in the mesh directory as directories
  set all_names [Util::getSubdirectoryList $mdir]
  set new_names ""
  set new_rows ""

  foreach an $all_names {

    if { -1 == [lsearch $MeshDefine(meshNames) $an] } {

      set elem_info [MeshDefine::readNofNodesAndElements $an]
      set nof_nodes [lindex $elem_info 0]
      set nof_elems [lindex $elem_info 1]
      
      # Mesh does not exists!
      #
      if { $nof_nodes == "" || $nof_elems == "" } {
        continue
      }

      set nm $an
      if { [string length $an] < 20 } {
        append nm [string repeat " " 20]
        set nm [string range $nm 0 20]
      }
      set nds [format "%-6d" $nof_nodes ]
      set elm $nof_elems
      set nr "$nm (Nd:$nds  El:$elm)"

      lappend new_names $an
      lappend new_rows $nr
    }
  }

  set w .selectMeshesW
  set wgeom  -100+100

  toplevel $w 
  focus $w

  #--Window properties
  wm title $w "Select meshes to be added"
  wm geometry $w $wgeom
  
  set flags [CheckBoxList::create $w $new_rows "" 50 15]

  # If nothing selected
  #
  if { -1 == [lsearch $flags 1] } {
    return
  }

  # Add selected meshes
  #
  set max_idx [llength $MeshDefine(meshNames)]
  
  foreach flag $flags nm $new_names {
    if { $flag } {
      lappend MeshDefine(meshNames) $nm

      lappend MeshDefine(meshBgMeshActives) 0
      lappend MeshDefine(meshBgMeshControls) 0
      lappend MeshDefine(meshBgMeshFiles) {}
      lappend MeshDefine(meshFs) -1.0
      lappend MeshDefine(meshHs) -1.0

      lappend MeshDefine(addedMeshIndices) $max_idx
      incr max_idx
    }   
  }

  MeshDefine::updateMeshNamesLB

  Panel::panelDataChanged 1 MeshDefine
  set MeshDefine(meshListChanged) 1
}


# Remove (detach) an existing mesh from model
#
# NOTE: Mesh direcotry or files are NOT deleted!
#
proc MeshDefine::meshRemove {} {
  global Info MeshDefine Solver

  set mn $MeshDefine(MESH_NAME)
  set set_solver_mesh 0

  if { [MeshDefine::checkSolverMeshUsage $mn] } {
 
    set msg [list "***WARNING: Mesh is used by solver(s)!\n\n" \
                  "Are you sure to remove the mesh!:\n$mn\n\n" \
                  $Info(continueOk) ]
 
    if { ![Message::verifyAction $Info(powerUser) $msg] } {
      return 
    }
    
    if { 2 == [llength $MeshDefine(meshNames)] } {

      set msg [list "***NOTE: Only one mesh will be left in the model!\n" \
                    "Should all solvers be set to use the remaining mesh?\n\n" \
                    $Info(continueOk) ]
 
      if { [Message::verifyAction $Info(powerUser) $msg] } {
        set set_solver_mesh 1
      }
    }

  }

  set MeshDefine(edit_mode) "delete"
  set remove_only 1
  MeshDefine::applyMeshDelete $mn $remove_only

  if { $set_solver_mesh } {
    MeshDefine::setSolverMesh [lindex $MeshDefine(meshNames) 0]
  }
}


# Display (load) the selected mesh
#
proc MeshDefine::meshDisplay {} {
  global MeshDefine

  MeshDefine::meshSelected

  set MeshDefine(edit_mode) "display"

  MeshDefine::panelSave

  MeshDefine::updateMeshNamesLB
}


# Use selected mesh in all solvers
#
proc MeshDefine::meshUse {} {
  global Info MeshDefine Solver

  set mn $MeshDefine(MESH_NAME)

  set msg [list "***NOTE: All solvers are set to use mesh:\n" \
                "$mn\n\n" \
                $Info(continueOk) ]

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return 
  }

  MeshDefine::setSolverMesh $mn
}


# Set data proc
#
proc MeshDefine::setCurrentH {} {
  global MeshDefine ModelProperty

  set MeshDefine(MESH_H) [lindex $MeshDefine(meshHs) $MeshDefine(MESH_INDEX)]

  Panel::panelDataChanged 0 MeshDefine
  set MeshDefine(MESH_H,mod) 0
  Widget::setFieldStatus MeshDefine MESH_H
}


# Set data proc
#
proc MeshDefine::setDefaultH {} {
  global MeshDefine ModelProperty

  set MeshDefine(MESH_H) $MeshDefine(MESH_H,default)

  Panel::panelDataChanged 1 MeshDefine
  set MeshDefine(MESH_H,mod) 1
  Widget::setFieldStatus MeshDefine MESH_H
}


# Set data proc
#
proc MeshDefine::setCurrentF {} {
  global MeshDefine ModelProperty

  set MeshDefine(MESH_F) [lindex $MeshDefine(meshFs) $MeshDefine(MESH_INDEX)]

  Panel::panelDataChanged 0 MeshDefine
  set MeshDefine(MESH_F,mod) 0
  Widget::setFieldStatus MeshDefine MESH_F
}


# Set data proc
#
proc MeshDefine::setDefaultF {} {
  global MeshDefine ModelProperty

  set MeshDefine(MESH_F) $MeshDefine(MESH_F,default)

  Panel::panelDataChanged 1 MeshDefine
  set MeshDefine(MESH_F,mod) 1
  Widget::setFieldStatus MeshDefine MESH_F
}



# Open MeshStructure/MeshDensity panels
#
proc MeshDefine::callPanel {globArray} {
  global MeshDefine
  upvar #0 $globArray theArray

  if { $globArray == "GridH" } {
    set MeshDefine(command_mode) "mesh_density"
  } else {
    set MeshDefine(command_mode) "mesh_structure"
  }

  #---Error in data
  if { ![MeshDefine::checkPanelData] } {
    return
  }
  
  # Set mesh index nad name, so that parameterId update proc
  # can check if there is any vali paramter for the current mesh
  #
  set theArray(MESH_INDEX) $MeshDefine(MESH_INDEX)
  set theArray(MESH_NAME) $MeshDefine(MESH_NAME)
  execNsProc $globArray updateObjectParameterId
  
  # NOTE: This seem to reset index nad name fields, so we have to
  # set them again after this call!
  #
  StdPanelCreate::openPanel $globArray

  set theArray(MESH_INDEX) $MeshDefine(MESH_INDEX)
  set theArray(MESH_NAME) $MeshDefine(MESH_NAME)

  set MeshDefine(command_mode) ""
}


# -------------------------------
# GENERATE NEW MESH (call Mesh2D)
# -------------------------------
#
# NOTE: Grid parameters are attached to the body layer only
# after this command is used!
#
# This is because we want to avoid confusion about what parameters
# are used for existing meshes
#
proc MeshDefine::meshGenerate { {call_mesh_generator 1} } {
  global Info MeshDefine Model ModelFlags ModelProperty ObjectTable

  #--If no active meshing bodies
  if { $call_mesh_generator } {

    set has_some_body 0
    foreach id $ObjectTable(ids) {

      if { "BL" != [Object::getType $id] } {
        continue
      }

      if { !$ObjectTable($id,excldMsh) } {
        set has_some_body 1
        break
      }
    }

    if { !$has_some_body } {
      set msg "None of the bodies is included in meshing!"
      set Info(messageIcon) error
      Message::dispOkMessage $msg "$Info(FRONT_NAME) error message" $MeshDefine(winName)
      return -1
    }
  }

  set MeshDefine(command_mode) "generate"

  if { ![MeshDefine::checkPanelData] } {
    return -1
  }

  #--Add possible new mesh name being entered
  #  to the list
  if { $MeshDefine(edit_mode) == "new" } {

    if { ![MeshDefine::applyMeshNew] } {
      return
    }
    
    # Set new mesh index
    set MeshDefine(MESH_INDEX) [expr [llength $MeshDefine(meshNames)] - 1]

  }
  
  #--If mesh definition panels were modified, apply modifications!
  #
  if { $MeshDefine(edited,GridParameter) } {
    GridParameter::updateObjectParameterIds
    StdPanelExec::panelSave GridParameter
  }

  if { $MeshDefine(edited,GridH) } {
    GridH::updateObjectParameterIds
    StdPanelExec::panelSave GridH
  }

  # Save data to cpp-side
  MeshDefine::panelSave

  #-Create mesh directory if it does not exist!
  #
  Util::checkAndCreateDirectory2B $ModelProperty(CURRENT_MESH_DIRECTORY,absolute) ""

  set Model(Mesh,needsUpdate) 0

  # Call ElmerMesh
  # ==============
  if { $call_mesh_generator } {

    # Set possible backround mesh file info
    #
    set Model(currentMeshBgMeshFile) ""
    set Model(currentMeshBgMeshControlFile) ""
    
    if { $MeshDefine(MESH_BG_MESH_ACTIVE) } {

      if { !$MeshDefine(MESH_BG_MESH_CONTROL) } {
        set Model(currentMeshBgMeshFile) $MeshDefine(MESH_BG_MESH_FILE)
      } else {
        set Model(currentMeshBgMeshControlFile) $MeshDefine(MESH_BG_MESH_FILE)
      }
    }

    # Elmer has a 2D generator  
    if { $Model(GEOMETRY_DIMENSION) == "2D" } {
      MenuExec::prepare_elmer_exec Mesh2D

    # Elmer does not have a 3D generator!
    } else {
      #MenuExec::prepare_elmer_exec Mesh3D
    }
  }

  #--Update mesh names lists (check that generated mesh is marked with *)
  MeshDefine::updateMeshNamesLB
  
  #--Update nof nodes and elements fields
  MeshDefine::updateNofNodesAndElements

  #--If no active meshes, make this the active mesh
  if { $Model(nofActiveMeshes) == 0 &&
       $MeshDefine(MESH_NAME) != ""
     } {
    MeshDefine::setSolverMesh $MeshDefine(MESH_NAME)
  }

  Solver::updateActiveSolvingOrderList

  $MeshDefine(displayButton) configure -state normal

  set MeshDefine(edit_mode) ""
  set MeshDefine(command_mode) ""
}


#=========================#
#   Button helper procs   #
#=========================#

proc MeshDefine::updateMeshNamesLB {} {
  global MeshDefine
  
  set mlist ""
  
  set idx 0
  foreach mn $MeshDefine(meshNames) {
    if { $idx == $MeshDefine(MESH_INDEX) } {
      lappend mlist "*$mn"
    } else {
      lappend mlist $mn
    }

    incr idx
  }

  ListBox::fill $MeshDefine(meshNamesLB) $mlist
  
  $MeshDefine(meshNamesLB) selection set $MeshDefine(MESH_INDEX)
  $MeshDefine(meshNamesLB) see $MeshDefine(MESH_INDEX)
}


# A mesh selected in the listbox
#
proc MeshDefine::meshSelected {} {
  global MeshDefine Info Model ModelProperty

  set index [$MeshDefine(meshNamesLB) curselection]

  if { $index == "" } {
    return
  }

  set rc 0

  if { $MeshDefine(edited,GridParameter) } {
    set msg [list "NOTE: Mesh structure parameters were modified, but the mesh was not\n" \
                  "regenerated. If you continue, modifications will be lost!\n\n"  \
                   "$Info(anywayOk)" ]

    if { ![Message::verifyAction $Info(advancedUser) $msg] } {
      set rc 1
    }

  } elseif { $MeshDefine(edited,GridH) } {
    set msg [list "NOTE: Mesh local density parameters were modified, but the mesh was not\n" \
                  "regenerated. If you continue, modifications will be lost!\n\n"  \
                   "$Info(anywayOk)" ]

    if { ![Message::verifyAction $Info(advancedUser) $msg] } {
      set rc 1
    }

  } elseif { $MeshDefine(dataChanged) } {
    set msg [list "NOTE: Mesh definition was changed, but the mesh was not\n" \
                  "(re)generated. If you continue, modifications will be lost!\n\n"  \
                   "$Info(anywayOk)" ]

    if { ![Message::verifyAction $Info(noviceUser) $msg] } {
      set rc 1
    }
  }

  if { $rc } {
    set idx $MeshDefine(MESH_INDEX,old)
    $MeshDefine(meshNamesLB) selection clear 0 end
    $MeshDefine(meshNamesLB) selection set $idx $idx
    return
  }

  # Reset panel
  # ===========
  $MeshDefine(nameEntry) configure -state disabled -bg $Info(nonActiveBg)
  set MeshDefine(command_mode) ""
  set MeshDefine(edit_mode) ""

  set MeshDefine(edited,GridParameter) 0
  set MeshDefine(edited,GridH) 0

  Panel::panelDataChanged 0 MeshDefine
  Panel::panelDataModified 0 MeshDefine

  # Update data
  # ===========
  # Store current index
  set MeshDefine(MESH_INDEX,old) $MeshDefine(MESH_INDEX)
  set MeshDefine(MESH_NAME,old) $MeshDefine(MESH_NAME)
  
  set MeshDefine(MESH_INDEX) $index

  MeshDefine::updateFields $index

  #set MeshDefine(MESH_NAME,old) $MeshDefine(MESH_NAME)

  #--Update nof nodes and elements fields
  MeshDefine::updateNofNodesAndElements

  MeshDefine::setAllButtonsState normal

  MeshDefine(MESH_BG_MESH_ACTIVE,fillProc)
}


# When Enter is pressed in the name field
#
proc MeshDefine::meshNameModified {} {
  global Info MeshDefine

  if { $MeshDefine(edit_mode) != "new" &&
       $MeshDefine(edit_mode) != "rename"
     } {
    return
  }

  set mn [string trim $MeshDefine(MESH_NAME)]

  # No mesh name entered!
  if { $mn == "" } {
    return
  }

  MeshDefine::setCommandButtonsState normal

  $MeshDefine(nameEntry) configure -state disabled -bg $Info(nonActiveBg)

  if { $mn == $MeshDefine(MESH_NAME,old) } {
    return
  }

  switch $MeshDefine(edit_mode) {

    new {
    }

    rename {
    
      MeshDefine::applyMeshRename 
      set MeshDefine(MESH_NAME,old) $MeshDefine(MESH_NAME)
      MeshDefine::setMeshEditButtonsState normal
    }
  }

}


proc MeshDefine::applyMeshNew {} {
  global Info MeshDefine

  set is_new_mesh 1

  set index 0
  #-Check if name exists
  foreach mn $MeshDefine(meshNames) {

    if { [string equal $mn $MeshDefine(MESH_NAME)] } {

      # If an existing name was given, when a new mesh
      # was being define
      #
      if { $MeshDefine(edit_mode) == "new" } {

        set msg [list "NOTE: mesh name already exist!\n\n" \
                      $Info(continueOk) ]
      
        if { ![Message::verifyAction $Info(powerUser) $msg] } {
          return 0
        }
      }
     
      set is_new_mesh 0
      set MeshDefine(MESH_INDEX) $index

      break
    }

    incr index
  }

  #-Add new mesh to the list
  if { $is_new_mesh } {
    MeshDefine::updateLists "add"
  }

  #-Update list box
  MeshDefine::updateMeshNamesLB

  $MeshDefine(nameEntry) configure -state disabled -bg $Info(nonActiveBg)

  return 1
}


# Delete mesh (delete also mesh files and directory) or just
# remove mesh from the model
#
# NOTE: When deleting or removing a mesh we have to update and save
# the model right away, otherwise mesh references in parameters (Solvers and
# GridStructute, GridH) could go into terrible mesh!!!
#
proc MeshDefine::applyMeshDelete { {mesh_name ""} {remove_only 0} } {
  global Info MeshDefine Model ModelProperty

  if { $mesh_name == "" } {
    set deleted_name $MeshDefine(MESH_NAME)
  } else {
    set deleted_name $mesh_name
  }

  set deleted_index [lsearch $MeshDefine(meshNames) $mesh_name]

  if { $deleted_index == -1 } {
    return
  }  

  # Delete mesh files if in delete mode!!!
  #
  if { !$remove_only } {
    set mesh_dir [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                           $Info(meshDirectoryName) \
                           $deleted_name ]

    if { [catch {file delete -force $mesh_dir} msg ] } {
  MSG "ERROR in mesh delete because: $msg"
      return
    }
  }
  # If current mesh was deleted, we must unload it!
  #
  if { $deleted_index == $Model(currentMeshIndex) } {
    MenuBuild::configureFileLoadMeshMenus 0
    MenuBuild::configureButtonOption loadMesh state disabled
    set unload_mesh 1
  } else {
    set unload_mesh 0
  }


  # Update panel data
  # =================
  MeshDefine::updateLists "delete" $deleted_index

  set MeshDefine(MESH_INDEX) $Info(NO_INDEX)
  set MeshDefine(MESH_NAME) ""

  MeshDefine::updateMeshNamesLB

  # None of the meshes in the listbox is currently selected, so these
  # buttons must be disabled!
  #
  MeshDefine::setMeshEditButtonsState disabled
  MeshDefine::setMeshUseButtonsState disabled

  set save_model 1

  # If some of the just added meshes is removed, we do not save the model!
  #
  set idx [lsearch $MeshDefine(addedMeshIndices) $deleted_index]
  if { $idx != -1  } {
    set MeshDefine(addedMeshIndices) \
      [lreplace $MeshDefine(addedMeshIndices) $idx $idx]
    set save_model 0
  }

  # Update model data!
  # ------------------
  #
  if { $save_model } {
    MeshDefine::panelSave
    MeshDefine::checkParameterMeshData $deleted_name $deleted_index

    # Save model file!!!
    #
    MenuExec::saveModelFile
  }

  if { $unload_mesh } {
    MenuExec::unloadMesh
  }
}


# Apply mesh rename when the user has prssed the renamed button
#
# NOTE: When renaming a mesh we have to update and save the model
# right away, otherwise mesh references in parameters (Solvers in this case)
# could go into terrible mesh!!!

#
proc MeshDefine::applyMeshRename {} {
  global Info MeshDefine Model ModelProperty Solver

  set index $MeshDefine(MESH_INDEX)
  set old_name [lindex $Model(meshNames) $index]
  set new_name $MeshDefine(MESH_NAME)

  set old_dir [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                         $Info(meshDirectoryName) \
                         $old_name ]

  set new_dir [file join $ModelProperty(MODEL_DIRECTORY,absolute) \
                         $Info(meshDirectoryName) \
                         $new_name ]

  if { [catch {file rename -force $old_dir $new_dir} msg ] } {
MSG "ERROR in mesh rename because: $msg"
    return
  }

  # Update data
  # =================
  set MeshDefine(meshNames) [lreplace $MeshDefine(meshNames) $index $index $new_name]

  set solvers_modified 0

  # Check solver meshes
  foreach id $Solver(ids) {
    set mn [DataField::getFieldValue Solver $id MESH]
    if { [string equal $mn $old_name] } {
      DataField::setFieldValue Solver $id MESH $new_name
      set solvers_modified 1
    }
  }

  if { $solvers_modified } {
    Solver::panelSave
  }

  MeshDefine::panelSave

  MeshDefine::updateMeshNamesLB


  # Save model file!!!
  #
  MenuExec::saveModelFile
}


# Update nof nodes/elements fields
#
proc MeshDefine::updateNofNodesAndElements {} {
  global MeshDefine ModelProperty

  set MeshDefine(nofNodes) ""
  set MeshDefine(nofElements) ""

  set result [MeshDefine::readNofNodesAndElements $MeshDefine(MESH_NAME)]

  set MeshDefine(nofNodes) [lindex $result 0]
  set MeshDefine(nofElements) [lindex $result 1]
}


# Read Elmer mesh nof nodes etc. from mesh header
#
proc MeshDefine::readNofNodesAndElements {mesh_name} {
  global MeshDefine ModelProperty

  set nof_nodes ""
  set nof_elements ""

  set dir [ModelProperty::getMeshDirectory $mesh_name]

  set fname [file join $dir "mesh.header"]

  if { ![catch {set ch [open $fname]}] } {

    # NOTE: Remove multiple inner spaces, so that split
    # gives correct list!!!
    #
    set einfo [split [Util::stringTrim [gets $ch]]]
    set einfo [split $einfo]

    catch { close $ch }

    set nof_nodes [lindex $einfo 0]
    set nof_elements [lindex $einfo 1]
  }

  return [list $nof_nodes $nof_elements]
}


#=================================#
#        Button state procs       #
#=================================#

proc MeshDefine::setAllButtonsState {state} {
  global MeshDefine 

  MeshDefine::setMeshAddButtonsState $state
  MeshDefine::setMeshEditButtonsState $state
  MeshDefine::setMeshUseButtonsState $state
  MeshDefine::setCommandButtonsState $state
}


proc MeshDefine::setMeshAddButtonsState {state} {
  global MeshDefine 

  if { $state != "normal" && $state != "disabled" } {
    return
  }

  if { $MeshDefine(acceptNewMesh) } {
    $MeshDefine(newButton) configure -state $state
  } else {
    $MeshDefine(newButton) configure -state disabled
  }

  $MeshDefine(addButton) configure -state $state
}


proc MeshDefine::setMeshEditButtonsState {state} {
  global MeshDefine 

  if { $state != "normal" && $state != "disabled" } {
    return
  }

  $MeshDefine(deleteButton) configure -state $state
  $MeshDefine(renameButton) configure -state $state

  $MeshDefine(removeButton) configure -state $state
  $MeshDefine(displayButton) configure -state $state
}


proc MeshDefine::setMeshUseButtonsState {state} {
  global MeshDefine 

  if { $state != "normal" && $state != "disabled" } {
    return
  }

  $MeshDefine(useButton) configure -state $state

  $MeshDefine(useButton) configure -state $state
}



proc MeshDefine::setCommandButtonsState {state} {
  global MeshDefine ModelFlags

  if { $state != "normal" && $state != "disabled" } {
    return
  }

  #--These are strictly geometry type dependent
  #
  # Mesh geometry states
  if { !$ModelFlags(GEOMETRY_TYPE_CAD) } {
    $MeshDefine(structureButton) configure -state disabled
    $MeshDefine(densityButton) configure -state disabled
    $MeshDefine(generateButton) configure -state disabled

  # Cad geometry states 
  } else {
    $MeshDefine(structureButton) configure -state $state
    $MeshDefine(densityButton) configure -state $state
    $MeshDefine(generateButton) configure -state $state
  }
}



#===========================================================#
#              Panel Save, Ok, Cancel, Apply procs          #
#===========================================================#

# ----------
# SAVE PANEL
# ----------
#
proc MeshDefine::panelSave { {inform_front 1} } {
  global Equation Info MeshDefine Model ModelProperty SolverSystem

  # Store previous mesh index
  #
  set old_index $Model(currentMeshIndex)

  # Update current field values to lists
  if { $MeshDefine(MESH_INDEX) >= 0 } {
    MeshDefine::updateLists "update" $MeshDefine(MESH_INDEX)
  }

  # Copy MeshDefine-lists to Model-lists
  MeshDefine::setModelLists

  # Compress MeshDefine-bg-mesh-lists
  # NOTE: We store Bg-mesh list as indexed and not as sparse
  # lists as they are used in the panel!
  #
  MeshDefine::makeIndexedBgMeshLists

  # Generated mesh will be always the current mesh!
  #
  set Model(currentMeshIndex) $MeshDefine(MESH_INDEX)
  set Model(currentMeshName) $MeshDefine(MESH_NAME)

  # Update mesh directory
  ModelProperty::setCurrentMeshDirectory
  Util::cpp_exec modelPropertyPanelOk

  # Check if mesh parameter has been changed
  # CORRECT THIS !!!###!!!
  if { 1 == 0 } {
    set Model(Mesh,needsUpdate) 1

    if {$inform_front} {
      set Model(Database,needsUpdate) 1
      set Model(Viewfactors,needsUpdate,mesh) 1
      set Model(Front,needsUpdate) 1
      set Model(Solver,needsUpdate) 1
    }
  }
  
  # Save data in cpp-side
  #
  Util::cpp_exec meshDefinePanelOk

  # NOTE: We have to copy back sparse Model-array list
  # to MeshDefine-list, because bg-mesh-list were
  # compressed before saving to cpp-side!!!
  #
  MeshDefine::setMeshDefineLists

  # Check if a new mesh must be loaded
  #
  if { $MeshDefine(edit_mode) == "display" } {
    MeshDefine::checkMeshUpdate $old_index
  }

  Panel::panelDataChanged 0 MeshDefine
  Panel::panelDataModified 0 MeshDefine

  set MeshDefine(meshListChanged) 0
  set MeshDefine(addedMeshIndices) ""

  set MeshDefine(MESH_NAME,err) 0
  set MeshDefine(MESH_NAME,mod) 0
  Widget::setFieldStatus MeshDefine MESH_NAME
  $MeshDefine(nameEntry) configure -state disabled -bg $Info(nonActiveBg)

  set MeshDefine(MESH_H,err) 0
  set MeshDefine(MESH_H,mod) 0
  Widget::setFieldStatus MeshDefine MESH_H

  set MeshDefine(MESH_F,err) 0
  set MeshDefine(MESH_F,mod) 0
  Widget::setFieldStatus MeshDefine MESH_F

  set MeshDefine(edit_mode) ""
  set MeshDefine(command_mode) ""

  set MeshDefine(edited,GridH) 0
  set MeshDefine(edited,GridParameter) 0
}


proc MeshDefine::panelOk {w} {
  global MeshDefine Model ModelProperty

  #---No changes
  if { !$MeshDefine(dataChanged) &&
       !$MeshDefine(meshListChanged)
     } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![MeshDefine::checkPanelData] } {
    return 0
  }

  #---Save data
  MeshDefine::panelSave

  Panel::cancel $w

  return 1
}


proc MeshDefine::panelApply {} {
  global MeshDefine

  #---No changes
  if { !$MeshDefine(dataChanged) &&
       !$MeshDefine(meshListChanged)
     } {
    return 0
  }

  if { ![MeshDefine::checkPanelData] } {
    return -1
  }

  MeshDefine::panelSave

  return 1
}


proc MeshDefine::panelCancel {w} {
  global ModelProperty

  if { ![Panel::verifyCancel MeshDefine] } {
    return
  }

  Panel::cancel $w

  return "cancel"
}


# Return 1 = ok, 0 = error
# 
proc MeshDefine::checkPanelData {} {
  global Info MeshDefine ModelFlags

  #-Check mesh name
  if { ( $MeshDefine(edit_mode) == "new" ||
         $MeshDefine(edit_mode) == "rename"
       )  && 
       ![MeshDefine::checkMeshName]
     } {
    return 0
  }

  #-These are relevant only for actual meshing
  #
  if { $MeshDefine(command_mode) == "mesh_structure" ||
       $MeshDefine(command_mode) == "mesh_density"   ||
       $MeshDefine(command_mode) == "generate"
     } {
      if { ![MeshDefine::checkH] ||
           ![MeshDefine::checkF] 
        } {
      return 0
    }
  }

  # Data ok
  return 1
}


#====================#
# Check helper procs #
#====================#

# Set mesh for all solvers
#
proc MeshDefine::setSolverMesh { mesh_name } {
  global MeshDefine Model Solver
    
  foreach id $Solver(ids) {
     DataField::setFieldValue Solver $id MESH $mesh_name
  }

  Solver::updateActiveMeshInfo

  Solver::panelSave
}


# Return 1 = ok, 0 = error
#
proc MeshDefine::checkF {} {
  global Info MeshDefine Model

  set f $MeshDefine(MESH_F)

  if { $f <= 0 } {
    set msg "Incorrect Mesh scaling factor:\nShould be > 0 !"

    set Info(messageIcon) error
    Message::dispOkMessage $msg "$Info(FRONT_NAME) error message" $MeshDefine(winName)
    return 0
  }

  return 1
}
  

# Return 1 = ok, 0 = error
#
proc MeshDefine::checkH {} {
  global Info MeshDefine Model

  set h $MeshDefine(MESH_H)

  if { $h <= 0 } {
    set msg "Incorrect Model mesh H:\nShould be > 0 !"

    set Info(messageIcon) error
    Message::dispOkMessage $msg "$Info(FRONT_NAME) error message" $MeshDefine(winName)
    return 0
  }

  return 1
}


# Return 1 = ok, 0 = error
#
proc MeshDefine::checkMeshName {} {
  global Info MeshDefine ModelFlags ModelProperty

  # If new mesh being added, but no generate command is applied
  if { $MeshDefine(edit_mode) == "new" &&
       $MeshDefine(command_mode) == "" 
     } {

    #--Cad geometry
    if { $ModelFlags(GEOMETRY_TYPE_CAD) } {


      set Info(messageIcon) error
      set msg [list "***WARNING: To store a new mesh, you must first generate the mesh!\n\n" \
                    "$Info(anywayOk)" ]

      if { ![Message::verifyAction $Info(powerUser) $msg] } {
        return 0
      }

      #set call_mesh_generator 0
      #MeshDefine::meshGenerate $call_mesh_generator

    #--Mesh geometry
    # NOTE: This is a "fake" generate, just to add
    # the mesh name to the model!
    } else {
      set call_mesh_generator 0
      MeshDefine::meshGenerate $call_mesh_generator
    }
  }
  
  set MeshDefine(MESH_NAME) [string trim $MeshDefine(MESH_NAME)]

  if { $MeshDefine(MESH_NAME) == "" } {

    set msg "Mesh name missing!"

    set Info(messageIcon) error
    Message::dispOkMessage $msg "$Info(FRONT_NAME) error message" $MeshDefine(winName)
    return 0
  }
  
   if { ![ModelProperty::checkNameForFile $MeshDefine(MESH_NAME) "mesh"] } {
    return 0
  }

  return 1
}


# Check if a new mesh should be loaded
#
proc MeshDefine::checkMeshUpdate { mesh_index } {
  global Info MeshDefine Model ModelFlags ModelProperty UserSetting

  # Load new mesh if current name changed or no mesh exists
  if { $Model(currentMeshIndex) != $mesh_index ||
       !$ModelFlags(DRAW_SOURCE_MESH)
     } {

    # Model file needs update!
    set Model(Front,needsUpdate) 1

    # Remove possible currently loaded mesh
    if { $ModelFlags(DRAW_SOURCE_MESH) } {
      MenuExec::unloadMesh
    }

    # Load nem mesh
    if { $UserSetting(AUTO_LOAD_MESH) && 
         $Model(currentMeshIndex) != $Info(NO_INDEX)
       } {
      MenuExec::loadMesh
    }
  }
}


# Check if mesh is used by any solver
#
proc MeshDefine::checkSolverMeshUsage {mesh_name} {
  global MeshDefine Solver

  foreach id $Solver(ids) {

    set mn [DataField::getFieldValue Solver $id MESH]

    if { [string equal $mn $mesh_name] } {
      return 1
    }
  }

  return 0
}


# Check and update other parameters (Solver, GridStrctureMeshH) when mesh
# is deleted or removed
#
proc MeshDefine::checkParameterMeshData { removed_names removed_indices } {
  global MeshDefine Solver ObjectTable
  
  set solvers_modified 0

  foreach name $removed_names idx $removed_indices {

    #--Update solver paramters mesh name
    #
    foreach id $Solver(ids) {

      set mn [DataField::getFieldValue Solver $id MESH]

      if { [string equal $mn $name] } {
        DataField::setFieldValue Solver $id MESH ""
        set solvers_modified 1
      }
    }

    #-Solver mesh info modified
    if { $solvers_modified } {
      Solver::panelSave

    #-Just update active mesh info
    } else {
      #--Update active mesh info
      Solver::updateActiveMeshInfo
    }

    # Update paramter data
    # ====================
    set gr_updated 0
    set gh_updated 0

    # Looop all object and check possible mesh related indices
    # ========================================================
    foreach id $ObjectTable(ids) {

      # Mesh structure (GridParameter) parameters
      # =========================================
      #
      if { [info exists ObjectTable($id,grMshIndcs)] &&
           $ObjectTable($id,grMshIndcs) != ""
         } {

        #-Mesh structure indices
        set pos 0
        foreach index $ObjectTable($id,grMshIndcs) {
        
          # Remove deleted
          if { $index == $idx } {
            set ObjectTable($id,grMshIndcs) [lreplace $ObjectTable($id,grMshIndcs) $pos $pos]
            set ObjectTable($id,grIds) [lreplace $ObjectTable($id,grIds) $pos $pos]
            set gr_updated 1
            incr pos -1
        
          # Decrease indices larger than deleted!
          } elseif { $index > $idx } {
            set new_index [expr $index -  1]
            set ObjectTable($id,grMshIndcs) [lreplace $ObjectTable($id,grMshIndcs) $pos $pos $new_index]
          }

          incr pos
        }
      }

      #-Exluded mesh indices
      if { [info exists ObjectTable($id,excldMshIndcs)] &&
           $ObjectTable($id,excldMshIndcs) != ""
         } {

        set pos 0
        foreach index $ObjectTable($id,excldMshIndcs) {

          # Remove deleted
          if { $index == $idx } {
            set ObjectTable($id,excldMshIndcs) [lreplace $ObjectTable($id,excldMshIndcs) $pos $pos]
            set gr_updated 1
            incr pos -1
        
          # Decrease indices larger than deleted!
          } elseif { $index > $idx } {
            set new_index [expr $index -  1]
            set ObjectTable($id,excldMshIndcs) [lreplace $ObjectTable($id,excldMshIndcs) $pos $pos $new_index]
          }

          incr pos
        }

        if { 0 == [llength $ObjectTable($id,excldMshIndcs)] } {

          set ObjectTable($id,excldMsh) 0
        }

      }

      # Mesh density (GridH) parameters
      # ===============================
      #
      if { [info exists ObjectTable($id,ghMshIndcs)] &&
           $ObjectTable($id,ghMshIndcs) != ""
         } {

        set pos 0
        foreach index $ObjectTable($id,ghMshIndcs) {

          # Remove deleted
          if { $index == $idx } {
            set ObjectTable($id,ghMshIndcs) [lreplace $ObjectTable($id,ghMshIndcs) $pos $pos]
            set ObjectTable($id,ghIds) [lreplace $ObjectTable($id,ghIds) $pos $pos]
            set gh_updated 1
            incr pos -1

          # Decrease indices larger than deleted!
          } elseif { $index > $idx } {
            set new_index [expr $index -  1]
            set ObjectTable($id,ghMshIndcs) [lreplace $ObjectTable($id,ghMshIndcs) $pos $pos $new_index]
          }

          incr pos
        }
      }
    } ;# ObjectTable-loop

    #-Save possibly modified parameter data
    if { $gr_updated } {
      StdPanelExec::panelSave GridParameter 1
    }

    if { $gh_updated } {
      StdPanelExec::panelSave GridH 1
    }

  } ;# For each removed meshname

}



#=============================#
# List and field update procs #
#=============================#

# Make compressed (indexed) list from bg-mesh lists
# Note these lists are stored in 'compressed' mode in model file
# because they are normally 'empty' and they should be constructed
# before saving data!!!
#
# On the other hand, in the panel we use "sparse" lists for convenience, so
# keep the sparse-mode versions in the equivalent
# Model-array lists and copy these to MeshDefine-lists when
# opening the panel!!
#
proc MeshDefine::makeIndexedBgMeshLists {} {
  global MeshDefine Model

  set Model(meshBgMeshFileIndices) ""

  set MeshDefine(meshBgMeshFiles) ""
  set MeshDefine(meshBgMeshActives) ""
  set MeshDefine(meshBgMeshControls) ""

  set index 0
  set count [llength $Model(meshNames)]

  while { $index < $count } {

    # If bg-mesh file given for the mesh
    #
    if { "" != [lindex $Model(meshBgMeshFiles) $index] } {

      lappend Model(meshBgMeshFileIndices) $index

      lappend MeshDefine(meshBgMeshFiles) [lindex $Model(meshBgMeshFiles) $index]
      lappend MeshDefine(meshBgMeshActives) [lindex $Model(meshBgMeshActives) $index]
      lappend MeshDefine(meshBgMeshControls) [lindex $Model(meshBgMeshControls) $index]
    }

    incr index
  }
}


proc MeshDefine::updateFields {index} {
global MeshDefine

  foreach fld $MeshDefine(fieldNames) lst $MeshDefine(listNames) {
    set MeshDefine($fld) [lindex $MeshDefine($lst) $index]
  }
}

proc MeshDefine::updateLists { mode {index ""} } {
global MeshDefine

  foreach fld $MeshDefine(fieldNames) lst $MeshDefine(listNames) {
    
    # Add
    if { [string equal -nocase $mode "add"] } {
      lappend MeshDefine($lst) $MeshDefine($fld)

    # Delete
    } elseif { [string equal -nocase $mode "delete"] } {
      set MeshDefine($lst) [lreplace $MeshDefine($lst) $index $index]

    # Update
    } else {
      set MeshDefine($lst) [lreplace $MeshDefine($lst) $index $index $MeshDefine($fld)]
    }
  }
}


proc MeshDefine::setMeshDefineLists {} {
global MeshDefine Model

  foreach lst $MeshDefine(listNames) {
    set MeshDefine($lst) $Model($lst)
  }
}


proc MeshDefine::setModelLists {} {
global MeshDefine Model

  foreach lst $MeshDefine(listNames) {
    set Model($lst) $MeshDefine($lst)
  }
}


#=================================#
# Mesh Define Fill and Edit procs #
#=================================#

proc MeshDefine(MESH_BG_MESH_ACTIVE,fillProc) {} {

  MeshDefine(MESH_BG_MESH_ACTIVE,editProc)
}


proc MeshDefine(MESH_BG_MESH_ACTIVE,editProc) {} {
  global MeshDefine GridParameter

  if { $MeshDefine(MESH_BG_MESH_ACTIVE) } {
    set MeshDefine(MESH_BG_MESH_FILE,act) 1
  } else {
    set MeshDefine(MESH_BG_MESH_FILE,act) 0
  }

  Widget::configureField MeshDefine MESH_BG_MESH_FILE

  if { [info exists GridParameter(allWidgets,MESH_BG_MESH_FILE)] &&
       [winfo exists $GridParameter(allWidgets,MESH_BG_MESH_FILE)]
     } {
    GridParameter(MESH_BG_MESH,editProc)
  }
}


# end ecif_tk_meshDefinePanel.tcl
# ********************
