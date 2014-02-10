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
#Module:    ecif_tk_bodyProperties.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for setting the body properties (name, color)
#
#************************************************************************


#--Edit body names and colors
#
proc BodyProperty::openPanel {} {
  global BodyProperty Info Model ModelFlags ObjectTable

  set w $BodyProperty(winName)
  set wgeom $BodyProperty(winGeometry)

  #--Store windows-id in globArray
  set id [winfo atom $w]
  set BodyProperty(propertiesWinId) $id

  #--Store last created window info in Info
  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow BodyProperty $id $BodyProperty(winTitle) $wgeom] } {
    return
  }  

  set BodyProperty(dataChanged) 0
  set BodyProperty(dataModified) 0
  
  set this $w 
  toplevel $w 
  focus $w

  #--Window properties
  wm title $w $BodyProperty(winTitle)
  wm geometry $w $wgeom

  BodyProperty::createPanelData

  #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)

  #---Outer frames
  set f1   [frame $w.f1]      ;# List box and entry button common frame
  set f11  [frame $w.f1.f1]   ;# List box frame
  set f12  [frame $w.f1.f2]   ;# Name edit entry and color handling buttons frame
  set f121 [frame $w.f1.f2.f1];# Name edit entry and color select button frame
  set f122 [frame $w.f1.f2.f2];# Color apply and reset buttons frame
  set f123 [frame $w.f1.f2.f3];# Delete body button frame
  set f2   [frame $w.f2]      ;# Apply+Ok+cancel buttons frame
  
  set lb [listbox $f11.lb -selectmode browse -height 15 -width 30 \
      -yscrollcommand [list $f11.sb_y set]]
  set sb_y [scrollbar $f11.sb_y -orient vertical  -command [list $f11.lb yview] ]

  pack $lb $sb_y -side left -fill y

  BodyProperty::fillBodyListBox $lb
    
  # Name entry
  # ==========
  set ne [entry $f121.ne -width 30  -exportselection 0]
  set cb [button $f122.cb -text "Select color" -width 12 -command "BodyProperty::selectBodyColor $w" ]
  pack $ne $cb -side top -padx $fpx2 -pady $fpy2

  # Apply, reset color buttons
  # ==========================
  set ab [button $f122.ab -text "Apply color" -width 12 -command "BodyProperty::applyBodyColor" ]
  set rb [button $f122.rb -text "Reset color" -width 12 -command "BodyProperty::resetBodyColor" ]

  $ab configure -state disabled
  $rb configure -state disabled

  set BodyProperty(applyBodyColorButton) $ab
  set BodyProperty(resetBodyColorButton) $rb

  pack $ab $rb -side top -padx $fpx2 -pady $fpy2
    
  set BodyProperty(listBox) $lb
  set BodyProperty(nameEntry) $ne
  set BodyProperty(colorButton) $cb
  
  bind $lb <ButtonRelease-1> BodyProperty::setNameAndColorContents 

  # Set name entry bindings
  # =======================
  set wdg $ne
  
  #-Key release
  set bnd "Panel::panelDataChanged 1 BodyProperty $wdg {%A %K}"
  bind $wdg <KeyRelease> $bnd

  #-Enter key
  set bnd "BodyProperty::applyBodyNameEntryContents"
  append bnd ";Panel::panelDataModified 0 BodyProperty $wdg {%A %K}"
  bind $wdg <KeyPress-Return> $bnd

  #-Ctrl+Enter key
  set bnd "BodyProperty::applyBodyNameEntryContents"
  append bnd ";BodyProperty::selectNextBody "
  append bnd ";Panel::panelDataModified 0 BodyProperty $wdg {%A %K}"
  bind $wdg <Control-Return> $bnd

  #-Leave field
  set bnd "BodyProperty::applyBodyNameEntryContents"
  append bnd ";BodyProperty::selectNextBody "
  append bnd ";Panel::panelDataModified 0 BodyProperty $wdg {%A %K}"
  bind $wdg <FocusOut> $bnd

  # Delete body button
  # ==================
  set db [button $f123.rb -text "Delete body" -width 12 -command "BodyProperty::deleteBody" ]

  if { $ModelFlags(GEOMETRY_TYPE_MESH) ||
       $Model(nofBodiesWithEquation) > 0
     } {
    $db configure -state disabled

  } else {
    $db configure -state normal
  }

  pack $db -side top -padx $fpx2 -pady [expr 10 * $fpy2]

  # Ok etc buttons
  # ==============
  set ap_st $Info(defaultApplyState)
  set ca_st $Info(defaultCancelState)

  set ok_btn [button $f2.ok -text OK -command "BodyProperty::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "BodyProperty::panelCancel $this" \
                                -state $ca_st ]
  set ap_btn [button $f2.apply -text Apply -command "BodyProperty::panelApply" \
                               -state $ap_st ]

  focus $ok_btn
    
  set BodyProperty(applyButton)  $ap_btn
  set BodyProperty(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1
  
  pack $f1 $f2 -side top -expand 1 -padx $fpx2 -pady $fpy2
  pack $f11 $f12 -side left 
  pack $f121 $f122 $f123 -side top 

  # Init panel
  $BodyProperty(listBox) selection set 0 0
  BodyProperty::setNameAndColorContents

}


proc BodyProperty::createPanelData {} {
  global BodyProperty ObjectTable

  # Pick body object ids
  set BodyProperty(ids) ""

  foreach id $ObjectTable(ids) {

    if { $ObjectTable($id,tp) != "B" } {
      continue
    }

    lappend BodyProperty(ids) $id
  }

  foreach id $BodyProperty(ids) {
    set BodyProperty($id,name) $ObjectTable($id,nm)
    set BodyProperty($id,name,old) $ObjectTable($id,nm)
    set BodyProperty($id,color) $ObjectTable($id,clr)
    set BodyProperty($id,color,old) $ObjectTable($id,clr)
  }
}


proc BodyProperty::fillBodyListBox { lb_wdg } {
  global BodyProperty

  set names ""
  foreach id $BodyProperty(ids) {
    lappend names $BodyProperty($id,name)
  }

  ListBox::fill $lb_wdg $names
}


###################
# SAVE etc. procs #
###################

proc BodyProperty::panelSave {} {
  global gMW BodyProperty Info Model ObjectTable

  set Model(Front,needsUpdate) 1

  Panel::panelDataChanged 0 BodyProperty 
  Panel::panelDataModified 0 BodyProperty 

  foreach id $BodyProperty(ids) {
    
    # Store new data to object table
    set ObjectTable($id,nm) $BodyProperty($id,name)
    set ObjectTable($id,clr) $BodyProperty($id,color)
  }

  Object::setBodyPairNames

  #-Store new values into Model
  Util::cpp_exec bodyPropertyPanelNamesOk 
  Util::cpp_exec bodyPropertyPanelColorsOk

  #-Update bodies menu list
  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""
}


proc BodyProperty::panelOk {w} {
  global BodyProperty

  #---No changes
  if { !$BodyProperty(dataChanged) } {
    Panel::cancel $w; return
  }
  
  #---Error in data
  if { ![BodyProperty::panelCheck] } {
    return
  }

  BodyProperty::panelSave

  Panel::cancel $w
}


proc BodyProperty::panelApply {} {
  global BodyProperty

  #---No changes
  if { !$BodyProperty(dataChanged) } {
    return
  }

  #---Error in data
  if { ![BodyProperty::panelCheck] } {
    return
  }

  # Backup panel data
  foreach id $BodyProperty(ids) {
    set BodyProperty($id,name,old) $BodyProperty($id,name)
    set BodyProperty($id,color,old) $BodyProperty($id,color)
  }
  
  BodyProperty::panelSave

  $BodyProperty(applyBodyColorButton) configure -state disabled
  $BodyProperty(resetBodyColorButton) configure -state disabled
}


proc BodyProperty::panelCancel {w} {
  global BodyProperty gMW ObjectTable

  if { ![Panel::verifyCancel BodyProperty] } {
    return
  }

  foreach id $BodyProperty(ids) {
    set ObjectTable($id,nm) $BodyProperty($id,name,old)
    set ObjectTable($id,clr) $BodyProperty($id,color,old)
  }

  #-Store new values into Model
  Util::cpp_exec bodyPropertyPanelNamesOk 
  Util::cpp_exec bodyPropertyPanelColorsOk

  #-Update bodies menu list
  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""

  Panel::panelDataChanged 0 BodyProperty 
  Panel::panelDataModified 0 BodyProperty 

  destroy $w
}


proc BodyProperty::panelCheck {} {
  global BodyProperty Info ObjectTable

  # Check body names
  foreach id1 $BodyProperty(ids) {

    foreach id2 $BodyProperty(ids) {

      if { $id2 <= $id1 } {
        continue
      }

      if { [string equal -nocase $BodyProperty($id1,name) $BodyProperty($id2,name)] } {
        
        set tg1 $ObjectTable($id,tg)
        set tg2 $ObjectTable($id,tg)

        set msg [list "NOTE: Bodies $tg1 and $tg2 have similar names:\n" \
                      "$BodyProperty($id1,name) \n" \
                      "$BodyProperty($id2,name) \n\n" \
                      $Info(anywayOk)]

        if { ![Message::verifyAction $Info(powerUser) $msg cancel question] } {
          return 0
        }
      }
    } ;# for id2
  } ;# for id1

  return 1
}



###############
# Other procs #
###############


# Body names panel bind proc: reacts to list box selection
#
proc BodyProperty::setNameAndColorContents {} {
  global BodyProperty

  set index [$BodyProperty(listBox) curselection]
  set BodyProperty(objectId) [lindex $BodyProperty(ids) $index]

  BodyProperty::setNameContents
  BodyProperty::setColorContents

  $BodyProperty(nameEntry) icursor 0
  $BodyProperty(nameEntry) selection range 0 end
}


# Body names panel bind proc: reacts to list box selection
#
proc BodyProperty::setNameContents {} {
  global BodyProperty ObjectTable

  set id $BodyProperty(objectId)

  #-Name entry
  $BodyProperty(nameEntry) delete 0 end
  #$BodyProperty(nameEntry) insert 0 $ObjectTable($id,nm)
  $BodyProperty(nameEntry) insert 0 $BodyProperty($id,name)
}


# Body names panel bind proc: reacts to list box selection
#
proc BodyProperty::setColorContents {} {
  global BodyProperty ObjectTable

  set id $BodyProperty(objectId)

  #-Color button
  set color $BodyProperty($id,color)
  $BodyProperty(colorButton) configure -bg $color -state normal
}


# Body names entry field bind proc:
#
proc BodyProperty::applyBodyNameEntryContents {} {
  global BodyProperty

  set index [$BodyProperty(listBox) curselection]
  set id $BodyProperty(objectId)

  set new_name [$BodyProperty(nameEntry) get]

  if { $new_name == "" } {
	return
  }
  
  #-NOTE: Do not trim in the entry field, only when storing!
  set new_name [string trim $new_name]

  set BodyProperty($id,name) $new_name

  $BodyProperty(listBox) delete $index $index
  $BodyProperty(listBox) insert $index $new_name
  $BodyProperty(listBox) selection set $index 
}


proc BodyProperty::selectBodyColor { parent_win } {
  global BodyProperty

  set id $BodyProperty(objectId)

  set new_color_tk ""
  set old_color_tk $BodyProperty($id,color)

  set new_color_tk [tk_chooseColor -parent $parent_win \
                                   -initialcolor $old_color_tk \
                                   -title "Choose body color" ]

  if {$new_color_tk == ""} {
    return
  }

  #-Color value
  #set new_name $new_color
  set BodyProperty($id,color) $new_color_tk

  BodyProperty::setColorContents

  $BodyProperty(applyBodyColorButton) configure -state normal
  $BodyProperty(resetBodyColorButton) configure -state normal

  Panel::panelDataChanged 1 BodyProperty
}


proc BodyProperty::applyBodyColor {} {
  global BodyProperty gMW ObjectTable

  set id $BodyProperty(objectId)

  set ObjectTable($id,clr) $BodyProperty($id,color)

  Util::cpp_exec bodyPropertyPanelColorsOk

  $BodyProperty(applyBodyColorButton) configure -state disabled
  $BodyProperty(resetBodyColorButton) configure -state normal

  #-Update bodies menu list
  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""
}


proc BodyProperty::resetBodyColor {} {
  global BodyProperty gMW ObjectTable

  set id $BodyProperty(objectId)

  set BodyProperty($id,color) $BodyProperty($id,color,old)
  set ObjectTable($id,clr) $BodyProperty($id,color)

  BodyProperty::setColorContents
  Util::cpp_exec bodyPropertyPanelColorsOk

  $BodyProperty(applyBodyColorButton) configure -state disabled
  $BodyProperty(resetBodyColorButton) configure -state disabled

  #-Update bodies menu list
  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""
}


proc BodyProperty::deleteBody {} {
  global BodyProperty gMW Info

  set msg [list \
     "NOTE: This command removes the selected body from the model!\n\n" \
     "\nPress OK to remove the body, otherwise press Cancel" \
     ]

  if { ![Message::verifyAction $Info(powerUser) $msg] } {
    return
  }


  Util::cpp_exec bodyPropertyPanelDeleteBodyOk
  
  BodyProperty::createPanelData
  BodyProperty::fillBodyListBox $BodyProperty(listBox)

  $BodyProperty(listBox) selection set 0 0
  BodyProperty::setNameAndColorContents

  MenuBuild::createBodyListMenu $gMW(Model,bodyList) ""
}



proc BodyProperty::getColorName {color_hex_value} {
  global Body Body(colorHex2Name) Info

  set Body(colorHex2Name) ""
  Util::cpp_exec colorHex2Name $color_hex_value

  tkwait variable Body(colorHex2Name)

  return $Body(colorHex2Name)
}


proc BodyProperty::selectNextBody {} {
  global BodyProperty

  Util::doWait 150

  set index [$BodyProperty(listBox) curselection]

  set count [expr [llength $BodyProperty(ids)] - 1]

  if { $index == $count } {
    set index 0
  } else {
    incr index
  }

  $BodyProperty(listBox) selection clear 0 end
  $BodyProperty(listBox) selection set $index $index

  event generate $BodyProperty(listBox) <ButtonRelease-1>
}


proc BodyProperty::setBodyName {} {
  global BodyProperty

  set id $BodyProperty(objectId)

  set new_name ""
  set old_name BodyProperty($id,name)

  gloabl entryResult
  Widget::dispDataEntry entryResult $old_name "Set body name"
  set new_name $entryResult
  
  if {$new_name != ""} {
    set BodyProperty($id,name) $new_name
  }

  Panel::panelDataChanged 1 BodyProperty

}


#ecif_bodyPropertiesPanel.tcl

