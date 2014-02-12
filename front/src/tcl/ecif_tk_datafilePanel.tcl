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
#Module:    ecif_tk_datafileDefPanel.tcl
#Language:  Tcl
#Date:      17.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  A panel for defining input and output directories and files
#
#************************************************************************


proc Datafile::openPanel { } {
  ## This procedure displays the input and output files definition panel
  ## Global variables
  global Info Common Model ModelProperty Datafile 
  global platform

  set w $Datafile(winName)
  set wgeom $Datafile(winGeometry)

  set id [winfo atom $w]
  set Datafile(winId) $id

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow Datafile $id $Datafile(winTitle) $wgeom] } {
    return
  }  

  set Datafile(dataChanged) 0
  set Datafile(dataModified) 0

  set this $w
  toplevel $w
  focus $w

  wm title $w $Datafile(winTitle)
  wm geometry $w $wgeom 

  Panel::resetFields Datafile

  Panel::initFields Datafile

  set id $Datafile(parameterId)
  DataField::formDataFields Datafile $Datafile($id,data)

  Panel::backupFields Datafile

  set Info(datafilePanelApplied) 0

  # Results directory to be displayed in the panel
  #
  if { $ModelProperty(RESULTS_DIRECTORY) != "" } {
    set Datafile(results_directory) $ModelProperty(RESULTS_DIRECTORY,absolute)
  } else {
    set Datafile(results_directory) $ModelProperty(MODEL_DIRECTORY,absolute)
  }

 #----WIDGET DEFINITION AND PACKING
  #
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)

  #---Outer frames
  set f1 [frame $w.f1] ;#File and dir entries
  set f2 [frame $w.f2] ;#Apply+Ok+cancel buttons frame
  #
  set lbl_wid 20
  set ent_wid 50

  # Editable-flag, globArray, fld, label
  set fields {
    { 0 Model         EMF_SAVE_DIR        "Model directory:"}
    { 0 Model         EMF_FILE            "Model file:"}
    { 0 Datafile      SOLVER_INPUT_FILE   "Solver input file:"}
    { 0 Datafile      MESH_INPUT_FILE     "Mesh input file:"}
    { 0 Datafile      results_directory   "Results directory:"}
    { 1 Datafile      OUTPUT_FILE         "Result file:"}
    { 1 Datafile      POST_FILE           "Postprocessor file:"}
    { 1 Datafile      RESTART_FILE        "Restart file:"}
    { 1 Datafile      RESTART_POSITION    "Restart timestep position:"}
    { 1 Datafile      GEBHARDT_FACTORS    "Gebhardt factors file:"}
    { 1 Datafile      VIEW_FACTORS        "View factors file:"}
  }

  set i 1
  foreach fld $fields {
    
    set active [lindex $fld 0]
    set arr [lindex $fld 1]
    set var [lindex $fld 2]
    set lbl [lindex $fld 3]

    
    if {$active} {
      set state normal
      set bg white
    } else {
      set state disabled
      set bg $Info(nonActiveBg)
    }

    set f [frame $f1.f$i]

    label $f.l  -text $lbl -width $lbl_wid -anchor w
    set Datafile(allWidgets,label,$var) $f.l
    
    set txtvar $arr
    append txtvar "($var)"    
 
    set wdg [entry $f.e  -textvariable $txtvar -width $ent_wid \
                         -font $Info(entryFont) -state $state -bg $bg]

    set Datafile(allWidgets,$var) $wdg
    
    if {$active} {
      bind  $wdg <KeyRelease> "Panel::panelDataChanged 1 Datafile $wdg {%A %K}"
      Widget::setEntryBindings non_standard Datafile $var $wdg
    }

    set Datafile($var,err) 0
    set Datafile($var,mod) 0

    pack $f.l $f.e -side left -anchor w
    pack $f -side top -padx $fpx2 -pady $fpy2 -anchor w

    incr i
  }

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $f2.ok -text OK -command "Datafile::panelOk $this"]
  set cn_btn [button $f2.cancel -text Cancel -command "Datafile::panelCancel $this" \
                                -state $ca]
  set ap_btn [button $f2.apply -text Apply -command "Datafile::panelApply" \
                               -state $ap]

  focus $ok_btn
  set Datafile(applyButton)  $ap_btn
  set Datafile(cancelButton) $cn_btn

  pack $ok_btn $cn_btn $ap_btn -side left -expand 1 -padx $fpx1


  pack $f1 -side top -expand 1 -anchor w -padx $fpx2 -pady $fpy2
  pack $f2 -side top -expand 1 -padx $fpx2 -pady $fpy2
}


proc Datafile::getDefaultSolverInputFile {} {
  global Model

  if { $Model(CASE_NAME) == "" } {
    return ""
  } else {
    return $Model(CASE_NAME).sif
  }
}


proc Datafile::getDefaultMeshInputFile {} {
  global ModelProperty

  if { $ModelProperty(MODEL_NAME) == "" } {
    return ""
  } else {
    return $ModelProperty(MODEL_NAME).mif
  }
}


proc Datafile::getDefaultOutputFile {} {
  global Model

  if { $Model(CASE_NAME) == "" } {
    return ""
  } else {
    return $Model(CASE_NAME).dat
  }
}


proc Datafile::getDefaultPostFile {} {
  global Model

  if { $Model(CASE_NAME) == "" } {
    return ""
  } else {
    return $Model(CASE_NAME).ep
  }
}


proc Datafile::panelSave { {inform_front 1} } {
  global Info Datafile Model

  #-Store old values 
  Panel::backupFields Datafile

  #--Form parameter data
  set Datafile(ids) 1
  DataField::formNonStandardParameter Datafile 1 "Datafile1"

  #--Write data into model
  if {$inform_front} {
    set Model(Front,needsUpdate) 1
  }

  Panel::panelDataChanged 0 Datafile 
  Panel::panelDataModified 0 Datafile 

  set modified 0
  Panel::compressParameters Datafile modified

  Util::cpp_exec datafilePanelOk

  Panel::uncompressParameters Datafile

}


proc Datafile::panelOk {w} {
  global Datafile

  #---No changes
  if { !$Datafile(dataChanged) } {
    Panel::cancel $w; return
  }

  #---Error in data
  if { ![Datafile::panelCheck] } {
    return
  }

  Datafile::panelSave
  Panel::cancel $w
}


proc Datafile::panelApply {} {
  global Datafile

  #---No changes
  if { !$Datafile(dataChanged) } {
    return
  }

  #---Error in data
  if { ![Datafile::panelCheck] } {
    return
  }

  Datafile::panelSave
}


proc Datafile::panelCancel {w} {
  global Datafile

  if { ![Panel::verifyCancel Datafile] } {
    return
  }

  #-Restore old values
  Panel::restoreFields Datafile

  Panel::cancel $w
}


# Return 1 = ok, 0 = error
#
proc Datafile::panelCheck {} {
  global Datafile Model Info

  Panel::initFields Datafile

  #--Check path format (backslash --> slash etc.)
  PanelCheck::checkFields Datafile $Datafile(allFields)

  #--NOTE: We accept only simple filenames for the following files:
  #
  set files {OUTPUT_FILE POST_FILE RESTART_FILE GEBHARDT_FACTORS VIEW_FACTORS}
  set names {"Result" "Postprocessor" "Restart" "Gebhardt factors" "View factors"}
 
  set name_list ""

  foreach fname $files nm $names {
 
    set fn $Datafile($fname)
    
    #--Only simple filenames are ok
    if { [Util::isSimpleFilename $fn] } {
      continue
    }

    #--Error, add name to the list!
    append name_list "$nm\n"
  }
 
  if { $name_list != "" } {

    set Info(messageIcon) error
    set msg [list "ERROR: Only simple file names are allowed for the following files:\n\n" "$name_list \n\n"]
    Message::dispOkMessage  $msg \
                            "$Info(FRONT_NAME) error message" \
                            $Datafile(winName) 

    return 0
  }

  #--Check restart file name
  if { ![Datafile::checkRestartFile] } {
    return 0
  }

  #--Check restart file position
  if { ![Datafile::checkRestartPosition] } {
    return 0
  }

  # NO OTHER CHECKING!!!
  return 1;

  # Old stuff
  # ---------
  set files {OUTPUT_FILE POST_FILE GEBHARDT_FACTORS VIEW_FACTORS}
 
  # If directory is given in some of the filenames, check that directory
  # is well defined!
  #
  foreach filename $files {
 
    set fn [set Datafile($filename)]
    
    # Simple filenames are always ok
    # ==============================
    if { [Util::isSimpleFilename $fn] } {
      continue
    }
 
    # Volume relative
    # ===============
    #-In Windows names like /elmer/Models or e:elmer/Models
    # are not absolute, but "volumerelative"
    if { "volumerelative" == [file pathtype $fn] } {
      set fn [Util::makePathAbsolute $fn]
      set Datafile($filename) $fn
    }

    # NOTE: Existence of the directory is NOT currently checked
    # It would be a bit tricky to do (it can be
    # below model-directory or results-directory and
    # these can even be changed after being checked here etc.
    # ==> Solver must create the final directories if needed!!!

    continue ;# We do NOT check the existence of the directory!

    #-Check if directoy exists
    set dname [file dirname $fn]
    if { ![file isdirectory $dname] } {
     
      set msg [list "Directory for the filename: $fn\n" \
                    "does not exist! Is it Ok to crete the directory?" ]

      if { ![Message::verifyAction $Info(advancedUser) $msg ok question] } {
        return 0

      } else {
        # NOTE: This would NOT be correct!!! We need the proper path
        # (model-directory or results-directory) in front of "dname"
        Util::createDirectory $dname
      }
    }
  }

  #-Ok
  return 1
}


proc Datafile::checkRestartFile {} {
  global Datafile Info

  # Impossible file name!
  #
  if { $Datafile(RESTART_FILE) == $Datafile(SOLVER_INPUT_FILE) ||
       $Datafile(RESTART_FILE) == $Datafile(MESH_INPUT_FILE) ||
       $Datafile(RESTART_FILE) == $Datafile(POST_FILE) ||
       $Datafile(RESTART_FILE) == $Datafile(GEBHARDT_FACTORS) ||
       $Datafile(RESTART_FILE) == $Datafile(VIEW_FACTORS) 
     } {

      set Info(messageIcon) error
      Message::dispOkMessage {"Invalid (reserved) Restart file name!"} \
                    "$Info(FRONT_NAME) error message" \
                    $Datafile(winName) 

      return 0
    }

  # Output file name same as restart file name!
  #
  if { $Datafile(RESTART_FILE) == $Datafile(OUTPUT_FILE) } {

    set msg [list "NOTE: Restart file will be overwritten by the output file!\n\n" \
                 $Info(anywayOk) ]

    if { ![Message::verifyAction $Info(powerUser) $msg] } {
      return 0
    }
  }


  return 1
}


proc Datafile::checkRestartPosition {} {
  global Datafile Info

  # Position given, but not the file name!
  #
  if { $Datafile(RESTART_POSITION) != "" &&
       $Datafile(RESTART_FILE) == ""
     } {

    set Info(messageIcon) error
    Message::dispOkMessage {"Restart file name missing!"} \
                  "$Info(FRONT_NAME) error message" \
                  $Datafile(winName) 

    return 0
  }

  
  # Filename give ,but not the position!
  #
  if { $Datafile(RESTART_FILE) != "" &&
       ( $Datafile(RESTART_POSITION) == "" ||
         $Datafile(RESTART_POSITION) == 0
       )
     } {

    set msg [list "No restart file position given\n\n" \
                  "Last position will be used!\n\n\n"  \
                  $Info(continueOk) ]

    if { ![Message::verifyAction $Info(powerUser) $msg] } {
      return 0
    }
  }

  return 1
}


#End ecif_datafileDefPanel.tcl

