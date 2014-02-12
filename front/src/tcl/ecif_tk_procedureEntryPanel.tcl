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
#Module:    ecif_tk_procedureEntryPanel.tcl
#Language:  Tcl
#Date:      06.11.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Displays a panel with for a field's procedure entry
#
#************************************************************************


#==========================#
#   ProcedureEntryPanel    # 
#==========================#

# Procedure displays a dialog with entry fields for procedure names.
# Result is put into the namespace variable resultId, where Id is the
# id-number for the namespace
#
# When panel is **saved**, **callback-function** is called in the form:
# $callback $id $globArray $fld
#
# Result: delivered using variable result$id
#
proc ProcedureEntryPanel::create { globArray fld callback} {
  global Info Equation Model
  upvar #0 $globArray theArray

  # Increment (unique) id
  variable id
  incr id

  # Namespace variables with id specifier
  #
  variable applyButton$id
  variable callback$id
  variable cancelButton$id
  variable compiler$id
  variable compilerEntry$id
  variable dataChanged$id
  variable field$id
  variable fieldName$id
  variable functionName$id ""
  variable functionNameEntry$id
  variable globArray$id
  variable isMatcProc$id
  variable libraries$id
  variable librariesEntry$id
  variable libraryName$id ""
  variable libraryName$id,abs 0
  variable libraryNameEntry$id
  variable matcWidgetGroup$id
  variable result$id
  variable sourceName$id ""
  variable sourceName$id,abs 0
  variable sourceNameEntry$id
  variable updated$id 0
  variable variableName$id "none"
  variable variableNameButton$id
  variable variableNameMenu$id
  variable variableNamePrev$id
  variable window$id

  set globArray$id $globArray
  
  if { $theArray(currentProcedureDir) == "" } {
    set theArray(currentProcedureDir) $Info(lastProcedureDir)
  }

  set callback$id $callback
  set field$id $fld

  set fname [DataField::getFieldProperty $globArray $fld SifName]

  if { $fname == "" } {
    set fname [DataField::fieldNameGuiToSif $fld]
  }

  set coord_index [DataField::getFieldProperty $globArray $fld CoordinateIndex]
  if { $coord_index > 0 } {
    ;#append fname "_$coord"
  }

  set fieldName$id $fname

  set dataChanged$id 0

  if { [info exists theArray($fld,proc)] &&
       $theArray($fld,proc)
     } {
    set fieldValues [DataField::splitProcedureEntryValue [set theArray($fld)] ]
    set libraryName$id [lindex $fieldValues 0]
    set functionName$id [lindex $fieldValues 1]
    set variableName$id [lindex $fieldValues 2]
    set sourceName$id [lindex $fieldValues 3]

  } else {
    if { [info exists theArray($fld,variables)] &&
          $theArray($fld,variables) != ""
       } {
      set variableName$id $theArray($fld,variables)
    }
  }

  if { "" == [set variableName$id] } {
    set variableName$id "none"
  }

  set variableNamePrev$id [set variableName$id]

  # Do we have a Matc-proc
  #
  if { $Info(matcMarker) == [string index [set functionName$id] 0] } {
    set isMatcProc$id 1
    set functionName$id [ProcedureEntryPanel::getMatcFunctionInName $id]
  } else {
    set isMatcProc$id 0 
  }

  # Set generic compiler command
  # ----------------------------
  if { ![info exist compiler$id] || "" == [set compiler$id] } {
    set EH $Info(ELMER_HOME)
    set idir [file join $EH include]
    
    #-Win32
    if { $Info(platform) == "windows" } {
      regsub -all "/" $idir "\\" idirWIN32

      set compiler$id "f90.exe -dll -I $idirWIN32"

    #-Unix  
    } else {  
      # NOTE: verbode flag (-v) is used  cause SGI f90
      # compiler is otherwise complety silent!

      set compiler$id "f90 -v -I$idir -shared -o"
    }
  }

  # Set default libraries
  # ---------------------
  set libraries$id $Info(procedureLibraries)

  # If nothing entered yet
  if { ![info exist libraries$id] || "" == [set libraries$id] } {
    set EH $Info(ELMER_HOME)
    set ldir [file join $EH lib]
  
    #-Win32
    if { $Info(platform) == "windows" } {
      regsub -all "/" $ldir "\\" ldirWIN32
      set libraries$id "$ldirWIN32\\solver.lib"

    #-Unix  
    } else {  
      set libraries$id ""
    }
  }

  # NOTE: Only one window per panel/field can exists
  set window$id .wPEP$globArray$fld

  set w [set window$id]
  set ProcedureEntry($id,winName) $w

  # Set window title
  set fn [string toupper [set fieldName$id]]
  set title "Procedure definition for $fn"
  set ProcedureEntry($id,winTitle) $title

  set wid [winfo atom $w]
  set ProcedureEntry($id,winId) $wid

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow ProcedureEntry $wid $ProcedureEntry($id,winTitle) ""] } {
    return
  }  

  toplevel $w
  wm title $w $ProcedureEntry($id,winTitle)

  # We collect Matc-proc value dependent widgets into list list
  set matcWidgetGroup ""

  # Principal frames
  frame $w.f1 ;# Library file, function, source, variable
  frame $w.f2 ;# Libraries, compiler
  frame $w.f3 ;# Edit, Compile buttons
  frame $w.f4 ;# Ok, Cancel buttons

  # Library file etc. subframes
  frame $w.f11
  frame $w.f12
  frame $w.f13
  frame $w.f14

  # Libraries,compiler subframes
  frame $w.f21
  frame $w.f22

  set lw 15
  set ew 35

  #---Entry fields and their labels

  #--Library file name 
  label $w.lLibrary -text "Library file:" -width $lw 

  set e [entry $w.eLibrary -width $ew -textvariable ProcedureEntryPanel::libraryName$id \
                          -font $Info(entryFont)]
  set libraryNameEntry$id $e
  bind $e <KeyRelease> "ProcedureEntryPanel::updated $id $e"

  lappend matcWidgetGroup $e

  # File browser button
  set browse_wdg [button $w.bbLibrary -text "..." -width 1 -height 0 \
                    -command "Util::fillFileEntry $e \
                              -array $globArray \
                              -absVar ProcedureEntryPanel::libraryName$id,abs \
                              -keepExtVar ProcedureEntryPanel::libraryName$id,abs \
                              -initDirFld currentProcedureDir \
                              -resultAbsFld currentProcedureDir" ]

  lappend matcWidgetGroup $browse_wdg

  # Abs check box
  set cb_wdg [checkbutton $w.abLibrary -text "Abs" -variable ProcedureEntryPanel::libraryName$id,abs \
                -command ""]

  lappend matcWidgetGroup $cb_wdg

  #--Function name
  label $w.lFunction -text "Function name:" -width $lw
  set e [entry $w.eFunction -width $ew -textvariable ProcedureEntryPanel::functionName$id\
                            -font $Info(entryFont)]
  set functionNameEntry$id $e
  bind $e <KeyRelease> "ProcedureEntryPanel::updated $id $e"

  # Is Matc-proc  check box
  set cb_wdg [checkbutton $w.isMatcProc -text "Is Matc proc" -variable ProcedureEntryPanel::isMatcProc$id \
                -command "ProcedureEntryPanel::checkMatcProc $id"]



  #--Variable name (an option menu)
  #
  if { [info exists theArray(targetMask)] } {
    set targetMask $theArray(targetMask)
  } else {
    set targetMask $Equation(problemMask)
  }

  label $w.lVariable -text "Variable:" -width $lw
	set accept_array 1
  set val_list [lindex [Panel::getCurrentVariables $targetMask $accept_array] 1]
  set val_list [linsert $val_list 0 "none"]

  set om [Panel::createOptionMenuWidget $w.mVariable ProcedureEntryPanel::variableName$id $val_list]
  set variableNameMenu$id $om
  set variableNameButton$id $w.mVariable

  $w.mVariable configure -width 15 -indicatoron 0 -relief raised \
                         -activebackground $Info(optionMenuBg)

  # NOTE: This seems to be the only event which is handled by the *.menu widget
  # We have to use this event to get the selection, when the menu is pulldown
  # and selection is done after that (not with a single "mouse-stroke")
  bind $om <<MenuSelect>> "ProcedureEntryPanel::checkVariable $id"

  # NOTE: ..menu widget does not respond to this event!
  #bind $om <ButtonRelease-1> "ProcedureEntryPanel::checkVariable $id $w.e4"
  # NOTE: This is not needed, when we use the <<MenuSelect>> event
  #bind $om <ButtonRelease-1> "ProcedureEntryPanel::checkVariable $id"


  #--Source name
  label $w.lSource -text "Source file:" -width $lw
  set e  [entry $w.eSource -width $ew -textvariable ProcedureEntryPanel::sourceName$id\
                           -font $Info(entryFont)]
  set sourceNameEntry$id $e
  bind $e <KeyRelease> "ProcedureEntryPanel::updated $id $e"
  bind $e <KeyPress-Return> "+ProcedureEntryPanel::editSource $id"

  lappend matcWidgetGroup $e

  # File browser button
  set browse_wdg [button $w.bbSource -text "..." -width 1 -height 0 \
                    -command "Util::fillFileEntry $e \
                                -absVar ProcedureEntryPanel::sourceName$id,abs \
                                -keepExtVar ProcedureEntryPanel::sourceName$id,abs"]

  lappend matcWidgetGroup $browse_wdg

  # Abs check box
  set cb_wdg [checkbutton $w.abSource -text "Abs" -variable ProcedureEntryPanel::sourceName$id,abs \
                -command ""]

  lappend matcWidgetGroup $cb_wdg

  #--Compiler command
  label $w.lCompiler -text "Compiler:" -width $lw
  set e  [entry $w.eCompiler -width $ew -textvariable ProcedureEntryPanel::compiler$id\
                             -font $Info(entryFont)]
  set compilerEntry$id $e

  lappend matcWidgetGroup $e

  #--Libraries
  label $w.lLibraries -text "Libraries:" -width $lw
  set e  [entry $w.eLibraries -width $ew -textvariable ProcedureEntryPanel::libraries$id\
                              -font $Info(entryFont)]
  set librariesEntry$id $e
  bind $e <KeyRelease> "ProcedureEntryPanel::updated $id $e"

  lappend matcWidgetGroup $e

  #---Button and button commands
  #-Edit, compile buttons
  eval "button $w.compileSource -text \"Compile source\" -command {ProcedureEntryPanel::compileSource $id }"
  eval "button $w.editSource -text \"Edit source\" -command {ProcedureEntryPanel::editSource $id }"
  
  lappend matcWidgetGroup $w.compileSource
  lappend matcWidgetGroup $w.editSource

  #-Ok, Cancel buttons

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)
  
  set ok_btn [button $w.ok -text OK -command "ProcedureEntryPanel::ok $id"]
  set cn_btn [button $w.cancel -text Cancel -command "ProcedureEntryPanel::cancel $id" -state $ca]
  set ap_btn [button $w.apply -text Apply -command "ProcedureEntryPanel::apply $id" -state $ap]

  set applyButton$id $ap_btn
  set cancelButton$id $cn_btn

  set matcWidgetGroup$id  $matcWidgetGroup


  #---Pack widgets
  #
  pack $w.f1 $w.f2 $w.f3 $w.f4 -side top -fill y -fill x -expand 1 -anchor c -pady 10 -padx 5

  #pack $textW -side top 

  pack $w.f11 $w.f12 $w.f13 $w.f14 -side top -in $w.f1 -anchor w -fill x -padx 4 -pady 4
  pack $w.f21 $w.f22 -side top -in $w.f2 -anchor w -fill x -padx 4 -pady 4

  #-Frame 1 stuff
  pack $w.lLibrary -side left -in $w.f11 -anchor w
  pack $w.eLibrary -side left -in $w.f11 
  pack $w.bbLibrary -side left -in $w.f11 -padx 2
  pack $w.abLibrary -side left -in $w.f11

  pack $w.lFunction -side left -in $w.f12 -anchor w
  pack $w.eFunction -side left -in $w.f12 
  pack $w.isMatcProc -side left -in $w.f12 

  pack $w.lVariable -side left -in $w.f13 -anchor w
  #pack $w.eVariable -side left -in $w.f13 -fill x
  pack $w.mVariable -side left -in $w.f13 

  pack $w.lSource -side left -in $w.f14 -anchor w
  pack $w.eSource -side left -in $w.f14 
  pack $w.bbSource -side left -in $w.f14 -padx 2
  pack $w.abSource -side left -in $w.f14 

  #-Frame 2 stuff
  pack $w.lCompiler -side left -in $w.f21 -anchor w
  pack $w.eCompiler -side left -in $w.f21

  pack $w.lLibraries -side left -in $w.f22 -anchor w
  pack $w.eLibraries -side left -in $w.f22

  #-Frame 3 stuff
  pack $w.editSource -side left -in $w.f3 -expand 1 -anchor e -padx 30 
  pack $w.compileSource -side left -in $w.f3 -expand 1 -anchor w -padx 30 

  #pack $w.variable -side left -in $w.f3 -expand 1 -anchor c

  #-Frame 4 stuff
  pack $ok_btn -side left -in $w.f4 -expand 1 -anchor e
  pack $cn_btn -side left -in $w.f4 -expand 0 -anchor e -padx 30
  pack $ap_btn -side left -in $w.f4 -expand 1 -anchor w


  #--Check widget states
  #
  ProcedureEntryPanel::checkMatcProc $id

  # Return the handle to ProcedureEntrypanel instance
  focus $w
  return $id 
}


proc ProcedureEntryPanel::delete {id} {
}


# Proc parses each inputline when the procedure stab file
# is read in
#
proc ProcedureEntryPanel::lineParser {id datalines } {

  set lines [split $datalines "\n"]

  set ln [set ProcedureEntryPanel::libraryName$id]
  set fn [set ProcedureEntryPanel::functionName$id]
  
	set arr [set ProcedureEntryPanel::globArray$id]
	set var [DataField::fieldNameSifToGui [set ProcedureEntryPanel::variableName$id]]
	set var_sz [lindex [lindex [DataField::getFieldSize Variables $var] 0] end]

  set new_lines ""

  foreach line $lines {
  
    # Skip stub "comment" lines (= line starting with "!#")
    #
    if { [regexp "^!#" $line] } {
      continue
    }
  
    # Replace:
    # patterns #library-name# with the acual library-name
    # patterns #function-name# with the acual function-name
    #
    regsub -all "#library-name#" $line $ln line
    regsub -all "#function-name#" $line $fn line

    if { "none" != [string tolower [set ProcedureEntryPanel::variableName$id]] } {
      regsub -all "#,n#" $line ",n" line
      regsub -all "#,X#" $line ",X" line
			if { $var_sz > 1 } {
				regsub -all "#,XX#" $line ",X($var_sz)" line
			} else {
				regsub -all "#,XX#" $line ",X" line
			}
    } else {
      regsub -all "#,n#" $line "" line
      regsub -all "#,X#" $line "" line
      regsub -all "#,XX#" $line "" line
    }

    append new_lines $line
    append new_lines "\n"
  }

  return [string trim $new_lines "\n"]
}


proc ProcedureEntryPanel::save {id} {
  global Info
  upvar #0 [set ProcedureEntryPanel::globArray$id] theArray


  # Store results
  # =============
  set result ""
  set sep $Info(dataListSeparator)

  set is_matc [set ProcedureEntryPanel::isMatcProc$id]
  
  if { $is_matc } {
    append result ""
    append result $sep
    append result [ProcedureEntryPanel::getMatcFunctionOutName $id]
    append result $sep
    append result [set ProcedureEntryPanel::variableName$id]
    append result $sep
    append result ""

  } else {
    append result [set ProcedureEntryPanel::libraryName$id]
    append result $sep
    append result [set ProcedureEntryPanel::functionName$id]
    append result $sep
    append result [set ProcedureEntryPanel::variableName$id]
    append result $sep
    append result [set ProcedureEntryPanel::sourceName$id]
  }

  set ProcedureEntryPanel::dataChanged$id 0
  set ProcedureEntryPanel::updated$id 1

  set ProcedureEntryPanel::result$id $result

  [set ProcedureEntryPanel::libraryNameEntry$id]    configure $Info(en_mod_option) $Info(en_nmod_color)
  [set ProcedureEntryPanel::functionNameEntry$id]   configure $Info(en_mod_option) $Info(en_nmod_color)
  [set ProcedureEntryPanel::sourceNameEntry$id]     configure $Info(en_mod_option) $Info(en_nmod_color)
  [set ProcedureEntryPanel::librariesEntry$id]      configure $Info(en_mod_option) $Info(en_nmod_color)
  [set ProcedureEntryPanel::variableNameButton$id]  configure $Info(om_mod_option) $Info(om_nmod_color)

  [set ProcedureEntryPanel::applyButton$id] configure -state disabled
  [set ProcedureEntryPanel::cancelButton$id] configure -state disabled

  set Info(procedureLibraries) [set ProcedureEntryPanel::libraries$id]
  set Info(lastProcedureDir) $theArray(currentProcedureDir) 


  # Call CallBack function
  # ======================
  set cb [set ProcedureEntryPanel::callback$id] 
  set ga [set ProcedureEntryPanel::globArray$id]
  set fld [set ProcedureEntryPanel::field$id]
  $cb $id $ga $fld

}


proc ProcedureEntryPanel::ok { id  {do_exit 1} } {
  global Info

  # If no changes
  if { ![set ProcedureEntryPanel::dataChanged$id] } {

    if {$do_exit} {
      ProcedureEntryPanel::cancel $id
    }

    return
  }


  # Check data
  if { ![ProcedureEntryPanel::checkData $id] } {
    return
  }

  ProcedureEntryPanel::save $id

  if { $do_exit } {
    Panel::cancel [set ProcedureEntryPanel::window$id]
  }
}


proc ProcedureEntryPanel::apply {id} {

  ProcedureEntryPanel::ok $id 0
}


proc ProcedureEntryPanel::cancel {id} {

  if { ![info exists ProcedureEntryPanel::window$id] } {
    return
  }

  set ProcedureEntryPanel::updated$id 0
  set ProcedureEntryPanel::result$id ""

  set ga [set ProcedureEntryPanel::globArray$id]
  set fld [set ProcedureEntryPanel::field$id]
  upvar #0 $ga theArray

  if { [info exist theArray($fld)] &&
       $theArray($fld) == "" &&
       [info exist theArray($fld,proc)]
     } {
    set theArray($fld,proc) 0
    Panel::procedureCheckBoxCommandProc $ga
  }

  Panel::cancel [set ProcedureEntryPanel::window$id]

  ProcedureEntryPanel::delete $id
}


# Check data, return possible error message
#
proc ProcedureEntryPanel::checkData {id} {
  global Info ModelProperty

  set ln [Util::stringTRQ [set ProcedureEntryPanel::libraryName$id]]
  set ln_abs [set ProcedureEntryPanel::libraryName$id,abs]

  set fn [Util::stringTRQ [set ProcedureEntryPanel::functionName$id]]

  set sn [Util::stringTRQ [set ProcedureEntryPanel::sourceName$id]]
  set sn_abs [set ProcedureEntryPanel::sourceName$id,abs]

  set is_matc [set ProcedureEntryPanel::isMatcProc$id]
  
  if { $is_matc } {
    set ln ""
    set needs [list 0 1]
  } else {
    set needs [list 1 1]
  }

  set msg ""

  foreach var  {ln fn} \
          type {Library Function} \
          name {Library Function} \
          need $needs {

    set val [set $var]

    if { $val == "" && $need } {
        set msg "ERROR: $name name missing!"
        Message::dispOkMessage $msg
        return 0
    }

    set check_proc "check"
    append check_proc $type
    append check_proc "Name"

    DataField::$check_proc $val "" msg
  
    if { $msg != "" } {
      Message::dispOkMessage $msg
      return 0
    }
  }

  # WIN32 .dll extension!
  # ======================
  if { $Info(platform) == "windows" } {
    if { $ln != "" } {
      if { ".dll" != [file extension $ln] } {
        append ln ".dll"
      }
    }
  }

  # Add model directory to the library-name and source paths if they are 
  # still relative
  set mdir $ModelProperty(MODEL_DIRECTORY,absolute)

  if { $ln != "" && "relative" == [file pathtype $ln] && $ln_abs } {
    set ln [ file join $mdir $ln ]
  }

  if { $sn != "" &&  "relative" == [file pathtype $sn] && $sn_abs} {
    set sn [ file join $mdir $sn ]
  }

  set ProcedureEntryPanel::libraryName$id $ln
  set ProcedureEntryPanel::sourceName$id $sn

  return 1
}


# Procedure set widget states when Matc-proc checkbox is turned on/off
#
proc ProcedureEntryPanel::checkMatcProc {id} {
  global Info

  if { [set ProcedureEntryPanel::isMatcProc$id] } {
    set state disabled
    set bg $Info(nonActiveBg)
  } else {
    set state normal
    set bg white
  }

  foreach wdg [set ProcedureEntryPanel::matcWidgetGroup$id] {

    if { [string equal -nocase "Entry" [winfo class $wdg]] } {
      Widget::configureSB $wdg $state $bg
    } else {
      Widget::configureS $wdg $state
    }
  }
}


# Procedure checks if the variablename entryfield in the Procedure definition
# panel should be updated
#
proc ProcedureEntryPanel::checkVariable {id} {
  global Info

  # Let option menu widget to update
  Util::doWait 50

  # Pick previous and current Variable name values
  #
  set var [set ProcedureEntryPanel::variableName$id]
  set var_prev [set ProcedureEntryPanel::variableNamePrev$id]

  # NOTE: We should not clear the variable entry field if "User defined"
  # is selected twice, otherwise we could delete laborously entered
  # data!!!
  #
  if { $var_prev != "" && $var == $var_prev } {
    return
  }

  set wdg [set ProcedureEntryPanel::variableNameButton$id]
  $wdg configure $Info(om_mod_option) $Info(om_mod_color)

  [set ProcedureEntryPanel::applyButton$id] configure -state normal
  [set ProcedureEntryPanel::cancelButton$id] configure -state normal
  
  set ProcedureEntryPanel::dataChanged$id 1
  set ProcedureEntryPanel::updated$id 1

  # Update previous variable name to be current name
  set ProcedureEntryPanel::variableNamePrev$id $var
}


proc ProcedureEntryPanel::compileSource {id} {
  global Info Model ModelProperty ProcessTable UserSetting

  if { ![ProcedureEntryPanel::checkData $id] } {
    return
  }

  set sn [set ProcedureEntryPanel::sourceName$id]
  set ln [set ProcedureEntryPanel::libraryName$id]
  set fn [set ProcedureEntryPanel::functionName$id]

  if { $sn ==  "" } {
    set msg [list "No source file name given!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return
  }
  
  # Test source existence
  set msg ""
  set rc [catch {set sn_ch [open $sn "r"] } msg]

  if { $rc != 0} {
    set msg [list "Cannot find source file!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return

  } else {
    close $sn_ch
  }

  if { $ln ==  "" } {
    set msg [list "No library name given!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return
  }

  if { $fn ==  "" } {
    set msg [list "No function name given!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return
  }

  set compiler [set ProcedureEntryPanel::compiler$id]
  set libs [set ProcedureEntryPanel::libraries$id]

  #--System specific call name specifier
  if { $Info(platform) == "windows" } {
    set sys WIN32
  } else {
    set sys UNIX
  }

  set use_eval 1
  set priority $Info(bgPriorityLevel)
  set add_to_process_table 1
  set show_in_process_list 1
  set browse_mode [string tolower $UserSetting(BROWSE_MODE_PROCEDURE_COMPILER)]
  set browse_delay $Info(browser,smallDelay)
  set process_name $Info(F90,processName)
  set process_tag $Info(F90,processTag)
  set logfile [file join $ModelProperty(LOG_DIRECTORY) $Info(F90,logfile)]
  set exec_name $Info(F90,proc,$sys)

  set EH $Info(ELMER_HOME)

  #-Win32
  if { $sys == "WIN32" } {
    #regsub -all "/" $compiler "\\" compiler
    regsub -all "/" $libs "\\" libs
    regsub -all "/" $sn "\\" sn
    regsub -all "/" $ln "\\" ln

    set exec_cmd $compiler
    set arg "$sn $libs /link /out:$ln"

  #-Unix
  } else {  
    set exec_cmd $compiler
    set arg "$ln $sn $libs"
  }

  #--Changedir to model's directory (store current)
  set Info(currentDirectory) [pwd]

  Util::cpp_exec setCurrentDirectory $ModelProperty(MODEL_DIRECTORY,absolute)

  set runnumber [MenuExec::do_elmer_exec \
    $process_name \
    $use_eval $exec_cmd $arg \
    $priority \
    $add_to_process_table \
    $show_in_process_list \
    $browse_mode \
    $exec_name \
    $process_tag \
    $logfile]

  # If call succesful
  # =================
  if { $runnumber != -1 } {

    set ProcessTable($runnumber,browserDelay) $Info(browser,smallDelay)

    if { $browse_mode == "logfile" } {
      MenuExec::browseProcessLog $runnumber
    }
  }

  #--Restore current directoy
  Util::cpp_exec setCurrentDirectory $Info(currentDirectory)

}


proc ProcedureEntryPanel::editSource {id} {
  global Info Model

  if { ![ProcedureEntryPanel::checkData $id] } {
    return
  }

  set sn [set ProcedureEntryPanel::sourceName$id]

  if { $sn ==  "" } {
    set msg [list "Give source file name!"]

    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return
  }

	set can_delete_cancelled 0
  
  # If not exist, copy stub to the file  
  if { ![file exists $sn] } {

		set can_delete_cancelled 1

    set stub [file join $Info(ELMER_FRONT_HOME) lib "procedure.stub.txt"]

    if { [file exists $stub] } {
      file copy $stub $sn
		
    } else {
      set msg [list "Cannot find stub file:\n\n" \
                    "$stub \n\n" \
                    "for the procedure source file!\n\n" \
                    "Press Ok to create the file and continue"]

      if { "cancel" == [Message::dispOkCancelMessage $msg] } {
        return
      }
      
      # Create the source file
      set ch_id [open $sn "w"]
      close $ch_id 

    }
  }

  set must_exist 0
  set can_edit 1
  set line_parser_cmd [list "ProcedureEntryPanel::lineParser $id"]
  set process_nbr 0
  set can_use_external 0

  set channel_id 0
  set textW [FileBrowser::browse .procedureEditor$id $sn \
                        $must_exist $can_edit $can_delete_cancelled \
												$line_parser_cmd $channel_id \
                        $process_nbr $can_use_external ]

  if { $textW == "" } {
    return
  }

  #$textW configure -font $Info(tableFont)
  set im $ProcedureEntryPanel::insertMarker
  set ii [$textW search -backwards $im end]

  if { $ii != "" } {
    $textW mark set insert $ii
  }

  focus $textW
  $textW configure -state normal

}


proc ProcedureEntryPanel::updated {id wdg} {
  global Info

  $wdg configure $Info(en_mod_option) $Info(en_mod_color)

  [set ProcedureEntryPanel::applyButton$id] configure -state normal
  [set ProcedureEntryPanel::cancelButton$id] configure -state normal

  set ProcedureEntryPanel::dataChanged$id 1
  set ProcedureEntryPanel::updated$id 1
}


proc ProcedureEntryPanel::getMatcFunctionInName {id} {

  set fn [Util::stringTRQ [set ProcedureEntryPanel::functionName$id]]

  set fn [string map {"$MATC" "" "(tx)" ""} $fn]

  set fn [Util::stringTRQ $fn]

  return $fn
}


proc ProcedureEntryPanel::getMatcFunctionOutName {id} {

  set fn "\$MATC "
  append fn [Util::stringTRQ [set ProcedureEntryPanel::functionName$id]]
  append fn "(tx)"

  return $fn
}


# end ecif_tk_procedureEntryPanel.tcl
# ********************
