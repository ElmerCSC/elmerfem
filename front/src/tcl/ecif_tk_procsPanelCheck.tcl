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
#Module:    ecif_tk_procsPanelCheck.tcl
#Language:  Tcl
#Date:      05.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Panel checking procedures.
#
#************************************************************************


##############################
### Field values checking ####
##############################


# Check one field
#
# NOTE: This is a "stand alone" function, which gives
# its own error msg boxes etc.
# Do NOT call it from PanelCheck::checkFields !!!
#
proc PanelCheck::checkField {globArray fld} {
  global Info
  upvar #0 $globArray theArray


  if { !$theArray($fld,act) } {
    return ""
  }

  set msg [PanelCheck::checkFieldValue $globArray $fld]

  #-Possible Error messages
  if {$msg != ""} {
    set label [DataField::getFieldProperty $globArray $fld Label]
    lappend totmsg [string toupper $label]
    lappend totmsg "\n\n$msg"

    set Info(messageIcon) error
    Message::dispOkMessage $totmsg  "$Info(FRONT_NAME) message!" $theArray(winName)

    set theArray($fld,err) 1
    return $msg
  }

  set theArray($fld,err) 0

  return ""
}


# Starts data checking in the panel referred by globArray.
# Each field should have a field-type indicator like InitialCondition(TEMP,"S")
# for: a scalar field in the TEMP-variable field in the init cond panel.
#NOTE: Checking of vectors should be more general (length as param. etc)
#
proc PanelCheck::checkFields {globArray {fields  ""} } {
  global Common Info
  upvar #0 $globArray theArray

  set totmsg ""

  if { $fields == "" } {
    set fields $theArray(allFields)
  }

  #-Check each field
  foreach fld $fields {

    if { !$theArray($fld,act) } {
      continue
    }

    set msg [PanelCheck::checkFieldValue $globArray $fld]

    #-Field error messages are collected!
    if { $msg != "" } {

      set theArray($fld,err) 1

      set label [DataField::getFieldProperty $globArray $fld Label]
      lappend totmsg "\n"
      lappend totmsg [string toupper "$label"]
      lappend totmsg "\n$msg\n"
    
    #-Field Ok
    } else {
      set theArray($fld,err) 0
    }
  }

  #-Ok (no msg-variable was created!)
  if { $totmsg == "" } {
    return 1

  #-Error message is displayed
  } else {
    set totmsg [linsert $totmsg 0 $Info(fieldMsg)]
    set Info(messageIcon) error
    Message::dispOkMessage $totmsg  \
                 "$Info(FRONT_NAME) message" \
                  $theArray(winName)

    return 0
  }
}


# Check data value for one field
#
proc PanelCheck::checkFieldValue {globArray fld} {
  global Common Info
  upvar #0 $globArray theArray

  # If not screen field, just run possible check proc for the field
  #
  if { ![Panel::isCheckableDataField $globArray $fld] } {
    set msg [PanelCheck::evalFieldCheckProc $globArray $fld]
    return $msg
  }

  set theArray($fld) [string trim $theArray($fld)]

  # Update data type if scalar value
  if { ( ![info exists theArray($fld,proc)] || !$theArray($fld,proc) ) &&
       ( ![info exists theArray($fld,table)] || !$theArray($fld,table) )
     } {
    set theArray($fld,valueData) $theArray($fld)
  }

  set theArray($fld,matc) ""
	
  # Find field properties
  # =====================
  set fvtype [DataField::getFieldProperty $globArray $fld FieldValueType]
  set frmt [DataField::getFieldProperty $globArray $fld FieldFormat]

  set min_max_sizes ""
  lappend min_max_sizes [DataField::getFieldProperty $globArray $fld MinDataSize1]
  lappend min_max_sizes [DataField::getFieldProperty $globArray $fld MaxDataSize1]
  lappend min_max_sizes [DataField::getFieldProperty $globArray $fld MinDataSize2]
  lappend min_max_sizes [DataField::getFieldProperty $globArray $fld MaxDataSize2]

  set needed_level [Panel::getFieldNeededLevel $globArray $fld]

  set range ""
  set valueset ""
  set tmp [DataField::getFieldProperty $globArray $fld Limits]

  if { "Interval" == [lindex $tmp 0] } {
    set range [lindex $tmp 1]
  } elseif { "Set" == [lindex $tmp 0] } {
    set valueset [lindex $tmp 1]
  }

  set accept_scalar 1
  set fdtype [DataField::getFieldProperty $globArray $fld FieldDataType]
  if { $fdtype == "Array" } {
    set accept_scalar 0
  }

  # Default value type is numeric
  set value_type "nbr"

	#-A Matc-script
  if { "$" == [string index $theArray($fld) 0] } {
    
    #-If second character is $, send field to Sif as it is
    #
		if { "$" == [string index [string trim [string range $theArray($fld) 1 end] ] 0] } {
			set msg ""
			return

    #-If next four non-blank characters are MATC, send field to Sif as it is
    #
    } elseif { "MATC" == [string toupper [string range [string trim [string range $theArray($fld) 1 end] ] 0 3 ] ] } {
      set mc_val [string trim [string range $theArray($fld) 1 end]]
      set mc_val [string trim [string range $mc_val 4 end]]
      if { "\"" == [string index $mc_val 0] } {
        set mc_val [string replace $mc_val 0 0 "\'"]
      }
      if { "\"" == [string index $mc_val end] } {
        set mc_val [string replace $mc_val end end "\'"]
      }
      
      # Format $MATC-field
      set theArray($fld) "\$MATC "
      append theArray($fld) $mc_val
			set msg ""
			return
    
    #-Otherwise send field to be evaluated by Front
    #
		} else {
			set theArray($fld,matc) $theArray($fld)
			set theArray($fld)  [string range $theArray($fld) 1 end]
			set value_type "matc"
		}

  #-File name
  } elseif { $fvtype == "File" ||
             ( [info exists theArray($fld,file)] &&
               $theArray($fld,file)
             )
           } {
    set value_type "file"

  #-Directory name
  } elseif { $fvtype == "Directory" } {
    set value_type "dir"

  #-String
  } elseif { $fvtype == "String" } {
    set value_type "str"

  #-Procedure name
  } elseif { $fvtype == "Procedure" ||
             ( [info exists theArray($fld,proc)] &&
               $theArray($fld,proc)
             )
           } {
    set value_type "proc"

  #-Library name
  } elseif { $fvtype == "Library" } {
    set value_type "lib"

  #-Function name
  } elseif { $fvtype == "Function" } {
    set value_type "func"

  #-Variable name
  } elseif { $fvtype == "Variable" } {
    set value_type "var"
  }

  # If field need special handling to read the data
  if { [info exists Common($globArray,$fld,widget2DataProc)] } {
    Util::execArrayProc $globArray $fld,widget2DataProc
  }

  # Start checking
  # ==============
  set msg ""
  set value $theArray($fld)

  # Pick possible variable and size description
  set description [DataField::extractVariableAndSizeDescription cleaned $value]
  set theArray($fld) $cleaned

  set check 1
  set dtype ""
  set msg [DataField::setVariableAndSizeValues $globArray $fld $dtype $description $check]

  # If error in description
  if { $msg != "" } {
    set theArray($fld,variables) $description
    set theArray($fld,dataSize) ""
    lappend totmsg "$msg"
  }

  # Possible variable values size in one entry
  if { [info exist theArray($fld,variables)] } {
    set var_size [llength [split $theArray($fld,variables) $Info(dataListSeparator)]]
  } else {
    set var_size 0
  }

  # Data item size in one entry
  # NOTE: Use this call, because it converts DIM-like defs into
  # numeric values
  set data_size [DataField::getFieldSize $globArray $fld]

  # Call check proc for the field
  # =============================
  set res_value ""
  set msg [DataField::checkValue res_value res_size $theArray($fld) $frmt $value_type    \
                      $needed_level $range $valueset $var_size $data_size $accept_scalar \
                      $min_max_sizes]

  set res_value [string trim $res_value]

  # If ok
  # =====
  if {$msg == ""} {

    set msg [PanelCheck::evalFieldCheckProc $globArray $fld]

    #--Finally data is Ok and can be stored to the field!!!
    #
    if { $msg == "" } {
      set theArray($fld) $res_value

			#set s1 [lindex [lindex $res_size 0] end]
			#set s2 [lindex [lindex $res_size 1] end]
      #set theArray($fld,dataSize) [list $s1 $s2]
			set theArray($fld,dataSize) [lindex [lindex $res_size 0] end]

      set theArray($fld,nofEntries) [lindex $res_size end]
    }
  }

  return $msg
}


# Run possible field-specific check-proc
#
proc PanelCheck::evalFieldCheckProc {globArray fld} {

  set pn $globArray
  append pn "($fld,checkProc)"

  if { "" != [namespace eval :: info procs $pn] } {
    set msg [eval $pn]
  }

  return ""
}


#-Procedure checks vector-type group field data
#
proc Common(VectorGroup,checkProc) {globArray var} {
  upvar #0 $globArray theArray
  global Common Info

  set result ""
  set glist $Common($globArray,$var,group)
  set frmt [DataField::getFieldProperty $globArray $var FieldFormat]
  set ndl [DataField::getFieldProperty $globArray $var NeededLevel]

  set procedure 0
  if { [info exist theArray($var,proc)] } {
      set procedure $theArray($var,proc)
  }

  #--Data is for a procedure name 
  #  It is put to the first field(index 0), like BoundaryCondition(VELOCITY,VELX)
  if {$procedure} {
    set fn_fld [lindex $glist 0]

    set msg [DataField::checkValue res_val res_size $theArray($fn_fld) "" "proc" $ndl]

    if {$msg == ""} {
      set theArray($fn_fld) $res_val
      set theArray($var,proc) 1
    }

    set result $msg

  #--All velocity fields are read as numbers
  #  NOTE: currently fields are supposed to be scalar-valued
  } else {
    foreach v $glist {
      set ndl [DataField::getFieldProperty $globArray $var NeededLevel]
      set msg [DataField::checkValue res_val res_size $theArray($v) $frmt "nbr" $ndl]
      if {$msg == ""} {
        set theArray($v) $res_val
      } else {
        lappend result [append [append $v ": "] $msg]
      }
    }
  }

  return $result
}



# Runs field specific check procs (for cross checking, widget states etc.)
# which are called when a field is filled with data

# NOTE: These are not tied to any button, but they are "panel wide"
# Check procs are called in the format:
# BoundaryCondition(RADIATION,fillProc)
# These proc are mostly in StdPanelCheck-file
#
proc PanelCheck::execPanelFillProcs {globArray} {
  global Common
  upvar #0 $globArray theArray

  foreach fld $theArray(allFields) {  

    # Only active fields for the case must be
    # checked!
    if { ![info exists theArray($fld,trg)] || !$theArray($fld,trg) } {
      continue
    }

    PanelCheck::execFillProc $globArray $fld
  }
}


# Call check procs for group fields
#
proc PanelCheck::execPanelGroupEditProcs {globArray} {
  global Common  $globArray
  upvar #0 $globArray theArray

  foreach var $theArray(allFields) {  

    # Only problem fields are interesting
    if { ![info exists theArray($fld,prb)] || !$theArray($fld,prb) } {
      continue
    }

    if {"" != [array names $globArray "$var,groupEditProc"] } {

      set proc_name $globArray
      append proc_name "($var,editProc)"
      $proc_name $var $var
    }
  }
}


#####################################
### Panel widget check procedures ###
#####################################


# ====================================#
# Generic widget status check procs   #
#=====================================#

# A generic proc to  respond to option-menu selection
# NOTE: Define "check_proc_call" method if you want to do some special
# processing for your option menu
#
proc PanelCheck::checkOptionMenu { om wdg globArray fld check_proc_call event_info } {
  global Info
  upvar #0 $globArray theArray

  # If mouse is just moved, but not yet released!
  if { $event_info >= 256 } {
    return 0
  }

  # Let option menu widget to update
  Util::doWait 50

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,old)] } {

    if { $theArray($fld) != $theArray($fld,old) } {

      set theArray($fld,mod) 1
      Panel::panelDataModified 1 $globArray

    } else {
      set theArray($fld,mod) 0
      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }
    }

    Widget::setOptionMenuStatus $globArray $fld $wdg

  }

  eval "$check_proc_call"

  return 1
}


# A generic proc to  respond to a selection-box selection
# NOTE: Define "check_proc_call" method if you want to do some special
# processing for your selection box
# event_mode: mouse action like Control-Button1
# local_x local_y: local widget coordinates for selection
#
proc PanelCheck::checkSelectionBox { wdg globArray fld check_proc_call event_mode {lb_x ""} {lb_y ""} } {
  global Info
  upvar #0 $globArray theArray

  set arguments {[list $event_mode $lb_x $lb_y]}

  eval "$check_proc_call $arguments"
}


# A generic proc to  respond to checkbox selection
#
proc PanelCheck::checkCheckBoxFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,old)] } {

    set cur_val $theArray($fld)
    set old_val $theArray($fld,old)

    if {$cur_val == "" } {
      set cur_val "False"
    }

    if {$old_val == "" } {
      set old_val "False"
    }

    if { $cur_val != $old_val } {
      set theArray($fld,mod) 1
      Panel::panelDataModified 1 $globArray

    } else {
      set theArray($fld,mod) 0

      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }
    }

    Widget::setCheckBoxStatus $globArray $fld $theArray(allWidgets,$fld)
  }
}


# A generic proc to  respond to checkbutton (yes/no button) selection
#
proc PanelCheck::checkCheckButtonFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray

  set cur_value $theArray($fld)

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,old)] } {

    if { $cur_value != $theArray($fld,old) } {
      set theArray($fld,mod) 1
      Panel::panelDataModified 1 $globArray

    } else {
      set theArray($fld,mod) 0

      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }

    }
  }
}


# A generic proc to  respond to radibutton selection
#
proc PanelCheck::checkRadioButtonFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray

  set cur_value $theArray($fld)

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,old)] } {

    if { $cur_value != $theArray($fld,old) } {
      set theArray($fld,mod) 1
      Panel::panelDataModified 1 $globArray

    } else {
      set theArray($fld,mod) 0

      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }

    }
  }
}


# A generic proc to  respond to Abs checkbox selection
#
proc PanelCheck::checkAbsFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray
}


# A generic proc to  respond to Use checkbox selection
#
proc PanelCheck::checkUseFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray

  set cur_value $theArray($fld,Use)

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,Use,old)] } {

    if { $cur_value != $theArray($fld,Use,old) } {
      set theArray($fld,Use,mod) 1

    } else {
      set theArray($fld,Use,mod) 0

      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }
    }
  }

  Widget::setFieldStatus $globArray $fld,Use "CheckBox"
}


# A generic proc to  respond to File checkbox selection
#
proc PanelCheck::checkFileFieldStatus { globArray fld } {
  global Info
  upvar #0 $globArray theArray

  set cur_value $theArray($fld,file)

  # Compare current value to "panel in" value for
  # modification marking
  if { [info exists theArray($fld,file,old)] } {

    if { $cur_value != $theArray($fld,file,old) } {
      set theArray($fld,file,mod) 1

    } else {
      set theArray($fld,file,mod) 0

      if { ![Panel::panelHasModifiedStatus $globArray] } {
        Panel::panelDataModified 0 $globArray
      }
    }
  }

  Widget::setFieldStatus $globArray $fld,file "CheckBox"
}


#===============================#
# Generic widget type fill procs
#===============================#

proc Common(CheckBox,fillProc) {arg} {

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  if { [info exists theArray(targetMask)] } {
    DataField::checkActivitySlaves $globArray $fld $theArray(targetMask)
  } else {
    DataField::checkActivitySlaves $globArray $fld ""
  }
}


proc Common(BrowsableFile,fillProc) {arg} {

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  if {1 == [DataField::getFieldProperty $globArray $fld HasUseCheckBox "" 0]} {
    if { $theArray($fld,act) && 
         $theArray($fld) != "" 
       } {
      set theArray($fld,Use) 1
    } else {
      set theArray($fld,Use) 0
    }
  }

  #-Set BrwosableFile-widget state
  Common(BrowsableFile,editProc) $arg
}


proc Common(IncludeFile,fillProc) {arg} {

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  if { $theArray(INCLUDE,act) && 
       $theArray(INCLUDE) != "" 
     } {
    set theArray(INCLUDE,Use) 1
  } else {
    set theArray(INCLUDE,Use) 0
  }

  #-Set INCLUDE-widget state
  Common(IncludeFile,editProc) $arg
}


proc Common(OptionMenu,fillProc) {arg} {

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  if { [info exists theArray(targetMask)] } {
    DataField::checkActivitySlaves $globArray $fld $theArray(targetMask)

  } else {
    DataField::checkActivitySlaves $globArray $fld ""
  }
}



#===============================#
# Generic widget type edit procs
#===============================#
#
# General proc activated from widget
#

proc Common(BrowsableFile,editProc) {arg} {
  global Info

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  # Clean trailing quotes and white space
  set theArray($fld) [string trim $theArray($fld) " \"\t\n"]

  # If we have a use-checkbox which controls entry fields state
  #
  if {1 == [DataField::getFieldProperty $globArray $fld HasUseCheckBox "" 0]} {

    if { !$theArray($fld,Use) } {
      set theArray($fld,act) 0
    } else {
      set theArray($fld,act) 1
    }    

    Widget::configureField $globArray $fld
  }
 
  return
}


proc Common(IncludeFile,editProc) {arg} {
  global Info

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]
  set key [lindex $arg 2]

  upvar #0 $globArray theArray

  # Clean trailing quotes and white space
  set theArray($fld) [string trim $theArray($fld) " \"\t\n"]

  # If we have a use-checkbox which controls entry fields state
  #
  if { [info exists theArray($fld,Use)] } {

    if { !$theArray($fld,Use) } {
      set theArray($fld,act) 0
    } else {
      set theArray($fld,act) 1
    }    

    Widget::configureField $globArray $fld
  }

  if { $key == "Double-1" } {
    set fn $theArray($fld)

    if { $fn != "" && [file isfile $fn] } {
      FileBrowser::browse .winEdit$globArray$fld $fn 1 1 
    }
  }

  return
}


proc Common(OptionMenu,editProc) {arg} {

  set globArray [lindex $arg 0]
  set fld [lindex $arg 1]

  upvar #0 $globArray theArray

  if { [info exists theArray(targetMask)] } {
    DataField::checkActivitySlaves $globArray $fld $theArray(targetMask)

  } else {
    DataField::checkActivitySlaves $globArray $fld ""
  }
}




# ========================================================== # 
# Check procs attached to the panel widgets at creation time #
# ========================================================== # 


proc PanelCheck::execFillProc {globArray fld } {

  # NOTE: Call first generic fill proc for the field type
  #
  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  Util::execArrayProc Common "$wtype,fillProc" [list $globArray $fld]

  # Then call possible array-field specific fill proc
  #
  if { 1 == [DataField::getFieldProperty $globArray $fld IsUserDefinedField "" 0] } {
    Util::execArrayProc $globArray "UserDefinedField,fillProc" $fld

  } else {
    Util::execArrayProc $globArray "$fld,fillProc"
  }

}


proc PanelCheck::execEditProc {globArray fld {arg ""}} {

  # NOTE: Call first generic edit proc for the field type
  #
  set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  Util::execArrayProc Common "$wtype,editProc" [list $globArray $fld $arg]

  # Then call possible array-field specific edit proc
  #
  if { 1 == [DataField::getFieldProperty $globArray $fld IsUserDefinedField "" 0] } {
    Util::execArrayProc $globArray "UserDefinedField,editProc" $fld

  } else {
    Util::execArrayProc $globArray "$fld,editProc" $arg
  }

}


# =====================
# Widget specific calls
# =====================

proc PanelCheck::singleButton {globArray fld } {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::groupButton {globArray group fld } {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleCheckButton {globArray fld } {
  PanelCheck::checkCheckButtonFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::groupCheckButton {globArray group fld } {
  PanelCheck::checkCheckButtonFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleCheckBox {globArray fld} {
  upvar #0 $globArray theArray

  if { [info exists theArray(targetMask)] } {
    DataField::checkActivitySlaves $globArray $fld $theArray(targetMask)
  } else {
    DataField::checkActivitySlaves $globArray $fld ""
  }

  PanelCheck::checkCheckBoxFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleListBox {globArray fld arg} {
  PanelCheck::execEditProc $globArray $fld $arg
}


proc PanelCheck::singleOptionMenu {globArray fld } {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleSelectionBox {globArray fld arg} {
  PanelCheck::execEditProc $globArray $fld $arg
}


proc PanelCheck::singleText {globArray fld } {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::groupCheckBox {globArray group fld} {
  PanelCheck::checkCheckBoxFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleRadioButton {globArray fld } {
  PanelCheck::checkRadioButtonFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::groupRadioButton {globArray group fld} {
  PanelCheck::checkRadioButtonFieldStatus $globArray $group
  if { ![Util::execArrayProc $globArray "$group,groupEditProc" $fld] } {
    PanelCheck::singleRadioButton $globArray $fld
  }
}


proc PanelCheck::singleFileButton {globArray fld} {
  PanelCheck::checkFileFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleProcedureButton {globArray fld} {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleBrowsableDirectory {globArray fld} {
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleBrowsableFile {globArray fld} {
  PanelCheck::execEditProc $globArray $fld
}

proc PanelCheck::singleAbsButton {globArray fld} {
  PanelCheck::checkAbsFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::singleIncludeFile {globArray fld arg} {
  PanelCheck::execEditProc $globArray $fld $arg
}

proc PanelCheck::singleUseButton {globArray fld} {
  PanelCheck::checkUseFieldStatus $globArray $fld
  PanelCheck::execEditProc $globArray $fld
}


proc PanelCheck::groupProcedureButton {globArray fld} {
  global Common 
  upvar #0 $globArray theArray

  PanelCheck::checkCheckBoxFieldStatus $globArray $fld

  set panel $theArray(parameterType)

  foreach FV [lrange $Common($panel,$fld,group) 1 end] {
    set m $entry_widget
    append m $FV
    if {1 == $theArray($fld,file)} {
      Widget::configureS $m disabled
    } else {
      Widget::configureS $m normal
    }
  }
}


# end ecif_tk_panelCheck.tcl
# ********************


