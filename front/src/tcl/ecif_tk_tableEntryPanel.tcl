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
#Module:    ecif_tk_tableEntryPanel.tcl
#Language:  Tcl
#Date:      06.11.99
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Displays a panel with for a field's table entry
#
#************************************************************************

#=====================#
#  TableEntryPanel    # 
#=====================#

# Procedure displays a data entry dialog for array data.
#
proc TableEntryPanel::create { globArray fld callback {sorted ""} {accept_multi_x ""} } {
  global Info Equation
  upvar #0 $globArray theArray

  # Increment (unique) id
  variable id
  incr id

  variable sortMarker

  # Namespace variables with id
  #
  variable acceptMultiX$id
  variable acceptMultiXInputArg$id
  variable applyButton$id
  variable callback$id
  variable cancelButton$id
  variable dataChanged$id
  variable dataSize1$id
  variable dataSize1Button$id
  variable dataSize1Menu$id
  variable dataSize1Prev$id
  variable dataSize2$id
  variable dataSize2Button$id
  variable dataSize2Menu$id
  variable dataSize2Prev$id
  variable dataType$id
  variable errorIndex$id -1
  variable errorMsg$id ""
  variable field$id
  variable fieldName$id
  variable globArray$id
  variable nofEntries$id
  variable result$id
  variable sort$id
  variable sortInputArg$id
  variable tableEntryId$id
  variable updated$id 0
  variable variableName$id
  variable variableNamePrev$id
  variable variableNameMenu$id
  variable variableNameButton$id
  variable window$id

  set globArray$id $globArray
  set callback$id $callback
  set field$id $fld
  set dataChanged$id 0

  set sortInputArg$id $sorted
  set acceptMultiXInputArg$id $accept_multi_x

  set fname [DataField::getFieldProperty $globArray $fld SifName]

  if { $fname == "" } {
    set fname [DataField::fieldNameGuiToSif $fld]
  }

  set coord_index [DataField::getFieldProperty $globArray $fld CoordinateIndex]
  if { $coord_index > 0 } {
    ;#append fname "_$coord_index"
  }

  set fieldName$id $fname

  # Init nofEntries
  if { [info exist theArray($fld,nofEntries)] && $theArray($fld,nofEntries) != "" } {
    set nofEntries$id $theArray($fld,nofEntries)
  } else {
    set nofEntries$id 0
  }

  # Data type: Scalar Array
  set dataType$id [DataField::getFieldProperty $globArray $fld FieldDataType]

  # Init size1 and size2
  if { [info exist theArray($fld,tableData,dataSize)] &&
       $theArray($fld,tableData,dataSize) != "" } {
    set size $theArray($fld,tableData,dataSize)
  } else {
    set tmp [DataField::getFieldSize $globArray $fld]
    lappend size [lindex [lindex $tmp 0] 0] ;#First size1
    lappend size [lindex [lindex $tmp 1] 0] ;#First size2
  }

  set dataSize1$id [lindex $size 0]
  set dataSize2$id [lindex $size 1]

  set dataSize1Prev$id [set dataSize1$id]
  set dataSize2Prev$id [set dataSize2$id]

  # Init variableName
  if { [info exist theArray($fld,variables)] && $theArray($fld,variables) != "" } {
    set variableName$id $theArray($fld,variables)
  } else {
    set variableName$id "none"
  }

  set variableNamePrev$id [set variableName$id]

  if { $sorted == "" } {
    if { [set variableName$id] != "none" } {
      set sorted 1
    } else {
      set sorted 0
    }
  }

  if { $accept_multi_x == "" } {

    if { [set variableName$id] != "none" } {
      set accept_multi_x 0
    } else {
      set accept_multi_x 1
    }
  }

  set sort$id $sorted
  set acceptMultiX$id $accept_multi_x

  # NOTE: Only one window per panel/field can exists
  set window$id .wTEP$globArray$fld

  set w [set window$id]
  set TableEntry($id,winName) $w

  # Set window title
  set fn [string toupper [set TableEntryPanel::fieldName$id]]
  set title "Table entry for $fn"
  set TableEntry($id,winTitle) $title

  set wid [winfo atom $w]
  set TableEntry($id,winId) $wid

  set Info(thisWindow) $w

  if { 1 == [Util::checkPanelWindow TableEntry $wid $TableEntry($id,winTitle) ""] } {
    return
  }  

  toplevel $w
  wm title $w $TableEntry($id,winTitle)

  frame $w.a
  frame $w.b
  frame $w.c

  # Principal frames
  frame $w.f1
  frame $w.f2
  frame $w.f3

  #---Entry fields stuff
  
  #-frames for the labels and entries
  frame $w.f11 ;# Variable name field
  frame $w.f12 ;# Size1 field
  frame $w.f13 ;# Size2 field
  
  #-labels and entries
  label $w.lVariable -text " Variable: "
  label $w.lSize1    -text " Size 1: "
  label $w.lSize2    -text " Size 2: "

  set values [DataField::getFieldDataAsList $globArray $fld]
  set sort [set TableEntryPanel::sort$id]
  set accept_multi_x [set TableEntryPanel::acceptMultiX$id]

  set tableEntryId$id [TableEntry::create $id \
                                          $w.f2 \
                                          $values \
                                          $sort $sortMarker \
                                          $accept_multi_x \
                                          55 15 \
                                          TableEntryPanel::checkEntry \
                                          TableEntryPanel::preCheckDataList \
                                          TableEntryPanel::postCheckDataList \
                                          $globArray]


  set omOptions "-indicatoron 0 -relief raised -activebackground $Info(optionMenuBg)"

  # Variable name option menu 
  #
  if { [info exists theArray(targetMask)] } {
    set targetMask $theArray(targetMask)
  } else {
    set targetMask $Equation(problemMask)
  }

	set accept_array 0
  set val_list [lindex [Panel::getCurrentVariables $targetMask $accept_array] 1]
  set val_list [linsert $val_list 0 "none"]
  set om [Panel::createOptionMenuWidget $w.mVariable TableEntryPanel::variableName$id $val_list]

  set TableEntryPanel::variableNameMenu$id $om
  set TableEntryPanel::variableNameButton$id $w.mVariable

  eval "$w.mVariable configure -width 15 $omOptions"
  bind $om <<MenuSelect>> "TableEntryPanel::checkVariable $id ;TableEntryPanel::postCheckDataList $id"

  # If field is simple array (no argument variables)
  #
  if { ![DataField::getFieldProperty $globArray $fld Variabled] } {
    $w.mVariable configure -state disabled
  }

  set size [DataField::getFieldSize $globArray $fld]

  # Size1 option menu 
  #
  set val_list [lindex $size 0]
  
  if { $val_list == "N" } {
    set dataSize1$id "N"
  }
  
  set om [Panel::createOptionMenuWidget $w.mSize1 TableEntryPanel::dataSize1$id $val_list]

  set TableEntryPanel::dataSize1Menu$id $om
  set TableEntryPanel::dataSize1Button$id $w.mSize1

  eval "$w.mSize1 configure -width 3 $omOptions"
  bind $om <<MenuSelect>> "TableEntryPanel::checkDataSize $id 1 {%s}"

  # Size2 option menu 
  #
  set val_list [lindex $size 1]

  if { $val_list == "N" } {
    set dataSize2$id "N"
  }

  set om [Panel::createOptionMenuWidget $w.mSize2 TableEntryPanel::dataSize2$id $val_list]

  set TableEntryPanel::dataSize2Menu$id $om
  set TableEntryPanel::dataSize2Button$id $w.mSize2

  eval "$w.mSize2 configure -width 3 $omOptions"
  bind $om <<MenuSelect>> "TableEntryPanel::checkDataSize $id 2 {%s}"


  #---Buttons
  #

  set ap $Info(defaultApplyState)
  set ca $Info(defaultCancelState)

  set ok_btn [button $w.ok -text OK -command "TableEntryPanel::ok $id"]
  set cn_btn [button $w.cancel -text Cancel -command "TableEntryPanel::cancel $id" -state $ca]
  set ap_btn [button $w.apply -text Apply -command "TableEntryPanel::apply $id" -state $ap]

  set applyButton$id $ap_btn
  set cancelButton$id $cn_btn

  #---Pack widgets
  #

  pack $w.a $w.b $w.c -side top -fill y -fill x -expand 1 -anchor c -pady 5

  pack $w.f1 -side top -in $w.a -fill y -fill x -expand 1 -anchor c -pady 5
  pack $w.f2 -side top -in $w.b -fill x -expand 1 -anchor c -pady 5
  pack $w.f3 -side top -in $w.c -fill y -fill x -expand 1 -anchor c -pady 5
 
  #pack $textW -side top 

  # Variable name, size frames
  pack $w.f11 $w.f12 $w.f13 -side left -in $w.f1 -anchor w 

  pack $w.lVariable $w.mVariable -side left -in $w.f11 -expand 1 -anchor w
  pack $w.lSize1 $w.mSize1 -side left -in $w.f12 -anchor w
  pack $w.lSize2 $w.mSize2 -side left -in $w.f13 -anchor w

  pack $ok_btn -side left -in $w.f3 -expand 1 -anchor e
  pack $cn_btn -side left -in $w.f3 -expand 0 -anchor e -padx 10 
  pack $ap_btn -side left -in $w.f3 -expand 1 -anchor w


  # Return the handle to TableEntryPanel instance
  focus $w
  return $id 
}


proc TableEntryPanel::delete {id} {
}


proc TableEntryPanel::checkVariable {id} {
  global Info

  # Let option menu widget to update
  Util::doWait 50

  set var_cur [set TableEntryPanel::variableName$id]
  set var_prev [set TableEntryPanel::variableNamePrev$id]

  if { $var_cur == $var_prev } {
    return
  }

  [set TableEntryPanel::applyButton$id] configure -state normal
  [set TableEntryPanel::cancelButton$id] configure -state normal

  # Mark widget modified
  set wdg [set TableEntryPanel::variableNameButton$id]
  $wdg configure $Info(om_mod_option) $Info(om_mod_color)
  
  # Mark table modified
  set TableEntryPanel::dataChanged$id 1

  # Set sort variable in the table entry widget
  #
  set tid [set TableEntryPanel::tableEntryId$id]

  #-If some argument variable, always sorted
  if { $var_cur != "none" } {
    set TableEntry::sort$tid 1

  #-Otherwise, use the original input value
  } else {
    set sort [set TableEntryPanel::sortInputArg$id]

    if { $sort == "" } {
      set TableEntry::sort$tid 0
    } else {
      set TableEntry::sort$tid $sort
    }
  }
}


# Check data consistency when one of the sizes
# has been changed. This is needed for Scalar/Matrix
# type data, where both sizes normally must be the same
#
proc TableEntryPanel::checkDataSize { id size_index {event_info ""} } {
  global Info

  # If mouse is just moved, but not yet released!
  if { $event_info >= 256 } {
    return
  }

  # Let option menu widget to update
  Util::doWait 50

  set sz_cur [set TableEntryPanel::dataSize$size_index$id]

  set prev_var TableEntryPanel::dataSize$size_index

  append prev_var Prev$id
  set sz_prev [set $prev_var]

  # Mark panel modified
  if { $sz_cur != $sz_prev } {

    [set TableEntryPanel::applyButton$id] configure -state normal
    [set TableEntryPanel::cancelButton$id] configure -state normal

    set wdg_name TableEntryPanel::dataSize$size_index
    append wdg_name "Button$id"
    set wdg [set $wdg_name]
    $wdg configure $Info(om_mod_option) $Info(om_mod_color)
  }

  TableEntryPanel::postCheckDataList $id
}


# Check that data list in the listbox matches current
# variable and size definitions already before a new entry is
# being added/deleted/inserted.
#
proc TableEntryPanel::preCheckDataList { id {mode ""} } {
  global Info

  set TableEntryPanel::dataChanged$id 1

  [set TableEntryPanel::applyButton$id] configure -state normal
  [set TableEntryPanel::cancelButton$id] configure -state normal

  # Only these modes need checking
  if { $mode != "add"    &&
       $mode != "insert" &&
       $mode != "delete" && 
       $mode != "update"
     } {
    return 1
  }

  set tid [set TableEntryPanel::tableEntryId$id]
  set data_entry [[set TableEntry::entry$tid] get]
  set data_list [TableEntry::getDataList $tid]
  set row_index [[set TableEntry::listbox$tid] curselection]

  set accept_multi_x [set TableEntryPanel::acceptMultiX$id]

  #--Check possible duplicates
  if { !$accept_multi_x && $mode != "delete" } {

    set entry_x [lindex $data_entry 0]
    set index 1
    set rindex 0

    foreach row $data_list {

      if { "" == [TableEntryPanel::checkEntry $id "" $row $index] } {
        set TableEntryPanel::errorIndex$id $index
        break
      }
    
      set row_x [lindex $row 0]

      # If multiple rows found
      #
      if { !$accept_multi_x &&
           $row_x == $entry_x &&
           !($mode == "update" && $rindex == $row_index)
         } {

        set msg [list "NOTE: Duplicate x-value with the entry at index $index ($entry_x)\n\n" \
                      $Info(noMoreNotifying) ]  

        set Info(messageIcon) question
        set result [Message::dispCancelYesNoMessage $msg]

        if { $result == "cancel" } {  
          return 0
        }

        if { $result == "yes" } {
          set accept_multi_x 1
          set TableEntryPanel::acceptMultiX$id  1
        }
        
        # Checking done for this time anyway!
        break
      }

      incr index
      incr rindex
    }
  }

  #--Check datalist size
  #
  set arr [set TableEntryPanel::globArray$id]
  set fld [set TableEntryPanel::field$id]

  set min_sz [DataField::getFieldProperty $arr $fld MinDataSize1]
  set max_sz [DataField::getFieldProperty $arr $fld MaxDataSize1]

  set sz1 [llength [TableEntry::getDataList $tid]]

  if { $mode == "delete" && $min_sz != "" && $sz1 == $min_sz } {
    set msg "NOTE: Number of rows must be >= $min_sz\n"
    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return 0
  }

  if { ($mode == "add" || $mode == "insert") && $max_sz != "" && $sz1 == $max_sz } {
    set msg "NOTE: Number of rows must be <= $max_sz\n"
    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return 0
  }

  return 1
}


# Check that data list in the listbox matches current
# variable and size definitions when a new (checked) entry
# is being added/deleted etc.
#
proc TableEntryPanel::postCheckDataList { id {mode ""} } {
  global Info

  # If rows deleted and no existing error, no
  # reason to check datalist
  if { $mode == "delete" &&
       [set TableEntryPanel::errorIndex$id] == -1
     } {
    return 1
  }

  # Let option menu widge to update
  Util::doWait 50

  set vn       [set TableEntryPanel::variableName$id]
  set vn_prev  [set TableEntryPanel::variableNamePrev$id]

  set sz1      [set TableEntryPanel::dataSize1$id]
  set sz1_prev [set TableEntryPanel::dataSize1Prev$id]

  set sz2      [set TableEntryPanel::dataSize2$id]
  set sz2_prev [set TableEntryPanel::dataSize2Prev$id]

  # If no changes in data criteria and no
  # delete operation, datalist is not rechecked
  # although there were a previous error.
  # NOTE: possible previous error is anyway checked
  # when Ok is pressed
  # NOTE: This prevents ugly checking "loops" with optionMenu control
  #
  if { $vn   == $vn_prev  &&
       $sz1  == $sz1_prev && 
       $sz2  == $sz2_prev && 
       $mode != "delete"  &&
       $mode != "update"
     } {
    return 1
  }

  set TableEntryPanel::variableNamePrev$id $vn
  set TableEntryPanel::dataSize1Prev$id $sz1
  set TableEntryPanel::dataSize2Prev$id $sz2

  set TableEntryPanel::errorIndex$id -1
  set TableEntryPanel::errorMsg$id ""

  return 1
}


# Check entry data, return evaluated data if ok, otherwise blank
#
# NOTE: Matc-expression are accepted in the entry, but they should
# entered quoted if they contain blanks:
# 1000 "$ 3 * 5"  
# "$1000 + 1" "$ 4 * 5" etc.
# $y=22.2  defines a new matc-variable value!!!
#
proc TableEntryPanel::checkEntry {id mode data {index -1} } {
  global Info

  set data [Util::stringTrim $data]
  set msg ""

  set result ""
  
  # Evaluate each entry component add build result string
  #  
  foreach val $data {
    
    # If the item is a Matc-expression
    #
    if { "\$" == [string index [string trim $val] 0] } {
      
      set is_var_def 0

      set to_do [string range [string trim $val] 1 end]
      
      # Check if the user is defining a matc-variable (ex. $y=...)
      # Verify action
      #
      set tmp [string map {" " ""} $to_do]
      set idx [string first "=" $tmp]
      if { $idx != -1 && "=" != [string index $tmp [expr 1 + $idx]] } {
        set msg "Are you sure to define a MATC variable value?"
        if { ![Message::verifyAction $Info(advancedUser) $msg] } {
          return ""
        }
      }

      set val [string trim [Util::doMatc $to_do]]

      if { "" == [string trim $val] } {
        return ""
      }
    }

    if { [catch {set val [expr $val]}] } {
      set msg "Invalid entry!"

      set Info(messageIcon) error
      Message::dispOkMessage $msg

      return ""
    }

    append result "$val "
  }    

  set tid [set TableEntryPanel::tableEntryId$id]
  set len [llength $data]

  set sz1 [set TableEntryPanel::dataSize1$id]
  set sz2 [set TableEntryPanel::dataSize2$id]

  if { $sz2 == 1 } {
    set sz 1
  } else {
    set sz $sz2
  }

  set vn  [set TableEntryPanel::variableName$id]
  
  # Variable's contribution to entry size
  if { $vn != "none" } {
    set vz 1
  } else {
    set vz 0
  }

  # If size is free
  set size1_free 0
  set size2_free 0

  if { $sz1 == "N" } {
    set size1_free 1
  }

  if { $sz2 == "N" } {
    set size2_free 1
  }

  if { $sz1 == "N" } {
    set asz1 1
  } else {
    set asz1 $sz1
  }

  if { $sz2 == "N" } {
    set asz2 1
  } else {
    set asz2 $sz2
  }

  set nof_rows [set TableEntry::nofEntries$tid]
  set nof_values [expr $nof_rows * $sz2]

  # If add or update mode
  # =====================
  #
  if { $mode == "add" || $mode == "insert" || $mode == "update" } {

    # No variable
    # -----------
    if { $vz == 0 } {

      set upd_sz 0
      
      # Scalar or vector (size2=1)
      #
      if { $sz2 == 1 && !$size1_free } {
        
        set upd_sz 1

        # Only scalar accepted
        #
        if { $sz1 == 1 && $len > 1 } {
          set msg "Only scalar value accepted (no variable defined)!"
        
        # And vector (array) accepted
        #
        } elseif { [expr $len + $nof_rows - $upd_sz] > $sz1 } {
          set msg "Max number of values is $sz1 !"
        }

      # Table
      #
      } else {

        if { !$size1_free && $nof_rows == $sz1 && $mode != "update"} {
          set msg "Number of data rows cannot exceed $sz1 !"
        }

        if { !$size2_free && $len != $sz2} {
          set msg "Incorrect data entry: given ($len) values,  expecting ($sz2) !"
        }
      }
      
    # Variable
    # --------
    } else {

      if { !$size1_free && !$size1_free } {
  
        set max_len [expr 1 + $sz1 * $sz2]

        # Each entry is a variable + one scalar
        #
        if { $sz1 == 1 && $sz2 == 1 && $len > 2 } {
          set msg "Only a variable value plus one scalar value is accepted!"
        
        # Each entry is a variable + an array or table 
        #
        } elseif { $len != $max_len } {
          set msg "Number of values should be $max_len (1+Size1*Size2)!"
        }

      }
    }
  }

  # If ERROR
  # ========
  #
  if { $msg != "" } {
    set fn [string toupper [set TableEntryPanel::fieldName$id]]

    # We are checking entry field data
    if { $index == -1 } {
      set msg [list "Data entry for $fn:\n\n" $msg]

    # We are checking listbox rows
    } else {
      incr index
      set TableEntryPanel::errorMsg$id $msg
      set msg [list "Data entry for $fn:\nData in row ($index)\n\n" $msg]
    }

    set Info(messageIcon) error
    Message::dispOkMessage $msg

    return ""
  }

  return [string trim $result]
}


proc TableEntryPanel::save {id} {
  global Info

  #-Create formatted result variable
  #
  set tid [set TableEntryPanel::tableEntryId$id]
  set data_list [TableEntry::getDataList $tid]

  set TableEntryPanel::result$id ""

  if { "none" == [set TableEntryPanel::variableName$id] &&
       1 == [set TableEntryPanel::dataSize2$id] 
     } {
    set sep " "
  } else {
    set sep "$Info(arrayDataSeparator)"
  }

  set counter 0
  foreach row $data_list {

    if { $counter > 0 } {
      append TableEntryPanel::result$id $sep
    }

    append TableEntryPanel::result$id $row
    incr counter
  }

  set TableEntryPanel::nofEntries$id $counter

  [set TableEntryPanel::applyButton$id] configure -state disabled
  [set TableEntryPanel::cancelButton$id] configure -state disabled

  [set TableEntryPanel::variableNameButton$id] configure $Info(om_mod_option) $Info(om_nmod_color)
  [set TableEntryPanel::dataSize1Button$id]  configure $Info(om_mod_option) $Info(om_nmod_color)
  [set TableEntryPanel::dataSize2Button$id]  configure $Info(om_mod_option) $Info(om_nmod_color)

  set TableEntryPanel::dataChanged$id 0
  set TableEntryPanel::updated$id 1

  # CallBack call
  set cb [set TableEntryPanel::callback$id] 
  set ga [set TableEntryPanel::globArray$id]
  set fld [set TableEntryPanel::field$id]
  set call [list $cb $ga $fld]
  $cb $id $ga $fld

} 


proc TableEntryPanel::ok {id {do_exit 1 } } {
  global Info

  # If no changes
  if { ![set TableEntryPanel::dataChanged$id] } {

    if {$do_exit} {
      TableEntryPanel::cancel $id
    }

    return
  }

  # If error
  #
  set index [set TableEntryPanel::errorIndex$id]
  set fn [string toupper [set TableEntryPanel::fieldName$id]]

  if { $index != -1 } {
    set msg "Data entry for $fn:\nERROR in data!\nData in row ($index)\n\n"
    set msg [list $msg $Info(anywayOk)]

    set Info(messageIcon) error

    # Continue with table
    if { "cancel" == [Message::dispOkCancelMessage $msg] } {
      return

    # Close table, but do not update result!
    } else {
     TableEntryPanel::cancel $id
     return
    }
  }

  # Check table size
  #
  set tid [set TableEntryPanel::tableEntryId$id]
  set data_list [TableEntry::getDataList $tid]

  set sz1 [set TableEntryPanel::dataSize1$id]
  set sz2 [set TableEntryPanel::dataSize2$id]

  set nof [llength $data_list]
  set vn [set TableEntryPanel::variableName$id]

  if { $vn == "none" } {
    if { $sz1 != "N" && $sz1 != $nof } {
      set msg "Number of rows does not match $sz1 (Size1) !!\n\n"
   
      set Info(messageIcon) error

      Message::dispOkMessage $msg
      return
    }
  
  } else {
    set row_sz [expr (1 + $sz1 * $sz2)]

    if { [llength [join $data_list]] != [expr $nof * $row_sz] } {
      set msg "All data rows must contain same number of values ($row_sz)!\n\n"
      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return      
    }
  }

  # Save data
  TableEntryPanel::save $id

  # Close window
  if { $do_exit } {
    set tid [set TableEntryPanel::tableEntryId$id]
    TableEntry::delete $tid
    Panel::cancel [set TableEntryPanel::window$id]
  }
} 


proc TableEntryPanel::apply {id} {

  TableEntryPanel::ok $id 0
}


proc TableEntryPanel::cancel {id} {

  if { ![info exists TableEntryPanel::window$id] } {
    return
  }

  set tid [set TableEntryPanel::tableEntryId$id]
  TableEntry::delete $tid

  set ga [set TableEntryPanel::globArray$id]
  set fld [set TableEntryPanel::field$id]
  upvar #0 $ga theArray

  if { [info exist theArray($fld)] &&
       $theArray($fld) == "" &&
       [info exist theArray($fld,table)]
     } {
    set theArray($fld,table) 0
    Panel::tableCheckBoxCommandProc $ga
  }
    
  set TableEntryPanel::updated$id 0
  #set TableEntryPanel::result$id "!cancel!"
  set TableEntryPanel::result$id ""

  Panel::cancel [set TableEntryPanel::window$id]
  
  TableEntryPanel::delete $id
}


# end ecif_tk_tableEntryPanel.tcl
# ********************
