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
#Module:    ecif_tk_procsWidget.tcl
#Language:  Tcl
#Date:      04.02.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Generic event handlers and procs for combibed widgets
#
#************************************************************************

#=========#
# Generic #
#=========#


# Check if a widget exist
#
proc Widget::exists { wdg } {
  upvar $wdg theWdg

  if { ![info exist $wdg] ||
       ![winfo exist $theWdg] 
     } {
    return 0
  } else {
    return 1
  }
}


# Check if an array widget exist
#
proc Widget::arrayWidgetExists { globArray fld } {
  upvar #0 $globArray theArray

  if { ![info exists theArray(allWidgets,$fld)] } {
    return 0
  }

  # First in the (possible) list must exists and also enough!
  #
  if { ![winfo exists [lindex $theArray(allWidgets,$fld) 0] ] } {
    return 0
  } else {
    return 1
  }
}

  
# Configure one widget and one option
#
proc Widget::configure1 { wdg option value } {

  set msg ""
  catch { $wdg configure $option $value} msg

  if { $msg != "" } {
    Message::showMessage "Illegal widget configure option ($option) in Widget::configure1 because: $msg"
  }
}


# Configure many widgets and one option
#
proc Widget::configureN { widgets option value } {

  set wdgs [join $widgets]

  foreach wdg $wdgs {

    set msg ""
    catch { $wdg configure $option $value} msg

    if { $msg != "" } {
      Message::showMessage "Illegal widget configure option ($option) in Widget::configure1 because: $msg"
      break
    }
  }
}


# Set widgets' state normal
#
proc Widget::activate {widgets} {
  set wdgs [join $widgets]
  foreach wdg $wdgs {
    set msg ""
    catch { $wdg configure -state normal} msg
    if { $msg != "" } {
      Message::showMessage "Illegal widget configure option (state) in Widget::deactivate because: $msg"
      break
    }
  }
}


# Set widgets's state disabled
#
proc Widget::deactivate {widgets} {
  set wdgs [join $widgets]
  foreach wdg $wdgs {
    set msg ""
    catch { $wdg configure -state disabled} msg
    if { $msg != "" } {
      Message::showMessage "Illegal widget configure option (state) in Widget::deactivate because: $msg"
      break
    }
  }
}


# Configure widgets State
#
proc Widget::configureS { wdg state } {

  set msg ""

  catch { $wdg configure -state $state} msg

  if { $msg != "" } {
    Message::showMessage "Illegal widget configure option in Widget::configureS because: $msg"
  }
}


# Configure widgets State, Background
#
proc Widget::configureSB { wdg state bg } {

  set msg ""

  catch { $wdg configure -state $state -bg $bg } msg

  if { $msg != "" } {
    Message::showMessage "Illegal widget configure option in Widget::configureSB because: $msg"
  }
}


# Configure widgets State, Background and Relief
#
proc Widget::configureSBR { wdg state bg relief } {

  set msg ""

  catch { $wdg configure -state $state -bg $bg -relief $relief } msg

  if { $msg != "" } {
    Message::showMessage "Illegal widget configure option in Widget::configureSBR because: $msg"
  }
}


# Configure widgets State, Background, Foreground and Relief
#
proc Widget::configureSBFR { wdg state bg fg relief } {

  set msg ""

  catch { $wdg configure -state $state -bg $bg -fg $fg -relief $relief } msg

  if { $msg != "" } {
    Message::showMessage "Illegal widget configure option in Widget::configureSBFR because: $msg"
  }
}


# ====================
# ARRAY field handling
# ====================

# Configure: (state, background etc.)
# ===================================

# Find widget for the field "fld" in the array "globArray"
# configure ist state (disabled, background etc)
#
proc Widget::configureField {globArray fld {wtype ""} } {
  upvar #0 $globArray theArray

  # Widget does not exist
  # =====================
  if { ![Widget::arrayWidgetExists $globArray $fld] } {
    return "Widget for the field=$fld does not exist!"
  }

  if { $theArray($fld,act) } {
    set mode "normal"
  } else {
    set mode "disabled"
  }

  if { $wtype == "" } {
    set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  }
  
  set wdg $theArray(allWidgets,$fld)

  # Call proper handler by widget type
  # If field not identified, return 0
  #
  switch $wtype {
    BrowsableDirectory    { Widget::configureBrowsableDirectory $globArray $fld $wdg $mode}
    BrowsableFile         { Widget::configureBrowsableFile $globArray $fld $wdg $mode}
    Entry                 { Widget::configureEntry $wdg $mode}
    Button                { Widget::configureButton $wdg $mode}
    IncludeFile           { Widget::configureIncludeFile $globArray $fld $wdg $mode}
    ListBox               { Widget::configureListBox $wdg $mode}
    OptionMenu            { Widget::configureOptionMenu $wdg $mode}
    CheckBox              { Widget::configureCheckBox $wdg $mode}
    RadioButton           { Widget::configureRadioButton $wdg $mode}
    RadioButtonOptionMenu { Widget::configureRadioButtonOptionMenu $wdg $mode}
    SelectionBox          { Widget::configureSelectionBox $wdg $mode}
    Text                  { Widget::configureText $wdg $mode}
    default               { return 0}
  }

  # Ok, field found and configured
  # ==============================
  return ""
}


# Configures one button widget according to
# the "mode" argument
#
proc Widget::configureButton {wdg mode} {
  global Info

  #--Set button foreground
  #-Normal mode
  if {$mode == "Normal" || $mode == "normal"} {
    Widget::configure1 $wdg -activeforeground black
  #-Disabled mode
  } elseif {$mode == "Disabled" || $mode == "disabled"} {
    Widget::configure1 $wdg -disabledforeground $Info(nonActiveFg)
  }

  Widget::configureS $wdg $mode

}


# Configures one browsable directory widget state according to
# the "mode" argument
#
proc Widget::configureBrowsableDirectory {globArray fld wdg mode} {
  global Info
  upvar #0 $globArray theArray

  # Directory entry
  # ----------
  Widget::configureEntry $wdg $mode

  # Browse button
  # -------------
  if { [info exists theArray(allWidgets,$fld,Browse)] } {
    Widget::configureButton $theArray(allWidgets,$fld,Browse) $mode
  }

  # Abs checkbox
  # ------------
  if { [info exists theArray(allWidgets,$fld,Abs)] } {
    Widget::configureButton $theArray(allWidgets,$fld,Abs) $mode
  }
}


# Configures one browsable file widget state according to
# the "mode" argument
#
proc Widget::configureBrowsableFile {globArray fld wdg mode} {
  global Info
  upvar #0 $globArray theArray

  # File entry
  # ----------
  Widget::configureEntry $wdg $mode

  # Browse button
  # -------------
  if { [info exists theArray(allWidgets,$fld,Browse)] } {
    Widget::configureButton $theArray(allWidgets,$fld,Browse) $mode
  }

  # Abs checkbox
  # ------------
  if { [info exists theArray(allWidgets,$fld,Abs)] } {
    Widget::configureButton $theArray(allWidgets,$fld,Abs) $mode
  }
}


# Configures one checkbox widget according to
# the "mode" argument
#
proc Widget::configureCheckBox {wdg mode} {
  global Info

  #-Normal mode
  if {$mode == "Normal" || $mode == "normal"} {
    Widget::configure1 $wdg -state normal
    Widget::configure1 $wdg -fg black

  #-Disabled mode
  } elseif {$mode == "Disabled" || $mode == "disabled"} {
    Widget::configure1 $wdg -state disabled
    Widget::configure1 $wdg -fg gray
    #Widget::configure1 $wdg -fg $Info(nonActiveFg)
  }
}


# Configures one entry widget state according to
# the "mode" argument
#
proc Widget::configureEntry {wdg mode} {
  global Info

  #-Normal mode
  if {$mode == "Normal" || $mode == "normal"} {
    Widget::configureSBFR $wdg normal white black groove 

  #-Disabled mode
  } elseif {$mode == "Disabled" || $mode == "disabled" } {
    Widget::configureSBFR $wdg disabled $Info(nonActiveBg) $Info(nonActiveFg) groove 
  }
}


# Configures one include file widget state according to
# the "mode" argument
#
proc Widget::configureIncludeFile {globArray fld wdg mode} {
  global Info
  upvar #0 $globArray theArray

  # NOTE: use-checkbox is not controilled here, that
  # would block the whole widget!!!
  # Use $fld,Use instead (like INCLUDE,Use)
  
  # Browse button
  # -------------
  if { [info exists theArray(allWidgets,$fld,Browse)] } {
    Widget::configureButton $theArray(allWidgets,$fld,Browse) $mode
  }

  # File entry
  # ----------
  #
  # Use-checkbox controls entry state
  if { [info exists theArray(allWidgets,$fld,Use)] } {

    # If normal mode and checkbox is clicked --> entry to normal state
    if { $mode == "normal" && $theArray($fld,Use) && $theArray($fld,Use,act) } {
      Widget::configureEntry $wdg "normal"

    # Otherwise entry to disabled state
    } else {
      Widget::configureEntry $wdg "disabled"
    }

  # No use-checkbox, set entry state by "mode"
  } else {
    Widget::configureEntry $wdg $mode
  }
}


# Configures one listBox widget according to
# the "mode" argument
#
proc Widget::configureListBox {wdg mode} {
  global Info

  # Listbox does not have a -state option!
  #Widget::configureS $wdg $mode
}


# Configures one optionMenu widget according to
# the "mode" argument
#
proc Widget::configureOptionMenu {wdg mode} {
  global Info
  Widget::configureS $wdg $mode
}


# Configures one radiobutton widget according to
# the "mode" argument
#
proc Widget::configureRadioButton {wdg mode} {
  global Info
  Widget::configureS $wdg $mode
}


# Configures one radioButtonOptionMenu widget according to
# the "mode" argument
# NOTE: Only OptionMenu part configured
#
proc Widget::configureRadioButtonOptionMenu {wdg mode} {
  global Info

  set rb_wdg [lindex $wdg 0]
  set om_wdg [lindex $wdg 1]

  Widget::configureOptionMenu $om_wdg $mode
}


# Configures one selectionBox widget according to
# the "mode" argument
#
# NOTE: A listbox does not have a -state option!
#
proc Widget::configureSelectionBox {wdg mode} {
  global Info

  #-Normal mode
  if {$mode == "Normal" || $mode == "normal"} {
    Widget::configure1 $wdg -fg black

  #-Disabled mode
  } elseif {$mode == "Disabled" || $mode == "disabled" } {
    Widget::configure1 $wdg -fg $Info(nonActiveFg)
  }
}


# Configures one text widget according to
# the "mode" argument
#
proc Widget::configureText {wdg mode} {
  global Info

  #-Normal mode
  if {$mode == "Normal" || $mode == "normal"} {
    Widget::configureSBFR $wdg normal white black groove 

  #-Disabled mode
  } elseif {$mode == "Disabled" || $mode == "disabled" } {
    Widget::configureSBFR $wdg disabled $Info(nonActiveBg) $Info(nonActiveFg) groove 
  }
}


# Field status: modified, error color etc.
# ========================================

# Find widget for the field "fld" in the array "globArray"
# and set its status (error/mofified etc.)
#
proc Widget::setFieldStatus {globArray fld {wtype ""} } {
  upvar #0 $globArray theArray

  # Widget does not exist
  # =====================
  if { ![Widget::arrayWidgetExists $globArray $fld] } {
    return -1
  }

  if { $wtype == "" } {
    set wtype [DataField::getFieldProperty $globArray $fld WidgetType]
  }
  
  set wdg $theArray(allWidgets,$fld)

  # Call proper handler by widget type
  # Return 0 if widget not identified
  #
  switch $wtype {
    BrowsableDirectory    { Widget::setEntryStatus $globArray $fld $wdg }
    BrowsableFile         { Widget::setEntryStatus $globArray $fld $wdg }
    Entry                 { Widget::setEntryStatus $globArray $fld $wdg }
    IncludeFile           { Widget::setEntryStatus $globArray $fld $wdg }
    OptionMenu            { Widget::setOptionMenuStatus $globArray $fld $wdg }
    RadioButtonOptionMenu { Widget::setRadioButtonOptionMenuStatus $globArray $fld $wdg }
    CheckBox              { Widget::setCheckBoxStatus $globArray $fld $wdg }
    RadioButton           { Widget::setRadioButtonStatus $globArray $fld $wdg }
    default               { return 0 }
  }

  # Ok, field found anf handled
  # ===========================
  return 1
}


# Entry widget status
#
proc Widget::setEntryStatus {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  if { [catch {
    #---Error status
    if { $theArray($fld,err) } {
      Widget::configure1 $wdg $Info(en_err_option) $Info(en_err_color) 

    #---Modified status
    } elseif { $theArray($fld,mod) } {
      Widget::configure1 $wdg $Info(en_mod_option) $Info(en_mod_color) 

    #---Normal status
    } else {
      Widget::configure1 $wdg $Info(en_err_option) $Info(en_nerr_color)
      Widget::configure1 $wdg $Info(en_mod_option) $Info(en_nmod_color)
    }
  } msg ] } {
    tk_messageBox -message $msg
  }
  
}


# CheckBox widget status
#
proc Widget::setCheckBoxStatus {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  if { [catch {
    #---Error status
    if { $theArray($fld,err) } {
      Widget::configure1 $wdg $Info(cb_err_option) $Info(cb_err_color) 

    #---Modified status
    } elseif { $theArray($fld,mod) } {
      Widget::configure1 $wdg $Info(cb_mod_option) $Info(cb_mod_color) 

    #---Normal status
    } else {
      Widget::configure1 $wdg $Info(cb_err_option) $Info(cb_nerr_color)
      Widget::configure1 $wdg $Info(cb_mod_option) $Info(cb_nmod_color)
    }
  } msg ] } {
    tk_messageBox -message $msg
  }
}


# OptionMenu widget status
#
proc Widget::setOptionMenuStatus {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  if { [catch {
    #---Error status
    if { $theArray($fld,err) } {
      Widget::configure1 $wdg $Info(om_err_option) $Info(om_err_color) 

    #---Modified status
    } elseif { $theArray($fld,mod) } {
      Widget::configure1 $wdg $Info(om_mod_option) $Info(om_mod_color) 

    #---Normal status
    } else {
      Widget::configure1 $wdg $Info(om_err_option) $Info(om_nerr_color)
      Widget::configure1 $wdg $Info(om_mod_option) $Info(om_nmod_color)
    }
  } msg ] } {
    tk_messageBox -message $msg
  }
  
}


# RadioButton widget status
#
proc Widget::setRadioButtonStatus {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  if { [catch {
    #---Error status
    if { $theArray($fld,err) } {
      Widget::configure1 $wdg $Info(rb_err_option) $Info(rb_err_color) 

    #---Modified status
    } elseif { $theArray($fld,mod) } {
      Widget::configure1 $wdg $Info(rb_mod_option) $Info(rb_mod_color) 

    #---Normal status
    } else {
      Widget::configure1 $wdg $Info(rb_err_option) $Info(rb_nerr_color)
      Widget::configure1 $wdg $Info(rb_mod_option) $Info(rb_nmod_color)
    }
  } msg ] } {
    tk_messageBox -message $msg
  }
  
}


# OptionMenu widget status
#
proc Widget::setRadioButtonOptionMenuStatus {globArray fld wdg} {
  global Info
  upvar #0 $globArray theArray

  set rb_wdg [lindex $wdg 0]
  set om_wdg [lindex $wdg 1]

  if { [catch {
    #---Error status
    if { $theArray($fld,err) } {
      Widget::configure1 $om_wdg $Info(om_err_option) $Info(om_err_color) 

    #---Modified status
    } elseif { $theArray($fld,mod) } {
      Widget::configure1 $om_wdg $Info(om_mod_option) $Info(om_mod_color) 

    #---Normal status
    } else {
      Widget::configure1 $om_wdg $Info(om_err_option) $Info(om_nerr_color)
      Widget::configure1 $om_wdg $Info(om_mod_option) $Info(om_nmod_color)
    }
  } msg ] } {
    tk_messageBox -message $msg
  }
}


################
### BINDINGS ###
################


proc Widget::setPanelBindings {globArray win} {
  upvar #0 $globArray theArray

  set events {
    Activate Deactivate
    FocusIn FocusOut
  }

  # Bind events to procs (named according to the event symbol)
  foreach event $events {
    bind $win <$event> "+Widget::panel$event $globArray"
  }
}


proc Widget::setLabelBindings {globArray {fields ""} } {
  upvar #0 $globArray theArray

  if { $fields == "" && [info exists theArray(allFields)] } {
    set fields $theArray(allFields)
  }

  foreach fld $fields {

    if { [info exists theArray(allWidgets,label,$fld)] &&
         [winfo exists $theArray(allWidgets,label,$fld)]
       } {
      set wdg $theArray(allWidgets,label,$fld)
      
      set lbl [string trim [$wdg cget -text]]

      # If field has a visible label, set bindings (for help)
      if { $lbl != "" } {
        bind $wdg <Button-3> "+Widget::genericButton-3 $globArray $fld $wdg"
      }
    }
  }
}



proc Widget::setEntryBindings { binding_grp globArray fld wdg } {

  # Select events by binding-group
  #
  switch $binding_grp {

    standard {
      set events {
        Button-1 Double-1
        Button-3 KeyPress-F1
        KeyPress-Return KeyPress-Down KeyPress-Up
        KeyPress-Escape Control-i Control-m
        KeyPress KeyRelease 
        FocusIn FocusOut
        Enter Leave 
      }
    }

    non_standard {
      set events {
        Button-3 KeyPress-F1
        KeyRelease
        KeyPress-Return KeyPress-Down KeyPress-Up KeyPress-Escape
      }
    }

    non_standard_single {
      set events {
        Button-3 KeyPress-F1
        KeyRelease
        KeyPress-Escape
      }
    }

    default {
      set events ""
      Message::showMessage "Unknown event group: $binding_grp"
    }

  }

  # Bind events to procs (named according to the event symbol)
  foreach event $events {
    bind $wdg <$event> "+Widget::entryEvent $event $globArray $fld $wdg {%A %K}"
  }
}


proc Widget::setCheckBoxBindings { binding_grp globArray fld wdg } {

  # Select events by binding-group
  #
  switch $binding_grp {

    standard {
      set events {
        Button-1 KeyPress-space KeyPress-Return
      }
    }

    non_standard {
      set events {
        Button-1 KeyPress-space KeyPress-Return
      }
    }


  }

  # Bind events to procs (named according to event symbol)
  foreach event $events {
    bind $wdg <$event> "+Widget::checkBox_generic $globArray $fld $wdg {%A %K}"
  }
}


proc Widget::setRadioButtonBindings { binding_grp globArray fld wdg } {

  # Select events by binding-group
  #
  switch $binding_grp {

    standard {
      set events {
        Button-1 KeyPress-space KeyPress-Return
      }
    }

    non_standard {
      set events {
        Button-1 KeyPress-space KeyPress-Return
      }
    }


  }

  # Bind events to procs (named according to event symbol)
  foreach event $events {
    bind $wdg <$event> "+Widget::radioButton_generic $globArray $fld $wdg {%A %K}"
  }
}



#########################
### EVENT classifiers ###
#########################

proc Widget::isEditingEvent {event_info} {

  set ascii [lindex $event_info 0]
  set symb  [lindex $event_info 1]

#-If uncommented, one way to see keycodes!
#print "ascii=$ascii symb=$symb"

  # If not an ascii-key or Delete, BackSpace, Ctrl+v (paste)
  if { $ascii == "" &&
       $symb != "Delete" &&
       $symb != "BackSpace" &&
       $symb != "v"
     } {
    return 0
  }
   
  # Tab or Enter are not editing events
  switch $symb {
    Return -
    Tab {
      return 0
    }
  }

  # Otherwise editing
  return 1
}


proc Widget::bindEditingEvent {action event_info} {

  if { [Widget::isEditingEvent $event_info] } {
    eval $action
  }
}


################### 
### EVENT procs ###
###################


# =================
# Panel event procs
# =================

proc Widget::panelFocusIn {globArray} {
  upvar #0 $globArray theArray

  set theArray(isActive) 1
}


proc Widget::panelFocusOut {globArray} {
  upvar #0 $globArray theArray

  set theArray(isActive) 0
}

proc Widget::panelFocusIn {globArray} {
  upvar #0 $globArray theArray

  set theArray(hasFocus) 1
}


proc Widget::panelFocusOut {globArray} {
  upvar #0 $globArray theArray

  set theArray(hasFocus) 0
}


# ===================
# Generic event procs
# ===================

# Help for a field
#
proc Widget::genericButton-3 { globArray fld wdg  {event_info ""} } {

  set parent [winfo parent $wdg]

  Message::helpMessage $globArray $fld $parent
}


# =================
# Entry event procs
# =================

proc Widget::entryEvent {event_name globArray fld wdg {event_info ""} } {

  set msg ""
  if { [catch {Widget::entry$event_name $globArray $fld $wdg $event_info } msg] } {
    Message::showMessage "$msg fld=$fld event=$event_name"
  }
}


# KeyPress in an entry-widget
#
proc Widget::entryKeyPress { globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray
  
  if { $event_info == "" } {
    return
  }

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }
  
  # A procedure or table-mode entry cannot be edited
  # (or of course can be, but only via Proc and Table widgets!)
  #
  if { $wdg != "" &&
       ( ([info exist theArray($fld,proc)] && $theArray($fld,proc)) ||
         ([info exist theArray($fld,table)] && $theArray($fld,table))
       )
     } {
    Widget::configureS $wdg disabled
  }

}


# KeyRelease in an entry-widget
#
proc Widget::entryKeyRelease { globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray
  
  if { $event_info == "" } {
    return
  }

  if { $event_info != "" &&
       ![Widget::isEditingEvent $event_info]
     } {
    return
  }

  # Show in normal model, although procedure or table-mode entry cannot
  # be edited
  #
  if { $wdg != "" &&
       ( ([info exist theArray($fld,proc)] && $theArray($fld,proc)) ||
         ([info exist theArray($fld,table)] && $theArray($fld,table))
       )
     } {
    Widget::configureS $wdg normal
    return
  }

  if { $theArray($fld,err) } {
    return
  }

  set theArray($fld,mod) 1

  Widget::setEntryStatus $globArray $fld $wdg
}


# ButtonClick event for an entry-widget
#
# NOTE: This is a trick to activate Pro/Table buttons
# when they are inactivated (eg. by selecting a body or parameter)
# and cursor is still in the entry whose FocusIn activated them!
#
proc Widget::entryButton-1 { globArray fld wdg {event_info ""} } {
  Widget::entryFocusIn $globArray $fld $wdg $event_info 
}


# Help for a field
# ----------------
#
proc Widget::entryButton-3 {globArray fld wdg  {event_info ""} } {

  set parent [winfo parent $wdg]
  Message::helpMessage $globArray $fld $parent
}

proc Widget::entryKeyPress-F1 {globArray fld wdg {event_info ""} } {

  set parent [winfo parent $wdg]
  Message::helpMessage $globArray $fld $parent
}


# Ctrl-i: Set initial value for the entry field
#
proc Widget::entryControl-i {globArray fld wdg {event_info ""} } {
  global Common Solver
  upvar #0 $globArray theArray

  set system_index ""

  if { $globArray == $Common(panelArray,SL) &&
       [info exist Solver(equationIndex)]
     } {
    set system_index $Solver(equationIndex)
  }

  set iv [DataField::getInitialValue $globArray $fld $system_index]

  if { $iv != "" } {
    set theArray($fld) $iv
  }

}

# Ctrl-m: Toggle with a possibl Matc-script and field value
#
proc Widget::entryControl-m {globArray fld wdg {event_info ""} } {
  global Solver
  upvar #0 $globArray theArray

	if { ![info exists theArray($fld,matc)] || $theArray($fld,matc) == "" } {
		return
	}
	
	# Swap values
	set tmp $theArray($fld)
	set theArray($fld) $theArray($fld,matc)
	set theArray($fld,matc) $tmp
}


# DoubleClick event for an entry-widget
#
proc Widget::entryDouble-1 { globArray fld wdg {event_info ""} } {
  global Info

  #-If entry is not active, do not open any
  # edit box
  if { "disabled" == [$wdg cget -state ] } {
    return
  }

  upvar #0 $globArray theArray

  #-Procedure editor
  if { [info exist theArray($fld,proc)]  &&
       $theArray($fld,proc)
     } {
    Panel::editProcedure $globArray $fld $wdg

  #-Table editor
  } elseif { [info exist theArray($fld,table)]  &&
       $theArray($fld,table)
     } {
    Panel::editTable $globArray $fld $wdg

  #-Normal (scalar) field, do nothing!
  } else {
    return
  }

}


# Enter (= Return) key pressed in an entry-widget
#
proc Widget::entryKeyPress-Return { globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  if { "" != [PanelCheck::checkField $globArray $fld] } {
    Widget::configure1 $wdg $Info(en_err_option) $Info(en_err_color) 
    return
  }

  Widget::configure1 $wdg $Info(en_err_option) $Info(en_nerr_color) 

  if { $theArray($fld) != $theArray($fld,prev) } {
    Widget::configure1 $wdg $Info(en_mod_option) $Info(en_mod_color) 
  } else {
    Widget::configure1 $wdg $Info(en_mod_option) $Info(en_nmod_color)
  } 

}


# Escape (Esc) key pressed in an entry-widget
#
proc Widget::entryKeyPress-Escape { globArray fld wdg {event_info ""} } {
  upvar #0 $globArray theArray

  if { ![info exists theArray($fld,prev)] } {
    return
  }

  set theArray($fld) $theArray($fld,prev)  
  set theArray($fld,mod) 0
  Widget::setEntryStatus $globArray $fld $wdg

  if { ![Panel::panelHasModifiedStatus $globArray] } {
    Panel::panelDataModified 0 $globArray
  }
}



# Arrow Down key pressed in an entry-widget
#
proc Widget::entryKeyPress-Down { globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  event generate $wdg <KeyPress-Tab> 
}


# Arrow Up key pressed in an entry-widget
#
proc Widget::entryKeyPress-Up { globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  event generate $wdg <Shift-KeyPress-Tab> 
}



# FocusIn event for an entry-widget
#
proc Widget::entryFocusIn {globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  if { $event_info != ""  && $wdg != ""} {

    #--Arrow cursor to mark a non-editable entry for proc/table modes
    if { ([info exist theArray($fld,proc)] && $theArray($fld,proc)) ||
         ([info exist theArray($fld,table)] && $theArray($fld,table))
       } {
     Widget::configure1 $wdg -cursor $Info(arrowCursor)

    #--Otherwise normal insert cursor
    } else {
     Widget::configure1 $wdg -cursor $Info(insertCursor)
    }
  }

  set theArray(currentEntryField) $fld

  # Activate proc and Table buttons, if needed
  #
  if { $theArray(hasProcButton) } {

    if { [DataField::getFieldProperty $globArray $fld Procedurable] } {
      set theArray(PROCEDURE_BUTTON,act) 1
      Widget::configureField $globArray PROCEDURE_BUTTON

      set f_var $globArray; append f_var "($fld,proc)"
      Widget::configure1 $theArray(allWidgets,PROCEDURE_BUTTON) -variable $f_var

    } else {
      set theArray(PROCEDURE_BUTTON,act) 0
      Widget::configureField $globArray PROCEDURE_BUTTON
    }
  }

  if { $theArray(hasTableButton) } {

    if { [DataField::getFieldProperty $globArray $fld Tableable] } {
      set theArray(TABLE_BUTTON,act) 1
      Widget::configureField $globArray TABLE_BUTTON

      set f_var $globArray; append f_var "($fld,table)"
      Widget::configure1 $theArray(allWidgets,TABLE_BUTTON) -variable $f_var

    } else {
      set theArray(TABLE_BUTTON,act) 0
      Widget::configureField $globArray TABLE_BUTTON
    }
  }
}


# FocusOut event for an entry-widget
#
proc Widget::entryFocusOut {globArray fld wdg {event_info ""} } {
  global Info
  upvar #0 $globArray theArray

  # Nothing currently!
  # ==================
  return


  if { [info exist theArray(hasFocus)] && !$theArray(hasFocus) } {
    return
  }

  if { [info exist theArray(isActive)] && !$theArray(isActive) } {
    return
  }
  
  if { $theArray($fld,err) } {
    return
  }

  set msg [PanelCheck::checkField $globArray $fld]

  if { $msg != "" } {
    Widget::configure1 $wdg $Info(en_err_option) $Info(en_err_color) 
    return
  }

  Widget::configure1 $wdg $Info(en_err_option) $Info(en_nerr_color) 

  #if { $theArray($fld) != $theArray($fld,prev) } {
  #  Widget::configure1 $wdg $Info(en_mod_option) $Info(en_mod_color) 
  #}
}


# Enter (into the widget) event for an entry-widget
#
proc Widget::entryEnter {globArray fld wdg {event_info ""} } {
  upvar #0 $globArray theArray

  if { ![info exists theArray(entryInfoLabel)] } {
    return
  }

  if { ![info exist theArray($fld,dataSize)] } {
    return
  }
  
  set vr ""
  set ne ""
  set sz ""
  set blanks "                    "

  if { [info exist theArray($fld,variables)] &&
    $theArray($fld,variables) != ""
  } {
    set vr $theArray($fld,variables)
    regsub -all "_" $vr " " vr
    append vr $blanks
    set vr [string range $vr 0 19]
  }

  if { [info exist theArray($fld,nofEntries)] &&
       $theArray($fld,nofEntries) != "" &&
       $theArray($fld,nofEntries) > 1
     } {
    set ne $theArray($fld,nofEntries)
    append ne $blanks
    set ne [string range $ne 0 19]
  }

  set size $theArray($fld,dataSize)
#MSG "sixe=$size"
  if { $size != "" &&
       1 != [expr [join $size *]]
     } {
    set sz $size
    append sz $blanks
    set sz [string range $sz 0 19]
  }

  set txt ""

  if { $vr != "" } {
    append txt "Var: $vr\n"
  }

  if { $ne != "" } {
    append txt "Nof entries: $ne\n"
  }
 
  if { $sz != "" } {
    append txt "Size: $sz"
  }

  Widget::configure1 $theArray(entryInfoLabel) -text $txt
}

 
# Leave event for an entry-widget
#
proc Widget::entryLeave {globArray fld wdg {event_info ""} } {
  upvar #0 $globArray theArray

  if { ![info exist theArray(entryInfoLabel)] } {
    return
  }
  
  Widget::configure1 $theArray(entryInfoLabel) -text ""
}


# ============= 
# Generic procs
# =============


# Generic for radio button
#
proc Widget::radioButton_generic { globArray fld wdg event_info} {
  global Info
  upvar #0 $globArray theArray

  if { $theArray($fld,err) } {
    return
  }

  # NOTE: Here current value is "lagging", so we do the equal
  # test to see the modification
  if { $theArray($fld) == $theArray($fld,old) } {
    set theArray($fld,mod) 1  
  } else {
    set theArray($fld,mod) 0
  }

  Widget::setRadioButtonStatus $globArray $fld $wdg
}


# Generic for check box
#
proc Widget::checkBox_generic { globArray fld wdg event_info} {
  global Info
  upvar #0 $globArray theArray

  if { $theArray($fld,err) } {
    return
  }

  # NOTE: Here current value is "lagging", so we do the equal
  # test to see the modification
  if { $theArray($fld) == $theArray($fld,old) } {
    set theArray($fld,mod) 1  
  } else {
    set theArray($fld,mod) 0
  }

  Widget::setCheckBoxStatus $globArray $fld $wdg
}



#=========================================================#
#=========================================================#
#              Combined widgets                           # 
#=========================================================#
#=========================================================#


# ========================================================#
#                      TableEntry
# ========================================================#

# Displays a table entry with unique widget id:
#
#-Entry widget for new rows
#-Listbox for all rows
#-Add,Insert,Update,Delete buttons
#-Ok,Cancel,Apply buttons
#
proc TableEntry::create { parent_id
                          frame
                          data_list
                          sort sort_marker
                          accept_multi_x
                          lb_width lb_height
                          entry_check_proc
                          datalist_pre_check_proc
                          datalist_post_check_proc
                          {globArray ""}
                        } {
  global Info

  variable id
  incr id

  # We need a unique id for these
  # variables, because multiple tables can exist
  #
  variable acceptMultiX$id
  variable dataListPreCheckProc$id
  variable dataListPostCheckProc$id
  variable entry$id
  variable entryCheckProc$id
  variable frame$id
  variable globArray$id
  variable listbox$id
  variable nofEntries$id
  variable parentId$id
  variable showIndex$id 0
  variable showLineNumbers$id 0
  variable sort$id
  variable sortMarker$id
  variable updated$id 0

  set acceptMultiX$id $accept_multi_x
  set dataListPreCheckProc$id $datalist_pre_check_proc
  set dataListPostCheckProc$id $datalist_post_check_proc
  set entryCheckProc$id $entry_check_proc
  set frame$id          $frame
  set globArray$id      $globArray
  set parentId$id       $parent_id
  set sort$id           $sort
  set sortMarker$id     $sort_marker

  set m [set frame$id]

  frame $m.ef   ;# Entry frame
  frame $m.bbf  ;# Box & buttons frame
  frame $m.lbf  ;# Listbox frame
  frame $m.btnf ;# Add etc. buttons frame

  # NOTE This entry is with table-font!
  entry $m.e -width $lb_width -font $Info(tableFont)
  set entry$id $m.e
  set wdg [set entry$id]
  bind $wdg <KeyPress-Return> "TableEntry::addListEntry $id"

  listbox $m.lb -width $lb_width -height $lb_height \
                -relief sunken -selectmode extended \
                -xscrollcommand [list $m.scrollbar_x set] \
                -yscrollcommand [list $m.scrollbar_y set] \
                -exportselection 0 -font $Info(tableFont)
  scrollbar $m.scrollbar_x -orient horizont -command [list $m.lb xview]
  scrollbar $m.scrollbar_y -orient vertical -command [list $m.lb yview]

  set listbox$id $m.lb
  set wdg [set listbox$id]
  bind $wdg <Double-1> "TableEntry::copyListEntry $id"

  button $m.add    -text "Add"     -width 8 \
                   -command "TableEntry::addListEntry $id"
  button $m.insert -text "Insert"  -width 8 \
                   -command "TableEntry::insertListEntry $id"
  button $m.update -text "Update"  -width 8 \
                   -command "TableEntry::updateListEntry $id"
  button $m.delete -text "Delete"  -width 8 \
                   -command "TableEntry::deleteListEntry $id"
  
  checkbutton $m.lines -text "Line nbrs" -command "TableEntry::toggleLineNumbers $id" \
                       -variable "TableEntry::showLineNumbers$id"

   
  #-----WIDGET PACKING
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-Frames
  set m $frame
  pack $m.ef $m.bbf -side top  -anchor w
  pack $m.lbf $m.btnf -in $m.bbf -side left -expand 1 -fill both
  
  #-Entry field
  pack $m.e -in $m.ef

  #-Listbox
  pack  $m.scrollbar_x -in $m.lbf -side bottom -fill x     
  pack  $m.lb $m.scrollbar_y -in $m.lbf -side left -expand 1 -fill y     

  #-Buttons
  pack  $m.lines -in $m.btnf -side bottom -anchor s -padx $fpx1 -pady $fpy1
  pack  $m.add $m.insert $m.update $m.delete \
                 -in $m.btnf -side top -expand 1 -anchor n -padx $fpx1 -pady $fpy1

  #-All
  pack $m -side top  -anchor nw -fill x -expand 1 -padx $fpx1 -pady $fpy1


  TableEntry::updateDataList $id $data_list 0

  return $id
}


# Delete TableEntry instance with id
#
proc TableEntry::delete {id} {
}


# Add a new row at the end listbox
#
proc TableEntry::addListEntry {id} {

  set pid [set TableEntry::parentId$id]

  #--Check entry
  #
  set wdg [set TableEntry::entry$id]

  set data_entry [$wdg get]

  $wdg icursor 0
  $wdg selection range 0 end

  set data_entry [[set TableEntry::entryCheckProc$id] $pid "add" $data_entry]

  if { $data_entry == "" } {
    return
  }

  #--Precheck data list
  #
  if { ![[set TableEntry::dataListPreCheckProc$id] $pid "add"] } {
    return
  } 

  #--Update datalist
  #
  set data_list [TableEntry::getDataList $id]

  set sort [set TableEntry::sort$id]
  set show_index end

  if {$sort} {
    set show_index -1
    append data_entry [set TableEntry::sortMarker$id]
  }
  
  lappend data_list $data_entry

  set TableEntry::showIndex$id $show_index

  TableEntry::updateDataList $id $data_list $sort

  #--Postcheck datalist
  #
  if { ![[set TableEntry::dataListPostCheckProc$id] $pid "add"] } {
    return
  }
}


proc TableEntry::copyListEntry {id} {
  global Timestep Info

  set index [[set TableEntry::listbox$id] curselection]

  if { $index == "" || [llength $index] > 1 } {

    set Info(messageIcon) info
    Message::dispOkMessage {"Select first one entry to be copied!"}                    
    return
  }

  set data_list [TableEntry::getDataList $id]

  #set new_row [string trim [lindex $data_list $index]]
  set new_row [lindex $data_list $index]

  set wdg [set TableEntry::entry$id]

  set old_row [$wdg get]
  set old_row [string trim $old_row]

  if { $old_row != "" } {
    set msg [list "NOTE: Current entry data will be replaced by the selected row!\n\n" \
                  $Info(continueOk) ]

    if { ![Message::verifyAction $Info(noviceUser) $msg] } {
      return
    }
  }

  $wdg delete 0 end
  $wdg insert 0 $new_row

}


# Delete current selection from the list
#
proc TableEntry::deleteListEntry {id} {
  global Info

  set pid [set TableEntry::parentId$id]

  set indices [[set TableEntry::listbox$id] curselection]

  if { $indices == "" } {

    set Info(messageIcon) info
    Message::dispOkMessage {"Select first the listrow(s) to be deleted!"} 
    return
  }

  #--Precheck datalist
  #
  if { ![[set TableEntry::dataListPreCheckProc$id] $pid "delete"] } {
    return
  } 

  #--Update datalist
  #
  set reversed_indices ""

  foreach index $indices {
    set reversed_indices [linsert $reversed_indices 0 $index]
  }

  set show_index [lindex $indices 0]


  set data_list [TableEntry::getDataList $id]

  foreach index $reversed_indices {
    set data_list [lreplace $data_list $index $index ]
  }

  set sort 0
  set TableEntry::showIndex$id $show_index

  TableEntry::updateDataList $id $data_list $sort

  #--Postcheck datalist
  #
  if { ![[set TableEntry::dataListPostCheckProc$id] $pid "delete"] } {
    return
  }

}


# Insert a new row into list above the current selection
#
proc TableEntry::insertListEntry {id} {
  global Info

  set pid [set TableEntry::parentId$id]

  set index [[set TableEntry::listbox$id] curselection]

  if {$index == ""} {
    set Info(messageIcon) info
    Message::dispOkMessage {"Select first the listrow above which data is to be inserted!"}
    return
  }
  
  #--Check entry
  #
  set data_entry [[set TableEntry::entry$id] get]
  set data_entry [[set TableEntry::entryCheckProc$id] $pid "insert" $data_entry]

  if { $data_entry == ""} {
    return
  }

  #--Precheck datalist
  #
  if { ![[set TableEntry::dataListPreCheckProc$id] $pid "insert"] } {
    return
  } 

  #--Update datalist
  #
  set show_index $index

  set sort [set TableEntry::sort$id]

  if {$sort} {
    set show_index -1
    append data_entry [set TableEntry::sortMarker$id]
  }

  set data_list [TableEntry::getDataList $id]
  set data_list [linsert $data_list $index $data_entry]

  set TableEntry::showIndex$id $show_index

  TableEntry::updateDataList $id $data_list $sort

  #--Postcheck datalist
  #
  if { ![[set TableEntry::dataListPostCheckProc$id] $pid "insert"] } {
    return
  }

}


# Updates current selection in the list using the value
# in the entry field
#
proc TableEntry::updateListEntry {id} {
  global Info

  set pid [set TableEntry::parentId$id]

  set index [[set TableEntry::listbox$id] curselection]

  if {$index == ""} {
    set Info(messageIcon) info
    Message::dispOkMessage {"Select first the listrow to be updated!"}                   
    return
  }
  
  #--Check entry
  #
  set data_entry [[set TableEntry::entry$id] get]
  set data_entry [[set TableEntry::entryCheckProc$id] $pid "update" $data_entry]

  if { $data_entry == ""} {
    return
  }

  #--Precheck datalist
  #
  if { ![[set TableEntry::dataListPreCheckProc$id] $pid "update"] } {
    return
  } 

  #--Update datalist
  #
  set show_index $index

  set sort [set TableEntry::sort$id]

  if {$sort} {
    set show_index -1
    append data_entry [set TableEntry::sortMarker$id]
  }
  
  set data_list [TableEntry::getDataList $id]
  set data_list [lreplace $data_list $index $index $data_entry]
  set TableEntry::showIndex$id $show_index

  TableEntry::updateDataList $id $data_list $sort

  #--Postcheck datalist
  #
  if { ![[set TableEntry::dataListPostCheckProc$id] $pid "update"] } {
    return
  }

}


# Update list box view
#
proc TableEntry::updateDataList { id data_list {sort 0} } {
  
  upvar #0 [set TableEntry::globArray$id] theArray

  set TableEntry::updated$id 1

  set marker [set TableEntry::sortMarker$id]

  if {$sort} {
    set data_list [lsort -index 0 -real $data_list]
  }

  set show_index [set TableEntry::showIndex$id]

  # If we should show a "marked" line
  #
  if {$show_index == -1 && $marker != ""} {

    set show_index [lsearch -regexp $data_list $marker]

    if {$show_index != -1} {
      set old_row [lindex $data_list $show_index]
      regsub -all $marker $old_row "" new_row

      set data_list [lreplace $data_list $show_index $show_index $new_row]
    }
  }

  # Add line numbers. if needed
  #
  if { 1 == [set TableEntry::showLineNumbers$id] } {
    set counter 0
    set new_list ""
    foreach row $data_list {
      incr counter
      lappend new_list "$counter  $row"
    }
    set data_list $new_list
  }

  set TableEntry::nofEntries$id [llength $data_list]

  ListBox::fill [set TableEntry::listbox$id]  $data_list

  [set TableEntry::listbox$id] see $show_index

}


# Get data list, so that possible line numbers
# are cleaned off
#
proc TableEntry::getDataList { id } {

  if { ![winfo exists [set TableEntry::listbox$id] ] } {
    return ""
  }

  set data_list [[set TableEntry::listbox$id] get 0 end]

  if { 0 == [set TableEntry::showLineNumbers$id] } {
    return $data_list
  }

  set new_list ""
  # Remove line numbers before returning "pure" data list
  foreach row $data_list {
    lappend new_list [lrange $row 1 end ]
  }

  return $new_list
}


# Toggles the display of line numbers
# Updates list view
#
proc TableEntry::toggleLineNumbers { id } {

  # NOTE: We have to temporarily change teh flag value
  # before we read the current data list, because
  # flag value was just changed when checkbutton was
  # clicked!
  #
  set sln [set TableEntry::showLineNumbers$id] 

  if { $sln == 1 } { 
    set TableEntry::showLineNumbers$id 0
  } else {
    set TableEntry::showLineNumbers$id 1
  }

  set data_list [TableEntry::getDataList $id]

  # Ok, now we have to restore to the current value!
  set TableEntry::showLineNumbers$id $sln

  TableEntry::updateDataList $id $data_list 0
}



# ========================================================#
#                        Data Entry
# ========================================================#

# Displays a dialog with one entry field.
# Result is put into variable referenced by argument
# 'resultvar'. It should be a global variable, because
# ok-button as a command-procedure works in global scope.
#
proc Widget::dispDataEntry  {
          {resultvar "entryResult"}
          {default ""} 
          {title "Data entry"}
          {label "Enter value: "}
                } {

  toplevel .w
  set w .w
  wm title $w $title
  frame $w.f1
  frame $w.f2
  label $w.l -text $label
  entry $w.e -width 50 
  $w.e insert 0 $default

  eval "button $w.ok -text OK -command {set " $resultvar "\[.w.e get\];destroy .w}"
  eval "button $w.cancel -text Cancel -command {set " $resultvar "\"!cancel!\";destroy .w}"

  pack $w.f1 $w.f2 -side top -fill y -expand 1 -padx 10 -pady 15
  pack $w.l -side left -in $w.f1
  pack $w.e -side left -in $w.f1 -fill x
  pack $w.ok $w.cancel -side left -in $w.f2 -padx 10
}



# ========================================================#
#                   Directory Entry
# ========================================================#

# Procedure displays a dialog with one entry field and
# directory browse button
# Result is put into variable referenced by argument
# 'resultvar'. It should be a global variable, because
# ok-button as a command-procedure works in global scope.
#
proc Widget::dispDirectoryEntry  {
          wname
          {resultvar "entryResult"}
          {default ""} 
          {title "Directory name entry"}
          {label "Enter name: "}
                } {
  global Info

  # Check that wname is a correct root window name
  #
  set wname [Util::stringTrim $wname]
  if { "." != [string index $wname 0] } {
    set wname ".$wname"
  }

  if { [winfo exists $wname] } {
    return
  }

  toplevel $wname
  set w $wname

  wm title $w $title
  frame $w.f1
  frame $w.f2

  label $w.l -text $label
  set wdg [entry $w.e -width 50 ]
  button $w.b -text "..." -width 1 \
              -command "Util::fillDirectoryEntry $wdg"

  $wdg insert 0 $default

  eval "button $w.ok -text OK -command {set " $resultvar "\[$w.e get\];destroy $w}"
  eval "button $w.cancel -text Cancel -command {set " $resultvar "\"!cancel!\";destroy $w}"

  pack $w.f1 $w.f2 -side top -fill y -expand 1 -padx 10 -pady 15
  pack $w.l -side left -in $w.f1
  pack $w.e -side left -in $w.f1 -fill x
  pack $w.b -side left -in $w.f1 -fill x -padx 10
  pack $w.ok $w.cancel -side left -in $w.f2 -padx 10
}



# ========================================================#
#                    File Entry
# ========================================================#

# Procedure displays a dialog with one entry field and
# file browse button
# Result is put into variable referenced by argument
# 'resultvar'. It should be a global variable, because
# ok-button as a command-procedure works in global scope.
#
proc Widget::dispFileEntry  {
          wname
          {resultvar "entryResult"}
          {default ""} 
          {title "File name entry"}
          {label "Enter name: "}
                } {
  global Info
  
  # Check that wname is a correct root window name
  #
  set wname [Util::stringTrim $wname]
  if { "." != [string index $wname 0] } {
    set wname ".$wname"
  }

  if { [winfo exists $wname] } {
    return
  }

  toplevel $wname
  set w $wname

  wm title $w $title
  frame $w.f1
  frame $w.f2

  label $w.l -text $label
  set wdg [entry $w.e -width 50 ]
  button $w.b -text "..." -width 1 \
              -command "Util::fillFileEntry $wdg"

  $wdg insert 0 $default

  eval "button $w.ok -text OK -command {set " $resultvar "\[$w.e get\];destroy $w}"
  eval "button $w.cancel -text Cancel -command {set " $resultvar "\"!cancel!\";destroy $w}"

  pack $w.f1 $w.f2 -side top -fill y -expand 1 -padx 10 -pady 10
  pack $w.l -side left -in $w.f1
  pack $w.e -side left -in $w.f1 -fill x
  pack $w.b -side left -in $w.f1 -fill x -padx 10
  pack $w.ok $w.cancel -side left -in $w.f2 -padx 10
}


# ========================================================#
#                  File Select Panel
# ========================================================#

# Procedure displays a dialog with one entry field and
# file browse button. In addition there is a list of radiobuttons
# where the user can spefify the input file type
#
# Result-file and selected-type-index (0,1,...) are put into global variables referenced
# by arguments: 'resultfile' and 'resulttype'.
# They should be global variables, because ok-button as a command-procedure works in global scope.
#
proc Widget::dispFileSelectPanel  {
          wname
          {wtitle "Select File"}
          {hdr "Select file type"}
          {resultpath "entryResult"}
          {resulttype "radioResult"}
          {filetypes "" }
          {extensions { {{Any} {*.* }} } }
          {default_type 0} 
          {default_dir ""} 
          {default_file ""} 
          {default_idx ""} 
                } {
  global Info $resultpath $resulttype
  upvar #0 $resultpath respath $resulttype restype
  
  # Check that wname is a correct root window name
  #
  set wname [Util::stringTrim $wname]
  if { "." != [string index $wname 0] } {
    set wname ".$wname"
  }

  if { [winfo exists $wname] } {
    return
  }

  toplevel $wname
  set w $wname
  
  wm title $w $wtitle
  frame $w.f1
  frame $w.f2
  frame $w.f3
  frame $w.f4

  # Create panel widgets
  # ====================

  # Panel header
  #
  label $w.hdr -text $hdr

  # Radiobuttons frames
  #
  frame $w.f2.sf1
  frame $w.f2.sf2

  #--Radiobuttons
  #
  # NOTE: Radiobuttons variable is the global tcl-variable name given in 'resultype'
  # Its value will be the index of the button selected (0,1,...)
  #
  set rbuttons ""
  set idx -1
  foreach ftp $filetypes {
    incr idx
    set wdg [radiobutton $w.rb$idx -anchor w \
              -indicatoron 1 -state normal -text $ftp \
              -selectcolor $Info(select_color) \
              -variable $resulttype -value $idx \
            ]
    lappend rbuttons $wdg

    if { $idx == $default_idx } {
      $wdg select
    }
  }

  #--File entry field
  #
  label $w.l -text "Enter file name  "
  set wdg [entry $w.e -width 50 ]

  $wdg insert 0 $default_file

  #--Buttons
  
  # File Browser button (...)
  # -------------------------
  #
  # Well, this is tricky, but couldn't find better way to build an argument which would be
  # evaluated correctly and at proper time!!!
  #
  # NOTE: The proper extensions list is selected at run-time by the index value in the
  # global reference variable 'resulttype'
  #
  set extensions [list $extensions]
  set cmd "Util::fillFileEntry $wdg \[list "
  append cmd  "-extensions \[lindex $extensions \$"
  append cmd "$resulttype \] "
  append cmd "-initDir \"$default_dir\" "
  append cmd "-initFile \"$default_file\" "
  append cmd "\]"

  # NOTE: Browser button command defined here!
  #
  #  -command "Widget::deactivate $rbList; wm focusmodel $w passive; $cmd; Util::doWait 50; wm focusmodel $w active; Widget::activate $rbList"
  set rbList [list $rbuttons]
  button $w.b -text "..." -width 1 \
    -command "Widget::deactivate $rbList; $cmd; Util::doWait 50; Widget::activate $rbList"
  
  # These button commands are not 'built', but handled with eval!
  #
  # NOTE: The global reference variable 'resultpath' will store the path name
  # when these buttons are pressed!
  #
  # NOTE: 'resultpath' should be set only by these buttons, because client must poll the
  # change in this variable (via 'tkwait varaible') to conclude that panel is closed/applied!!!
  #
  eval "button $w.ok -text OK -command {set " $resultpath "\[$w.e get\];destroy $w}"
  eval "button $w.cancel -text Cancel -command {set " $resultpath "\"!cancel!\";destroy $w}"

  # Pack panel
  # ==========
  
  pack $w.f1 $w.f2 $w.f3 -side top -fill both -expand 1 -padx 10 -pady 10
  pack $w.f4 -side top -fill y -expand 0 -padx 10 -pady 10
 
  # Pack radiobutton frames and radiobuttons into proper subframe
  #
  pack $w.f2.sf1 $w.f2.sf2 -side left -fill both -expand 1 -in $w.f2
  
  set hcount [expr round([llength $rbuttons] / 2)]
  set idx 0

  foreach rbtn $rbuttons {
    incr idx
    if { 1 == [expr $idx % 2] } {
      set f $w.f2.sf1
    } else {
      set f $w.f2.sf2
    }
    pack $rbtn -side top -anchor nw -expand 0 -in $f
  }

  pack $w.hdr -side left -in $w.f1
  pack $w.l -side left -in $w.f3
  pack $w.e -side left -in $w.f3 -fill x
  pack $w.b -side left -in $w.f3 -fill x -padx 10
  pack $w.ok $w.cancel -side left -in $w.f4 -padx 10

}


# ========================================================#
#                    CheckBoxList
# ========================================================#

# Displays a listbox where each row has a 'checkbox'
#
#-Listbox for all rows
#-All, None buttons
#-Ok,Cancel,Apply buttons
#
# Return: selections  (0/1 flags for each row)
#
proc CheckBoxList::create { win
                            datalist
                            in_selections
                            lb_width lb_height
                          } {
  global Info

  set selections ""

  if { $in_selections == "" } {
    foreach item $datalist {
      lappend selections 0
    }
  } else {
    set selections $in_selections
  }

  variable id
  incr id

  # We need a unique id for these
  # variables, because multiple boxes can exist
  #
  variable win$id
  variable listbox$id
  variable datalist$id
  variable selections$id
  variable updated$id
  variable closed$id
  
  set win$id $win
  set updated$id 0
  set closed$id 0
  set datalist$id $datalist
  set selections$id $selections

  set m [frame $win.f]

  frame $m.bbf  ;# Box & All,None buttons frame
  frame $m.lbf  ;# Listbox frame
  frame $m.btnf1 ;# All, None buttons frame
  frame $m.btnf2 ;# Ok etc. buttons frame

  listbox $m.lb -width $lb_width -height $lb_height \
                -relief sunken -selectmode browse \
                -xscrollcommand [list $m.scrollbar_x set] \
                -yscrollcommand [list $m.scrollbar_y set] \
                -exportselection 0 -font $Info(tableFont)
  scrollbar $m.scrollbar_x -orient horizont -command [list $m.lb xview]
  scrollbar $m.scrollbar_y -orient vertical -command [list $m.lb yview]

  set listbox$id $m.lb

  set wdg [set listbox$id]
  bind $wdg <Button-1> "CheckBoxList::setSelection $id Button-1 %x %y"
  
 
  button $m.btnf1.all -text "All" -width 3 \
                  -command "CheckBoxList::all $id"

  button $m.btnf1.none -text "None" -width 3 \
                   -command "CheckBoxList::none $id"

  button $m.btnf2.ok -text "Ok" -width 5 \
                 -command "CheckBoxList::ok $id"

  button $m.btnf2.cancel -text "Cancel" -width 5 \
                     -command "CheckBoxList::cancel $id"

  #-----WIDGET PACKING
  set fpx1 $Info(framePadX1)
  set fpy1 $Info(framePadY1)
  set fpx2 $Info(framePadX2)
  set fpy2 $Info(framePadY2)
  set fpx3 $Info(framePadX3)
  set fpy3 $Info(framePadY3)

  #-Frames
  pack $m.bbf -side top -anchor w -expand 1 -fill both
  pack $m.lbf -in $m.bbf -side left -expand 1 -fill both
  pack $m.btnf1 -in $m.bbf -side left -expand 0 -fill x
  pack $m.btnf2 -side top -expand 0 -fill y
  
  #-Listbox
  pack  $m.scrollbar_x -in $m.lbf -side bottom -expand 0 -fill x     
  pack  $m.lb -in $m.lbf -side left -expand 1 -fill both     
  pack  $m.scrollbar_y -in $m.lbf -side left -expand 0 -fill y     

  #-Buttons
  pack  $m.btnf1.all $m.btnf1.none -side top -expand 0 -padx $fpx2 -pady $fpy3

  pack  $m.btnf2.ok $m.btnf2.cancel -side left -expand 0 -padx $fpx3 -pady $fpy2

  #-All
  pack $m -side top -anchor nw -fill both -expand 1 -padx $fpx1 -pady $fpy1

  # Fill list box
  #
  CheckBoxList::fillListBox $id 

  # Wait until window is closed
  #
  vwait CheckBoxList::closed$id
  
  if { 1 == [set CheckBoxList::updated$id] } {
    return [set CheckBoxList::selections$id]
  } else {
    return $in_selections
  }
}


# Fill listbox with current/updated contents
#
proc CheckBoxList::fillListBox {id} {
  global Info

  set wdg [set CheckBoxList::listbox$id]
  set dlist [set CheckBoxList::datalist$id]
  set sels [set CheckBoxList::selections$id]
  
  set rows ""

  set index 0
  foreach item $dlist dmode $sels {

    if { $dmode == 1 } {
      set row $Info(selectionBoxTrueMarker)
    } else {
      set row $Info(selectionBoxFalseMarker)
    }

    append row " $item"
    
    lappend rows $row

    incr index
  }

  ListBox::fill $wdg $rows
}


# NOTE: Tcl special arg-name 'args' <--> all aruments in one list!
#
proc CheckBoxList::setSelection {args} {

  set id [lindex $args 0] ;# Widget id

  set event_mode [lindex $args 1] ;# Button and keys pressed
  set lb_x [lindex $args 2] ;# Local lb x-coordinate for selection
  set lb_y [lindex $args 3] ;# Local lb y-coordinate for selection 

  set wdg [set CheckBoxList::listbox$id]
  set sels [set CheckBoxList::selections$id]

  #-Selected ("clicked") row
  set idx [$wdg index @$lb_x,$lb_y]

  if { $idx == "" || $idx < 0 } {
    return
  }

  set dmode [lindex $sels $idx]

  # Toggle selection
  if { $dmode == 1 } {
    set sels [lreplace $sels $idx $idx 0]
  } else {
    set sels [lreplace $sels $idx $idx 1]
  }

  set CheckBoxList::selections$id $sels

  CheckBoxList::fillListBox $id
}


# Set all selected
#
proc CheckBoxList::all {id} {

  set sels [set CheckBoxList::selections$id]
  
  set new_sels ""
  foreach dmode $sels {
    lappend new_sels 1
  }

  set CheckBoxList::selections$id $new_sels

  CheckBoxList::fillListBox $id

  set wdg [set CheckBoxList::listbox$id]
  $wdg selection clear 0 end
}


# Set none selected
#
proc CheckBoxList::none {id} {

  set sels [set CheckBoxList::selections$id]
  
  set new_sels ""
  foreach dmode $sels {
    lappend new_sels 0
  }

  set CheckBoxList::selections$id $new_sels

  CheckBoxList::fillListBox $id

  set wdg [set CheckBoxList::listbox$id]
  $wdg selection clear 0 end
}


# Save proc
#
proc CheckBoxList::save {id} {
  
  set CheckBoxList::updated$id 1
}


# Close proc
#
proc CheckBoxList::close {id} {
    
    set wdg [set CheckBoxList::win$id]
    destroy $wdg
    set CheckBoxList::closed$id 1
}


# Ok proc
#
proc CheckBoxList::ok {id} {

  CheckBoxList::save $id
  CheckBoxList::cancel $id
}


# Cancel proc
#
proc CheckBoxList::cancel {id} {
  
  CheckBoxList::close $id
}



# ========================================================#
#                     File Browser
# ========================================================#

# Procedure creates a file browser (or editor)
#
proc FileBrowser::browse {windowName 
                          filename
                          {file_must_exist 1}
                          {canEdit 0}
                          {canDeleteCancelled 0}
                          {line_parsing_call none}
                          {channel_id 0}
                          {process_nbr 0}
                          {canUseExternal 1}
                          {width 60}
                          {height 30}
                        } {

  global Info UserSetting ProcessTable

  #--If target file should exist
  if { $file_must_exist && $filename != "" } {

    if { ![file exists $filename] } {
      set msg "File:\n $filename\nnot found!\n\n"
      set Info(messageIcon) error
      Message::dispOkMessage $msg
      return ""
    }
  }

  set source_ok 0

  #--Try to open the file

  #-Try first the channel id
  if { $channel_id != "" &&
       $channel_id != 0 &&
       ![catch {eof $channel_id}]
     } {
    set in $channel_id
    set source_ok 1

  #-Next try the filename
  } elseif { ![catch {set channel_id [open $filename] } ] } {
    set in $channel_id
    set source_ok 1

  #-Next try source array
  } elseif { 2 == [llength $filename] } {

    set globSourceArray [lindex $filename 0]
    set theSourceVar  [lindex $filename 1]

    upvar #0 $globSourceArray theSourceArray

    if { [info exists theSourceArray($theSourceVar)] } {
      set source_ok 1
      set filename ""
      set in 0
    }
  }

  #-Ouch, no data!
  if { !$source_ok } {
    set msg "Cannot open file:\nFile: $filename\n or\n ChannelId: $channel_id"
    set Info(messageIcon) error
    Message::dispOkMessage $msg
    return ""
  }

  #--We check if we should and can use an external program
  set external_command ""
  set external_ok 1
  set external_msg ""

  if { $canUseExternal } {

    if { $canEdit } {
      set external_command $UserSetting(EDITOR_COMMAND)

    } else {
      set external_command $UserSetting(BROWSER_COMMAND)
    }
  }

  #--If external command given, try it
  if { $external_command != "" } {

    #-Reformat to correct platform format
    set fname [file nativename $filename]
    
    #-Create first proper execuable name (we try separate all flags!)
    set dir_part [file dirname $external_command]
    set exe_part [lindex [file tail $external_command] 0]
 
    if { $dir_part == "." } { 
      set dir_part ""
    }

    set external_exec_command [file join $dir_part $exe_part]

    set flags [join [lrange [file tail $external_command] 1 end] ]

#MSG "cmd=$external_exec_command flags=$flags"


    #-Create first proper executable name (we try separate all flags!)
    #set external_exec_command [lindex $external_command 0]
    #set flags [lindex $external_command 1]

    if { $flags != "" } {
      set code [ catch { exec $external_exec_command $flags $fname & } msg ]
    } else {
      set code [ catch { exec $external_exec_command $fname & } msg ]
    }

    #-if something went wrong
    if {$code != 0} {
      set external_ok 0
      set external_msg "External editor/browser was not run properly, because: $msg"
    }
  }


  #--If external command was given
  if { $external_command != "" } {
    
    #-If succes, we can stop here  
    if {$external_ok} {
      close $in
      return ""

    #-Otherwise we give message and continue with 
    # internal browser/editor
    } else {
      Message::showMessage $external_msg
    }
  }
  

  #--Process item: use the possible earlier browser id
  #
  if { $process_nbr > 0 } {

    #-Check if this reopened browser for an existing process
    #
    if { [info exists ProcessTable($process_nbr,bid)] &&
         $ProcessTable($process_nbr,bid) > 0
       } {
      set bid $ProcessTable($process_nbr,bid)

      # Cancel old delayed scripts
      after cancel $Info(browser,$bid,monitorId)
      after cancel $Info(browser,$bid,readHandlerId)

      # Start reading from the beginning of the logfile!!!
      seek $channel_id 0 start

      if { $line_parsing_call == "" } {
        set line_parsing_call $Info(browser,$bid,parser)
      }

    #-Create new unique browser id for process
    } else {
      set bid [incr Info(browser,id)]
      set ProcessTable($process_nbr,bid) $bid
    }
  
  #--Non-process item: always an unique broser id
  #
  } else {
    set bid [incr Info(browser,id)]
  }

  #--Data filling line at a time using fileevent
  # NOTE: Fileevent seems to be still blocking in NT!!!
  #fconfigure $in -blocking 0 -buffering line
  set Info(browser,$bid,filename) $filename
  set Info(browser,$bid,infile) $in
  set Info(browser,$bid,act) 1
  set Info(browser,$bid,processNbr) $process_nbr
  set Info(browser,$bid,parser) $line_parsing_call
  set Info(browser,$bid,readHandlerId) -1
  set Info(browser,$bid,monitorId) -1
  set Info(browser,$bid,update) 1
  set Info(browser,$bid,modified) 0
  set Info(browser,$bid,cancelled) 0
  set Info(browser,$bid,canDeleteCancelled) $canDeleteCancelled
  set Info(browser,$bid,charCount) 0

  set Info(browser,$bid,mtime) 0

  if { ![catch { set mtime [file mtime $filename] } ] } {
    set Info(browser,$bid,mtime) $mtime
  }
  
  set Info(browser,$bid,readFromStart) 1
  set Info(browser,$bid,writeFromStart) 1

  #--What is our window like?
  set w $windowName
  set wclass ""

  if {[winfo exists $w]} {
    set wclass [winfo class $w]
  }

  #---If our widget is a toplevel window
  if { $wclass != "Frame" } {

    #--Run number is used also as window id !!!
    set id $Info(browser,$bid,processNbr)
    set w $windowName$id

    if {[winfo exists $w]} {
      wm deiconify $w
      focus $w
      return $w.text
    }

    toplevel $w
    wm title $w $filename

    # Check also 'violent' closes!
    wm protocol $w WM_DELETE_WINDOW "FileBrowser::closed $w $bid"

    # Set process table related accelerators
    if { $Info(browser,$bid,processNbr) > 0 } {
      bind $w <Control-KeyPress-t> "Util::execProc ProcessTable::openPanel"
      bind $w <Control-KeyPress-T> "Util::execProc ProcessTable::openPanel"
      bind $w <Control-KeyPress-k> "Util::execProc ProcessTable::kill $Info(browser,$bid,processNbr); focus $w; focus $w.f2.fB.ok"
      bind $w <Control-KeyPress-K> "Util::execProc ProcessTable::kill $Info(browser,$bid,processNbr); focus $w; focus $w.f2.fB.ok"
    }

    set xpos [winfo x $Info(mainWindow)]
    incr xpos -40
    set ypos [winfo y $Info(mainWindow)]
    incr ypos 20
    wm geometry $w +$xpos+$ypos
  }

  set Info(browser,$bid,window) $w

  #---Frames, widgets and packing
  frame $w.f1
  frame $w.f2
  frame $w.f2.fB

  pack $w.f2 -side bottom -fill x -expand 0 -anchor c
  pack $w.f1 -side top -fill both -expand 1

  #--Text frame
  set m $w.f1
  text $m.text -width $width -height $height -font $Info(browserFont) \
    -setgrid true -wrap none \
    -xscrollcommand [list $m.xscroll set] \
    -yscrollcommand [list $m.yscroll set]
  scrollbar $m.xscroll -orient horizontal \
    -command [list $m.text xview]
  scrollbar $m.yscroll -orient vertical  \
    -command [list $m.text yview]
  pack $m -side top -fill both -expand true
  pack $m.xscroll -side bottom -fill x
  pack $m.yscroll -side right -fill y
  pack $m.text -side left -fill both -expand true
  
  #--Buttons frame
  pack $w.f2.fB
  set m $w.f2.fB

  button $m.ok -text OK -command "FileBrowser::closed $w $bid"
  button $m.refresh -text Refresh -command "FileBrowser::refresh $bid $w.f1.text $filename"
  button $m.updateMode -text Freeze -width 8 \
                       -command "FileBrowser::toggleUpdateMode $bid $m.updateMode"
  button $m.kill -text Kill -width 5 \
                       -command "Util::execProc ProcessTable::kill $Info(browser,$bid,processNbr); focus $w; focus $w.f2.fB.ok"

  set Info(browser,$bid,panelRefreshButton) $m.refresh
  $Info(browser,$bid,panelRefreshButton) configure -state disabled

  focus $m.ok

  set mT $w.f1
  set mB $w.f2.fB

  #--In browse mode
  if {!$canEdit} {
    if { $Info(browser,$bid,processNbr) <= 0 } {
      pack $mB.refresh $mB.ok -side left -padx 10
    } else {
     pack $mB.kill $mB.updateMode $mB.ok -side left -padx 10
    }

  #--In edit mode
  } else {
    button $mB.editOk -text Ok -command "FileBrowser::save $mT.text $bid $filename $in; FileBrowser::closed $w $bid"
    button $mB.editApply -text Apply -command "FileBrowser::save $mT.text $bid $filename $in"
    button $mB.cancel -text Cancel -command "set Info(browser,$bid,cancelled) 1 ;FileBrowser::closed $w $bid"

    set Info(browser,$bid,panelApplyButton) $mB.editApply
    $mB.editApply configure -state disabled

    bind $mT.text <KeyPress> "FileBrowser::markBrowserModified 1 $bid"

    bind $mT.text <Control-s> "FileBrowser::save $mT.text $bid $filename $in"
    bind $mT.text <Control-S> "FileBrowser::save $mT.text $bid $filename $in"

    pack $mB.refresh $mB.editOk $mB.cancel $mB.editApply  -side left -padx 10

  }

  # If we are reading some process output, we have to use
  # some delay time in after call
  #
  if { $Info(browser,$bid,processNbr) > 0 } {
    set delay [set ProcessTable($Info(browser,$bid,processNbr),browserDelay)]
  } else {
    set delay 0
  }

  #--Read data from the file using a read handler
  #   
  if { $in > 0 } {
    
    FileBrowser::readHandler $bid $delay $in $mT.text $line_parsing_call
    FileBrowser::monitor $bid
  
    #-Wait until browse job is done
    #tkwait variable Info(browser,$bid,act)

    #-Remove event handlers
    #fileevent [set Info(browser,$bid,infile)] readable ""
  
    #catch { close [set Info(browser,$bid,infile)] }

  #--Read data from an array
  #
  } else {
    
    Widget::configureS $mT.text normal
 
    foreach line $theSourceArray($theSourceVar) {
      $mT.text insert end $line
      $mT.text insert end \n
    }
 
    Widget::configureS $mT.text disabled
    Widget::configureS $Info(browser,$bid,panelRefreshButton) disabled
  }
          

  if {$canEdit} {
    Widget::configureS $mT.text normal
  }

  #--Return the text widget pointer
  return $mT.text
}


# ---------------
# BROWSER MONITOR
# ---------------
#
# Monitors when a browser can be detached from
# the file channel
#
proc FileBrowser::monitor { bid } {
  global Info

  if { ![info exists Info(browser,$bid,act)] } {
    return
  }

# Uncomment this to see 'lost' browser!
#MSG "Browser  bid=$bid  process=$Info(browser,$bid,processNbr) running!"

  FileBrowser::checkRefreshStatus $bid

  #--Inactive process, remove event handlers and close file
  if { -1 == $Info(browser,$bid,act) } {

    catch { fileevent [set Info(browser,$bid,infile)] readable "" }

    set msg ""
    set ch_id $Info(browser,$bid,infile)

		if { $Info(browser,$bid,cancelled) && $Info(browser,$bid,canDeleteCancelled) } {
			catch { [file delete -force $Info(browser,$bid,filename)] msg}
		}

    catch { [close $ch_id] msg}

# Uncomment this to see closing of a browser
#MSG "Browser bid=$bid  nbr=$Info(browser,$bid,processNbr) closed and unset!"

    # Remove browser info variables
    # =============================
    unset Info(browser,$bid,act)
    unset Info(browser,$bid,processNbr)
    unset Info(browser,$bid,parser)
    unset Info(browser,$bid,readHandlerId)
    unset Info(browser,$bid,monitorId)
    unset Info(browser,$bid,filename) 
    unset Info(browser,$bid,infile) 
    unset Info(browser,$bid,mtime) 
    unset Info(browser,$bid,modified) 
    unset Info(browser,$bid,update) 
    unset Info(browser,$bid,cancelled) 
    unset Info(browser,$bid,canDeleteCancelled) 
    unset Info(browser,$bid,window)
    unset Info(browser,$bid,charCount)

    catch { unset Info(browser,$bid,panelApplyButton) }
    catch { unset Info(browser,$bid,panelRefreshButton) }
    catch { unset Info(browser,$bid,readFromStart) }
    catch { unset Info(browser,$bid,writeFromStart) }


  # Otherwise continue monitoring
  } else {
    set Info(browser,$bid,monitorId) [after 500 "FileBrowser::monitor $bid"]
  }
}


proc FileBrowser::checkRefreshStatus {bid} {
  global Info
  
  set mtime0 $Info(browser,$bid,mtime)

  if { ![catch { set mtime [file mtime $Info(browser,$bid,filename)] } ] } {
    set Info(browser,$bid,mtime) $mtime
  }

  set mtime1 $Info(browser,$bid,mtime)

  if { $mtime0 != $mtime1 } {

    if { [info exists Info(browser,$bid,panelRefreshButton) ] &&
         [winfo exists $Info(browser,$bid,panelRefreshButton)]
       } {
      $Info(browser,$bid,panelRefreshButton) configure -state normal
      $Info(browser,$bid,panelRefreshButton) configure -fg red
    }
  }
}


# --------------------
# BROWSER READ HANDLER
# --------------------
#
# This proc is the browser update handler
# bid : browser id
#
proc FileBrowser::readHandler { bid delay infile text_widget line_parsing_call } {
  global Info ProcessTable

  # If browser unset or stop from outside
  if { ![info exist Info(browser,$bid,act)] ||
       -1 == $Info(browser,$bid,act)
     } {
    return
  }

#MSG "parser=$line_parsing_call"

  # If browser freezed, continue polling, but do not
  # update browser
  if { !$Info(browser,$bid,update) } {

    # STRANGE: This is needed to maintain list format
    # of the line-parser argument. It must be a list, otherwise
    # its possible (predefined) arguments would be arguments for
    # browseReadhandler!
    #
    if { [llength $line_parsing_call] > 1 } {
      set line_parsing_call [list "$line_parsing_call"]
    }

    after $delay "FileBrowser::readHandler $bid $delay $infile $text_widget $line_parsing_call"
    return
  }

  # If browser widget is closed
  #if { ![winfo exists $text_widget] } {
  #  Message::showMessage "FileBrowser::readHandler: $text_widget browser window closed!"
  #  set Info(browser,$bid,act) 0
  #  return
  #}

  # If process is dead, we can read the rest of the file
  # as quickly as possible!
  if { $Info(browser,$bid,processNbr) > 0 &&
       ![Process::exists $Info(browser,$bid,processNbr)]
     } {
    set delay 0
  }

  set count 0

  #-Try to (block) read data
  if { ![catch {set dataline [read -nonewline $infile]}] } {
    
    if { [catch {set count [llength $dataline]} ] } {
      set dataline [string map { \" \'} $dataline]
      catch {set count [llength $dataline]} 
    }
  }

  # Data available
  # ==============
  if { $count > 0 } {
    # Add new line(s) into the browser

    # Do we apply some line parsing proc
    if { $line_parsing_call != "none" } {
      set dataline [eval [join $line_parsing_call] [list $dataline]]
    }

    # If something left to insert into widget
    if { $dataline != "" && 
         [info exists text_widget] &&
         [winfo exists $text_widget]
       } {

      Widget::configureS $text_widget normal
      
      # Always start reading from the beginning of the data (which apprarenty is not
      # then output line by line)
      #
      if { $Info(browser,$bid,writeFromStart) } {
        $text_widget delete 0.0 end
        set Info(browser,$bid,writeFromStart) 0

      # If process browser's char-count limit is exceeded
      #
      } elseif { $Info(browser,$bid,processNbr) > 0 &&
                 $Info(browser,$bid,charCount) >= $Info(browser,maxCharCount)
               } {

        set del_count [expr round(0.1 * $Info(browser,maxCharCount))]
        set idx1 "0.0"
        set idx2 "0.0 +$del_count chars"

        $text_widget delete $idx1 $idx2
        incr Info(browser,$bid,charCount) -$del_count

        set msg "----- NOTE: Log truncated -----\n"
        $text_widget insert 0.0 $msg
        incr Info(browser,$bid,charCount) [string length $msg]
      }


      #--Insert the new line into browser data
      #
      $text_widget insert end $dataline
      $text_widget insert end \n
      incr Info(browser,$bid,charCount) [string length $dataline]
      incr Info(browser,$bid,charCount) 1
      #Widget::configureS $text_widget disabled

      if { $Info(browser,$bid,processNbr) > 0 } {
        $text_widget see end
      }
    }

  # No data available
  # =================
  } else {

    set pnbr $Info(browser,$bid,processNbr)

    #--A process started
    if { $pnbr > 0 } {

      #--Process ended
      if { ![Process::exists $pnbr] } {
        
        if { $Info(displayExecMessages) } {
          Message::showMessage "Process $pnbr at END!"
        }
        
        if { ![info exists ProcessTable($pnbr,pid)] } {
          set Info(browser,$bid,act) -1
        } else {
          set Info(browser,$bid,act) 0
        }

        return

      #--File at EOF state
      } elseif { $ProcessTable($pnbr,state) == "running" &&
                 $Info(browser,$bid,readFromStart)
               } {
        catch {close $infile}
        set infile [open $ProcessTable($pnbr,logfile)]
        set Info(browser,$bid,infile) $infile
        set Info(browser,$bid,readFromStart) 0
        #fileevent $infile readable "FileBrowser::readHandler $bid $delay $infile $text_widget $line_parsing_call"
      }

    #--Simple file read (no process started)
    } else {
      if { 1 == [eof $infile] } {
        catch { close $infile }
        set Info(browser,$bid,act) 0
        return
      }
    }
  }

  # Continue data polling
  # =====================
  
  # STRANGE: This is needed to maintain list format
  # of the line-parser argument. It must be a list, otherwise
  # its possible (predefined) arguments would be arguments for
  # browseReadhandler!
  #
  if { [llength $line_parsing_call] > 1 } {
    set line_parsing_call [list "$line_parsing_call"]
  }

  set Info(browser,$bid,readHandlerId) \
    [after $delay "FileBrowser::readHandler $bid $delay $infile $text_widget $line_parsing_call"]

  return
}


proc FileBrowser::closed { w bid } {
  global Info

  if { [info exists Info(browser,$bid,processNbr)] } {

    # NOTE: This will set the browser data to be deleted!!!
    #
    if { $Info(browser,$bid,processNbr) == 0 } {
      set Info(browser,$bid,act) -1

    # NOTE: This keeps browser data active and available to the process, although
    # it is not visible!!!
    #
    } else {
      set Info(browser,$bid,act) 0
    }
  }

  destroy $w
}


proc FileBrowser::save {text_widget bid filename channel_id} {
  global Info
	
  set Info(browser,$bid,canDeleteCancelled) 0

  if { !$Info(browser,$bid,modified) } {
    return
  }

  set success 0
  set counter 0

  # Try to open a temporary file for saving
  #
  while { !$success && $counter < 64 } {
    set tmp "ElmerFrontEdit_tmp.$counter"

    if { ![catch {set out [open $tmp "w+"] } msg] } {
      set success 1
    }

    incr counter
  }

  if { !$success } {
    Message::showMessage "Can not save the file because could not open a temporary file"
  }

  puts $out [$text_widget get 0.0 end]
  close $out

  set msg ""
  catch { close $Info(browser,$bid,infile) } msg

  if { $msg != "" } {
    Message::showMessage "Can not save the file because: $msg"
  }

  # Rename the temporary file
  set msg ""
  catch { file rename -force $tmp $filename } msg

  if { $msg != "" } {
    Message::dispMessage [list "Can not save the file because: $msg"]
    return
  }

  FileBrowser::markBrowserModified 0 $bid
}


proc FileBrowser::toggleUpdateMode {bid button_widget} {
  global Info

  if { ![info exists Info(browser,$bid,act)] ||
       -1 == $Info(browser,$bid,act) 
     } {
    Widget::configure1 $button_widget -state disabled
    return
  }

  set mode [set Info(browser,$bid,update)]

  # Toggle mode
  if { $mode == 1 } {
    set mode 0
    Widget::configure1 $button_widget -text "Unfreeze"
    $button_widget configure -fg red
  } else {
    set mode 1
    Widget::configure1 $button_widget -text "Freeze"
    $button_widget configure -fg black
  }

  set Info(browser,$bid,update) $mode
}


proc FileBrowser::refresh {bid text_widget filename {append 0} } {
  global Info

  if { ![file exists $filename] } {
    return
  }

  if { $Info(browser,$bid,modified) } {
    set msg [list "NOTE: Data has been modified in the editor!\n\nAre you sure to replace the modified data?\n\n"  \
             $Info(anywayOk) ]

    if { ![ Message::verifyAction $Info(advancedUser) $msg] } {
      return
    }
  }


  set in [open $filename ]
  #set in $channel

  set tw $text_widget

  set old_state [$tw cget -state]

  Widget::configureS $tw normal

  #--Pick the current top of the window x,y-position
  set xpos [lindex [$tw xview] 0]
  set ypos [lindex [$tw yview] 0]

  #--Clear old data
  if {!$append} {
    $text_widget delete 0.0 end
  }

  #---Data filling
  while {![eof $in]} {

    set count [gets $in dataline]

    if {$count != -1 } {
      $text_widget insert end $dataline
      $text_widget insert end \n
    }
  }
  
  #--Note: The previous insert command resets the window to the beginning
  #  of the text, so we try to restore the position...
  $tw xview moveto $xpos
  $tw yview moveto $ypos

  #$text_widget insert 0.0 [read $in]

  Widget::configureS $text_widget $old_state

  if { ![catch { set mtime [file mtime $filename] } ] } {
    set Info(browser,$bid,mtime) $mtime
  }

  $Info(browser,$bid,panelRefreshButton) configure -fg black
  $Info(browser,$bid,panelRefreshButton) configure -state disabled

  close $in

  FileBrowser::markBrowserModified 0 $bid
}


proc FileBrowser::markBrowserModified { is_modified bid } {
  global Info

    if { ($is_modified && $Info(browser,$bid,modified)) ||
         (!$is_modified && !$Info(browser,$bid,modified))
       } {
      return
    }

    # Apply widget state
    if { [info exists Info(browser,$bid,panelApplyButton)] &&
         [winfo exist $Info(browser,$bid,panelApplyButton)]
       } {
      
      if { $is_modified } {
        $Info(browser,$bid,panelApplyButton) configure -state normal
      } else {
        $Info(browser,$bid,panelApplyButton) configure -state disabled
      }
    }
    
    # Windows title *-marking
    set wtitle [string trimright [wm title $Info(browser,$bid,window)] "*"]

    if {$is_modified} {
      append wtitle "*"
    }

    wm title $Info(browser,$bid,window) $wtitle

    set Info(browser,$bid,modified) $is_modified

}



# ========================================================#
#                      Text Browser
# ========================================================#

# Procedure creates a text browser (or editor)
#
proc TextBrowser::browse {windowName 
                          {title ""} 
                          {data ""}
                          {canEdit 0}
                          {scrollToEnd 0}
                          {width 60}
                          {height 30}
                          {selected_butt_call ""}
                          {selected_butt_text "Show"}
                        } {

  global Info

  #--Unique broser session id
  set bid [incr Info(browser,id)]

  set Info(browser,$bid,act) 1
  set Info(browser,$bid,update) 1
  set Info(browser,$bid,modified) 0
  set Info(browser,$bid,cancelled) 0

  #--What is our window like?
  set w $windowName
  set wclass ""

  if {[winfo exists $w]} {
    set wclass [winfo class $w]
  }

  #---If our widget is a toplevel window
  if { $wclass != "Frame" } {

    #--Bid number is used also as window id !!!
    set id $bid
    set w $windowName$id

    if {[winfo exists $w]} {
      wm deiconify $w
      focus $w
      return $w.text
    }

    toplevel $w
    wm title $w $title

    set xpos [winfo x $Info(mainWindow)]
    incr xpos -40
    set ypos [winfo y $Info(mainWindow)]
    incr ypos 20
    wm geometry $w +$xpos+$ypos
  }

  set Info(browser,$bid,window) $w

  #---Frames, widgets and packing
  frame $w.f1
  frame $w.f2

  pack $w.f1 -side top -fill both -expand 1
  pack $w.f2 -side bottom -fill x -expand 0 -anchor c

  #--Text frame
  set m $w.f1
  text $m.text -width $width -height $height -font $Info(browserFont) \
    -setgrid true -wrap word \
    -xscrollcommand [list $m.xscroll set] \
    -yscrollcommand [list $m.yscroll set]
  scrollbar $m.xscroll -orient horizontal \
    -command [list $m.text xview]
  scrollbar $m.yscroll -orient vertical  \
    -command [list $m.text yview]
  pack $m -side top -fill both -expand 1
  pack $m.xscroll -side bottom -fill x
  pack $m.yscroll -side right -fill y
  pack $m.text -side left -fill both -expand 1

  set mT $w.f1
  set mB $w.f2
  
  set $mB.selected ""

  if { $selected_butt_call != "" } {
    button $mB.selected -text $selected_butt_text -command "TextBrowser::selectedTextCall $mT.text $selected_butt_call"
  }

  #--In browse mode
  if {!$canEdit} {
    button $mB.ok -text OK -command "TextBrowser::closed $w $mT.text $bid"

    if { [winfo exists $mB.selected] } {
      pack $mB.ok $mB.selected -side left -expand 1 -padx 10
    } else {
      pack $mB.ok -side left -expand 1 -padx 10
    }

    focus $mB.ok

  #--In edit mode
  } else {
    button $mB.editOk -text Ok -command "TextBrowser::save $mT.text $bid ; TextBrowser::closed $w $mT.text $bid"
    button $mB.editApply -text Apply -command "TextBrowser::save $mT.text $bid"
    button $mB.cancel -text Cancel -command "set Info(browser,$bid,cancelled) 1 ;TextBrowser::closed $w $mT.text $bid"

    set Info(browser,$bid,panelApplyButton) $mB.editApply
    $mB.editApply configure -state disabled

    bind $mT.text <KeyPress> "TextBrowser::markBrowserModified 1 $bid"

    bind $mT.text <Control-s> "TextBrowser::save $mT.text $bid"
    bind $mT.text <Control-S> "TextBrowser::save $mT.text $bid"

    if { [winfo exists $mB.selected] } {
      pack $mB.editOk $mB.cancel $mB.editApply $mB.selected -side left -expand 1 -padx 10
    } else {
      pack $mB.editOk $mB.cancel $mB.editApply -side left -expand 1 -padx 10
    }

  }


  #--Insert data 
  #
  Widget::configureS $mT.text normal
  $mT.text insert 0.0 $data

  if { $scrollToEnd } {
    $mT.text see end
  }
  
  if {$canEdit} {
    Widget::configureS $mT.text disabled
  } else {
    Widget::configureS $mT.text normal
  }

  #--Return the text widget pointer
  return $mT.text
}


proc TextBrowser::closed { w text_widget bid } {
  global Info

  set Info(browser,$bid,act) 0

  set data [$text_widget get 0.0 end]

  destroy $w

  return $data
}


proc TextBrowser::selectedTextCall {text_widget func} {
  global Info

  if { [catch {set data [$text_widget get sel.first sel.last]}] } {
    set data ""
  }

  $func $data
}


proc TextBrowser::save {text_widget bid}  {
  global Info
	
  set Info(browser,$bid,canDeleteCancelled) 0

  if { !$Info(browser,$bid,modified) } {
    return
  }

  TextBrowser::markBrowserModified 0 $bid
}


proc TextBrowser::markBrowserModified { is_modified bid } {
  global Info

  if { ($is_modified && $Info(browser,$bid,modified)) ||
       (!$is_modified && !$Info(browser,$bid,modified))
     } {
    return
  }

  # Apply widget state
  if { [info exists Info(browser,$bid,panelApplyButton)] &&
       [winfo exist $Info(browser,$bid,panelApplyButton)]
     } {
    
    if { $is_modified } {
      $Info(browser,$bid,panelApplyButton) configure -state normal
    } else {
      $Info(browser,$bid,panelApplyButton) configure -state disabled
    }
  }
  
  # Windows title *-marking
  set wtitle [string trimright [wm title $Info(browser,$bid,window)] "*"]

  if {$is_modified} {
    append wtitle "*"
  }

  wm title $Info(browser,$bid,window) $wtitle

  set Info(browser,$bid,modified) $is_modified

}


# end ecif_tk_procsWidget.tcl
# ********************
