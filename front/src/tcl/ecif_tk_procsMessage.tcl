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
#Module:    ecif_tk_procsMessage.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Message box handling procedures
#
#************************************************************************


# For debugging!
# ==============
# Displays debug message in a Tk-message box
# Arguments {msg header} with default values
#
proc print { {msg "Debugging!"} {title "ELMER_FRONT_DEBUG_TCL_PRINT"} } {
  global Info

  if { !$Info(ELMER_FRONT_DEBUG_TCL_PRINT) } {
    return
  }

  #set msg [join $msg ]
  #tk_dialog .mBox $title $msg error 0 "OK"

  tk_messageBox -title $title -message $msg -default ok
}

proc MESSU0 { {msg ""} } {
  global Info

  set Info(messuCounter) 0
  
  MSG "MESSU-0 $msg  "
}


proc MESSU { {msg ""} } {
  global Info

  incr Info(messuCounter)
  
  MSG "MESSU-$Info(messuCounter) $msg  "
}


# Show debug message in the main window info area
#
proc MSG {messu {color ""} {nof_line_feeds 0} {append 1}} {
  global Info

  if { !$Info(ELMER_FRONT_DEBUG_TCL_PRINT) } {
    return
  }
  
  Message::showMessage $messu $color $nof_line_feeds $append
}


# Show debug message in the main window info area
# for a matching name
#
proc MSGN { name wanted_name messu {color ""} {nof_line_feeds 0} {append 1}} {
  global Info

  if { !$Info(ELMER_FRONT_DEBUG_TCL_PRINT) } {
    return
  }
  
  if { $name != $wanted_name } {
    return
  }

  Message::showMessage $messu $color $nof_line_feeds $append
}



#########################
### MESSAGE BOX procs ###
#########################


# Show message in the main window info area (a listbox)
#
proc Message::showMessage {messu {color ""} {nof_line_feeds 0} {append 1}} {
  global Info

  set last_row [$Info(messageWindow,listBox) size]

  # Limit the size of the listbox
  if { $last_row > $Info(messageWindow,maxSize) } {
    $Info(messageWindow,listBox) delete 0 0
  }

  # Replace last line (delete first last line)
  if { $append == 0 } {
    incr last_row -1
    $Info(messageWindow,listBox) delete $last_row end
  }

  # Insert new line at the end
  $Info(messageWindow,listBox) insert end $messu

  if { $color != "" } {
    $Info(messageWindow,listBox) itemconfigure end -fg $color
  }

  # Possible extra ( > 0) linefeeds
  set x 0
  while {$x < $nof_line_feeds} {
    $Info(messageWindow,listBox) insert end ""
    incr x
  }

  # Make last line visible
  $Info(messageWindow,listBox) see end
}


# Clear all messages in messagae area
#
proc Message::clearMessageArea {} {
  global Info
  
  set msg [list "All messages in the log area will be deleted!\n\n" \
                $Info(continueOk)
          ]

  if { ![Message::verifyAction $Info(noviceUser) $msg] } {
    return
  }

  $Info(messageWindow,listBox) delete 0 end
}


# Displays a Tk-message box
# NOTE: Global variable Info(messageIcon) is used for icon!
#
proc Message::dispMessage { {msg ""} {title ""} {parent ""} {type ok} {default ok} } {
  global Info

  if { $msg == "" } {
    set msg "Error!"
  }

  if { $title == "" } {
    set title "$Info(FRONT_NAME) message"
  }

  #--Set defaults and check
  if { $parent == "" && [info exists Info(thisWindow)] } {
    set parent $Info(thisWindow)
  } 

  if { $parent == "" || ![winfo exists $parent] } {
    set parent .
  }

  # Icon possibly set outside
  if { [info exists Info(messageIcon)] } {
    set icon $Info(messageIcon)  
  } else {
    set icon info
  }

  set msg [join $msg ]

  return [tk_messageBox -message $msg     \
                        -title $title     \
                        -type $type       \
                        -default $default \
                        -icon $icon       \
                        -parent $parent]

  # Back to default icon
  set Info(messageIcon) info

}


# Displays an Ok message in a Tk-message box
# Arguments {msg title} with default values
proc Message::dispOkMessage { {msg ""}
                     {title ""}
                     {parent "" }
                   } {
  return [Message::dispMessage $msg $title $parent ok ok]
}


# Displays an OkCancel verify box
#
# NOTE: verify_level is the highest user level which IS asked
# If it is 'powerUser', then all are asked!
#
# 
proc Message::verifyAction {verify_level {msg ""} {default ok} {icon warning} {title ""} {parent "" }} {
  global Info

  if { $verify_level < $Info(userLevel) } {
    return 1
  }

  set Info(messageIcon) $icon

  if { "cancel" == [Message::dispMessage $msg $title $parent okcancel $default] } {
    return 0
  } else {
    return 1
  }
}


# Displays an OkCancel message in a Tk-message box, ok as default
# Arguments {msg title} with default values
#
proc Message::dispOkCancelMessage {{msg ""} {title ""} {parent "" }} {
  return [Message::dispMessage $msg $title $parent okcancel ok]
}


# Displays an OkCancel message in a Tk-message box, cancel as default
# Arguments {msg title} with default values
#
proc Message::dispCancelOkMessage {{msg ""} {title ""} {parent "" }} {
  return [Message::dispMessage $msg $title $parent okcancel cancel]
}


# Displays an YesNoCancel message in a Tk-message box, Yes as default
# Arguments {msg title} with default values
#
proc Message::dispYesNoCancelMessage {{msg ""} {title ""} {parent "" } } {
  return [Message::dispMessage $msg $title $parent yesnocancel yes]
}


# Displays an YesNoCancel message in a Tk-message box, Cancel as default
# Arguments {msg title} with default values
#
proc Message::dispCancelYesNoMessage { {msg ""} {title ""} {parent "" }} {
  return [Message::dispMessage $msg $title $parent yesnocancel cancel]
}

# Displays an YesNo message in a Tk-message box, No as default
# Arguments {msg title} with default values
#
proc Message::dispYesNoMessage { {msg ""} {title ""} {parent "" }} {
  return [Message::dispMessage $msg $title $parent yesno no]
}


# Displays a yes-no dialog box
# Arguments {msg header} with default values
#
proc Message::dispYesNo {{msg "Error!"} {header ""}} {
  global Info

  set msg_str [join $msg ]

  if { $header == "" } {
    set header "$Info(FRONT_NAME) message!"
  }

  set result [tk_dialog .mBox $header $msg_str error 0 "Cancel" "OK"]

  switch $result {
    0 {return "cancel"}
    1 {return "ok"}
  } 
}


# Displays field help message
#
proc Message::helpMessage {globArray fld parent} {
  global Info
  upvar #0 $globArray theArray 

  # If field help defined
  #
  if { [info exist theArray(help,$fld)] } {
    set msg $theArray(help,$fld)

  # If no help for the field
  # NOTE: currently nothing is displayed
  } else {
    set msg "No help available for this field!"

    if { !$Info(ELMER_FRONT_DEBUG_TCL) } {
      return
    }
  }

  set arr [Panel::panelNameGuiToMenu $globArray]
  set var [DataField::fieldNameGuiToSif $fld]

  #set title "$arr:    $var"
  #tk_messageBox -message $msg -title $title -parent $parent -icon question

  set title "$var"
  tk_dialog  $parent.help $title $msg "" 0 Ok

}


# end ecif_tk_procsMessage.tcl
# ****************************
