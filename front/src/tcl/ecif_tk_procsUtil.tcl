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
#Module:    ecif_tk_procsUtil.tcl
#Language:  Tcl
#Date:      16.11.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Miscellanous utility procedures
#
#************************************************************************


##################
### Exec procs ###
##################

# NOTE: Not in use!!###!!
#
proc StdPanelExec::execStandardPanel {panel} {
	global Info

  set msg "All ok!"
  set code [catch {standardPanelCreate $panel} msg]
  if {0 != $code} {
    print $msg
  }
}


# Exec given tcl-proc
#
proc Util::execProc {proc_name {arguments ""} } {
	global Info

  set msg "All ok!"

  if { $arguments == "" } {
    set code [catch {$proc_name} msg]

  } else {
    set code [catch {$proc_name $arguments} msg]
  }

  if {0 != $code} {
    print $msg
  }
}


# Call interface to cpp-proc cpp_exec_
#
proc Util::cpp_exec {cmd {arg ""} } {
	global Info
	
	set do_exec 1

	set cmsg ""
	set code [catch { if { [string match -nocase "*PanelOk" $cmd] } {
						        	if { !$Info(SAVE_PANEL_DATA) } {
							          set do_exec 0
						          }
					          }
				          } cmsg ]
  
	if { $code != 0 } {
		print $cmsg
	}
	
	if { !$do_exec } {
		print "Not done: cmd=$cmd"	
		return
	}

  set emsg [cpp_exec_ $cmd $arg]

  if {$emsg != ""} {
    print "cpp-exec msg=$emsg"
  }
}


# Executes an 'array'proc like: InitialCondition(objectSelected)
# Returns 1 if proc found, otherwise returns 0
#
proc Util::execArrayProc {arrayName procName {arguments ""}} {
	global Info

  set array_proc_name $arrayName
  append array_proc_name "($procName)"

#MSG "proc=$array_proc_name"

  #-If proc exists
  if { "" != [info commands $array_proc_name] } {

    if { $arguments == "" } {
      $array_proc_name 

    } else {
      $array_proc_name $arguments
    }
#MSG "proc=$array_proc_name exits"

    return 1

  #-Otherwise
  } else {
    return 0
  }
}


# Call Elmer MATC mtc_domath routine via cpp
#
# NOTE: If precision given, Matc is set to default precision
# after evaluation!
# 
proc Util::doMatc { arg {precision ""} } {
  global Info 
  
  if { $precision == "" } {
    set to_do $arg

  } else {
    set frmt "format($precision,\"rowform\");"
    set def_frmt ";format($Info(defaultMatcPrecision),\"rowform\")"
    set to_do $frmt$arg$def_frmt
  }

	set Info(MATC_result) ""
  
  Util::cpp_exec doMatc $to_do

	while { $Info(MATC_result) == "" } {}

	return $Info(MATC_result)
}


# Get MATC variables
# 
proc Util::getMatcVariables { {all 0} } {
  global Info 

	set Info(MATC_result) ""

  Util::cpp_exec doMatc "who"

	while { $Info(MATC_result) == "" } {}
  
  # Remove info text
  #
  regsub -nocase "Constants" $Info(MATC_result) "^" allNames
  regsub -nocase "currently defined VARIABLES" $allNames "^" allNames

  #-Include constants
  if { $all } {
    set names [lrange [split $allNames "^"] 1 end]

  #-Drop constants
  } else {
    set names [lrange [split $allNames "^"] 2 end]
  }
  
  # Make a flat ws-separated list
  set names [split [join $names]]

  set result ""

  foreach nm $names {
    set nm [string trim $nm]
    if { $nm == "" } continue
    if {!$all && $nm == "ans" } continue
    lappend result $nm
  }

  return $result
}


# Get MATC function names
# 
proc Util::getMatcFunctions { {all 0} } {
  global Info 

	set Info(MATC_result) ""

  Util::cpp_exec doMatc "help"

	while { $Info(MATC_result) == "" } {}
  
  # Remove info text
  #
  regsub -nocase "Builtin Functions" $Info(MATC_result) "^" allNames
  regsub -nocase "User Functions" $allNames "^" allNames

  #-Include internal functions
  if { $all } {
    set names [lrange [split $allNames "^"] 1 end]

  #-Drop internal functions
  } else {
    set names [lrange [split $allNames "^"] 2 end]
  }
  
  # Make a flat ws-separated list
  set names [split [join $names]]

  set result ""

  foreach nm $names {
    set nm [string trim $nm]
    if { $nm == "" } continue
    lappend result $nm
  }

  return $result
}


# Return info concerning a Matc defnition
# Resul list: name type (F/V)
#
proc Util::getMatcDefInfo { def } {
  global Info
  
  set name ""
  set type ""

  #-Function name, remove function kwrd and argument part from the name
  #
  if { [string equal -nocase -length 8 "function" $def] } {

    set name [string range $def 8 end]
    set idx1 [string first "\{" $name]
    set idx2 [string first "(" $name]

    if { $idx1 < 0 } {
      set idx $idx2
    } elseif { $idx2 < 0 } {
      set idx $idx1
    } elseif { $idx1 < $idx2 } {
      set idx $idx1
    } else {
      set idx $idx2
    }

    incr idx -1
    if { $idx >= 0 } {
      set type "F"
      set name [string trim [string range $name 0 $idx]]
    }

  #-A variable name
  #
  } else {
    set type "V"
    set name [string trim [lindex [split $def "="] 0]]
  }

  return [list $name $type]
}


proc Util::generateEvent {} {
  global Info

  #event generate $Info(clbWinId) <FocusIn>
  #event generate $Info(clbWinId) <FocusIn>
  update
}


##########################
#### General routines ####
##########################


# Read one line from an input file
#
# NOTE: Skip comments and empty lines, strip end of
# line comment
# Returns: one line
#
proc Util::readInputFileLine { infile {skip_comments 1} {continue_lines 1} {previous_line ""}} {

  #-Read until a non-empty line or end of file
  while { 1 } {

    set count -1

    #--Read one line
    catch { set count [gets $infile line] }

    #--If no succes!
    if { $count == -1 } {
      return $previous_line
    }

    set line [string trim $line]

    #--Empty line
    if { $line == "" } {
      continue
    }

    #--Find possible comment character positions
    set idx [string first "!"  $line]
    set idx2 [string first "#" $line]

    if { $idx < 0 } {
      set idx $idx2
    }
    if { $idx2 < 0 } {
      set idx2 $idx
    }
    if { $idx2 < $idx } {
      set idx $idx2
    }

    #--Strip possible comments
    if { $skip_comments && $idx != -1 } {
      incr idx -1
      set line [string range $line 0 $idx]
    }
    
    #--If empty line (after stripping comment)  
    if { $line == "" } {
      continue
    }

    #--Check possible line continuation mark ("\" or "&")
    #  and append next line to the result if necessary
    #
    set ch [string index $line end]

    # Skip the marker and add next line(s) to the result
    #
    if { $continue_lines && ($ch == "\\" || $ch == "&") } {
      set line [string range $line 0 end-1]
      return [Util::readInputFileLine $infile $skip_comments 1 $line]
    }

    #--Build result
    #
    set result $previous_line
    append result $line

    return $result
  }
}


# Try to get the value for an global array variable
# If variable does not exists, return blank
#
proc Util::getArrayValue { globArray var} {
  upvar #0 $globArray theArray

  if { ![info exists theArray($var)] } {
    return ""
  } else {
    return $theArray($var)
  }
}


# Remove (unset) all array variables
#
proc Util::unsetArray {globArray} {

  Util::unsetArrayVariables $globArray *
}


# Remove (unset) array variables (like Material(DENSITY,err) etc.)
#
proc Util::unsetArrayVariables { globArray {vars} } {
  upvar #0 $globArray theArray
  
  if { ![array exists theArray] } {
#MSG "Array $globArray does NOT exists!"
    return
  }

  foreach vn $vars {

    #-Form variable pattern (like DENSITY,*)
    append pat $vn

    #-Pick all array names by pattern
    global $globArray
    set cnames [array names $globArray $pat]

#MSG "Unsetting vars=$vars pat=$pat cnames=$cnames"

    #-Unset all picked names
    foreach cname $cnames {
      set aname $globArray
      append aname ($cname)
      catch {unset $aname}
    }
  }
}


# Remove (unset) array id-variables (like Equation(1,mask) etc.)
# NOTE: If no vars given, delete all vars! (pattern is "*")
#
proc Util::unsetArrayIdVariables { globArray {ids ""} {flds "*"} } {
  global $globArray
  upvar  #0 $globArray theArray

  if { $ids == "" } {
    if { ![info exists theArray(ids)] } {
      return
    }

    set ids $theArray(ids)
  }

  foreach id $ids {

    foreach vn $flds {

      #-Variable pattern
      set pat $id
      append pat ","
      append pat $vn

      #-Pick all array names by pattern
      set cnames [array names $globArray $pat]

      #-Unset all picked names
      foreach cname $cnames {
        set aname $globArray
        append aname ($cname)
        catch {unset $aname}
      }
    }
  }
}


# Convert user entered panel name to Gui format
#
# NOTE: Accepted name pairs (input-output) are given in
# argument 'name_table'
#
# NOTE: simple upvar here, NOT like upvar #0 etc.
#
proc Util::convertPanelNameToGuiFormat { name_table panel gui_panel} {
  upvar $gui_panel gui_pnl

  foreach pair $name_table {

    # Skip table comment lines
    set tmp [string trim $pair]
    if {"#" == [string index $tmp 0] } {
      continue
    }

    set in_file [lindex $pair 0]
    set in_gui  [lindex $pair 1]
    
    # Remove quotes and all spaces
    set panel [Util::stringRemoveQuotes $panel]
    set panel [string map {" " ""} $panel]

    if { [string equal -nocase $in_file $panel] } {
      set gui_pnl $in_gui
      return 1
    }
  }

  # Error!
  return 0
}


# make a panel array name from a string
# Like:  Boundary condition --> BoundaryCondition
#
proc Util::makeArrayName {name} {
  
  set arr_name ""

  foreach nm $name {
    append arr_name [string totitle [string trim $nm]]
  }
  
  return $arr_name
}

# Remove surrounding spaces and quotes from a string
#
proc Util::stringTRQ {name {always_remove 1} {remove_ws 1} } {
  set name [Util::stringTrim $name]
  set name [Util::stringRemoveQuotes $name $always_remove $remove_ws]
  return $name
}


# Remove surrounding quotes from a string
#
proc Util::stringRemoveQuotes {name {always_remove 1} {remove_ws 1} } {

  set tmp [string trim $name " \t"]

  if { -1 != [string first " " $tmp] } {
    set has_blanks 1
  } else {
    set has_blanks 0
  }

  # If not force-flag, we do not remove quotes if the
  # string contains internal blanks!
  #
  if { $has_blanks && !$always_remove } {
    return $name
  }

  # Remove quotes
  set name [string trim $name "\""]
  set name [string trim $name "\'"]

  # Remove also white space possibly created after removing the quotes
  if { $remove_ws } {
    set name [string trim $name " \t"]
  }

  return $name
}


#--Remove surrounding spaces
#--Change multiple inner spaces to single spaces
#
proc Util::stringTrim {name} {

  set name [string trim $name " \t"]

  while {1} {
    set val [string map {"  " " "} $name]

    if { 0 == [string compare $val $name] } {
      return $name
    }

    set name $val
  }
}


# String non-case sensitive equal (s1=s2)
#
proc Util::nce { s1 s2 } {
  return [string equal -nocase $s1 $s2]
}


# Format number
#
proc Util::format { frmt nbr } {

  if { [catch { set result [ ::format $frmt $nbr] } msg ] } {
    set test [ ::format $frmt 0 ]
    set result [string repeat "-" [string length $test]]
  }

  return $result
}
    

proc Util::doWait { {time 50} } {
  # NOTE: This wait is needed to make the tk_optionMenu variable
  # to update!!! (strange indeed!!!)
  set done 0
  after $time "set done 1"
  vwait done
}


# Set new title for the Front main window
# If argument is empty, default (No model) title is set
# Add "*" at the end of the name if model data modified
#
proc Util::setMainWindowTitle { {title ""} } {
  global Model

  #-Title given
  if { $title != "" } {
    set wtitle $title

  #-Form default title
  } else {

    set wtitle "ELMER Front - "

    if { $Model(CASE_NAME) != "" } {
      
      append wtitle $Model(CASE_NAME)

    } else {
      append wtitle "No model name"
    }
  }

  Util::updateMainWindowTitle $wtitle
}


# Update title for the gui main window
# Add "*" at the end of the name if model data modified
#
proc Util::updateMainWindowTitle { {wtitle ""}} {
  global Info Model

  # If no title given, use current title
  if { $wtitle == "" } {
    set wtitle [wm title $Info(mainWindow)]
  }

  set wtitle [string trimright $wtitle "*"]

  # If model data modified
  #
  if { $Info(markMainWindowTitle) && 
       $Model(Front,needsUpdate)
     } {
    append wtitle "*"
  }

  wm title $Info(mainWindow) $wtitle
}


proc Util::setTearOffTitle {window_varname menu window} {
  global Model Info
  
  # Store windo name  
  set Info($window_varname) $window

  set title "$Model(CASE_NAME) "
  append title [string trim [wm title $window] "."]

  wm title $window $title
}


# Form timestring of current time and store
# it also into global Info(currentTimestamp)
#
proc Util::setCurrentTimestamp {} {
  global Info

  set Info(currentTimestamp) [clock format [clock seconds] -format %c]

  return $Info(currentTimestamp)
}


# Format a number to a minimum value
# of decimals (add trailing zeroes if needed)
#
proc Util::formatNumber {value width nof_decimals} {

  #--Trim spaces
  set tvalue [string trim $value]

  #-trim leading zeroes
  set tvalue [string trimleft $value "0"]
  
  set tmp [split $tvalue "e"]
  set tvalue1 [lindex $tmp 0]
  set tvalue2 [lindex $tmp 1]

  set width1 $width

  if { $tvalue2 != "" } {
    set width1 [expr $width1 - 1 - [string length $tvalue2]]
  }

  #-first zero is needed
  if { [string index $tvalue1 0] == "." } {
    set tmp "0"
    append tmp $tvalue1
    set tvalue1 $tmp
  }

  #--Field is already "full"
  if { $width1 == [string length $tvalue1] } {
    set result $tvalue1
    if { $tvalue2 != "" } {
      append result "e"
      append result $tvalue2
    }
    return $result
  }

  set dpos [string first "." $tvalue1]

  #--Count current decimals
  # If there are no decimals
  if { -1 == $dpos } {
    set current_nof_decimals 0
    append tvalue1 "."

  # Count current nof decimals
  } else {
    set tlen [string length $tvalue1]
    set current_nof_decimals [ expr $tlen - $dpos - 1]
  }

  #--Append "enough" trailing zeroes
  set nof_to_add [expr $nof_decimals - $current_nof_decimals]

  set count 0
  while { $count < $nof_to_add } {
    incr count
    append tvalue1 "0"
  }

  #--Append "enough" leading spaces
  set result ""
  set tlen [string length $tvalue1]
  set nof_to_add [expr $width1 - $tlen]

  set count 0
  while { $count < $nof_to_add } {
    incr count
    append result " "
  }

  append result $tvalue1

  if { $tvalue2 != "" } {
    append result "e"
    append result $tvalue2
  }

  return $result

}


# Converts time string to clock second
#
proc Util::timeStr2Clock {time_str} {
  return [clock scan $time_str]
}


# Converts clock seconds to time string
#
proc Util::timeClock2Str {time} {
  
  if { $time == "" || $time <= 0 } {
    return ""
  }

  return [clock format $time -format %c]
}


# Converts absolute seconds to time string
#
proc Util::timeSec2Str {time} {

  if { $time == "" } {
    return ""
  }

  set t $time

  # Days
  set d [expr $t / 86400]
  set t [expr $t - $d * 86400]

  # Hours
  set h [expr $t / 3600]
  set t [expr $t - $h * 3600]

  # Minutes
  set m [expr $t / 60]

  # Seconds
  set s [expr $t - $m * 60]

  if { $d != 0 } {
    set d [format "% 3" $d]
  } else {
    set d "   "
  }
  set h [format "%02u" $h]
  set m [format "%02u" $m]
  set s [format "%02u" $s]

  return "$d $h:$m:$s"
}


# Converts absolute seconds to time H:M:S string
#
proc Util::timeSec2HourStr {time} {

  set t $time

  # Hours
  set h [expr $t / 3600]
  set t [expr $t - $h * 3600]

  # Minutes
  set m [expr $t / 60]

  # Seconds
  set s [expr $t - $m * 60]

  set h [format "%02u" $h]
  set m [format "%02u" $m]
  set s [format "%02u" $s]

  return "$h:$m:$s"
}


# Converts absolute seconds exactly to time H:M:S.xxx string
#
proc Util::timeSec2HourStrExact {time} {

  set t round($time)

  # Parts of Second
  set p [expr $time - $t]

  # Hours
  set h [expr $t / 3600]
  set t [expr $t - $h * 3600]

  # Minutes
  set m [expr $t / 60]

  # Seconds
  set s [expr $t - $m * 60]


  set h [format "%02u" $h]
  set m [format "%02u" $m]
  set s [format "%02u" $s]

  if { $p != 0 } {
    set p [lindex [split $p "."] 1]
    return "$h:$m:$s.$p"
  } else {
    return "$h:$m:$s"
  }
}


# Compare two times in numeric format
# Return:
# time1 <  time2  --> -1
# time1 == time2  -->  0
# time1 >  time2  -->  1
#
proc Util::compareTimestamps {$time1 $time2} {

  return 0
}

proc Util::normalizeVector {x y z} {
  set norm [expr sqrt([expr $x*$x + $y*$y + $z*$z])]

  if { $norm != 0 } {
    set x_ [expr $x / $norm]
    set y_ [expr $y / $norm]
    set z_ [expr $z / $norm]
    return [list $x_ $y_ $z_]
  } else {
    return [list $x $y $z]
  }
}


# Combine two masks to a unique list of mask items
#
proc Util::combineMasks {m1 m2} {

  set masks [lsort "$m1 $m2"]

  set prev_m ""
  set result ""

  # Pick unique items
  foreach m $masks {
    if { $prev_m != $m } {
      append result $m
      append result " "
      set prev_m $m
    }
  }

  return [string trim $result]
}


# Remove "mask_to_remove" from "mask"
#
proc Util::removeMask {mask mask_to_remove} {

  set result ""

  # Pick unique items
  foreach m $mask {
    if { $m != $mask_to_remove } {
      append result $m
      append result " "
    }
  }

  return [string trim $result]
}


# If a variable mask (vmask) matches the panel mask (pmask)
#
proc Util::masksAreMatching {vmask pmask} {

  return [Util::patternMatchesString $vmask $pmask]
}



proc Util::masksAreEqual { mask1  mask2 } {

  if { [llength $mask1] != [llength $mask2] } {
    return 0
  }

  foreach m1 $mask1 m2 $mask2 {

    if { $m1 != $m2 } {
      return 0
    }
  }

  return 1
}


proc Util::patternMatchesString_org {pattern str} {
  global Info

  if { $pattern == "" } {
    return 1
  }

  set len [string length $pattern]
  set i 0
  set in_exp 0

  # Parse RE items like (HT1) in ((FT)&(HT)&!(HT1))
  # ===============================================
  while { $i < $len } {

    # Read next character
    # ===================
    set c [string index $pattern $i]

    #--Check if expression could start
    if {$c == "("} {
      set in_exp 1
    }

    #--Check if expression ends
    if {$c == ")"} {
      set in_exp 0
    }

    # Pick subexpression
    # ==================
    # If a letter (mask) group starts after "("
    if { $in_exp && [regexp "\[A-Z\]" $c] } {
      set sub $c
      incr i

      #--Read characters untill the next ")"
      while { $i < $len } {
        set c [string index $pattern $i]
        if { $c != ")" } {
          append sub $c
        } else {
          set in_exp 0
          break
        }
        incr i
      }

      #--Replace all sub_expressions in the "pattern" with "½"(1) or "§"(0)  
      if { [regexp $sub $str] } {
        regsub -all $sub $pattern $Info(maskTrueMark) pattern
      } else {
        regsub -all $sub $pattern $Info(maskFalseMark) pattern
      }

      #--Start from beginning again!
      set i -1
      set len [string length $pattern]
    }

    incr i
  }

  regsub -all $Info(maskTrueMark) $pattern 1 pattern
  regsub -all $Info(maskFalseMark) $pattern 0 pattern

  set result [expr $pattern]

  return $result
}


# pattern: typically a field mask like (HT1)
# str : typically a problem mask like ((FT)&(HT)&!(HT1))
#
proc Util::patternMatchesString {pattern str {empty_str_matches_all 1} } {
  global Info

  if { $pattern == "" || ($str == "" && $empty_str_matches_all) } {
    return 1
  }

  set len [string length $pattern]
  set i 0

  # Parse RE items like (HT1) in ((FT)&(HT)&!(HT1))
  # ===============================================
  while { $i < $len } {
  
    set new_pos $i
    set sub [Util::findNextSubPattern $pattern $i "new_pos"]


    if { $sub == "" } {
      break
    }

    #--Replace all sub_expressions in the "pattern" with "½"(1) or "§"(0)  

    if { [regexp $sub $str] } {
      regsub -all $sub $pattern $Info(maskTrueMark) pattern

    } else {
      regsub -all $sub $pattern $Info(maskFalseMark) pattern
    }

    #--Start from beginning again!
    set i -1
    set len [string length $pattern]

    incr i
  }

  regsub -all $Info(maskTrueMark) $pattern 1 pattern
  regsub -all $Info(maskFalseMark) $pattern 0 pattern

  set result [expr $pattern]

  return $result
}


proc Util::findNextSubPattern {pattern current_pos result_pos  } {
  global Info
  upvar 1 $result_pos pos

  set pos $current_pos

  set len [string length $pattern]

  set sub ""
  set in_exp 0

  while { $pos < $len } {

    # Read next character
    # ===================
    set c [string index $pattern $pos]

    #--Check if expression could start
    while { $c == "(" } {
      set in_exp 1
      incr pos
      set c [string index $pattern $pos]
    }
    
    # Pick subexpression
    # ==================
    # NOTE: We use the special character "¤" to mark
    # the beginning of the mask!!!
    if { $in_exp && [regexp $Info(maskStart) $c] } {
      set sub $c
      incr pos
    
      #--Read characters untill the next ")"
      while { $pos < $len } {
        set c [string index $pattern $pos]
        if { $c != ")" } {
          append sub $c
        } else {
          return $sub
        }
        incr pos
      }

    } else {
      set in_exp 0
    }
 
    incr pos
  }

  return $sub
} 



####################
### Window procs ###
####################


# Check if panel window exists
#
proc Util::panelExists {globArray} {
  upvar #0 $globArray theArray

  if { [winfo exists $theArray(winName) ] } {
    return 1
  } else {
    return 0
  }
}


# Check if window already exists (before opening it)
# NOTE: wid: the atom name for window
#       wname: the widget name for window
#       wgeom: can be empty ("") 
#
proc Util::checkPanelWindow {globArray wid wname wgeom {update_window_list 1} } {
  global Info gMW
  upvar #0 $globArray theArray

  if { $wid == -1 } {
    return 0
  }

  # Check if window is open
  set w [winfo atomname $wid]

  #-If window is distroyed
  #
  if {[winfo exists $w]} {
    set exists 1
  } else {
    set exists 0
  }

  if {$update_window_list} {
    set tmp ""
    
    # Skip this window (wid), but add all other as they are
    foreach row $Info(windowList) {
      if { $wid != [lindex $row 1] } {
        lappend tmp $row
      }
    }

    # Insert this (wid) in the beginning of the list

    if { ( [info exist theArray(dataChanged)] && $theArray(dataChanged) )  ||
         ( [info exist theArray(dataModified)] && $theArray(dataModified) )
       } {
      set Info(windowList) [linsert $tmp 0 [list $globArray $wid "*$wname"]]

    } else {
      set Info(windowList) [linsert $tmp 0 [list $globArray $wid $wname]]
    }

    MenuBuild::createWindowListMenu $gMW(Window,windowList)
  }


  #-If not exists
  if { !$exists } {
    return 0

  #-Windo exists, activate it
  } else {
    wm deiconify $w

    if {$wgeom != ""} {
      wm geometry $w $wgeom
    }

    focus $w
    raise $w

    return 1
  }
}


# Form a full path
#
# NOTE: Check ./path and ../path formats
#
proc Util::formFullPath {path } {

  set path [string trim $path]

  #--Check that path is legal!!!
  if { [catch { file pathtype $path } ] } {
    return $path
  }


  #--If abslute path, ok
  if { "absolute" == [file pathtype $path] } {
    return $path
  }

  #--If relative path
  if { "relative" == [file pathtype $path] } {

    if { ".." == [string range $path 0 1] } {
      set pth [string range $path 2 end]
      set pth [string trim $pth]
      set pth [string range $pth 1 end]
      set pth [string trim $pth]
      set cur_pwd [pwd]     ;# Store current pwd
      cd ".."               ;# Back to parent
      set bs [pwd]          ;# Get parent cd
      cd $cur_pwd           ;# back to current pwd

    } elseif { "." == [string range $path 0 0] } {
      set pth [string range $path 1 end]
      set pth [string trim $pth]
      set pth [string range $pth 1 end]
      set pth [string trim $pth]
      set bs [pwd]

    } else {
      set pth $path
      set bs ""
    }

    set dir  [file join $bs $pth]
    return  $dir
  }

  return $path
}



# Make a path absolute
# relative_base: a possible base for relative paths
# modified_flag: a possible reference variable to used by caller
#
proc Util::makePathAbsolute { inpath {base ""} {modified_flag ""} } {

  set path [Util::formFullPath $inpath]
  set base [string trim $base]

  if { $modified_flag != "" } {
    upvar  $modified_flag modified
  } else {
    set modified 0
  }

  #--Check that path is legal!!!
  if { [catch { file pathtype $path } ] } {
    return $path
  }

  # If no base given for relative values, we use
  # current working directory as the base
  if { $base == "" } {
    set base [pwd]
  }

  #--If abslute path, ok
  if { "absolute" == [file pathtype $path] } {
    set modified 0
    return $path
  }

  #--If relative path, add base
  if { "relative" == [file pathtype $path] } {
    set modified 1
    return  [file join $base $path]
  }

  #--Volumerelative is the only version left which we 
  #  know how to handle!
  if { "volumerelative" != [file pathtype $path] } {
    set modified 0
    return $path
  }

  #--Ok we have volumerelative!
  set vols [file volume]

  set path_parts [file split $path]

  set path_vol [string toupper [lindex $path_parts 0]]
  set path_vol [string trim $path_vol "/"]
  append path_vol "/"

  set path_start ""
  set path_rest [lrange $path_parts 1 end]

  # If no volume-part in the path, we have to add current volume
  if { -1 == [lsearch $vols $path_vol] } {
    set path_start [join [lindex [file split [pwd]] 0]]

  # Pick current directory in the path's volume
  } else {
    set cur_pwd [pwd]     ;# Store current pwd
    cd $path_vol          ;# Cd to volume
    set path_start [pwd]  ;# Get pwd in the volume
    cd $cur_pwd           ;# Restore to curent pwd
  }

  # Build absolute path
  set new_path $path_start

  foreach pp $path_rest {
    append new_path " $pp"
  }

  set new_path [eval "file join $new_path"]

  set modified 1
  return $new_path
}


proc Util::makeTclFilename { fname } {
  global Info

  #-Windows
  if { $Info(platform) == "windows" } {
    set fn [string map {"\\" "/"} $fname]
  } else {
    set fn $fname
  }

  return [string trim $fn]
}


# Checks that directory exists and ask the user
# the permission to create it if not existing
# Two Button (Ok/Cancel) version which does not allow
# continuing without creating the directory
#
proc Util::checkAndCreateDirectory2B { dir_name
                                     name_msg
                                     {ask_to_create 0}
                                     {create_msg ""}
                                     {verify_existing 0}
                                     {verify_msg ""}
                                   } {
  global Info

  # If directory exists verify using it.
  # ===================================
  if { [file isdirectory $dir_name] } {

    if {$verify_existing} {
      set msg [list "Directory:\n\n" \
                    "$dir_name\n\n" \
                    "$name_msg already exists.\n\n" \
                    "Use it anyway?" ]

      if { $verify_msg != "" } {
        lappend msg "\n\n$verify_msg"
      }

      if { ![Message::verifyAction $Info(noviceUser) $msg ok question] } {
        return "cancel"
      } else {
        return "ok"
      }

    } else {
      return "ok"
    }
  }
  
  # Does not exist, ask if it can be created
  # ========================================

  # We ask novices always!
  if {$ask_to_create} {
    set level $Info(powerUser)
  } else {
    set level $Info(noviceUser)
  }

  set msg [list "Directory:\n\n" \
                "$dir_name\n\n" \
                "$name_msg does not exist.\n\n" \
                "Create the directory?" ]

  if { $create_msg != "" } {
     lappend msg "\n\n$create_msg"
  }

  if { ![Message::verifyAction $level $msg ok question] } {
    return "cancel"
  }

  #-Ok, try to create the directory
  if { ![Util::createDirectory $dir_name] } {
    return "cancel"
  }

  return "ok"
}


# Checks that directory exists and ask the user
# the permission to create it if not existing
# Three Button (Yes/No/Cancel) version which allows
# continue without creating the directory
#
proc Util::checkAndCreateDirectory3B { dir_name
                                     name_msg
                                     {ask_to_create 0}
                                     {create_msg ""}
                                     {verify_existing 0}
                                     {verify_msg ""}
                                   } {
  global Info

  # If directory exists verify using it.
  # ===================================
  if { [file isdirectory $dir_name] } {

    if {$verify_existing} {
      set msg [list "Directory:\n\n" \
                    "$dir_name\n\n" \
                    "$name_msg already exists.\n\n" \
                    "Use it anyway?" ]

      if { $verify_msg != "" } {
        lappend msg "\n\n$verify_msg"
      }

      if { ![Message::verifyAction $Info(noviceUser) $msg ok question] } {
        return "cancel"
      } else {
        return "ok"
      }

    } else {
      return "ok"
    }
  }
  
  # Does not exist, ask if it can be created
  # =========================================

  set rc "ok"

  # We ask novices always, power users never!
  #
  if { ($ask_to_create && ( $Info(userLevel) < $Info(powerUser) )) ||
        $Info(userLevel) == $Info(noviceUser)
     } {

    set msg [list "Directory:\n\n" \
                  "$dir_name\n\n" \
                  "$name_msg does not exist.\n\n" \
                  "Create the directory?" ]

    if { $create_msg != "" } {
       lappend msg "\n\n$create_msg"
    }
    
    # Ask
    set Info(messageIcon) question
    set rc [Message::dispYesNoCancelMessage $msg]
  }


  # Check answer
  # ============

  # Should not continue
  if { $rc == "cancel" } {
    return "cancel"
  }
  
  # Continue, but do not create directory
  if { $rc == "no" } {
    return "no"
  }

  # Try to create the directory
  # ===========================
  if { ![Util::createDirectory $dir_name] } {
    return "cancel"
  }

  return "ok"
}



proc Util::createDirectory {dir_name} {
  global Info

  #-NOTE: Max directory depth 16!
  set path_levels {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15}

  set path_list [file split $dir_name]
 
  set count [llength $path_list]
  set counter 0

  while { $counter < $count } {

    set levels [lrange $path_levels 0 $counter]

    set path_to_create ""

    foreach i $levels {

      set path_to_create [file join $path_to_create [ lindex $path_list $i]]
      set msg ""
      set code [catch { file mkdir $path_to_create} msg]

      if { $code != 0 } {
        set msg [list "Directory name $path_to_create :\n" \
                      "Cannot create the directory because:\n" \
                      $msg ]

        set Info(messageIcon) error
        Message::dispOkMessage $msg
        return 0
      }
    }
    incr counter
  }

  return 1
}


# Return subdirectoies under the parent as a list
#
proc Util::getFileList {parent_dir fpattern} {
  global Info ModelProperty

  if { ![file isdirectory $parent_dir] } {
    return ""
  }

  set tmp_file [file join $ModelProperty(MODEL_DIRECTORY,absolute) "dirlist.txt"]

  set pdir [file join $parent_dir $fpattern]

  set pdir [file nativename $pdir]


  set msg ""

  # Windows (note check this for Win2000)
  # =======
  if { $Info(platform) == "windows" } {
    
    # WinNT and Win2000
    if { $Info(os) == $Info(os,WindowsNT) ||
         $Info(os) == $Info(os,Windows2K)
       } {
      catch {set rc [exec cmd.exe /c dir /B /A-D $pdir > $tmp_file] } msg
    
    # Win98 and W95
    } else {
      catch {set rc [exec command.com /c dir /B /A-D $pdir > $tmp_file] } msg
    }

  # Unix
  # ====
  } else {
    catch { set rc [exec find $pdir -type f > $tmp_file] } msg
  }

  # If some error catched!
  #
  if { $msg != "" } {
    set msg [ list "Cannot find files: $pdir\n\n" \
                   "because: \n\n" \
                   "$msg\n\n"]

    ;#Message::dispOkMessage $msg

    return ""
  }

  set ch [open $tmp_file]

  set files ""
  
  while { ![eof $ch] } {
 
    set nm [gets $ch]

    # Pick filenames
    if { [file isfile [file join $parent_dir $nm]] } {
      lappend files [file tail $nm]
    }
  }
  
  close $ch

  file delete $tmp_file
 
  return $files
}


# Return subdirectoies under the parent as a list
#
proc Util::getSubdirectoryList {parent_dir} {
  global Info ModelProperty

  if { ![file isdirectory $parent_dir] } {
    return ""
  }

  set tmp_file [file join $ModelProperty(MODEL_DIRECTORY,absolute) "dirlist.txt"]

  set pdir [file nativename $parent_dir]

  set msg ""

  # Windows (note check this for Win2000)
  # =======
  if { $Info(platform) == "windows" } {
    
    # WinNT and Win2000
    if { $Info(os) == $Info(os,WindowsNT) ||
         $Info(os) == $Info(os,Windows2K)
       } {
      catch {set rc [exec cmd.exe /c dir /B /AD $pdir > $tmp_file] } msg
    
    # Win98 and W95
    } else {
      catch {set rc [exec command.com /c dir /B /AD $pdir > $tmp_file] } msg
    }

  # Unix
  # ====
  } else {
    catch { set rc [exec find $pdir -type d > $tmp_file] } msg
  }

  # If some error catched!
  #
  if { $msg != "" } {
    set msg [ list "Cannot find directories under: $pdir\n\n" \
                   "because: \n\n" \
                   "$msg\n\n"]

    Message::dispOkMessage $msg

    return ""
  }

  set ch [open $tmp_file]

  set dirs ""
  
  while { ![eof $ch] } {
 
    set nm [gets $ch]

    # Skip non-directories (and . ..)
    if { $nm != ""    &&
         $nm != "."   &&
         $nm != ".."  &&
         $nm != $pdir &&
         [file isdirectory [file join $parent_dir $nm]]
       } {

      lappend dirs [file tail $nm]
    }
  }
  
  close $ch

  file delete $tmp_file
 
  set sub_dirs ""

  foreach nm $dirs {

    set nm [string trim $nm]

    if { $nm != "" } {
      lappend sub_dirs $nm
    }
  }

  return $sub_dirs
}


# Check that "path" is writable file
# NOTE: If file does not exist, it is writable!
#
proc Util::isWritableFile {path} {

  # Ok, clients can (possibly) create the file
  if { ![file exists $path] } {
    return 1
  }

  # File is not writable!
  if { ![file writable $path] } {
    return 0
  }

  # NOTE Take care that file is closed after testing it!!!

  # File is not write openable!
  if { [catch { set ch_id [open $path r+] }] } {
    catch { close $ch_id }
    return 0
  }

  catch { close $ch_id }

  return 1
}


proc Util::isSimpleFilename {name} {
  if { "." == [file dirname $name] } {
    return 1
  } else {
    return 0
  }
}


# Open a directory browser to set a directory
# the entry widget "wdg"
# NOTE: uses a Tcl-extension which
# seems to be a bit slow
#
proc Util::fillDirectoryEntry { wdg args } {
  global Info

  if {![winfo exists $wdg]} {
    return
  }

  set initial_dir ""
  set append 0
  set panel ""
  set field ""

  set args [split [join $args] "-"]

  foreach arg $args {

    set kwd [lindex $arg 0]
    set val [lindex $arg 1]

    switch $kwd {
      append     { set append $val }
      initDir    { set initial_dir $val }
      panel      { set panel $val }
      field      { set field $val }
    }
  }

  set arguments "-title \"Select directory\""

  # Clean curent value
  set cur_val [$wdg get]
  set cur_val [string trim [string trim $cur_val] ";"]

  if { $initial_dir != "" } {
    append arguments "  -initialdir \"$initial_dir\""

  } elseif { !$append &&
       $cur_val != ""
       && [file isdirectory $cur_val]
     } {
    append arguments "  -initialdir \"$cur_val\""
  }

  # Open browser
  set path [tk_getDirectory $arguments ]

  #-If pressed "cancel"
  if { $path == "" } {
    return
  }

  #-Append mode
  if {$append} {

    set new_val $cur_val

    if { $new_val != "" } {
      append new_val "; "
    }

    append new_val $path

    # Add separator for convenience!
    append new_val ";"

  #-Normal mode (replace)
  } else {
    set new_val $path
  }

  $wdg delete 0 end
  $wdg insert 0 $new_val

  if { $panel != "" } {
    Panel::panelDataChanged 1 $panel 
    Panel::panelDataModified 1 $panel

    if { $field != "" } {
      upvar #0 $panel theArr
      set theArr($field,mod) 1
      Widget::setFieldStatus $panel $field
    }
  }

}


# Open a file browser to set a file name
# into the entry widget "wdg"
#
proc Util::fillFileEntry { wdg args } {
  global Info ModelProperty

  if { ![winfo exists $wdg] } {
    return
  }

  set init_dir ""
  set init_file ""
  set absolute 1
  set append 0
  set drop_inc_path 0
  set keep_extension 1

  set panel ""
  set field ""
  set array ""

  set absolute_fld ""
  set append_fld ""
  set drop_inc_path_fld ""
  set init_dir_fld ""
  set init_file_fld ""
  set keep_extension_fld ""

  set absolute_fld ""
  set absolute_var ""
  set drop_inc_path_var ""
  set file_types ""
  set init_dir_var ""
  set init_file_var ""
  set keep_extension_var ""

  set result_fld ""
  set result_abs_fld ""

  set result_var ""
  set result_abs_var ""

  set args [split [join $args] "-"]

  foreach arg $args {

    set kwd [lindex $arg 0]
    set val [lindex $arg 1]

    switch $kwd {
      array       { set array $val }
      field       { set field $val }
      panel       { set panel $val }
      abs         { set absolute $val }
      absVar      { set absolute_var $val }
      absFld      { set absolute_fld $val }
      append      { set append $val }
      appendVar   { set append_var $val }
      appendFld   { set append_fld $val }
      dropInc     { set drop_inc_path $val }
      dropIncVar  { set drop_inc_path_var $val }
      dropIncFld  { set drop_inc_path_fld $val }
      extensions  { set file_types $val }
      initDir     { set init_dir $val }
      initDirVar  { set init_dir_var $val }
      initDirFld  { set init_dir_fld $val }
      initFile    { set init_file $val }
      initFileVar { set init_file_var $val }
      initFileFld { set init_file_fld $val }
      keepExt     { set keep_extension $val }
      keepExtVar  { set keep_extension_var $val }
      keepExtFld  { set keep_extension_fld $val }
      resultVar     { set result_var $val }
      resultFld     { set result_fld $val }
      resultAbsVar  { set result_abs_var $val }
      resultAbsFld  { set result_abs_fld $val }
    }
  }
  
  #--If array name given
  #
  if { $array != "" } {
    upvar #0 $array theArr
  }

  #--If result arguments given
  #
  if { $result_var != "" } {
    upvar $result_var result
  }

  if { $result_abs_var != "" } {
    upvar $result_abs_var result_abs
  }

  #--Set control variables (possibly from array)
  #  
  foreach vn {absolute drop_inc_path keep_extension init_dir init_file} {

    set fld_var [set $vn\_fld]
    set var_var [set $vn\_var]

    #-Control in an array variable
    if { $array != "" && $fld_var != "" } {
      set $vn $theArr($fld_var)

    #-Control in a variable
    } elseif { $var_var != "" } {
        set $vn [set $var_var]
    }
  }

  set title "Select file"
  set idir ""
  set ifile ""

  set inc_dirs [string trim $ModelProperty(INCLUDE_PATH) ";"]
  set inc_dirs [split $inc_dirs ";"]
  set inc_dir [lindex $inc_dirs 0]

  #--Clean current entry value
  set cur_val [$wdg get]
  set cur_val [string trim [string trim $cur_val] ";"]

  if { $cur_val == "" ||
       [Util::isSimpleFilename $cur_val]
     } {
    set cur_dir $inc_dir

  } else {
    set cur_dir [file dirname $cur_val]
  }

  set cur_val [file tail $cur_val]

  if { $cur_dir != "" &&
       [file isdirectory $cur_dir]
     } {
    set idir $cur_dir
  }

  if { $init_dir != "" } {
    set idir $init_dir
  }

  if { $init_file != "" } {
    set ifile $init_file
  }

  #--Open file browser
  set path [tk_getOpenFile -title $title -initialdir $idir -initialfile $ifile -parent $wdg -filetypes $file_types]

  #--If pressed "cancel"
  if { $path == "" } {
    return
  }

  set new_dir  [file dirname $path]
  set new_file [file tail $path]
  
  #--Pick results for possible result arguments
  #-Result to array
  #
  if { $array != "" } {
   
    if { $result_abs_fld != "" } {
      set theArr($result_abs_fld) $new_dir
    }
   
    if { $result_fld != "" } {
      set theArr($result_fld) $new_file
    }
  }

  #-Result to (possible) upvar-variable
  #
  set result_abs $new_dir
  set result $new_file

  #--If file extension should be dropped
  if { !$keep_extension } {
    set new_file [file rootname $new_file]
  }

  #--This is the default result
  set new_path [file join $new_dir $new_file]
  set new_val $new_path

  #-If file-path belong to include paths
  # we do not keep the path
  if { $drop_inc_path } {

    # If new-dir is not among include paths,
    # we use the full path
    #
    if { -1 == [List::fileListSearch $inc_dirs $new_dir] } {
      set new_val $new_path

    # Otherwise new value is just the file name
    #
    } else {
      set new_val $new_file
    }
  }
  
  #-If path should be dropped
  if { !$absolute } {
    set new_val $new_file
  }

  $wdg delete 0 end
  $wdg insert 0 $new_val

  if { $panel != "" } {
    Panel::panelDataChanged 1 $panel 
    Panel::panelDataModified 1 $panel 

    if { $field != "" } {
      upvar #0 $panel theArr
      set theArr($field,mod) 1
      Widget::setFieldStatus $panel $field
    }
  }
}


# Select proper coordinate label (X R P etc.)
# for the coordinate (1 2 3)
#
proc Util::getCoordinateLabel {coordinate} {
  global Coordinate

  if { $coordinate == "" || $coordinate < 1 || $coordinate > 3 } {
    return ""
  }

	if { [info exists Coordinate(COORDINATE_MAPPING)] } {
		set coord_index [lindex $Coordinate(COORDINATE_MAPPING) [incr coordinate -1] ]
	} else {
		set coord_index $coordinate
	}

  incr coord_index -1

  set labels [lindex $Coordinate(coordinateLabels) \
                     $Coordinate(coordinateLabelIndex)]

  return [lindex $labels $coord_index]
}



# Set currently valid flag names for each group
#
proc Util::updateFlagNameValues {} {
  global Info ModelFlags

  # Groups
  set draw_target { DRAW_TARGET_BODIES
                    DRAW_TARGET_SURFACES
                    DRAW_TARGET_EDGES }

  set select_method { SELECT_METHOD_SINGLE
                      SELECT_METHOD_ALL
                      SELECT_METHOD_BY_NEIGHBOR
                      SELECT_METHOD_BY_NORMAL
                      SELECT_METHOD_BY_PLANE
                      SELECT_METHOD_BY_BOX
                      SELECT_METHOD_BY_RECTANGLE }

  set select_mode { SELECT_MODE_TOGGLE
                    SELECT_MODE_EXTEND
                    SELECT_MODE_REDUCE }

  set grp_flag_vars  { draw_target select_method select_mode }
  set grp_names      { DRAW_TARGET SELECT_METHOD SELECT_MODE }

  #---Loop all groups
  foreach grp_flags $grp_flag_vars  grp_name $grp_names {

    set ModelFlags($grp_name) ""

    #--Loop all flags within the group
    foreach flag [set $grp_flags] {

      #-Find first flag values which is ON, and
      # make that flag the current
      #
      if { [info exists ModelFlags($flag)] &&
           $ModelFlags($flag)
         } {
        set ModelFlags($grp_name) $flag
        break
      }
    }
  }
 
#print "draw_target=$ModelFlags(DRAW_TARGET)"
#print "select_method=$ModelFlags(SELECT_METHOD)"
#print "select_mode=$ModelFlags(SELECT_MODE)"
}





# Set flag value in the Gui and in the model (cpp)
#
proc Util::setFlagValue { flag_group flag_name {value ""} } {
  global Info ModelFlags
  
  set old_value ModelFlags($flag_name,oldValue)

  #---If new value is not given, read it from
  #   the flag array where it is already set
  if {$value == ""} {
    set new_value $ModelFlags($flag_name)
  } else {
    set new_value $value
  }

  #---Reset first exlusive type of flags
  switch $flag_group {
    DRAW_TARGET {
      set ModelFlags(DRAW_TARGET) $flag_name
      set ModelFlags(DRAW_TARGET_BODIES) 0
      set ModelFlags(DRAW_TARGET_SURFACES) 0
      set ModelFlags(DRAW_TARGET_EDGES) 0
    }
    SELECT_METHOD {
      set ModelFlags(SELECT_METHOD) $flag_name
      set ModelFlags(SELECT_METHOD_SINGLE) 0
      set ModelFlags(SELECT_METHOD_ALL) 0
      set ModelFlags(SELECT_METHOD_BY_NEIGHBOR) 0
      set ModelFlags(SELECT_METHOD_BY_NORMAL) 0
      set ModelFlags(SELECT_METHOD_BY_PLANE) 0
      set ModelFlags(SELECT_METHOD_BY_BOX) 0
      set ModelFlags(SELECT_METHOD_BY_RECTANGLE) 0
   }
    SELECT_MODE {
      # You can't turn these completely off
      if {$new_value == $old_value} {
        return
      }
      set ModelFlags(SELECT_MODE) $flag_name
      set ModelFlags(SELECT_MODE_TOGGLE) 0
      set ModelFlags(SELECT_MODE_EXTEND) 0
      set ModelFlags(SELECT_MODE_REDUCE) 0
    }
  }

  #---Then set the updated value
  set ModelFlags($flag_name) $new_value
  set ModelFlags($flag_name,oldValue) $new_value

  #---Send to the model
  set sp $Info(cmdSeparator)
  set data "$flag_group $flag_name $new_value"

  Util::cpp_exec "setFlagValue" $data
}


proc Util::restoreFlagValue {globArray flag_grp {flag_name ""} } {
  global ModelFlags
  upvar #0 $globArray theArray 

  # Was flag value stored at group level or flag-name level
  if { $flag_name == "" } {
    set array_var $flag_grp
  } else {
    set array_var $flag_name
  }

  if { ![info exists theArray($array_var,cur)] ||
       ![info exists theArray($array_var,prev)]
     } {
    return
  }

  # Flag name given is stored array variable
  if { $flag_name == "" } {

    set flag_name $theArray($array_var,cur)

    # If user has turned off the flag which was set when panel was
    # open, do nothing!
    if { !$ModelFlags($flag_name) } {
      return
    }

    # Turn on the previous flag
    set flag_name $theArray($flag_grp,prev)

    Util::setFlagValue $flag_grp $flag_name 1

  # Flag name given as an argument and flag value stored in array variable
  } else {

    if { $ModelFlags($flag_name) != $theArray($array_var,cur) } {
      return
    }

    Util::setFlagValue $flag_grp $flag_name $theArray($array_var,prev)
  }

}


proc Util::updateGui { {only_idle_tasks 1} } {
  global Info

  if { $only_idle_tasks } {
    update idletasks

  } else {
    if { $Info(platform) == "windows" } {
      update
    } else {
      update idletasks
    }
  }
}


# Add a file path to list
#
proc Util::addToFileList { listname filename {only_latest 1} } {
  global Info
  upvar #0 $listname lname
  
  # Pick new file modified time
  if { [file exists $filename] } {
    set mt [file mtime $filename]
  } else {
    set mt ""
  }

  # Drop earlier instances if only the lates should be stored
  #  
  if { $only_latest } {
    
    set new_list ""
    
    # Loop old list
    foreach itm $lname {

      # Skip entry if names match and time is earlier
      #
      set n [lindex $itm 0]
      set t [lindex $itm 1]
      
      #-In Windows filenames are not case-sensitive
      if { $Info(platform) == "windows" } {
        if { [string equal -nocase $filename $n] && $mt >= $t } {
          continue
        }

      #-In Unix filenames are case-sensitive
      } else {
        if { [string equal $filename $n] && $mt >= $t } {
          continue
        }
      }
      
      # Old entry is stored normally
      #
      lappend new_list $itm
    }
    
    # This is now our list to add to
    #
    set lname $new_list
  }

  # Store new filename,mtime 
  #
  lappend lname [list $filename $mt]
}

# end ecif_tk_procsUtil.tcl
# *************************
