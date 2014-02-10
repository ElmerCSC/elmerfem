
proc math_CommandReturn {w l t} {
    global math_str math_file

    if { $math_str == "" } return

    if { $l != "" } {
        $l insert end $math_str
        history add "math \{ $math_str \}"

        set n [$l size]

        set k [@ $n-5]
        if { $k < 0 } { set k 0}
        $l yview $k

        $l select clear 0 end
        $l select set end end

        if { $t != "" } {
            $t insert end [math $math_str]
            $t insert end ">\n"

            set n [split [$t index end] .]
            set row [lindex $n 0]

            set h [$t cget -height]
            set row [@ $row-$h]

            $t yview $row
        }
    } else {
        uplevel 1 "math $math_str\n"
    }

    $w delete 0 end
}

proc math_CommandUpKey {w l} {
   if { $l != "" } {
       set n [$l curselection]
       if { $n == "" } {
           $l select clear 0 end
           $l select set end end
           set n [$l curselection]
       }

       set n [@ $n-1]
       if { $n < 0 } { set n 0 }

        $l select clear 0 end
        $l select set $n $n

        set k [@ $n-5]
        if { $k < 0 } { set k 0 }
        $l yview $k
   }

   $w delete 0 end
   if { $l != "" } { $w insert 0 [$l get $n] }
}

proc math_CommandDownKey {w l} {
   if { $l != "" } {
       set n [$l curselection]
       if { $n == "" } {
           $l select clear 0 end
           $l select set end end
           $l yview  end
           set n [$l curselection]
       }
       set n [@ $n+1]
       if { $n >= [$l size] } { set $n [@ [$l size]-1] }

       $l select clear 0 end
       $l select set $n $n

       set k [@ $n-5]
       if { $k < 0 } { set k 0 }
       $l yview $k
   }

   $w delete 0 end
   if { $l != "" } { $w insert 0 [$l get $n] }
}

proc math_CommandLeftKey {w} {
   set t [$w index insert]
   set t [@ $t-1]

   if { $t < 0 } { set t 0 }

   $w icursor $t
}

proc math_CommandRightKey {w} {
   set t [$w index insert]
   set t [@ $t+1]

   if {$t >= [$w index end]} {set t [$w index end]}

   $w icursor $t
}

proc math_CommandToEnd {w} { $w icursor end }

proc math_CommandToBegin {w} { $w icursor 0 }

proc math_CommandDelete {w} { $w delete [$w index insert] }

proc math_CommandDeleteEnd {w} { $w delete [$w index insert] end }

proc math_TextUpKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set row [@ $row-1]

   $w yview -pickplace $row
   $w mark set insert $row.$col
}

proc math_TextDownKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set row [@ $row+1]

   $w yview -pickplace $row
   $w mark set insert $row.$col
}

proc math_TextLeftKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set col [@ $col-1]
   $w mark set insert $row.$col
}

proc math_TextRightKey {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col [lindex $t 1]

   set col [@ $col+1]
   $w mark set insert $row.$col
}

proc math_TextHomeKey {w} {
    $w yview -pickplace 0
    $w mark set insert 0.0
}

proc math_TextEndKey {w} {
    $w yview -pickplace end
    $w mark set insert end
}

proc math_TextPageDownKey {w} {
    set h [$w cget -height]
    set h [@ $h-1]

    set t [split [$w index insert] .]

    set row [lindex $t 0]
    set col [lindex $t 1]

    set row [@ $row+$h]

    $w yview -pickplace $row
    $w mark set insert $row.$col
}

proc math_TextPageUpKey {w} {
    set h [$w cget -height]
    set h [@ $h-1]

    set t [split [$w index insert] .]

    set row [lindex $t 0]
    set col [lindex $t 1]

    set row [@ $row-$h]

    $w yview -pickplace $row
    $w mark set insert $row.$col
}

proc math_TextToBegin {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col 0

   $w mark set insert $row.$col
}

proc math_TextToEnd {w} {
   set t [split [$w index insert] .]

   set row [lindex $t 0]
   set col 0
   $w mark set here $row.$col

   set r [@ $row+1]
   $w mark set insert $r.$col

   set str [$w get here insert]

   set col [@ [string length $str]-1]

   $w mark set insert $row.$col
   $w mark unset here
}

proc math_TextDelete {w} { $w delete insert }

proc math_Nil {} {}

proc math { {str ""} args } {
   if { $str == "" } { return [math_win] }

   return [c_MathCommand "$str $args"]
}

proc math_win {} {
    global math_str math_cmd math_log math_out math_file

    if { [winfo exists .math] } {
        wm iconify .math
        wm deiconify .math
        return
    }

    toplevel .math
    place_window .math

    frame .math.lbox -bd 5 -bg lightblue
    frame .math.tbox -bd 5 -bg lightblue
    frame .math.cbox
    label .math.cbox.lab -text "math: "
    entry .math.cbox.cmd -relief sunken -textvariable math_str

    text .math.tbox.text -yscroll ".math.tbox.scroll set" -wrap none
    scrollbar .math.tbox.scroll -orient vertical -command ".math.tbox.text yview"
    pack .math.tbox.scroll -side left -fill y
    pack .math.tbox.text -side left -fill both -expand 1 

    listbox .math.lbox.list -yscroll ".math.lbox.scroll set"
    scrollbar .math.lbox.scroll -orient vertical -command ".math.lbox.list yview"
    pack .math.lbox.scroll -side left -fill y
    pack .math.lbox.list -side left -fill both -expand 1 

    set math_out .math.tbox.text
    set math_log .math.lbox.list
    set math_cmd .math.cbox.cmd

    bind $math_log <Button-1> {
        $math_cmd delete 0 end
        set cur [$math_log curselection]
        if { $cur != "" } {
            $math_cmd insert 0 [$math_log get $cur]
        }
    }

    bind $math_log <Double-Button-1> {
        $math_cmd delete 0 end

        set cur [$math_log curselection]
        if { $cur != "" } {
            $math_cmd insert 0 [$math_log get $cur]
            math_CommandReturn $math_cmd $math_log $math_out
        }
    }

    bind $math_cmd <Return>    {math_CommandReturn %W $math_log $math_out}

    bind $math_cmd <Up>        {math_CommandUpKey %W $math_log}
    bind $math_cmd <Down>      {math_CommandDownKey %W $math_log}

    bind $math_cmd <Control-p> {math_CommandUpKey %W $math_log}
    bind $math_cmd <Control-n> {math_CommandDownKey %W $math_log}

    bind Text <Escape> {Nil}
    bind Text <Escape>v     {math_TextPageUpKey %W}
    bind Text <Control-v>   {math_TextPageDownKey %W}

    pack .math.lbox -side top -fill both -expand 1
    pack .math.tbox -fill both -expand 1

    pack .math.cbox.lab -side left
    pack .math.cbox.cmd -side right -expand 1 -fill both
    pack .math.cbox -fill x -side bottom

    $math_out configure -setgrid 1 -width 40 -height 10
    $math_log configure -setgrid 1 -width 40 -height 10 
    $math_cmd configure -width 40

    $math_out insert end "\n\nMATC COMMAND WINDOW\n"
    $math_out insert end "\n"

    set math_str ""

    focus $math_cmd
}

