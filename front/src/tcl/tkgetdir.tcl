#########################################################
# Directory Selector TCL version 1.2
#
# Daniel Roche, <dan@lectra.com>
# thanks to Cyrille Artho <cartho@netlink.ch> for the 'saving pwd fix'
#
#########################################################
 

#########################################################
# 
# tk_getDirectory [option value ...]
#
#  options are :
#   [-initialdir dir]     display in dir
#   [-title string]       make string title of dialog window
#   [-ok string]          make string the label of OK button
#   [-open string]        make string the label of OPEN button
#   [-cancel string]      make string the label of CANCEL button
#   [-msg1 string]        make string the label of the first directory message
#   [-msg2 string]        make string the label of the second directory message
#
#########################################################

proc tk_getDirectory {args} {
    variable fini
    global tcl_platform drives

    #
    # arguments
    #
    set _titre "Directory Selector"
    set _ldir Directory:
    set _ldnam "Directory Name:"
    set _open Ok
    set _expand Open
    set _cancel Cancel
    
	set args [join $args]
	
    set ind 0
    set max [llength $args]
    set pwd [pwd]
    while { $ind < $max } {
	switch -exact -- [lindex $args $ind] {
	    "-initialdir" {
		incr ind
		cd [lindex $args $ind]
		incr ind
	    }
	    "-title" {
		incr ind
		set _titre [lindex $args $ind]
		incr ind
	    }
	    "-ok" {
		incr ind
		set _open [lindex $args $ind]
		incr ind
	    }
	    "-open" {
		incr ind
		set _expand [lindex $args $ind]
		incr ind
	    }
	    "-cancel" {
		incr ind
		set _cancel [lindex $args $ind]
		incr ind
	    }
	    "-msg1" {
		incr ind
		set _ldir [lindex $args $ind]
		incr ind
	    }
	    "-msg2" {
		incr ind
		set _ldnam [lindex $args $ind]
		incr ind
	    }
	    default {
		#puts "unknown option [lindex $args $ind]"
		set opt [lindex $args $ind]
		print "Unknown option:\n args=$args \n len=$max \n ind=$ind \n opt=$opt"
		return ""
	    }
	}
    }
    
    #
    # variables et data
    #
    set fini 0
    
    image create bitmap b_up -data "
    #define up_width 31
    #define up_height 23
    static unsigned char up_bits[] = {
	0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
	0x00, 0x00, 0x00, 0x80, 0x00, 0x3f, 0x00, 0x80, 0x80, 0x40, 0x00, 0x80,
	0x40, 0x80, 0x00, 0x80, 0xe0, 0xff, 0xff, 0x83, 0x20, 0x00, 0x00, 0x82,
	0x20, 0x04, 0x00, 0x82, 0x20, 0x0e, 0x00, 0x82, 0x20, 0x1f, 0x00, 0x82,
	0x20, 0x04, 0x00, 0x82, 0x20, 0x04, 0x00, 0x82, 0x20, 0x04, 0x00, 0x82,
	0x20, 0xfc, 0x0f, 0x82, 0x20, 0x00, 0x00, 0x82, 0x20, 0x00, 0x00, 0x82,
	0xe0, 0xff, 0xff, 0x83, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
	0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80};"

    image create bitmap b_dir -background #ffff80 -data "
    #define dir_width 17
    #define dir_height 16
    static unsigned char dir_bits[] = {
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe0, 0x01, 0x00, 0x10, 0x02, 0x00,
	0x08, 0x04, 0x00, 0xfc, 0x7f, 0x00, 0x04, 0x40, 0x00, 0x04, 0x40, 0x00,
	0x04, 0x40, 0x00, 0x04, 0x40, 0x00, 0x04, 0x40, 0x00, 0x04, 0x40, 0x00,
	0x04, 0x40, 0x00, 0xfc, 0x7f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};" \
		-maskdata "
    #define dirm_width 17
    #define dirm_height 16
    static unsigned char dirm_bits[] = {
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xe0, 0x01, 0x00, 0xf0, 0x03, 0x00,
	0xf8, 0x07, 0x00, 0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00,
	0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00,
	0xfc, 0x7f, 0x00, 0xfc, 0x7f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};"
		
    switch -exact $tcl_platform(platform) {
	unix {
	    font create myfont -family lucida -size 12 -weight bold
	}
	windows {
	    font create myfont -family courier -size 10
	}
    }

    #
    # widgets
    #
    toplevel .dirsel
    grab set .dirsel
    wm geometry .dirsel 450x400
    wm title .dirsel $_titre

    frame .dirsel.f1 -relief flat -borderwidth 0
    frame .dirsel.f2 -relief sunken -borderwidth 2 
    frame .dirsel.f3 -relief flat -borderwidth 0
    frame .dirsel.f4 -relief flat -borderwidth 0
    
    pack .dirsel.f1 -fill x
    pack .dirsel.f2 -fill both -expand 1 -padx 6 -pady 6
    pack .dirsel.f3 -fill x
    pack .dirsel.f4 -fill x
    
    label .dirsel.f1.lab -text $_ldir
    menubutton .dirsel.f1.dir -relief raised -indicatoron 1 -anchor w \
	    -menu .dirsel.f1.dir.m
    menu .dirsel.f1.dir.m -tearoff 0
    button .dirsel.f1.up -image b_up -command UpDir
    
    pack .dirsel.f1.up -side right -padx 4 -pady 4
    pack .dirsel.f1.lab -side left -padx 4 -pady 4
    pack .dirsel.f1.dir -side right -padx 4 -pady 4 -fill x -expand 1
    
    canvas .dirsel.f2.cv -borderwidth 0 -yscrollcommand ".dirsel.f2.sb set"
    if ![string compare $tcl_platform(platform) windows] {
	    .dirsel.f2.cv configure -background white
    }

    scrollbar .dirsel.f2.sb -command ".dirsel.f2.cv yview"
    set scw 16
    place .dirsel.f2.cv -x 0 -relwidth 1.0 -width [expr -$scw ] -y 0 \
	    -relheight 1.0
    place .dirsel.f2.sb -relx 1.0 -x [expr -$scw ] -width $scw -y 0 \
	    -relheight 1.0
    unset scw
    
    .dirsel.f2.cv bind TXT <Any-Enter> EnterItem
    .dirsel.f2.cv bind TXT <Any-Leave> LeaveItem
    .dirsel.f2.cv bind TXT <Any-Button> ClickItem
    .dirsel.f2.cv bind TXT <Double-Button> DoubleClickItem
    .dirsel.f2.cv bind IMG <Any-Enter> EnterItem
    .dirsel.f2.cv bind IMG <Any-Leave> LeaveItem
    .dirsel.f2.cv bind IMG <Any-Button> ClickItem
    .dirsel.f2.cv bind IMG <Double-Button> DoubleClickItem
    
    label .dirsel.f3.lnam -text $_ldnam
    entry .dirsel.f3.chosen -takefocus 0
    pack .dirsel.f3.lnam -side left -padx 4 -pady 4
    pack .dirsel.f3.chosen -side right -fill x -expand 1 -padx 4 -pady 4
    
    button .dirsel.f4.open -text $_open -command { 
	set tmp [.dirsel.f3.chosen get]
	if [ string length $tmp ] {
	    set fini 1 
	}
    }
    button .dirsel.f4.expand -text $_expand -command DownDir
    button .dirsel.f4.cancel -text $_cancel -command { 
	set fini -1 
    }
    
    pack .dirsel.f4.open .dirsel.f4.expand -side left -padx 10 -pady 4
    pack .dirsel.f4.cancel -side right -padx 10 -pady 4
    
    #
    # realwork
    #
    ShowDir [pwd]
    
    #
    # wait user
    #
    tkwait variable fini

    if { $fini == 1 } {
	set curdir [.dirsel.f1.dir cget -text]
	set nnam [.dirsel.f3.chosen get]
	set retval [ file join $curdir $nnam ]
    } else {
	set retval ""
    }
    
    font delete myfont 
    destroy .dirsel
    unset drives fini
    cd $pwd
    return $retval
}

proc ShowDir {curdir} {

    global tcl_platform
    variable drives 
    
    cd $curdir
    .dirsel.f1.dir configure -text $curdir
    
    set hi1 [font metrics myfont -linespace]
    set hi2 [image height b_dir]
    if { $hi1 > $hi2 } {
	set hi $hi1
    } else {
	set hi $hi2
    }
    set wi1 [image width b_dir]
    incr wi1 4
    set wi2 [winfo width .dirsel.f2.cv]
    
    set lidir [list]
    foreach file [ glob -nocomplain * ] {
	if [ file isdirectory [string trim $file "~"] ] { 
	    lappend lidir $file
	}
    }
    set sldir [lsort $lidir]
    
    .dirsel.f2.cv delete all
    set ind 0
    foreach file $sldir {
	if [ file isdirectory $file ] { 
	    .dirsel.f2.cv create image 2 [expr $ind * $hi] \
		    -anchor nw -image b_dir -tags IMG
	    .dirsel.f2.cv create text $wi1 [expr $ind * $hi] \
		    -anchor nw -text $file -font myfont -tags TXT
	    set ind [ expr $ind + 1 ]
	}
    }

    set ha [expr $ind * $hi]
    .dirsel.f2.cv configure -scrollregion [list 0 0 $wi2 $ha]
    
    set curlst [file split $curdir]
    set nbr [llength $curlst]
    
    .dirsel.f1.dir.m delete 0 last
    incr nbr -2
    for {set ind $nbr} {$ind >= 0} {incr ind -1} {
	set tmplst [ lrange $curlst 0 $ind] 
	set tmpdir [ eval file join $tmplst] 
	.dirsel.f1.dir.m add command -label $tmpdir -command "ShowDir {$tmpdir}"
    }
    if {[info exist drives] == 0} {
	set drives [file volume]
    }
    if ![string compare $tcl_platform(platform) windows] {
	foreach drive $drives {
	    .dirsel.f1.dir.m add command -label $drive -command "ShowDir {$drive}"
	}
    }
    
}

proc UpDir {} {
    set curdir [.dirsel.f1.dir cget -text]
    set curlst [file split $curdir]
    
    set nbr [llength $curlst]
    if { $nbr < 2 } {
	return
    }
    set tmp [expr $nbr - 2]
    
    set newlst [ lrange $curlst 0 $tmp ]
    set newdir [ eval file join $newlst ]
    
    .dirsel.f3.chosen delete 0 end
    ShowDir $newdir
}

proc DownDir {} {
 set curdir [.dirsel.f1.dir cget -text]
 set nnam [.dirsel.f3.chosen get]

 set newdir [ file join $curdir $nnam ]

 .dirsel.f3.chosen delete 0 end
 ShowDir $newdir
}

proc EnterItem {} {
 global tcl_platform

 set id [.dirsel.f2.cv find withtag current]
 set wt [.dirsel.f2.cv itemcget $id -tags]
 if {[lsearch -exact $wt IMG] >= 0} {
  set id [.dirsel.f2.cv find above $id]
 }
 if [string compare $tcl_platform(platform) windows] {
   set cocol #00FF00
 } else {
   set cocol #0000FF
 }
 .dirsel.f2.cv itemconfigure $id -fill $cocol
}

proc LeaveItem {} {
 set id [.dirsel.f2.cv find withtag current]
 set wt [.dirsel.f2.cv itemcget $id -tags]
 if {[lsearch -exact $wt IMG] >= 0} {
  set id [.dirsel.f2.cv find above $id]
 }
 .dirsel.f2.cv itemconfigure $id -fill black
}

proc ClickItem {} {
 .dirsel.f2.cv delete BOX
 set id [.dirsel.f2.cv find withtag current]
 set wt [.dirsel.f2.cv itemcget $id -tags]
 if {[lsearch -exact $wt IMG] >= 0} {
  set id [.dirsel.f2.cv find above $id]
 }
 set bxr [.dirsel.f2.cv bbox $id]
 eval .dirsel.f2.cv create rectangle $bxr -fill #a2a2ff -outline #a2a2ff -tags BOX
 .dirsel.f2.cv lower BOX
 set nam [.dirsel.f2.cv itemcget $id -text]
 .dirsel.f3.chosen delete 0 end
 .dirsel.f3.chosen insert 0 $nam
}

proc DoubleClickItem {} {
 set id [.dirsel.f2.cv find withtag current]
 DownDir
}


