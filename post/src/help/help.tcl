#
# here is a sample html viewer to demonstrate the library usage
# Copyright (c) 1995 by Sun Microsystems
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#
# This REQUIRES Tk4.0 -- make sure "wish" on the next line is a 4.0 version
# The next line is a TK comment, but a shell command
#

#if {$tk_version < 4.0 || [regexp {b[123]} $tk_patchLevel] } {
#	puts stderr "This library requires TK4.0, this is only $tk_version, \
#			patchlevel $tk_patchLevel"
#	exit 1
#}
#
#if {[catch {array get env *}]} {
#	puts stderr "This library requires tcl7.4, this version is too old!"
#	exit 1
#}
#puts stderr "Starting sample HTML viewer..."

source $ELMER_POST_HOME/help/http.tcl
source $ELMER_POST_HOME/help/html_library-0.3/html_library.tcl

array set helpTypeAction {
   {}       helpRender
  .html     helpRender
}
# .mpg     "exec mpeg_play"
# .gif     "exec xv"
# .jpg     "exec xv"
  


# Sample hypertext link callback routine - should be replaced by app
# This proc is called once for each <A> tag.
# Applications can overwrite this procedure, as required, or
# replace the HMevents array
#   win:   The name of the text widget to render into
#   href:  The HREF link for this <a> tag.

array set HMevents {
	Enter	{-borderwidth 2 -relief raised}
	Leave	{-borderwidth 2 -relief flat}
	1		{-borderwidth 2 -relief sunken}
	ButtonRelease-1	{-borderwidth 2 -relief raised}
}

# We need to escape any %'s in the href tag name so the bind command
# doesn't try to substitute them.

proc HMlink_setup {win href} {
	global HMevents helpMessage
	regsub -all {%} $href {%%} href2
	foreach i [array names HMevents] {
             if { $i == "Enter" } {
		eval {$win tag bind  L:$href <$i> "set helpMessage $href; $win tag configure \{L:$href2\} $HMevents($i)"}
             } else {
                if { $i == "Leave" } { 
  	          eval {$win tag bind  L:$href <$i> "set helpMessage \{\}; $win tag configure \{L:$href2\} $HMevents($i)"}
                } else {
           	  eval {$win tag bind  L:$href <$i> "$win tag configure \{L:$href2\} $HMevents($i)"}
                }
             }
	}
}


# construct a simple user interface

proc helpSetup {} {

	frame .helpWindow.frame

	menubutton .helpWindow.file -relief raised -bd 2 -text File...    -menu .helpWindow.file.menu
	menubutton .helpWindow.options -relief raised -bd 2 -text Options... -menu .helpWindow.options.menu

#	button .helpWindow.quit  -command { destroy .helpWindow }  -text "Quit"
#        bind .helpWindow <Control-q> "destroy .helpWindow"

	label .helpWindow.status -textvariable helpRunning -width 6 -relief ridge \
			-bd 2 -padx 9 -pady 3

	label .helpWindow.msg -textvariable helpMessage -relief sunken

	frame .helpWindow.dspl
	scrollbar .helpWindow.dspl.scrollbar  -command ".helpWindow.dspl.text yview"  -orient v
	option add *Text.height 40 startup
	option add *Text.width 80 startup
	text .helpWindow.dspl.text  -yscrollcommand ".helpWindow.dspl.scrollbar set" -padx 3 -pady 3 -takefocus 0

	button .helpWindow.back -text "Back" -command \
	      { 
		set helpCurrent [expr $helpCurrent-1]; if { $helpCurrent < 0 } { incr helpCurrent };
		helpRender $helpHistory($helpCurrent);
		set helpCurrent [expr $helpCurrent-1]; if { $helpCurrent < 0 } { incr helpCurrent };
	      }

	pack .helpWindow.frame -side top -expand 1 -fill x
	pack .helpWindow.file .helpWindow.options -in .helpWindow.frame -side left
	pack .helpWindow.status .helpWindow.back -in .helpWindow.frame -side right

	frame .helpWindow.url

	label .helpWindow.url.label  -text "Url: "
	entry .helpWindow.url.entry -textvariable helpUrl

	pack .helpWindow.url.label -side left
	pack .helpWindow.url.entry -side left -expand 1 -fill x
	pack .helpWindow.url -expand 1 -side left -side top -fill x

	pack .helpWindow.msg -side top -expand 1 -fill x

	pack .helpWindow.dspl.scrollbar -side left -expand 0 -fill y
	pack .helpWindow.dspl.text -side left -fill both -expand 1
	pack .helpWindow.dspl -side top -fill both -expand 1

	# set up some sample keyboard bindings for the text widget
	bind .helpWindow.url.entry <Return> { helpRender $helpUrl }
	bind all <End> {.helpWindow.dspl.text yview end}
	bind all <Home> {.helpWindow.dspl.text yview 0.0}
	bind all <Next> {.helpWindow.dspl.text yview scroll 1 page}
	bind all <Prior> {.helpWindow.dspl.text yview scroll -1 page}

	# I'm constantly being criticized for never using menus.
	# so here's a menu.  So there.
	menu .helpWindow.file.menu
	.helpWindow.file.menu add command -label "Quit...Ctrl+Q" -command { destroy .helpWindow }
        bind .helpWindow <Control-q> { destroy .helpWindow }

	menu .helpWindow.options.menu
	.helpWindow.options.menu add command -label "Options menu"
	.helpWindow.options.menu add separator
	.helpWindow.options.menu add command -label "Font size" -foreground red 

	.helpWindow.options.menu add radiobutton -label small -value 0   -variable helpFontSize \
		-command {HMset_state .helpWindow.dspl.text -size $helpFontSize; helpRender $helpUrl}

	.helpWindow.options.menu add radiobutton -label medium -value 4  -variable helpFontSize \
		-command {HMset_state .helpWindow.dspl.text -size $helpFontSize; helpRender $helpUrl}

	.helpWindow.options.menu add radiobutton -label large -value 12  -variable helpFontSize \
		-command {HMset_state .helpWindow.dspl.text -size $helpFontSize; helpRender $helpUrl}

	.helpWindow.options.menu add separator
	.helpWindow.options.menu add command -label "Indent level" -foreground red

	.helpWindow.options.menu add radiobutton -label small -value 0.6 -variable helpIndent \
		-command {HMset_indent .helpWindow.dspl.text $helpIndent}

	.helpWindow.options.menu add radiobutton -label medium -value 1.2 -variable helpIndent \
		-command {HMset_indent .helpWindow.dspl.text $helpIndent}

	.helpWindow.options.menu add radiobutton -label large -value 2.4 -variable helpIndent \
		-command {HMset_indent .helpWindow.dspl.text $helpIndent}

}

# Go helpRender a page.  We have to make sure we don't render one page while
# still helpRendering the previous one.  If we get here from a recursive 
# invocation of the event loop, cancel whatever we were helpRendering when
# we were called.
# If we have a fragment name, try to go there.

proc helpRender { file } {
	global HM.text helpUrl
	global helpRunning helpMessage helpHistory helpCurrent

	incr helpCurrent
	set helpHistory($helpCurrent) $file

	set fragment ""
	regexp {([^#]*)#(.+)} $file dummy file fragment
	if {$file == "" && $fragment != ""} {
		HMgoto .helpWindow.dspl.text $fragment
		return
	}
	HMreset_win .helpWindow.dspl.text
	set helpRunning Busy
	set helpMessage "Displaying $file"
	update idletasks
	set helpUrl $file
	if { $fragment != "" } {
		HMgoto .helpWindow.dspl.text $fragment
	}
	HMparse_html [helpGetHtml $file] {HMrender .helpWindow.dspl.text}
	set helpRunning Ready
	HMset_state .helpWindow.dspl.text -stop 1       ;# stop helpRendering previous page if busy
	set helpMessage ""
}

#
# if html is at the end of http-request this will be called by http_get.
#
proc httpCallback { token } {
    upvar #0 $token state
    global helpText

    set helpText $state(body)
}

# given a file name, return its html, or invent some html if the file can't
# be opened.

proc helpGetHtml { file } {
      global helpHome helpText

      if { [string match http:* $file] } {
	 set token [http_get $file -command httpCallback] 
	 http_wait $token
	 return $helpText
      } else {
	if {[catch {set fd [open $file]} msg]} {
		return "<title>Bad file $file</title>
			<h1>Error reading $file</h1><p>
			$msg<hr>
			<a href=$helpHome>Back</a>"
	}
	set result [read $fd]
	close $fd
	return $result
     }
}

proc helpGetImage { file } {
     set img_file hlp_tmp_img

     set img [open $img_file w]
     fconfigure $img -translation binary
     set token [http_get $file -command httpCallback -channel $img]
     http_wait $token
     close $img

     return $img_file
}


# Override the library link-callback routine for the sample app.
# It only handles the simple cases.

proc HMlink_callback {win href} {
	global helpUrl helpHome helpTypeAction

	if { [string match #* $href] } {
		helpRender $href
		return
	}

        if { [string last ? $href] > 0 } {
          set href [string range $href 0 [expr [string last ? $href]-1]]
        }

	if { [string match /* $href] } {
	   if { [string match http://* $helpUrl] } {
		set tmp [string range $helpUrl 7 end]
		set helpUrl http:/
		set helpUrl $helpUrl/[string range $tmp 0 [string first / $tmp]]
		set helpUrl $helpUrl/$href
	   } else {
		set helpUrl [file dirname $helpHome]/$href
	   }
	} else {
	   if { [string match *:* $href] } {
	       set helpUrl $href
	   } else {
	       set helpUrl [string range $helpUrl 0 [string last / $helpUrl]]$href
	   }
	}

        set extension [file extension $helpUrl]

        if { ![string match *#* $helpUrl] } {
           if { ![string match "*$extension*" [join [array names helpTypeAction] " "]] } {

              set name [tk_getSaveFile -parent .helpWindow -defaultextension [file extension $helpUrl]]

              if { $name != "" } {
                 set file [open $name w]
                 fconfigure $file -translation binary
                 set token [http_get $helpUrl -command httpCallback -channel $file]
                 http_wait $token
                 close $file
              }
              return
           }
        }

        update

	if { $extension == "" } {
	   if { [string last / $helpUrl] != [expr [string length $helpUrl]-1] } { set helpUrl $helpUrl/ }
   	   helpRender $helpUrl
	} else {
#           eval "$helpTypeAction($extension) $helpUrl"
           helpRender $helpUrl
        }

}

# Supply an image callback function
# Read in an image if we don't already have one
# callback to library for display

proc HMset_image {win handle src} {
	global helpUrl helpMessage helpHome

	if {[string match /* $src]} {
	      if { [string match http://* $helpUrl] } {
		set tmp [string range $helpUrl 7 end]
		set image http:/
		set image $image/[string range $tmp 0 [string first / $tmp]]
		set image $image/$src
	      } else {
		set image [file dirname $helpHome]/$src
	      }
	} else {
	   if {[string match *:* $src] } {
		   set image $src
	   } else {
		  set image [string range $helpUrl 0 [string last "/" $helpUrl]]$src
	   }
	}

	set helpMessage "Fetching image: $image"
	update

	if { [string first " $image " " [image names] "] >= 0 } {
	   HMgot_image $handle $image
	} else {
		set type photo

		if { [string match http:* $image] } {
		    set file [helpGetImage $image]
		} else {
		    set file $image
		}

	       if { [file extension $image] == ".bmp" } {set type bitmap}

	       image create $type $image -file $file
	       HMgot_image $handle $image
	}
}

# Handle base tags.  This breaks if more than 1 base tag is in the document

proc HMtag_base {win param text} {
	global helpUrl
	upvar #0 HM$win var
	HMextract_param $param href helpUrl
}

# downloading fonts can take a long time.  We'll override the default
# font-setting routine to permit better user feedback on fonts.  We'll
# keep our own list of installed fonts on the side, to guess when delays
# are likely

proc HMset_font {win tag font} {
	global helpMessage Fonts
	if {![info exists Fonts($font)]} {
		set Fonts($font) 1
		.helpWindow.msg configure -fg blue
		set helpMessage "downloading font $font"
		update
	}
	.helpWindow.msg configure -fg black
	set helpMessage ""
	catch {$win tag configure $tag -font $font} helpMessage
}

# Lets invent a new HTML tag, just for fun.
# Change the color of the text. Use html tags of the form:
# <color value=blue> ... </color>
# We can invent a new tag for the display stack.  If it starts with "T"
# it will automatically get mapped directly to a text widget tag.

proc HMtag_color {win param text} {
	upvar #0 HM$win var
	set value bad_color
	HMextract_param $param value
	$win tag configure $value -foreground $value
	HMstack $win "" "Tcolor $value"
}

proc HMtag_/color {win param text} {
	upvar #0 HM$win var
	HMstack $win / "Tcolor {}"
}

# Add a font size manipulation primitive, so we can use this sample program
# for on-line presentations.  sizes prefixed with + or - are relative.
#  <font size=[+-]3>  ..... </font>.  Note that this is not the same as
# Netscape's <font> tag.

proc HMtag_font {win param text} {
	upvar #0 HM$win var
	set size 0; set sign ""
	HMextract_param $param size
	regexp {([+-])? *([0-9]+)} $size dummy sign size
	if {$sign != ""} {
		set size [expr [lindex $var(size) end] $sign $size]
	}
	HMstack $win {} "size $size"
}

# This version is closer to what Netscape does

proc HMtag_font {win param text} {
	upvar #0 HM$win var
	set size 0; set sign ""
	HMextract_param $param size
	regexp {([+-])? *([0-9]+)} $size dummy sign size
	if {$sign != ""} {
		set size [expr [lindex $var(size) end] $sign  $size*2]
		HMstack $win {} "size $size"
	} else {
		HMstack $win {} "size [expr 10 + 2 * $size]"
	}
}

proc HMtag_/font {win param text} {
	upvar #0 HM$win var
	HMstack $win / "size {}"
}


proc help { { address "index.html" } } {
 
     global ELMER_POST_HOME helpHome helpHistory helpCurrent

     if { [winfo exists .helpWindow] } {
	wm iconify .helpWindow
	wm deiconify .helpWindow
     } else {
	toplevel .helpWindow
	wm title .helpWindow "ELMER POST HELP"
    }


  # set initial values
  set helpFontSize 0                            ;# font size adjustment
  set helpIndent 1.2                            ;# tab spacing (cm)
  if { [string match http:* $address] } {
     set helpHome $address
  } else {
     set helpHome $ELMER_POST_HOME/help/$address        ;# home document
  }
  set helpUrl $helpHome                         ;# current file
  set helpRunning Busy                          ;# page status
  set helpMessage ""                            ;# message line

  # make the interface and render the home page

  set helpCurrent -1

   catch helpSetup                              ;# the catch lets us re-source this file

  HMinit_win   .helpWindow.dspl.text
  HMset_state  .helpWindow.dspl.text -size $helpFontSize
  HMset_indent .helpWindow.dspl.text $helpIndent

  helpRender $helpHome
}
