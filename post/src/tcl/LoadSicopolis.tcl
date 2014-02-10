proc LoadSicopolis { } {
    wm geometry . -200+200 
    toplevel .sico
    bind .sico <Return> {if {[file exists "$SicoID.log" ]} { if {[file exists "$SicoID$sliceno.erg" ]} {SicoRead $SicoID $sliceno $writeASCII} else {echo "File $SicoID$sliceno.erg does not exist"}} else { echo "File $SicoID.log does not exist"} }


  label .sico.top -text "\nSICOPOLIS input:\n" -font "Helvetica-Bold 12"
  pack  .sico.top -side top

  frame .sico.frame1 -relief sunken
    label  .sico.frame1.runnamel -text "run id (5 characters)"
    entry .sico.frame1.runname -width 5 -textvariable SicoID
    label  .sico.frame1.slicel -text "slice file number (2 digits)"
    entry .sico.frame1.slice -width 2 -textvariable  sliceno
    pack .sico.frame1.runnamel .sico.frame1.runname -side top
    pack .sico.frame1.slicel .sico.frame1.slice -side top
  pack .sico.frame1 -side top
  
  frame .sico.frame2 -relief sunken
    checkbutton .sico.frame1.ascii -text "Output of results in ASCII-files" \
      -variable writeASCII
    pack .sico.frame1.ascii  -side top -expand 1 -fill x
  pack .sico.frame2 -side bottom


    button .sico.execute -text "Read File" \
      -command {if {[file exists "$SicoID.log" ]} { if {[file exists "$SicoID$sliceno.erg" ]} {SicoRead $SicoID $sliceno $writeASCII}\
							else {echo "File $SicoID$sliceno.erg does not exist"}} \
							 else { echo "File $SicoID.log does not exist"} }

  pack .sico.execute -side left
  button .sico.close -text "Close" -command "destroy .sico"
  pack .sico.close -side right
}

proc SicoRead { id slice wascii } {
    global ELMER_POST_HOME

    # open shell script file that runs sico2elmer-process
    if { [file exists ".RUNLOADSICO" ]} {
	file delete -force ".RUNLOADSICO"
    }    
    set  sicorunloadfile [open ".RUNLOADSICO" {WRONLY CREAT}]
    # write orders to shell script file that runs sico2elmer-process
    puts $sicorunloadfile "#!/bin/csh"
    puts $sicorunloadfile "$ELMER_POST_HOME/../../bin/sico2elmer << EOF >&! sico2elmer.log"
    puts $sicorunloadfile "$id"
    puts $sicorunloadfile "1"
    puts $sicorunloadfile "$slice"
    echo "Sico2elmer: Writing grid file $id$slice.ep";
    puts $sicorunloadfile "2"
    echo "Sico2elmer: Writing data file $id$slice.dat";
    puts $sicorunloadfile "3"
    if { $wascii } {
	echo "Sico2elmer: Writing ascii-data files $id$slice\_2d\.asc"
	echo "and $id$slice\_3d\.asc"
	puts $sicorunloadfile "5"
    }
    puts $sicorunloadfile "0" 
    puts $sicorunloadfile "EOF"
    # close shell script file that runs sico2elmer-process
    close $sicorunloadfile    
    file attributes ".RUNLOADSICO" -permissions 01755
    # run sico2elmer-process
    echo "Sico2elmer: running Sico2elmer; output in sico2elmer.log"
    $ ./.RUNLOADSICO; 
    # load grid
    readfile "$id$slice.ep"; 
    echo "Read grid file $id$slice.ep"; 
    # load data
    math source("$ELMER_POST_HOME/tcl/loadsingle");
    math sicoreadfile("$id",$slice,0);
    echo "Read data file $id$slice.dat"; 
}
