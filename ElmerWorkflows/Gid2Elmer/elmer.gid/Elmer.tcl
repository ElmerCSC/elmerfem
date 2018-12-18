proc InitGIDProject {dir} {

      ##set materials  [.central.s info materials]
      ##set conditions [.central.s info conditions ovpnt]

      set conditions [GiD_Info conditions ovpnt]
      set materials [GiD_Info materials]

      ## ----- lines for 2D
      ##GiD_Set ForceMeshEntities 2


      ## ----- Surfaces for 3D
      GiD_Set ForceMeshEntities 4
      set mesh_options [GiD_Set ForceMeshEntities]

      CreateWindow $dir $materials $conditions $mesh_options

}


proc CreateWindow {dir mat cond mes} {

    set w .gid.win_elmer

    InitWindow $w "ELMER3D.TCL" Elmer "" \
	    "" 1

    frame $w.top
    label $w.top.title_text -text " Problem type: Elmer_3D "
   
    frame $w.information  -relief ridge -bd 2 
    label $w.information.path        -text " Problem Type path: $dir "
    label $w.information.materials   -text " Avalaible materials: $mat"
    label $w.information.conditions  -text " Avalaible conditions: $cond"
    ## ----- lines for 2D
    ##label $w.information.mesh_options  -text " Mesh always by default option : lines($mes) set"
    ## ----- Surfaces for 3D
    label $w.information.mesh_options  -text " Mesh always by default option : surfaces($mes) set"

    frame $w.bottom
    button $w.bottom.start -text "CONTINUE" \
           -height 1 -width 14 -command "destroy $w"
 
    pack $w.top.title_text -pady 10
    pack $w.information.path $w.information.materials \
         $w.information.conditions $w.information.mesh_options -side top -anchor w

    pack $w.bottom.start  -side left -anchor center
    pack $w.top 
    pack $w.information -expand yes -fill both
    pack $w.bottom -side top -padx 6 -pady 10 -ipady 2
}

