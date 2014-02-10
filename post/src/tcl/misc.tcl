#/*****************************************************************************
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

#*******************************************************************************
#*
#* Misc utilities...
#*
#*******************************************************************************
#*
#*                     Author:       Juha Ruokolainen
#*
#*                    Address: CSC - IT Center for Science Ltd.
#*                                Keilaranta 14, P.O. BOX 405
#*                                  02101 Espoo, Finland
#*                                  Tel. +358 0 457 2723
#*                                Telefax: +358 0 457 2302
#*                              EMail: Juha.Ruokolainen@csc.fi
#*
#*                       Date: 26 Sep 1995
#*
#*                Modified by:
#*
#*       Date of modification:
#*
#*******************************************************************************
proc quit {} { Exit }

proc trnpriority { args } {
    set Usage "Usage: trnpriority value\n\nValue is one of (trs,tsr,rts,rst,str,srt).\n\n"

    set value(trs) 0
    set value(tsr) 1
    set value(rts) 2
    set value(rst) 3
    set value(str) 4
    set value(srt) 5

    check_args trnpriority $Usage 0 1 opt opt_val 1 1 arg_val $args

    if { ![info exists arg_val] } { return $Usage }

    set sid [array startsearch arg_val]

    while { [array anymore arg_val $sid] != 0 } {

        set n [array nextelement arg_val $sid]
        set val $arg_val($n)
        set n [string range $n 1 [string length $n]]
    }

    array donesearch arg_val $sid

    if { ![info exists value($val)] } { return $Usage }

    cTrnPriority $value($val)
}

proc rotpriority { args } {
    set Usage "Usage: rotpriority value\n\nValue is one of (xyz,xzy,yxz,yzx,zxy,zyx,local,global)."

    set value(xyz) 0
    set value(xzy) 1
    set value(yxz) 2
    set value(yzx) 3
    set value(zxy) 4
    set value(zyx) 5
    set value(local)  6
    set value(global) 7

    set value(aapo)   8

    check_args rotpriority $Usage 0 1 opt opt_val 1 1 arg_val $args

    if { ![info exists arg_val] } { return $Usage }

    set sid [array startsearch arg_val]

    while { [array anymore arg_val $sid] != 0 } {

        set n [array nextelement arg_val $sid]
        set val $arg_val($n)
        set n [string range $n 1 [string length $n]]
    }

    array donesearch arg_val $sid

    if { ![info exists value($val)] } { return $Usage }

    cRotPriority $value($val)
}

proc rotate { args } {
    set Usage "Usage: rotate \[-relative -x -y -z\] ang1 \[ang2\] \[ang3\]\n\n"

    set opt(1,name) -relative; set opt(1,args) 0
    set opt(2,name) -x; set opt(2,args) 0
    set opt(3,name) -y; set opt(3,args) 0
    set opt(4,name) -z; set opt(4,args) 0

    check_args rotate $Usage 4 0 opt opt_val 1 3 arg_val $args

    if { ![info exists arg_val] } { return $Usage }

    set x 0
    set y 0
    set z 0
    set rel 0

    if { [info exists opt_val] } {
        set sid [array startsearch opt_val]

        while { [array anymore opt_val $sid] != 0 } {

            set n [array nextelement opt_val $sid]
            set value $opt_val($n)
            set n [string range $n 1 [string length $n]]

            switch $n {
               relative { set rel 1 }
               x { set x 1 }
               y { set y 1 }
               z { set z 1 }
           }
        }

        array donesearch opt_val $sid
    }

    if { [@ $x+$y+$z] > 1 } {
        return -code error "rotate: Give only one of \[x,y,z\].\n\n$Usage"
    }

    if { $x } {
        set cmd0 "x"
    } elseif { $y } {
        set cmd0 "y"
    } elseif { $z } {
        set cmd0 "z"
    } else { set cmd0 "all" }

    if { $rel } { set cmd1 "relative" } else { set cmd1 "absolute" }

    set cmd2 $arg_val(0)
    if { [info exists arg_val(1)] } { set cmd2 "$cmd2 $arg_val(1)" }
    if { [info exists arg_val(2)] } { set cmd2 "$cmd2 $arg_val(2)" }

    cRotate $cmd0 $cmd1 $cmd2
}

proc scale { args } {
    set Usage "Usage: scale \[-relative -x -y -z -all\] scl1 \[scl2\] \[scl3\]\n\n"

    set opt(1,name) -relative; set opt(1,args) 0
    set opt(2,name) -x; set opt(2,args) 0
    set opt(3,name) -y; set opt(3,args) 0
    set opt(4,name) -z; set opt(4,args) 0
    set opt(5,name) -all; set opt(5,args) 0

    check_args scale $Usage 5 0 opt opt_val 1 3 arg_val $args

    if { ![info exists arg_val] } { return $Usage }

    set x 0
    set y 0
    set z 0
    set all 0
    set rel 0

    if { [info exists opt_val] } {
        set sid [array startsearch opt_val]
        while { [array anymore opt_val $sid] != 0 } {

            set n [array nextelement opt_val $sid]
            set value $opt_val($n)
            set n [string range $n 1 [string length $n]]

            switch $n {
               relative { set rel 1 }
               x { set x 1 }
               y { set y 1 }
               z { set z 1 }
               all { set all 1 }
            }
        }
        array donesearch opt_val $sid
    }

    if { [@ $x+$y+$z+$all] > 1 } {
        return -code error "scale: Give only one of \[x,y,z,all\].\n\n$Usage"
    }

    if { $x } {
        set cmd0 "x"
    } elseif { $y } {
        set cmd0 "y"
    } elseif { $z } {
        set cmd0 "z"
    } else { set cmd0 "all" }

    if { $rel } { set cmd1 "relative" } else { set cmd1 "absolute" }

    set cmd2 $arg_val(0)
    if { [info exists arg_val(1)] } { set cmd2 "$cmd2 $arg_val(1)" }
    if { [info exists arg_val(2)] } { set cmd2 "$cmd2 $arg_val(2)" }

    cScale $cmd0 $cmd1 $cmd2
}

proc translate { args } {
    set Usage "Usage: translate \[-relative -x -y -z\] trs1 \[trs2\] \[trs3\]\n\n"

    set opt(1,name) -relative; set opt(1,args) 0
    set opt(2,name) -x; set opt(2,args) 0
    set opt(3,name) -y; set opt(3,args) 0
    set opt(4,name) -z; set opt(4,args) 0

    check_args translate $Usage 4 0 opt opt_val 1 3 arg_val $args

    if { ![info exists arg_val] } { return $Usage }

    set x 0
    set y 0
    set z 0
    set rel 0

    if { [info exists opt_val] } {
        set sid [array startsearch opt_val]

        while { [array anymore opt_val $sid] != 0 } {
            set n [array nextelement opt_val $sid]
            set value $opt_val($n)
            set n [string range $n 1 [string length $n]]

            switch $n {
               relative { set rel 1 }
               x { set x 1 }
               y { set y 1 }
               z { set z 1 }
            }
        }
        array donesearch opt_val $sid
    }

    if { [@ $x+$y+$z] > 1 } {
        return -code error "scale: Give only one of \[x,y,z,all\].\n\n$Usage"
    }

    if { $x } {
        set cmd0 "x"
    } elseif { $y } {
        set cmd0 "y"
    } elseif { $z } {
        set cmd0 "z"
    } else { set cmd0 "all" }

    if { $rel } { set cmd1 "relative" } else { set cmd1 "absolute" }

    set cmd2 $arg_val(0)
    if { [info exists arg_val(1)] } { set cmd2 "$cmd2 $arg_val(1)" }
    if { [info exists arg_val(2)] } { set cmd2 "$cmd2 $arg_val(2)" }

    cTranslate $cmd0 $cmd1 $cmd2
}
