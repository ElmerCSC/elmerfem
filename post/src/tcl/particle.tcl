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
#* Particles display parameter settings
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
#
#

#
# 23 Apr 1996
#

set ParticleNofParticles      0
set ParticleLineStyle         0
set ParticleArrowStyle        0
set ParticleQuality           1
set ParticleRadius            1
set ParticleColor             "none"
set ParticleVelocity          "none"
set ParticleParticle          "none"
set ParticleOutDT             1.0E-1
set ParticleMaxDT             1.0E-3
set ParticleTolerance         1.0E-5
set ParticleIntegMethod       1
set ParticleIntegPolicy       1
set ParticleStyle             0

set ParticleVariableNames(0)  "none"
set NumberOfParticleVariables 0

set particle "none"

proc get_particle_variable { w } {
    global particle

    set particle [$w get [$w curselection]]
}

#
# List box containing particle variable names
#
proc make_particle_list { } {
   global ParticleVariableNames NumberOfParticleVariables particle saveparticle

   if { ![info exists ParticleVariableNames] } { return }

   set saveparticle $particle

   toplevel .vlist
   place_window .vlist

   frame .vlist.vari -relief sunken -bg lightblue
   listbox .vlist.vari.list -yscroll ".vlist.vari.scroll set"
   scrollbar .vlist.vari.scroll -command ".vlist.vari.list yview"

   pack .vlist.vari.list -side left -fill y
   pack .vlist.vari.scroll -side left -expand 1 -fill both
   pack .vlist.vari -side top -expand 1 -fill both

   bind .vlist.vari.list <Double-1> { set particle [get_particle_variable %W] }

   frame .vlist.equ
   label .vlist.equ.label -text "math: "
   entry .vlist.equ.entry -width 20 -relief sunken -textvariable mathcmd
   pack .vlist.equ -side top
   pack .vlist.equ.label -side left
   pack .vlist.equ.entry -side left -fill x

   bind .vlist.equ.entry <Return> { math $mathcmd; set particle "tryagain"; }

   frame .vlist.close
   button .vlist.close.ok  -text  "OK" -command { set particle [get_particle_variable .vlist.vari.list]  }
   button .vlist.close.cancel -text  "Cancel" -command { set particle $saveparticle; }

   pack .vlist.close -side top
   pack .vlist.close.ok -side right
   pack .vlist.close.cancel -side right

   set oldfocus [focus]
   set particle "tryagain"
   while { $particle == "tryagain" } {
       .vlist.vari.list delete 0 end

       .vlist.vari.list insert end none
       do i 0 [@ $NumberOfParticleVariables-1] {
           set val $ParticleVariableNames($i)
           .vlist.vari.list insert end $val
       }

       vwait particle
   }
   focus $oldfocus
   destroy .vlist

   return $particle
}


proc particle_edit { } {

    global ParticleLineStyle ParticleQuality ParticleRadius ParticleStyle
    global ParticleOutDT ParticleMaxDT  ParticleTolerance ParticleArrowStyle
    global ParticleColor ParticleVelocity ParticleParticle ParticleIntegMethod
    global ParticlePolicy

    if { [winfo exists .particle] } {
        wm iconify .particle
        wm deiconify .particle
        return
    }

    toplevel .particle
    place_window .particle

#
# style
#
    frame .particle.style
    label .particle.style.lab -text "Particle Style: "
    radiobutton .particle.style.vec -value 0 -variable ParticleStyle -text "Arrow"
    radiobutton .particle.style.sph -value 1 -variable ParticleStyle -text "Sphere"

    pack .particle.style -side top
    pack .particle.style.lab -side left
    pack .particle.style.vec -side left  -fill x
    pack .particle.style.sph -side left -fill x

#
# line style
#
    frame .particle.line
    label .particle.line.label -text "Display Style: "
    radiobutton .particle.line.line -value 0 -variable ParticleLineStyle -text "Line"
    radiobutton .particle.line.cyli -value 1 -variable ParticleLineStyle -text "Solid"

    pack .particle.line -side top
    pack .particle.line.label -side left
    pack .particle.line.line -side left -fill x
    pack .particle.line.cyli -side left  -fill x
#
# arrow style
#
    frame .particle.arrow
    label .particle.arrow.label -text "Arrow Style: "
    radiobutton .particle.arrow.stick -value 0 -variable ParticleArrowStyle -text "Stick"
    radiobutton .particle.arrow.arrow -value 1 -variable ParticleArrowStyle -text "Arrow"

    pack .particle.arrow -side top
    pack .particle.arrow.label -side left
    pack .particle.arrow.stick -side left  -fill x
    pack .particle.arrow.arrow -side left -fill x

#
# disp qual
#
    frame .particle.qual
    label .particle.qual.label -text "Display Quality: "
    entry .particle.qual.entry -relief sunken -width 5 -textvariable ParticleQuality

    pack .particle.qual -side top
    pack .particle.qual.label -side left
    pack .particle.qual.entry -side left -fill x

#
# Display scale
#
    frame .particle.radi
    label .particle.radi.label -text "Display Scale: "
    entry .particle.radi.entry -relief sunken -width 5 -textvariable ParticleRadius

    pack .particle.radi -side top
    pack .particle.radi.label -side left
    pack .particle.radi.entry -side left -fill x

#
# color
#
    frame .particle.color
    label .particle.color.label -text "Color Variable: "
    button .particle.color.but -textvariable ParticleColor \
           -command { set ParticleColor [make_scalar_list]; UpdateVariable "ParticleColor" }
    UpdateVariable "ParticleColor"

    pack .particle.color -side top
    pack .particle.color.label -side left
    pack .particle.color.but -side left -fill x

#
# velocity
#
    frame .particle.velocity
    label .particle.velocity.label -text "Velocity Variable: "
    button .particle.velocity.but -textvariable ParticleVelocity \
               -command { set ParticleVelocity [make_vector_list]; }

    pack .particle.velocity -side top
    pack .particle.velocity.label -side left
    pack .particle.velocity.but -side left -fill x

#
# particle
#
    frame .particle.particle
    label .particle.particle.label -text "Particle Variable: "
    button .particle.particle.but -textvariable ParticleParticle \
               -command { set ParticleParticle [make_particle_list]; }

    pack .particle.particle -side top
    pack .particle.particle.label -side left
    pack .particle.particle.but -side left -fill x
#
# space
#
    frame .particle.space0
    label .particle.space0.lab -text " "

    pack .particle.space0 -side top
    pack .particle.space0.lab -side left

#
# integ method
#
    frame .particle.integ
    label .particle.integ.label -text "Integration Method: "
    radiobutton .particle.integ.euler -value 0 -variable ParticleIntegMethod -text "Euler"
    radiobutton .particle.integ.runge -value 1 -variable ParticleIntegMethod -text "Runge Kutta"

    pack .particle.integ -side top
    pack .particle.integ.label -side left
    pack .particle.integ.euler -side left -fill x
    pack .particle.integ.runge -side left  -fill x
#
# integ policy
#
    frame .particle.policy
    label .particle.policy.label -text "Integration Policy: "
    radiobutton .particle.policy.fixed -value 0 -variable ParticleIntegPolicy -text "Fixed"
    radiobutton .particle.policy.adapt -value 1 -variable ParticleIntegPolicy -text "Adaptive"

    pack .particle.policy -side top
    pack .particle.policy.label -side left
    pack .particle.policy.fixed -side left -fill x
    pack .particle.policy.adapt -side left  -fill x
#
# output timestep
#
    frame .particle.time
    label .particle.time.timelab -text "Output timestep: "
    entry .particle.time.timeent -relief sunken -width 7 -textvariable ParticleOutDT

    pack .particle.time -side top
    pack .particle.time.timelab -side left
    pack .particle.time.timeent -side left -fill x

#
# max timestep
#
    label .particle.time.maxdtlab -text "Maximum timestep: "
    entry .particle.time.maxdtent -relief sunken -width 7 -textvariable ParticleMaxDT

    pack .particle.time.maxdtlab -side left
    pack .particle.time.maxdtent -side left -fill x

#
# tolerance
#
    frame .particle.err
    label .particle.err.tolelab -text "Absolute error: "
    entry .particle.err.toleent -relief sunken -width 7 -textvariable ParticleTolerance

    pack .particle.err -side top
    pack .particle.err.tolelab -side left
    pack .particle.err.toleent -side left -fill x
#
# space
#
    frame .particle.space1
    label .particle.space1.lab -text " "

    pack .particle.space1 -side top
    pack .particle.space1.lab -side left

#
# advance / close & apply
#
    frame .particle.buttons
    button .particle.buttons.advance -text "Advance" -command "pintegrate"
    button .particle.buttons.apply -text "Apply" -command "UpdateObject; play"
    button .particle.buttons.close -text "Close" -command "destroy .particle"

    pack .particle.buttons -side top
    pack .particle.buttons.advance -side left
    pack .particle.buttons.apply -side left
    pack .particle.buttons.close -side left -fill x
}

proc pintegrate { } {
    global ParticleAdvance;

    set ParticleAdvance 1
    UpdateObject

    play
}
