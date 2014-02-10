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
#Module:    ecif_tk_screenParams.tcl
#Language:  Tcl
#Date:      05.10.98
#Version:   1.00
#Author(s): Martti Verho
#Revisions: 
#
#Abstract:  Global screen parameter variables (for Win32 and Unix).
#
#************************************************************************


# NOTE: This is a big proc!!!
#
proc Screen::setParameters {} {
global Info

  # WIN32 font settings
  # ===================
  if { $Info(platform) == "windows" } {
    
    if { $Info(ELMER_FRONT_FONT_SIZE) != "" } {
      set fsize $Info(ELMER_FRONT_FONT_SIZE)
    } else {
      set fsize 8
    }


  # Unix font settings
  # ==================
  } else {

    if { $Info(ELMER_FRONT_FONT_SIZE) != "" } {
      set fsize $Info(ELMER_FRONT_FONT_SIZE)
    } else {
      set fsize 11
    }
  }

  #---Font types and sizes

  set Info(fontFamilyNormal) "Helvetica"
  set Info(fontFamilyFixed)  "Courier"

  set meFn [list "Helvetica"]; lappend meFn $fsize
  set buFn [list "Helvetica"]; lappend buFn $fsize
  set lbFn [list "Helvetica"]; lappend lbFn $fsize
  set msFn [list "Helvetica"]; lappend msFn $fsize
  set laFn [list "Helvetica"]; lappend laFn $fsize

  #option add *Menu.font $meFn
  #option add *Menubutton.font $buFn
  #option add *Button.font $buFn
  #option add *Checkbutton.font $buFn
  #option add *Radiobutton.font $buFn
  #option add *Listbox.font $lbFn
  #option add *Dialog.msg.font $msFn
  #option add *Label.font $laFn

  # Main window command buttons
  set Info(cmdBtnFont) [list "Helvetica"]; lappend Info(cmdBtnFont) $fsize

  set Info(fontSizeEntry) $fsize
  set Info(fontSizeTable) $fsize
  set Info(fontSizeMsg)   $fsize

  #---Frame padding
  set Info(framePadX1) [ expr int( 0.005 * $Info(maxWinSizeX) ) ]
  set Info(framePadX2) [ expr int( 0.010 * $Info(maxWinSizeX) ) ]
  set Info(framePadX3) [ expr int( 0.015 * $Info(maxWinSizeX) ) ]

  set Info(framePadY1) [ expr int( 0.002 * $Info(maxWinSizeY) ) ]
  set Info(framePadY2) [ expr int( 0.004 * $Info(maxWinSizeY) ) ]
  set Info(framePadY3) [ expr int( 0.006 * $Info(maxWinSizeY) ) ]

  #---Standard panel objects/elements/parameters list box sizes
  set Info(standardLbHeight) 8
  set Info(standardObjectLbWidth2) 24
  set Info(standardObjectLbWidth3) 13
  set Info(standardBoundaryLbWidth) 24
  set Info(standardParameterLbWidth) 13

  # Fonts
  set fam   $Info(fontFamilyNormal)
  set fam_f $Info(fontFamilyFixed)

  set fsz_e $Info(fontSizeEntry)
  set fsz_t $Info(fontSizeTable)
  set fsz_m $Info(fontSizeMsg)

  set entryFont     [font create -family $fam -size $fsz_e]
  set entryFontB    [font create -family $fam -size $fsz_e -weight bold]
  set entryFontI    [font create -family $fam -size $fsz_e -slant italic]
  set entryFontBI   [font create -family $fam -size $fsz_e -weight bold -slant italic]
  set entryFontFix  [font create -family $fam_f -size $fsz_e]

  set tableFont    [font create -family $fam -size $fsz_t]
  set tableFontB   [font create -family $fam -size $fsz_t -weight bold]
  set tableFontI   [font create -family $fam -size $fsz_t -slant italic]
  set tableFontBI  [font create -family $fam -size $fsz_t -weight bold -slant italic]
  set tableFontFix [font create -family $fam_f -size $fsz_t]

  set msgFont     [font create -family $fam -size $fsz_m]
  set msgFontB    [font create -family $fam -size $fsz_m -weight bold]
  set msgFontI    [font create -family $fam -size $fsz_m -slant italic]
  set msgFontBI   [font create -family $fam -size $fsz_m -weight bold -slant italic]
  set msgFontFix  [font create -family $fam_f -size $fsz_m]


  #############
  # Windows32 #
  #############

  if { $Info(platform) == "windows" } {

    set Info(arrowCursor) arrow
    set Info(insertCursor) xterm


    set select_color  white
    set option_menu_bg gray81
    set non_active_bg gray85
    set non_active_bg_light  gray90
    set non_active_fg gray65

    set window_bg_clr systemWindow

    # Win32 Error option and colors
    #
    set en_err_option -bg
    set om_err_option -bg
    set rb_err_option -bg
    set cb_err_option -bg

    set en_err_color yellow
    set om_err_color yellow
    set rb_err_color systemButtonFace
    set cb_err_color systemButtonFace

    set en_nerr_color systemWindow
    set om_nerr_color systemButtonFace
    set rb_nerr_color systemButtonFace
    set cb_nerr_color systemButtonFace


    # Win32 Modified option and colors
    #
    set en_mod_option -fg
    set om_mod_option -fg
    set rb_mod_option -bg
    set cb_mod_option -selectcolor

    set en_mod_color red
    set om_mod_color red
    set rb_mod_color systemButtonFace
    set cb_mod_color pink

    set en_nmod_color systemWindowText
    set om_nmod_color systemWindowText
    set rb_nmod_color systemButtonFace
    set cb_nmod_color systemWindow

    # Win32 widget fonts
    #
    set entryFont     $entryFontFix
    set tableFont     $tableFontFix 
    set msgFont       $msgFont

  # End Win32 specific


  ########
  # Unix #
  ########
  } else {

    set Info(arrowCursor) arrow
    set Info(insertCursor) xterm

    set select_color   gray66 ;# darkGray is  66.5
    set option_menu_bg gray81 ;# lightGray is 82.5
    set non_active_bg  gray85
    set non_active_bg_light  gray90
    set non_active_fg  gray65

    set window_bg_clr "#d9d9d9"

    # Unix Error option and colors
    #
    set en_err_option -bg
    set om_err_option -bg
    set rb_err_option -bg
    set cb_err_option -bg

    set en_err_color yellow
    set om_err_color yellow
    set rb_err_color $window_bg_clr
    set cb_err_color $window_bg_clr

    set en_nerr_color white
    set om_nerr_color $window_bg_clr
    set rb_nerr_color $window_bg_clr
    set cb_nerr_color $window_bg_clr


    # Unix Modified option and colors
    #
    set en_mod_option -fg
    set om_mod_option -fg
    set rb_mod_option -bg
    set cb_mod_option -selectcolor

    set en_mod_color red
    set om_mod_color red
    set rb_mod_color $window_bg_clr
    set cb_mod_color pink

    set en_nmod_color black
    set om_nmod_color black
    set rb_nmod_color $window_bg_clr
    set cb_nmod_color white

    # Unix widget fonts
    #
    set entryFont    $entryFontFix
    set tableFont    $tableFontFix 
    set msgFont      $msgFont

  }
  # End Unix specific


  ####################################
  # Screen outlook control variables #
  ####################################

  #---Fonts by object type
  set Info(entryFont)   $entryFont
  set Info(tableFont)   $tableFont
  set Info(msgAreaFont) $msgFont

  set Info(browserFont)      $Info(entryFont)
  #set Info(processTableFont) $Info(entryFont)
  set Info(processTableFont) $Info(tableFont)

  #---Colors
  # check box selection color
  set Info(select_color) $select_color
  set Info(optionMenuBg) $option_menu_bg
  set Info(nonActiveBg)  $non_active_bg
  set Info(nonActiveBgLight)  $non_active_bg_light
  set Info(nonActiveFg)  $non_active_fg

  #---Widget relief
  set Info(activeRelief) sunken
  set Info(nonActiveRelief) sunken

  set Info(en_err_option) $en_err_option
  set Info(om_err_option) $om_err_option
  set Info(rb_err_option) $rb_err_option
  set Info(cb_err_option) $cb_err_option

  set Info(en_err_color) $en_err_color
  set Info(om_err_color) $om_err_color
  set Info(rb_err_color) $rb_err_color
  set Info(cb_err_color) $cb_err_color

  set Info(en_nerr_color) $en_nerr_color
  set Info(om_nerr_color) $om_nerr_color
  set Info(rb_nerr_color) $rb_nerr_color
  set Info(cb_nerr_color) $cb_nerr_color


  set Info(en_mod_option) $en_mod_option
  set Info(om_mod_option) $om_mod_option
  set Info(rb_mod_option) $rb_mod_option
  set Info(cb_mod_option) $cb_mod_option

  set Info(en_mod_color) $en_mod_color
  set Info(om_mod_color) $om_mod_color
  set Info(rb_mod_color) $rb_mod_color
  set Info(cb_mod_color) $cb_mod_color

  set Info(en_nmod_color) $en_nmod_color
  set Info(om_nmod_color) $om_nmod_color
  set Info(rb_nmod_color) $rb_nmod_color
  set Info(cb_nmod_color) $cb_nmod_color

} ;# End proc Screen::setParameters 

# end ecif_tk_screenParams.tcl
# ********************
