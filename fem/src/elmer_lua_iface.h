/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/******************************************************************************
 *
 *  Authors: Juhani Kataja
 *  Email:   juhani.kataja@csc.fi
 *  Web:     http://www.csc.fi/elmer
 *  Address: CSC - IT Center for Science Ltd.
 *           Keilaranta 14
 *           02101 Espoo, Finland
 *
 *  Original Date: 08 Jun 1997
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <math.h>

typedef struct {
  double a;
  double b;
} history_t; 

typedef struct NumArray {
  int size;
  double values[1];
} NumArray;

int luaopen_array(lua_State *L);

static void stackDump(lua_State *L) ;

int getfield(lua_State* L, const char* key) ;

void load(lua_State* L, char * filename, int* width, int* height) ;

lua_State* lua_init();

void lua_runfile(lua_State* L, const char* filename);


