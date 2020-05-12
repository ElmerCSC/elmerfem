/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library (in file ../LGPL-2.1); if not, write 
 *  to the Free Software Foundation, Inc., 51 Franklin Street, 
 *  Fifth Floor, Boston, MA  02110-1301  USA
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

#include "elmer_lua_iface.h"


static void stackDump(lua_State *L) {
  int i;
  int top = lua_gettop(L);
  for (i = 1; i<= top; i++) {
    int t = lua_type(L,i);
    switch (t) {
      case LUA_TSTRING:
        printf("`%s'", lua_tostring(L, i));
        break;

      case LUA_TBOOLEAN:  
        printf(lua_toboolean(L, i) ? "true" : "false");
        break;

      case LUA_TNUMBER:  
        printf("%g", lua_tonumber(L, i));
        break;

      default: 
        printf("%s", lua_typename(L, t));
        break;
    }
    printf("  ");  
  }
  printf("\n"); 
}

lua_State* lua_init() {
  lua_State *L = lua_open();
  luaL_openlibs(L);
  luaopen_array(L);
  return L;
}

void lua_runfile(lua_State* L, const char* filename) {
  if (strlen(filename) > 0 && luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0)) {
    luaL_error(L, "cannot run configuration file: %s", lua_tostring(L, -1));
  }
}

void lua_pop_c(lua_State* L, int n) {
  lua_pop(L, n);
}

/* static methods and structs for handling the tx array */
static int newarray(lua_State *L) {
  int n = luaL_checkint(L, 1);
  size_t nbytes = sizeof(NumArray) + (n-1)*sizeof(double);
  NumArray *a = (NumArray *)lua_newuserdata(L, nbytes);

  luaL_getmetatable(L, "LuaBook.array");
  lua_setmetatable(L, -2);

  a->size = n;
  return 1;
}

static NumArray *checkarray(lua_State *L) {
  void* ud = luaL_checkudata(L, 1, "LuaBook.array");
  luaL_argcheck(L, ud != NULL, 1, "'array' expected");
  return (NumArray *) ud;
}

/* Get element pointer */ 
static double *getelem(lua_State *L) {
  NumArray *a = checkarray(L);
  /* int index = luaL_checkint(L, 2); */
  int index = luaL_checkint(L, 2);

  /* luaL_argcheck(L, 0 <= index && index < a->size, 2, */
  luaL_argcheck(L, 0 <= index && index < a->size, 2,
                "index out of range");
  return &(a->values[index]);
}

/* Set element (from lua stack) */
static int setarray(lua_State* L) {
  double newvalue = luaL_checknumber(L, 3);
  *getelem(L) = newvalue;
  return 0;
}

/* Get element (to lua stack) */
static int getarray(lua_State* L) {
  lua_pushnumber(L, *getelem(L));
  return 1;
}

/* Return size of array (to lua stack) */
static int getsize(lua_State* L) {
  NumArray *a = checkarray(L);
  lua_pushnumber(L, a->size);
  return 1;
}

/* arraylib method table (to be used with tx arrays)*/
static const struct luaL_reg arraylib [] = {
    {"new", newarray},
    {"set", setarray},
    {"get", getarray},
    {"size", getsize},
    {NULL, NULL}
};

#pragma omp threadprivate(arraylib)
/* end of tx array associated static methods */

/* turn arraylib instantes into metatables with set and get methods
 * obeying a[n] = x semantics. Indexing starts at 0 due to tradition!*/
int luaopen_array(lua_State *L) {
  luaL_newmetatable(L, "LuaBook.array");
  luaL_openlib(L, "array", arraylib, 0);
  lua_pushstring(L, "__index");
  lua_pushstring(L, "get");
  lua_gettable(L, 2);
  lua_settable(L, 1);
  lua_pushstring(L, "__newindex");
  lua_pushstring(L, "set");
  lua_gettable(L, 2);
  lua_settable(L, 1);
  return 0;
}

void get_userdataptr(lua_State *L, void* cp_raw, double** table, int* len) {
     NumArray* T = (NumArray*) cp_raw;
     *len = T -> size;
     *table = &(T -> values[0]);
}
