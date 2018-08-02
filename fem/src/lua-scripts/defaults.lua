-- /*****************************************************************************/
--  *
--  *  Elmer, A Finite Element Software for Multiphysical Problems
--  *
--  *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
--  *
--  *  This program is free software; you can redistribute it and/or
--  *  modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE version 2.1
--  *  as published by the Free Software Foundation. 
--  * 
--  *  This program is distributed in the hope that it will be useful,
--  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  *  GNU General Public License for more details.
--  *
--  *  You should have received a copy of the GNU Lesser General Public
--  *  License along with this library; if not, write to the Free Software
--  *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
--  *
--  *****************************************************************************/
-- 
-- /******************************************************************************
--  *
--  *  Authors: Juhani Kataja
--  *  Email:   juhani.kataja@csc.fi
--  *  Web:     http://www.csc.fi/elmer
--  *  Address: CSC - IT Center for Science Ltd.
--  *           Keilaranta 14
--  *           02101 Espoo, Finland
--  *
--  *  Original Date: 08 Jun 1997
--  *
--  *****************************************************************************/


-- Default Lua scripts for Elmer 

-- Define some mathematical constants and functions

pi = math.pi
sin = math.sin
cos = math.cos
tan = math.tan
asin = math.asin
acos = math.acos
sinh = math.sinh
cosh = math.cosh
tanh = math.tanh
log = math.log
log10 = math.log10


-- Reads filename and extracts strings between '!---LUA BEGIN' and '!---LUA END'
-- removes optional trailing '!' from the extracted strings.
-- Returns the extracted string

function readsif(fname)
  local f = assert(io.open(fname), 'r')

  local luadata = ""
  local luablock = false 
  local line = f:read()
  local linenum, beginline

  local function remove_trailing_char(line, char)
    local i,j = string.find(line, char)
    if i == 1 then
      return line:sub(2)
    end
    return line
  end

  linenum = 1
  beginline = linenum
  repeat
    if not(luablock) then
      local i, j = string.find(line, "!---LUA BEGIN")
      if i == 1 then
        luablock = true
        beginline = linenum
      end 
    else
      local i, j = string.find(line, "!---LUA END")
      if i == 1 then -- found end of block
        luablock = false
      else
        luadata = luadata .. remove_trailing_char(line, '!') .. "\n"
      end
    end
    line = f:read()
    linenum = linenum + 1
  until line == nil
  if luablock then
    error("unmatched '!---LUA BEGIN' at line ".. beginline)
  end 
  f:close()
  return luadata
end 


-- This (create_new_fun) will create unique names for functions with given
-- prefix and body used to transform `REAL LUA "..."` to function with body
-- matching `...`. Uses global variable ELMER_FUNCTION_COUNTER to ensure
-- uniqueness of the name
function create_new_fun(prefix,  body)
  -- Do naive string sanitization 
  local sane_prefix = string.gsub(prefix, "{", "_OCB_")
  local sane_prefix = string.gsub(sane_prefix, "}", "_CCB_")
  local sane_prefix = string.gsub(sane_prefix, " ", "_")
  local sane_prefix = string.gsub(sane_prefix, "%(", "_ORB_")
  local sane_prefix = string.gsub(sane_prefix, "%)", "_CRB_")
  local sane_prefix = string.gsub(sane_prefix, "%[", "_OSB_")
  local sane_prefix = string.gsub(sane_prefix, "%]", "_CSB_")

  if ELMER_FUNCTION_SUFFIX_TABLE == nil then
    ELMER_FUNCTION_SUFFIX_TABLE = {}
  end 

  if ELMER_FUNCTION_SUFFIX_TABLE[sane_prefix] == nil then
    ELMER_FUNCTION_SUFFIX_TABLE[sane_prefix] = 1
  else
    ELMER_FUNCTION_SUFFIX_TABLE[sane_prefix] = ELMER_FUNCTION_SUFFIX_TABLE[sane_prefix] + 1
  end
  local counter = ELMER_FUNCTION_SUFFIX_TABLE[sane_prefix]

  local underscored, num_space = string.gsub(sane_prefix, " ", "_")
  local fname = underscored .. "_" .. counter
  local codestr = "function " .. fname .. "() return " .. body .. " end"
  local code = loadstring(codestr)
  code()
  return fname
end  

-- print("pe, thread " .. ELMER_PARALLEL["pe"] .. ", " .. ELMER_PARALLEL["thread"])
