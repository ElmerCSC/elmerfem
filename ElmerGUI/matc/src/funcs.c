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

/*******************************************************************************
 *
 *     MATC user function utilities.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 30 May 1996
 *
 *                Modified by:
 *
 *       Date of modification:
 *
 ******************************************************************************/
/***********************************************************************
|
|  FUNCS.C - Last Edited 7. 8. 1988
|
***********************************************************************/

/*======================================================================
|Syntax of the manual pages:
|
|FUNCTION NAME(...) params ...
|
$  usage of the function and type of the parameters
?  explane the effects of the function
=  return value and the type of value if not of type int
@  globals effected directly by this routine
!  current known bugs or limitations
&  functions called by this function
~  these functions may interest you as an alternative function or
|  because they control this function somehow
^=====================================================================*/


/*
 * $Id: funcs.c,v 1.2 2005/05/27 12:26:20 vierinen Exp $ 
 *
 * $Log: funcs.c,v $
 * Revision 1.2  2005/05/27 12:26:20  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.3  2003/05/06 09:14:49  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:34:39  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

FUNCTION *fnc_check(name) char *name;
/*======================================================================
?  Look for specified user defined function from the FUNCTIONS list
|
=  NULL if not found, otherwise FUNCTION *fnc
&  lst_find()
^=====================================================================*/
{
  return (FUNCTION *)lst_find(FUNCTIONS, name);
}

VARIABLE *fnc_delete(ptr) VARIABLE *ptr;
/*======================================================================
?  Unlink given function definition from list FUNCTION *FUNC_HEAD,
|  and free associated memory. 
|
|  user command fdel("name")
|
@  FUNC_HEAD
&  FREEMEM, var_to_string(), fprintf(), fnc_free_entry(), fnc_check()
^=====================================================================*/
{
   FUNCTION *fnc;                  /* all these exist just because        */
   char *s;                        /* i can't get this done without them  */

   /*
       convert string from ptr
   */
   s = var_to_string(ptr);

   /*
       function exists. Unlink from list, and free memory.
   */
   if ((fnc = fnc_check(s)) != (FUNCTION *)NULL) {

     fnc_free_entry(fnc);

   }

   /*
        we did not found the function.
   */
   else {
     error("Function definition not found: %s.\n", s);
   }

   FREEMEM(s);

   return (VARIABLE *)NULL;
}

VARIABLE *fnc_list(ptr) VARIABLE *ptr;
/*======================================================================
?  Print given function definition from list FUNCTION *FUNC_HEAD,
|
|  user command flist("name")
|
&  FREEMEM, var_to_string(), printclause(), fnc_check()
^=====================================================================*/
{
   FUNCTION *fnc;                  /* all these exist just because    */
   char *s, *file;                 /* i can't get this done without   */
   int i;                          /* them.                           */

   FILE *fp = math_out;

   /*
       convert string from ptr 
   */
   s = var_to_string(ptr);
 
   /* 
       function exists. try listing the definition
   */
   if ((fnc = fnc_check(s)) != (FUNCTION *)NULL) {

     /*
         If file name given try opening it.
     */
     if (NEXT(ptr) != (VARIABLE *)NULL) {
       file = var_to_string(NEXT(ptr));
       if ((fp = fopen(file, "a")) == (FILE *)NULL) {
         error( "flist: can't open file: %s.",file );
       }
       FREEMEM(file);
     }

     /*
      *   print function header.
      */
     PrintOut( "function %s", NAME(fnc) );
     if ( fnc->parcount != 0 )
     {
         PrintOut( "(%s", fnc->parnames[0] );
         for( i = 1; i < fnc -> parcount; i++ ) 
             PrintOut( ",%s", fnc -> parnames[i] );
         PrintOut( ")" );
     }
     PrintOut( "\n" );

     /*
           and then the body
     */
     /*
           printclause(fnc->body, fp, 1); PrintOut( "end\n" );
      */
     if ( fp != math_out ) fclose(fp);
   }

   /*
        we did not found the function.
   */
   else {
     error( "Function definition not found: %s\n", s );
   }

   FREEMEM(s);

   return (VARIABLE *)NULL;
}


void fnc_free_entry(fnc) FUNCTION *fnc;
/*======================================================================
?  Free allocated memory from FUNCTION structure.
|
&  FREEMEM, free_clause(), lst_free()
^=====================================================================*/
{
  int i;

  free_clause(fnc->body);      /* function body */
  if (fnc -> parcount > 0) {
    for(i = 0; i < fnc -> parcount; i++) {
      FREEMEM(fnc -> parnames[i]);     /* parameter names, if any */
    }
    FREEMEM((char *)fnc -> parnames);  /* parameter name array    */
  }

  if (fnc -> imports) {
    for(i = 0; fnc->imports[i] != NULL; i++) {
      FREEMEM(fnc -> imports[i]);     /* imported variable names, if any */
    }
    FREEMEM((char *)fnc -> imports);  /* name array */
  }

  if (fnc -> exports) {
    for(i = 0; fnc->exports[i] != NULL; i++) {
      FREEMEM(fnc -> exports[i]);     /* exported variable names, if any */
    }
    FREEMEM((char *)fnc -> exports);  /* name array */
  }

  lst_free(FUNCTIONS, (LIST *)fnc);
}

void fnc_free()
/*======================================================================
?  Deallocate memory reserved for user defined functions 
| and unlink the list FUNCTION *FUNC_HEAD.
|
@  FUNCTION *FUNC_HEAD
&  free_clause(), FREEMEM
^=====================================================================*/
{
   FUNCTION *fnc, *fnc1;

   for(fnc = (FUNCTION *)FUNC_HEAD; fnc;)
   {
     fnc1 = NEXT(fnc);
     fnc_free_entry(fnc);   /* just plain and cold */
     fnc = fnc1;
   }

   FUNC_HEAD = (LIST *)NULL;     /* security */
}

VARIABLE *fnc_exec(fnc, par) FUNCTION *fnc; VARIABLE *par;
/*======================================================================
?  Execute function from parameter FUNCTION *fnc, with it's 
|  parameters in VARIABLE VARIABLE *par;
|
=  Return value is the executed function's value, which is 
|  given in VARIABLE _function_name, or if nonexeistent, 
|  the return value of the last executed statement in 
|  function body.
|
@  VAR_HEAD
&  ALLOCMEM, FREEMEM, STRCOPY, strcpy(), fprintf(),
|  lst_unlink, var_free(), evalclause()
^=====================================================================*/
{
   VARIABLE *ptr, *imp, *res, *headsave, *var;
   char *str;
   int i;

   /*
      we make new global VARIABLE list for this function,
      have to save the old one.
   */
   headsave = (VARIABLE *)VAR_HEAD;

   /*
    *    rename parameter from function header
    */
   for(i = 0, ptr = par; ptr; ptr = NEXT(ptr), i++)
   {
     if (ptr == NULL) break;
     if (i < fnc->parcount)
       NAME(ptr) = STRCOPY(fnc -> parnames[i]);
     else
       NAME(ptr) = ALLOCMEM(1);
   }

   /*
    * check for imported variables
    */
   if (fnc->imports != NULL)
     for(i = 0; fnc->imports[i] != NULL; i++)
      if ((ptr = var_check(fnc->imports[i])) != NULL)
      {
        VAR_HEAD = (LIST *)par;
        if (var_check(fnc->imports[i]) == NULL) 
        {
          ptr = var_temp_copy(ptr);
          NAME(ptr) = STRCOPY(fnc->imports[i]);
          lst_add(VARIABLES, (LIST *)ptr);
        }
        par = (VARIABLE *)VAR_HEAD;
        VAR_HEAD = (LIST *)headsave;
      }
      else 
        PrintOut( "WARNING: %s: imported variable [%s] doesn't exist\n",
                          NAME(fnc), fnc->imports[i]);

 
   /*
       parameters to functions own list of VARIABLES.
   */
   VAR_HEAD = (LIST *)par;

   /*
       initializations done, execute the function body.
   */
   res = evalclause(fnc->body);

   par = (VARIABLE *)VAR_HEAD;
   /*
    * check for exported variables
    */
   if (fnc->exports != NULL)
     for(i = 0; fnc->exports[i] != NULL; i++)
       if ((ptr = var_check(fnc->exports[i])) != NULL)
       { 
         VAR_HEAD = (LIST *)headsave;
#if 0
         ptr = var_temp_copy(ptr);
         NAME(ptr) = STRCOPY( fnc->exports[i] );
#else
         var = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
         var->this = ptr->this;
         REFCNT(ptr)++;
         NAME(var) = STRCOPY( fnc->exports[i] );
#endif
         var_delete( fnc->exports[i] );
         lst_add( VARIABLES, (LIST *)var );
         headsave = (VARIABLE *)VAR_HEAD;

         VAR_HEAD = (LIST *)par;
       }

   /*
       check for explicit return value from
       VARIABLE named "_function_name"
   */
   str = ALLOCMEM(strlen(NAME(fnc)) + 2);
   str[0] = '_'; strcat(str, NAME(fnc));

   if ((res = var_check(str)) != NULL)
   {
     lst_unlink(VARIABLES, (LIST *)res);
     FREEMEM(NAME(res));
     NEXT(res) = NULL;
   }
   else {
     var_delete_temp(res);
     res = NULL;
   }

   FREEMEM(str);

   /*
      rebuild the environment and return
   */
   var_free();
   VAR_HEAD = (LIST *)headsave;

   return res;
}


void fnc_com_init()
/*======================================================================
?  Initialize function handling commands.
|
&  com_init()
~  com_init()
^=====================================================================*/
{
  com_init(
             "funcdel",  FALSE, FALSE, fnc_delete, 1, 1,
             "funcdel(name)\nDelete function definition from parser.\n"
          );

  com_init(
             "funclist", FALSE, FALSE, fnc_list, 1, 2,
             "funclist(name)\nGive header of a given function.\n\nSEE ALSO: help.\n"
          );
}
