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
 *     MATC variable manipulation.
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
/*********************************************************************
|
|  VARIABLE.C - Last Edited 9. 8. 1988
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
 * $Id: variable.c,v 1.6 2007/05/11 07:53:32 jpr Exp $ 
 *
 * $Log: variable.c,v $
 * Revision 1.6  2007/05/11 07:53:32  jpr
 * *** empty log message ***
 *
 * Revision 1.5  2006/02/07 10:21:42  jpr
 * Changed visibility of some variables to local scope.
 *
 * Revision 1.4  2006/02/02 06:54:44  jpr
 * small formatting changes.
 *
 * Revision 1.2  2005/05/27 12:26:22  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:58  jpr
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

VARIABLE *const_new(name, type, nrow, ncol) int type, ncol, nrow; char *name;
/*======================================================================
?  return a new global VARIABLE given name, type, and matrix size.
|  VARIABLE is linked to CONSTANTS lists.
|
=  pointer to a new VARIABLE 
&  mat_new(), lst_add(), ALLOCMEM, FREEMEM, STRCOPY
^=====================================================================*/
{
  VARIABLE *ptr;
  
  /*
       Allocate the structure and link to global list of VARIABLES.
  */

  ptr =   (VARIABLE *)ALLOCMEM(VARIABLESIZE);    /* list entry          */
  ptr->this = mat_new(type, nrow, ncol);         /* allocate new MATRIX */
  REFCNT(ptr) = 1;                               /* one reference       */
  NAME(ptr) = STRCOPY(name);                     /* name as given       */

  lst_add(CONSTANTS, (LIST *)ptr);               /* add to list   */

  return ptr;
}

VARIABLE *var_new(name, type, nrow, ncol) int type, ncol, nrow; char *name;
/*======================================================================
?  return a new global VARIABLE given name, type, and matrix size.
|  VARIABLE is linked to VARIABLES list.
|
=  pointer to a new VARIABLE 
&  var_check(), lst_add(), ALLOCMEM, FREEMEM, STRCOPY
^=====================================================================*/
{
  VARIABLE *ptr;
  
  /*
   * Delete old definition of name if any...
   */
  var_delete(name);

  /*
   *    Allocate the structure and link to global list of VARIABLES.
   */
  ptr =   (VARIABLE *)ALLOCMEM(VARIABLESIZE);    /* list entry          */
  ptr->this = mat_new(type, nrow, ncol);         /* allocate new MATRIX */
  REFCNT(ptr) = 1;                               /* one reference       */
  NAME(ptr) = STRCOPY(name);                     /* name as given       */

  lst_addhead(VARIABLES, (LIST *)ptr);           /* add to list */

  return ptr;
}

void var_create_vector( char *name, int ntime, int ncol, double *data )
{
    VARIABLE *var = var_new( name,TYPE_DOUBLE, ntime, ncol );
    int i;

    FREEMEM( MATR(var) );
    MATR(var) = data;
}

VARIABLE *var_rename(ptr, str) VARIABLE *ptr; char *str;
{
  VARIABLE *res;

  if (ptr == (VARIABLE *)NULL) return NULL;

  res = (VARIABLE *)lst_find( VARIABLES, str );

  if (res == NULL && REFCNT(ptr) > 1)
  {
    res = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
    NAME(res) = STRCOPY(str);
    res->this = mat_copy(ptr->this);
    REFCNT(res) = 1;
    lst_addhead(VARIABLES, (LIST *)res);
  }
  else if (res == NULL)
  {
    res = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
    NAME(res) = STRCOPY(str);
    res->this = ptr->this;
    REFCNT(res)++;
    lst_addhead(VARIABLES, (LIST *)res);
  }
  else
  {
     if ( res != ptr )
     {
#if 1
         if ( NROW(res) == NROW(ptr) && NCOL(res) == NCOL(ptr) )
         {
             memcpy( MATR(res),MATR(ptr), NROW(res)*NCOL(res)*sizeof(double) );
         }
         else 
#endif
         {
            if (--REFCNT(res) == 0)
            {
                 FREEMEM( (char *)MATR(res) );
                 FREEMEM( (char *)res->this );
             }
             res->this = ptr->this;
             REFCNT(res)++;
         }
     }
  }

  if ( res != ptr ) var_delete_temp(ptr);

  return res;
}

static int var_pprec = 3,
    var_pinp = FALSE, var_rowintime = FALSE;

VARIABLE *var_format(var) VARIABLE *var;
{
  if (*MATR(var) > 0 && *MATR(var) < 20)
  {
    var_pprec = *MATR(var);
  }

  if (NEXT(var) != NULL)
  {
    char *frm = var_to_string(NEXT(var));

    if (strcmp(frm,"input") == 0)
    {
      var_pinp = TRUE;
    }
    else
    {
      var_pinp = FALSE;
      if ( strcmp(frm,"rowform") == 0)
         var_rowintime = TRUE;
      else
         var_rowintime = FALSE;
    }
    FREEMEM(frm);
  }

  return (VARIABLE *)NULL;
}

void var_print(ptr) VARIABLE *ptr;
{
  double maxp, minp, maxx;
  int i, j, k;
  char fmt[80];
  
  if (ptr == (VARIABLE *)NULL) return;
  
  if (TYPE(ptr) == TYPE_STRING)
  {
    if (var_pinp)
      PrintOut( "%d %d %% \"",NROW(ptr),NCOL(ptr) );

    for(i = 0; i < NROW(ptr); i++)
    {
      for(j = 0; j < NCOL(ptr); j++)
        PrintOut( "%c",  (char)M(ptr,i,j));
      if (var_pinp)
      {
        if (i < NROW(ptr)-1)
          PrintOut("\"\\");
        else
          PrintOut("\"");
      }
      PrintOut( "\n");
    }
    return;
  }

  k = 0;
  do
  {
    if (var_pinp)
      PrintOut("%d %d %% ", NROW(ptr), NCOL(ptr));
    else if (NCOL(ptr) > 8 && !var_rowintime ) 
      PrintOut( "\nColumns %d trough %d\n\n", 
              k, min(NCOL(ptr) - 1, k + 7));

    if (var_pinp || var_rowintime )
      sprintf(fmt, "%%.%dg",var_pprec );
    else 
      sprintf(fmt, "%% %d.%dg",var_pprec+7,var_pprec); 
     
    for(i = 0; i < NROW(ptr); i++)
    {
      if ( var_rowintime ) {
         for( j=0; j<NCOL(ptr); j++ ) {
           if ( j>0 ) PrintOut(" ");
           PrintOut( fmt, M(ptr,i,j));
         }
      } else {
         for(j = 0; j < 80/(var_pprec+7) && k + j < NCOL(ptr); j++)
           PrintOut( fmt, M(ptr,i,j+k));

         if (var_pinp)
           if (i < NROW(ptr)-1) PrintOut("\\");
       }

      PrintOut("\n");
    }

    k += j;
  } while(k < NCOL(ptr));
}

void var_delete(str) char *str;
{
    VARIABLE *ptr;

    ptr = var_check(str);

    if ( ptr != (VARIABLE *)NULL )
    {
        if ( --REFCNT(ptr) == 0 )
        {
            FREEMEM((char *)MATR(ptr));
            FREEMEM((char *)ptr->this);
        }
        lst_free(VARIABLES, (LIST *)ptr);
     }
  
     return;
}

VARIABLE *var_vdelete( var ) VARIABLE *var;
{
   var_delete( var_to_string( var ) );
   return (VARIABLE *)NULL;
}


void var_free()
{
    VARIABLE *ptr; 
  
    for( ptr = (VARIABLE *)VAR_HEAD; ptr; ptr = NEXT(ptr) )
    {
        if ( --REFCNT(ptr) == 0 )
        {
            FREEMEM((char *)MATR(ptr));
            FREEMEM((char *)ptr->this);
        }
     }

     lst_purge(VARIABLES);
  
     return;
}

void const_free()
{
    VARIABLE *ptr; 
  
    for( ptr = (VARIABLE *)CONST_HEAD; ptr; ptr = NEXT(ptr) )
    {
        if ( --REFCNT(ptr) == 0 )
        {
            FREEMEM((char *)MATR(ptr));
            FREEMEM((char *)ptr->this);
        }
    }

    lst_purge(CONSTANTS);
  
    return;
}

VARIABLE *var_varlist()
/*======================================================================
?  print a list of VARIABLES for the user
|
=  (VARIABLE *)NULL
&  lst_print()
^=====================================================================*/
{
    lst_print(CONSTANTS); lst_print(VARIABLES);

    return NULL;
}

VARIABLE *var_ccheck(var) VARIABLE *var;
/*======================================================================
?  look for a VARIABLE from the global list of VARIABLES and return
|  it or (VARIABLE *)NULL if not found.
|
=  VARIABLE *
&  var_check(), var_to_string()
^=====================================================================*/
{
    VARIABLE *res;
    char *str;
    int i, n;

    for(n = 0, res = var; res != NULL; n++, res=NEXT(res));
    res = var_temp_new(TYPE_DOUBLE, 1, n);

    for( i=0; i<n; i++, var=NEXT(var) )
    {
        str = var_to_string(var);

       if ( var_check(str) == NULL )
           M(res,0,i) = FALSE;
       else
           M(res,0,i) = TRUE;

       FREEMEM(str);
    }

    return res;
}

VARIABLE *var_check(str) char *str;
/*======================================================================
?  look for a VARIABLE from the global list of VARIABLES and return
|  it or (VARIABLE *)NULL if not found.
|
=  VARIABLE *
&  lst_find()
^=====================================================================*/
{
    VARIABLE *res;

    if ( (res = (VARIABLE *)lst_find(VARIABLES, str)) == NULL )
    {
        res = (VARIABLE *)lst_find(CONSTANTS, str);
    }

  return res;
}

VARIABLE *var_temp_copy(from) VARIABLE *from;
/*======================================================================
?  Make a temporary (not linked to global list of VARIABLES) 
|  copy of a VARIABLE *from and.
|
=  pointer to new VARIABLE
&  ALLOCMEM
^=====================================================================*/
{
    VARIABLE *to;
  
    /*
     *  if there's nothing to copy return.
     */
    if ( from == NULL ) return NULL;

    to = (VARIABLE *)ALLOCMEM(VARIABLESIZE);  /* list entry */
    to->this = mat_copy(from->this);
    REFCNT(to) = 1;

    return to;
}

VARIABLE *var_temp_new(type,nrow,ncol) int type, nrow, ncol;
/*======================================================================
?  Make a new temporary (not linked to global list of VARIABLES) 
|  VARIABLE, type and matrix dimensions from function parameters.
|
=  pointer to new VARIABLE entry
&  ALLOCMEM
^=====================================================================*/
{
    VARIABLE *ptr;
  
    ptr =   (VARIABLE *)ALLOCMEM(VARIABLESIZE); /* list entry */
    ptr->this = mat_new(type, nrow, ncol);
    REFCNT( ptr ) = 1;

    return ptr; 
}

void var_delete_temp_el( VARIABLE *ptr )
{  
    if ( ptr != NULL )
    {
        if ( --REFCNT(ptr) == 0 )
        {
           FREEMEM((char *)MATR(ptr));
           FREEMEM((char *)ptr->this);
        }
        FREEMEM((char *)ptr);
    }
    return;
}

void var_delete_temp( VARIABLE *head )
{
    VARIABLE *ptr, *ptr1;

    for( ptr = head; ptr; )
    {
        ptr1 = NEXT(ptr);
        var_delete_temp_el(ptr);
        ptr = ptr1;
    }
    return;
}

char *var_to_string(ptr) VARIABLE *ptr;
{
    char *str;
    int i;

    str = ALLOCMEM(NCOL(ptr)+1);

    for( i=0; i<NCOL(ptr); i++ )
    {
        str[i] = (char)M(ptr, 0, i);
    }

    return str;
}

void var_reset_status(char *name)
{
   VARIABLE *ptr = var_check(name);

   if ( ptr ) ptr->changed = 0;
}


int var_get_status(char *name)
{
   VARIABLE *ptr = var_check(name);

   if ( ptr )
      return ptr->changed;
   else
      return 0;
}

void var_com_init()
{
   static char *existsHelp =
   {
       "exists(name)\n"
       "Return TRUE if variable by given name exists otherwise return FALSE.\n"
   };

   static char *whoHelp =
   {
       "who\n"
       "Gives list of currently defined variables.\n"
   };

   static char *formatHelp =
   {
      "format(precision)\n"
      "Set number of digits used in printing values in MATC.\n\n"
   };

   static char *deleteHelp =
   {
      "delete(name)\n"
      "Delete a variable with given name.\n"
   };

   com_init( "exists",  FALSE, FALSE, var_ccheck , 1, 1000, existsHelp );
   com_init( "who"   ,  FALSE, FALSE, var_varlist, 0, 0,    whoHelp    );
   com_init( "format" , FALSE, FALSE, var_format, 1, 2,     formatHelp );
   com_init( "delete",  FALSE, FALSE, var_vdelete, 1, 1,    deleteHelp );
}
