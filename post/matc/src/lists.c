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
 *     List handling utilities.
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
|  LISTS.C - Last Edited 7. 8. 1988
|
|  Handling of global lists. Starts of global lists are hold in a
|  global array named listheaders and type of LIST. For each routine
|  managing one of lists you provide an index to this array. There
|  are definitions in MATC.H to make these indexes more abstract
|  entities.
|
|  List structure is currently as follows
|
|       struct list 
|       {
|         struct list *next;
|         char *name;
|       };
|
|  This is typedef'ed to LIST.
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
 * $Id: lists.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: lists.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:44  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

void lst_addtail(list, item) int list; LIST *item;
/*======================================================================
?  Add specified item to end of a list given.
^=====================================================================*/
{
   LIST *lst;

  /*
   *   check if the list exists, if not just make this item first in list
   */
   if (listheaders[list].next == (LIST *)NULL)
     listheaders[list].next = item;

  /*
   * else look for current last item in list
   */
   else
   { 
     for(lst = listheaders[list].next; NEXT(lst); lst = NEXT(lst));
 
     /*
      *  and make the link
      */
     NEXT(lst) = item;
   }
}

void lst_addhead(list, item) int list; LIST *item;
/*======================================================================
?  add specified item to start of list given.
^=====================================================================*/
{
  /*
   *  make the link.
   */
   NEXT(item) = listheaders[list].next;
   listheaders[list].next = item;
}

void lst_add(list, item) int list; LIST *item;
/*======================================================================
?  add item to lexically right place in the list.
|
&  strcmp()
^=====================================================================*/
{
   LIST *lst, *lstn;
   
   /* 
    *  if the list is empty make this item first and return.
    */
   if ((lst = listheaders[list].next) == (LIST *)NULL)
   {
     lst_addhead(list, item); return;
   }

   /* 
    *  if the name of the new item is lexically 
    *  smaller than first item in the list, add it
    *  to beginning of the list.
    */
   if (strcmp(NAME(lst), NAME(item)) > 0)
   {
     lst_addhead(list, item); return;
   }
 
   /*
    *  look for right place to add.
    */
   for(; NEXT(lst); lst = NEXT(lst)) 
     if (strcmp(NAME(NEXT(lst)), NAME(item)) > 0)
     {
       lstn = NEXT(lst); 
       NEXT(lst) = item; 
       NEXT(item) = lstn;
       return;
     }
 
   /*
    *  fell of the loop. item should be added to the tail of the list.
    */
   NEXT(lst) = item; 
}

void lst_unlink(list, item) int list; LIST *item;
/*======================================================================
?  unlink specified item from a list given.
^=====================================================================*/
{
  LIST *lst;

  /*
   *  if the list is empty return
   */
  if ((lst = listheaders[list].next) == (LIST *)NULL) return;
 
  /*
   *  it's not the header, look if it is in list at all  
   */
  if (lst != item)
  {

    for(; NEXT(lst); lst = NEXT(lst))
    {
      if (NEXT(lst) == item) break;
    }

    /*
     *  item was not found from the list. do nothing.
     */
    if (NEXT(lst) == (LIST *)NULL) return;

    /*
     *   found, unlink
     */
    NEXT(lst) = NEXT(item);
  }

  /*
   *  item was the header, unlink it    
   */
  else
    listheaders[list].next = NEXT(item);
}

void lst_free(list, item) int list; LIST *item;
/*======================================================================
?  Unlink item from list and free memory used by it.
|
&  lst_unlink(), FREEMEM
^=====================================================================*/
{
   lst_unlink(list, item);

   FREEMEM(NAME(item));
   FREEMEM((char *)item);
}

LIST *lst_find(list, name) int list; char *name;
/*======================================================================
?  Look for a named item from given list.
|
&  strcmp()
^=====================================================================*/
{
  LIST *lst;

  /*
   *   look for item
   */
  for( lst = listheaders[list].next; lst; lst = NEXT(lst) )
  {
      if ( NAME(lst) && strcmp(name, NAME(lst)) == 0 ) break;
  }

  return lst;
}

void lst_purge(list) int list;
/*======================================================================
?  Delete list and free memory allocated to it.
|
&  FREEMEM
^=====================================================================*/
{
  LIST *lst, *lstn;

  /*
   *  free memory allocated for this list
   */
  for(lst = listheaders[list].next; lst;)
  {
    lstn = NEXT(lst);
    FREEMEM(NAME(lst));
    FREEMEM((char *)lst);
    lst = lstn;
  }

  listheaders[list].next = (LIST *)NULL;  /* security */
}

VARIABLE *lst_print(list) int list;
/*======================================================================
?  Print list name and item names from given list
|
!  Output looks real ugly, should do something to that. The command
|  to get here is "help" but it really is not very helpful.
|
&  fprintf(), strlen()
^=====================================================================*/
{
   LIST *lst;      

   int i, spc, len;
  
   /*
    *  if empty list return.
    */
   if ( listheaders[list].next == (LIST *)NULL )
       return (VARIABLE *)NULL;

   /*
    *  name of the list
    */
    PrintOut( "\n%s\n\n", listheaders[list].name );

   /*
    * and finally list all item names one by one.
    * try printing as many names as possible to
    * each output line.
    */
    for(len = 0,lst = listheaders[list].next; lst; lst = NEXT(lst))
    {
      if ( NAME(lst) )
      {
        if (len >= 80)
        {
          PrintOut("\n");
          len = 0;
        } else len += 20;
        PrintOut("%-20s\t", NAME(lst));
		if ( strlen(NAME(lst)) >= 20 ) { PrintOut("%-20%s\t", " "); len+= 20; }
		/*
        spc = 20 - spc;
        for(i=0; i<spc && len<80; i++)
        {
          PrintOut(" "); len++;
        }
		*/
      }
    }
    PrintOut("\n");
 
    return (VARIABLE *)NULL;
}
