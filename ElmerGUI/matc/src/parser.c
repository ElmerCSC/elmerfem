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
 *     MATC language/expression parser.
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
|  PARSER.C - Last Edited 8. 8. 1988
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
 * $Id: parser.c,v 1.5 2006/11/22 10:57:14 jpr Exp $ 
 *
 * $Log: parser.c,v $
 * Revision 1.5  2006/11/22 10:57:14  jpr
 * *** empty log message ***
 *
 * Revision 1.4  2006/02/02 06:54:44  jpr
 * small formatting changes.
 *
 * Revision 1.2  2005/05/27 12:26:21  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:54  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

static SYMTYPE symbol, bendsym;
static char *str, csymbol[4096], buf[4096];

int char_in_list(int ch, char *list)
{
  char *p;
  
  for(p = list; *p != '\0'; p++)
    if (*p == ch) return TRUE;
  
  return FALSE;
}

void scan() 
{
  char *p, ch;
  int i;
  
  symbol = nullsym;
  if ( *str == '\0' ) return; 
  
  while( isspace(*str) ) str++;
  if (*str == '\0') return; 
  
  p = str;

  if (isdigit(*str) || (*str == '.' && isdigit(*(str+1))))
  {
    str++; while(isdigit(*str)) str++;

    if (*str == '.')
    {
      str++; 
      if (isdigit(*str))
      {
	while(isdigit(*str)) str++;
      }
      else if ( *str != '\0' && *str != 'e' && *str != 'E' && *str != 'd' && *str != 'D'  )
      {
	error("Badly formed number.\n");
      }
    }

    if ( *str == 'd' || *str == 'D' ) *str = 'e';

    if (*str == 'e' || *str=='E' )
    {
      str++; 
      if (isdigit(*str))
      {
        while(isdigit(*str)) str++;
      }
      else if (char_in_list(*str,"+-"))
      {
        str++;
        if (isdigit(*str))
        {
          while(isdigit(*str)) str++;
        }
        else
        {
          error("Badly formed number.\n");
        }
      }
      else
      {
        error("Badly formed number.\n");
      }
    }  
    symbol = number;
  }
  
  else if (isalpha(*str) || char_in_list(*str, symchars))
  {
    while(isalnum(*str) || char_in_list(*str, symchars)) str++;
    ch = *str;  *str = '\0';

    for(i = 0; reswords[i] != NULL; i++)
      if (strcmp(p, reswords[i]) == 0)
      {
        symbol = rsymbols[i]; break;
      }
    if (reswords[i] == NULL) symbol = name;

    *str = ch;
  }
  
  else if (*str == '"')
  {
    str++;
    while(*str != '"' && *str != '\0')
    {
      if (*str++ == '\\') str++;
    }

    if (*str == '\0')
    {
      error("String not terminated.\n");
    }
    str++; symbol = string;
  }
  
  else if (char_in_list(*str, csymbols))
  {  
    for(i = 0; *str != csymbols[i]; i++);
    symbol = ssymbols[i];
    
    str++;

    if (*str == '=')
      switch(symbol)
      {
      case assignsym:
        symbol = eq; str++; break;

      case lt:
        symbol = le; str++; break;

      case gt:
        symbol = ge; str++; break;

      case indclose: case rightpar:
      break;

      default:
        error("Syntax error.\n");
      }
    
    if (*str == '>')
      if (symbol == lt)
      {
        symbol = neq; str++;
      }
  }
  
  else
  {
    error("Syntax error.\n");
  }

  ch  =  *str;
  *str = '\0';

  strcpy( csymbol, p );
  *str = ch;

  return;
}

TREE *newtree()
{
  return (TREE *)ALLOCMEM(sizeof(TREE));
}

TREE *args(minp, maxp)
     int minp, maxp;
{
  TREE *treeptr, *root;
  int numgot = 0;
  
  root = treeptr = equation();
  numgot++;

  while(symbol == argsep)
  {
    scan();
    NEXT(treeptr) = equation();
    treeptr = NEXT(treeptr); 
    numgot++;
    if (numgot > maxp) error("Too many parameters.\n");
  }
  
  if (numgot < minp) error("Too few parameters.\n");
  
  return root;
}
      

TREE *nameorvar()
{
  TREE *root, *treeptr, *prevtree, *tp;

  SYMTYPE sym = nullsym;

  int i, slen;
 
  char *tstr;

  root = treeptr = prevtree = newtree();

  if (symbol == minus && !isspace(*str) &&
     (str-2<buf || isspace(*(str-2)) || char_in_list(*(str-2),"{};=[(\\<>&|+-*/^,")))
  {
    sym = minus; scan();
  }

  if (symbol != name   && symbol != number  && 
      symbol != string && symbol != leftpar)
  {
    error("Expecting identifier, constant or leftpar.\n");
  }

  while(symbol == name   || symbol == number  || 
        symbol == string || symbol == leftpar)
  {

    switch(symbol)
    {
      case name:
        SDATA(treeptr) = STRCOPY(csymbol);
        ETYPE(treeptr) = ETYPE_NAME;
        if (*str == '(' || *str == '[')
        {
          scan(); scan(); ARGS(treeptr) = args(0, 10000);
          if (symbol != rightpar && symbol != indclose)
          {
            error("Expecting closing parenthesis.\n");
          }
        }
      break;

    case string:
      tstr = csymbol + 1;
      tstr[strlen(tstr)-1] = '\0';
      slen = strlen(tstr);
      for(i = 0; i < strlen(tstr); i++)
        if (tstr[i] == '\\') 
          switch(tstr[++i])
          {
            case 'n': break;
            default: slen--;
            break;
          }
      SDATA(treeptr) = (char *)ALLOCMEM(slen+1);
      for(i = 0; *tstr != '\0'; i++, tstr++)
        if (*tstr == '\\')
          switch(*++tstr)
          {
            case 'n':
              SDATA(treeptr)[i++] = '\r'; 
              SDATA(treeptr)[i]   = '\n'; 
            break;

            case 't':
              SDATA(treeptr)[i]   = '\t'; 
            break;

            case 'v':
              SDATA(treeptr)[i]   = '\v'; 
            break;

            case 'b':
              SDATA(treeptr)[i]   = '\b'; 
            break;

            case 'r':
              SDATA(treeptr)[i]   = '\r'; 
            break;

            case 'f':
              SDATA(treeptr)[i]   = '\f'; 
            break;

            case 'e':
              SDATA(treeptr)[i]   = 27; 
            break;

            default:
              SDATA(treeptr)[i] = *tstr;
            break;
          }
        else
          SDATA(treeptr)[i] = *tstr;
      ETYPE(treeptr) = ETYPE_STRING;
      break;

    case number:
      DDATA(treeptr) = atof(csymbol);
      ETYPE(treeptr) = ETYPE_NUMBER;
      break;

    case leftpar:
      scan(); LEFT(treeptr) = equation();
      if (symbol != rightpar)
      {
        error("Right paranthesis missing.\n");
      }
      ETYPE(treeptr) = ETYPE_EQUAT;
      break;
    }

    if (*str == '[')
    {
      scan(); scan(); SUBS(treeptr) = args(1,2);
      if (symbol != rightpar && symbol != indclose)
      {
        error("Expecting closing parenthesis.\n");
      }
    }

    if (sym == minus)
    {
      tp = newtree();
      VDATA(tp) = opr_minus;
      ETYPE(tp) = ETYPE_OPER;
      LEFT(tp) = treeptr;
      if (root == treeptr)
        root = treeptr = tp;
      else 
        LINK(prevtree) = treeptr = tp;
    }

    sym = symbol;
    scan();

    if (symbol == minus && !isspace(*str) &&
         (str-2<buf || isspace(*(str-2)) || char_in_list(*(str-2),"{};=([\\<>&|+-*/^,")))
    {
      sym = minus;

      if (*str == '-' && !isspace(*(str + 1)))
      {
        break;
      }
      else if (*str == '-') 
        error("Syntax error.\n");

      scan();

      if (symbol != name   && symbol != number  && 
          symbol != string && symbol != leftpar)
      {
        error("Expecting identifier, constant or leftpar.\n");
      }
    }

    if (symbol == name   || symbol == number || 
        symbol == string || symbol == leftpar)
    {
      prevtree = treeptr; LINK(treeptr) = newtree(); treeptr = LINK(treeptr);
    }
  }

  return root;
}

TREE *par_apply(root)
	TREE *root;
{ 
  TREE *newroot;

  newroot = newtree();

  switch(symbol)
  { 
    case apply:
      VDATA(newroot) = opr_apply;
    break;

    case not:
      VDATA(newroot) = opr_not;
    break;
  }

  ETYPE(newroot) = ETYPE_OPER;
  scan();

  if (symbol == apply || symbol == not)
    LEFT(newroot) = par_apply(newroot);
  else
    LEFT(newroot) = nameorvar(); 

  return newroot;
}


TREE *par_trans(root)
	TREE *root;
{ 
  TREE *newroot;

  while(symbol == transpose)
  {
    newroot = newtree();
    LEFT(newroot) = root;
    VDATA(newroot) = opr_trans;
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;
    scan();
  }

  return newroot;
}

TREE *par_pow(root)
	TREE *root;
{
  TREE *newroot;

  while(symbol == power)
  {
    newroot = newtree();
    LEFT(newroot) = root; 
    VDATA(newroot) = opr_pow;
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;

    scan(); RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot));
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }

  return newroot;
}

TREE *par_timesdivide(root)
	TREE *root;
{
  TREE *newroot;

  while(symbol == times || symbol == ptimes || symbol == divide)
  {
    newroot = newtree();
    LEFT(newroot) = root;
    switch(symbol)
    {
      case times:
        VDATA(newroot) = opr_mul;
      break;

      case ptimes:
        VDATA(newroot) = opr_pmul;
      break;

      case divide:
        VDATA(newroot) = opr_div;
      break;
    }
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;    

    scan(); RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot));
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot));
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }
  
  return newroot;
}


TREE *par_plusminus(root)
	TREE *root;
{
  TREE *newroot;

  while(symbol == plus || symbol == minus)
  {
    newroot = newtree();
    LEFT(newroot) = root;

    switch(symbol)
    {
      case plus:
        VDATA(newroot) = opr_add; 
      break;

      case minus:
        VDATA(newroot) = opr_subs; 
      break;
    }
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;

    scan(); RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case times: case ptimes: case divide:
        RIGHT(newroot) = par_timesdivide(RIGHT(newroot));
      break;

      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot));
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot));
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }
  
  return newroot;
}

TREE *par_compare(root)
	TREE *root;
{
  TREE *newroot;
  
  while(symbol == eq  || symbol == neq || symbol == lt ||
        symbol == gt  || symbol ==  le || symbol == ge)
  {

    newroot = newtree();
    LEFT(newroot) = root;
    switch(symbol)
    {
      case eq:
        VDATA(newroot) = opr_eq; 
      break;

      case lt:
        VDATA(newroot) = opr_lt; 
      break;

      case gt:
        VDATA(newroot) = opr_gt; 
      break;

      case neq:
        VDATA(newroot) = opr_neq; 
      break;

      case le:
        VDATA(newroot) = opr_le; 
      break;

      case ge:
        VDATA(newroot) = opr_ge; 
      break;
    }
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;

    scan(); RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case plus: case minus:
        RIGHT(newroot) = par_plusminus(RIGHT(newroot));
      break;

      case times: case ptimes: case divide:
        RIGHT(newroot) = par_timesdivide(RIGHT(newroot)); 
      break;

      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot)); 
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot)); 
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }
  
  return newroot;
}

TREE *par_vector(root)
	TREE *root;
{
  TREE *newroot;
  
  while(symbol == vector)
  {
    newroot = newtree();
    LEFT(newroot) = root;
    VDATA(newroot) = opr_vector; 
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;
    scan();
    RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case eq: case neq: case lt: case gt: case le: case ge:
        RIGHT(newroot) = par_compare(RIGHT(newroot)); 
      break;

      case plus: case minus:
        RIGHT(newroot) = par_plusminus(RIGHT(newroot));
      break;

      case times: case ptimes: case divide:
        RIGHT(newroot) = par_timesdivide(RIGHT(newroot)); 
      break;

      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot)); 
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot)); 
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }
  
  return newroot;
}

TREE *par_logical(root)
	TREE *root;
{
  TREE *newroot;
  
  while(symbol == and  || symbol == or)
  {

    newroot = newtree();
    LEFT(newroot) = root;
    switch(symbol)
    {
      case and:
        VDATA(newroot) = opr_and; 
      break;

      case or:
        VDATA(newroot) = opr_or; 
      break;
    }
    ETYPE(newroot) = ETYPE_OPER;
    root = newroot;
    scan(); RIGHT(newroot) = nameorvar();

    switch(symbol)
    {
      case vector:
        RIGHT(newroot) = par_vector(RIGHT(newroot)); 
      break;

      case eq: case neq: case lt: case gt: case le: case ge:
        RIGHT(newroot) = par_compare(RIGHT(newroot)); 
      break;

      case plus: case minus:
        RIGHT(newroot) = par_plusminus(RIGHT(newroot));
      break;

      case times: case ptimes: case divide:
        RIGHT(newroot) = par_timesdivide(RIGHT(newroot)); 
      break;

      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot)); 
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot)); 
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    }
  }
  
  return newroot;
}

TREE *par_reduction(root)
	TREE *root;
{
  TREE *newroot;

  while(symbol == reduction)
  {
    newroot = newtree();
    VDATA(newroot) = opr_reduction;
    ETYPE(newroot) = ETYPE_OPER;
    scan(); RIGHT(newroot) = nameorvar();
    LEFT(newroot) = root;
    root = newroot;

    switch(symbol)
    {
      case and: case or:
        RIGHT(newroot) = par_logical(RIGHT(newroot));
      break;

      case vector:
        RIGHT(newroot) = par_vector(RIGHT(newroot)); 
      break;

      case eq: case neq: case lt: case gt: case le: case ge:
        RIGHT(newroot) = par_compare(RIGHT(newroot)); 
      break;

      case plus: case minus:
        RIGHT(newroot) = par_plusminus(RIGHT(newroot)); 
      break;

      case times: case ptimes: case divide:
        RIGHT(newroot) = par_timesdivide(RIGHT(newroot)); 
      break;

      case power:
        RIGHT(newroot) = par_pow(RIGHT(newroot)); 
      break;

      case transpose:
        RIGHT(newroot) = par_trans(RIGHT(newroot)); 
      break;

      case apply: case not:
        RIGHT(newroot) = par_apply(RIGHT(newroot));
      break;
    } 
  }
  
  return newroot;
}

TREE *par_resize(root)
	TREE *root;
{
  TREE *newroot;

  while(symbol == resize)
  {
    newroot = newtree();
    VDATA(newroot) = opr_resize;
    ETYPE(newroot) = ETYPE_OPER;
    scan(); LEFT(newroot) = nameorvar();
    RIGHT(newroot) = root;
    root = newroot;

    switch(symbol)
    {
      case reduction:
        LEFT(newroot) = par_reduction(LEFT(newroot));
      break;

      case and: case or:
        LEFT(newroot) = par_logical(LEFT(newroot));
      break;

      case vector:
        LEFT(newroot) = par_vector(LEFT(newroot)); 
      break;

      case eq: case neq: case lt: case gt: case le: case ge:
        LEFT(newroot) = par_compare(LEFT(newroot)); 
      break;

      case plus: case minus:
        LEFT(newroot) = par_plusminus(LEFT(newroot)); 
      break;

      case times: case ptimes: case divide:
        LEFT(newroot) = par_timesdivide(LEFT(newroot)); 
      break;

      case power:
        LEFT(newroot) = par_pow(LEFT(newroot)); break;

      case transpose:
        LEFT(newroot) = par_trans(LEFT(newroot)); 
      break;

      case apply: case not:
        LEFT(newroot) = par_apply(LEFT(newroot));
      break;
    } 
  }
  
  return newroot;
}

TREE *equation()
{
  TREE *treeptr;

  switch(symbol)
  {
    case apply: case not:
    break;

    default:
      treeptr = nameorvar();
    break;
  }
  
  while(TRUE)
  {
    switch(symbol)
    {
      case resize:
        treeptr = par_resize(treeptr);
      break;

      case reduction: 
        treeptr = par_reduction(treeptr);
      break;

      case and: case or:
        treeptr = par_logical(treeptr);
      break;

      case vector:
        treeptr = par_vector(treeptr);
      break;

      case eq: case neq: case lt: case gt: case le: case ge:
        treeptr = par_compare(treeptr); 
      break;

      case plus: case minus:
        treeptr = par_plusminus(treeptr);
      break;

      case times: case ptimes: case divide:
        treeptr = par_timesdivide(treeptr);
      break;

      case power: 
        treeptr = par_pow(treeptr);
      break;

      case transpose: 
        treeptr = par_trans(treeptr);
      break;

      case apply: case not:
        treeptr = par_apply(treeptr);
      break;

      default: 
        return treeptr;
    }    
  }
}
 
CLAUSE *commentparse()
{
  char *p = str;

  CLAUSE *root = NULL;
   
  while( *str!='\n' && *str!='\0' ) str++;
  scan();

  return root;
}

CLAUSE *scallparse()
{
  char *p = str;

  CLAUSE *root = NULL;
   
  while( *str!='\n' && *str != ';' && *str!='\0' ) str++;
  if ( *str ) *str++ = '\0';
 
  if ( *p )
  {
      root = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
      root->data = systemcall;

      root->this = newtree();
      SDATA(root->this) = STRCOPY( p );
      ETYPE(root->this) = ETYPE_STRING;
  }

  scan();

  return root;
}

CLAUSE *statement()
{
  char *csymbcopy, *p;

  CLAUSE *root = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));

  if (symbol == name)
  {
    p = str; 
    csymbcopy = STRCOPY(csymbol);

    do
    {
       scan();
    } while( symbol != assignsym && symbol != nullsym && symbol != statemend );

    strcpy(csymbol, csymbcopy);
    FREEMEM(csymbcopy);
    str = p;

    if (symbol == assignsym)
    {
      symbol = name; root -> this = nameorvar(); scan();
    }
    else 
      symbol = name;
  }

  LINK(root) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  LINK(root) -> this = equation();

  root->data = assignsym;

  return root;
}

CLAUSE *blockparse()
{
  CLAUSE *root, *ptr;

  root = (CLAUSE *)NULL;

  if (symbol != beginsym)
    error("if|while|function: missing block open symbol.\n");

  scan();

  if (symbol == nullsym)
  {
    dogets(str, PMODE_BLOCK);
    scan();
  }

  if (symbol != endsym)
  {
    root = ptr = parse();
    while(LINK(ptr) != NULL)
    {
      ptr = LINK(ptr);
    }
  }

  while(symbol != endsym && symbol != elsesym)
  {
    if (symbol == nullsym)
    {
      dogets(str, PMODE_BLOCK); scan();
    }
    if (symbol != endsym && symbol != elsesym)
    {
      LINK(ptr) = parse();
      while(LINK(ptr) != NULL)
      {
        ptr = LINK(ptr);
      }
    }
  }
         
  bendsym = symbol;
  scan();

  return root;
}

CLAUSE *funcparse()
{
  CLAUSE *root, *ptr;
  SYMTYPE sym;
  TREE *lptr, *rptr,*help;

  int ch,n;

  char *p = str;

  root = ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  ptr->data = funcsym;

  scan();
  ptr->this = nameorvar();

  help = SUBS(root->this) = newtree();
  SDATA( help ) = STRCOPY( p );
  p = str;

  while ( symbol == nullsym || symbol == comment )
  {
      dogets( str, PMODE_CONT );
      scan();

      if ( symbol == comment )
      {
          NEXT(help) = newtree();
          help = NEXT(help);

          while( *str != '\n' && *str != '\0' ) str++;
          ch = *str;
          if ( *str ) *++str = '\0';
          *str = ch;
          SDATA(help) = STRCOPY( p ); 

          p = str;
      }
  }

  while(symbol == import || symbol == export)
  {
    if (symbol == import) 
      lptr = LEFT(root->this);
    else
      lptr = RIGHT(root->this);

    sym = symbol;
    scan();
    rptr = args(1,1000);

    if (lptr == NULL)
    { 
      if (sym == import) 
        LEFT(root->this) = rptr;
      else
        RIGHT(root->this) = rptr;
    }
    else
    {
      while(NEXT(lptr)) lptr=NEXT(lptr);
      NEXT(lptr) = rptr;
    }

    if (symbol == nullsym)
    {
      dogets(str, PMODE_CONT);
      scan();
    }
  }

  if (symbol == beginsym)
  {
    LINK(ptr) = blockparse();
    if (bendsym != endsym)
      error("function: missing end.\n");
  }
  else
    LINK(ptr) = parse();

  return root;
}

CLAUSE *ifparse()
{
  CLAUSE *root, *ptr, *parse();
  int block = FALSE;

  scan();
  if (symbol != leftpar)
  {
    error("Missing leftpar.\n");
  }

  root = ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  ptr->data = ifsym; 

  scan();
  ptr -> this = equation();

  if (symbol != rightpar)
  {
    error("Missing rightpar.\n");
  }
  scan();

  if (symbol == thensym) scan();

  if (symbol == nullsym)
  {
    dogets(str, PMODE_CONT);
    scan();
  }

  if (symbol == beginsym) 
  {
    block = TRUE;
    LINK(ptr) = blockparse();
  }
  else
    LINK(ptr) = parse();

  while(LINK(ptr) != NULL)
  {
    ptr = LINK(ptr);
  }

  root->jmp = LINK(ptr) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  ptr = LINK(ptr); ptr->data = endsym;
  
  if (symbol == elsesym || bendsym == elsesym)
  { 
    root -> jmp = LINK(ptr) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
    ptr = LINK(ptr); ptr->data = elsesym; 

    if (symbol == elsesym) scan();  

    if (symbol == nullsym)
    {
      dogets(str, PMODE_CONT);
      scan();
    }

    if (symbol == beginsym)
    {
      LINK(ptr) = blockparse();
      if (block && bendsym != endsym)
        error("else: missing end.\n");
    }
    else
      LINK(ptr) = parse();

    while(LINK(ptr) != NULL)
    {
      ptr = LINK(ptr);
    }
    root->jmp->jmp = LINK(ptr) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
    LINK(ptr)->data = endsym;
  }

  return root;
}

CLAUSE *whileparse()
{
  CLAUSE *root, *ptr;

  scan();

  if (symbol != leftpar)
  {
    error("Missing leftpar.\n");
  }

  root = ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  ptr->data = whilesym; 

  scan();
  ptr->this = equation();

  if (symbol != rightpar)
  {
    error("Missing rightpar.\n");
  }
  scan();

  if (symbol == nullsym)
  {
    dogets(str, PMODE_CONT);
    scan();
  }

  if (symbol == beginsym)
  {
    LINK(ptr) = blockparse();
    if (bendsym != endsym) 
      error("while: missing end.\n");
  }
  else
    LINK(ptr) = parse();

  while(LINK(ptr) != NULL)
  {
    ptr = LINK(ptr);
  }

  root -> jmp = LINK(ptr) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  LINK(ptr)->data = endsym;

  return root;
}

CLAUSE *forparse()
{
  CLAUSE *root, *ptr;

  scan();

  if (symbol != leftpar)
  {
    error("for: missing leftpar.\n");
  }

  root = ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  ptr->data = forsym; 

  scan();
  ptr -> this = nameorvar();
  if (symbol != assignsym)
  {
     error("for: missing equalsign\n");
  }
  scan();

  LINK(ptr->this) = equation();

  if (symbol != rightpar)
  {
    error("Missing rightpar.\n");
  }
  scan();

  if (symbol == nullsym)
  {
    dogets(str, PMODE_CONT);
    scan();
  }

  if (symbol == beginsym)
  {
    LINK(ptr) = blockparse();
    if (bendsym != endsym) 
      error("for: missing end.\n");
  }
  else
    LINK(ptr) = parse();

  while(LINK(ptr) != NULL)
  {
    ptr = LINK(ptr);
  }

  root -> jmp = LINK(ptr) = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));
  LINK(ptr)->data = endsym;

  return root;
}

CLAUSE *parse()
{
  CLAUSE *ptr = (CLAUSE *)NULL;

  switch(symbol)
  {
    case funcsym:
      ptr = funcparse();
    break;

    case beginsym:
      ptr = blockparse();
      if (bendsym != endsym)
        error("begin: missing end.\n");
    break;

    case ifsym:
      ptr = ifparse();
    break;

    case whilesym:
      ptr = whileparse();
    break;

    case forsym:
      ptr = forparse();
    break;

    case systemcall:
      ptr = scallparse();
    break;

    case comment:
      ptr = commentparse();
    break;

    default:
      ptr = statement();
    break;
  }

  while( symbol == statemend ) scan();

  if (ptr == (CLAUSE *)NULL)
    ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));

  return ptr;
}

void free_treeentry(root)
   TREEENTRY *root;
{
   if (root == NULL) return;

   free_tree(root->args);

   free_tree(root->subs);
   if ( root->entrytype == ETYPE_STRING || root->entrytype == ETYPE_NAME )
        FREEMEM(root->entrydata.s_data);
   else if ( root->entrytype == ETYPE_CONST )
        var_delete_temp(root->entrydata.c_data);
}

void free_tree(root)
   TREE *root;
{
   if (root == NULL) return;

   free_tree(NEXT(root));
   free_tree(LINK(root));
   free_tree(LEFT(root));
   free_tree(RIGHT(root));
   free_treeentry(&root->tentry);
   FREEMEM((char *)root);
}

void free_clause(root)
    CLAUSE *root;
{
    if (root == NULL) return;

    free_clause(LINK(root));
    free_tree(root->this);
    FREEMEM((char *)root);
}

VARIABLE *doit(line)
	char *line;
{
  CLAUSE *ptr, *root;
  VARIABLE *res;

  str = buf;
  strcpy( str, line );

  root = ptr = (CLAUSE *)ALLOCMEM(sizeof(CLAUSE));

  scan();

  while(symbol != nullsym)
  {
    LINK(ptr) = parse();
    while(LINK(ptr) != NULL)
    {
      ptr = LINK(ptr);
    }
  }

/*  root = optimclause(root); */
/*  printclause(root, math_out, 0);   */
  res = evalclause(root);

  free_clause(root);

  return res;
}
