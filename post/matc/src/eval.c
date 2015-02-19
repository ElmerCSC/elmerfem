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
 *     Evaluate expression threes in format parsed by parser.c.
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
|  EVAL.C - Last Edited 6. 8. 1988
|
***********************************************************************/
#include "../config.h"
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
 * $Id: eval.c,v 1.5 2006/02/02 06:51:16 jpr Exp $ 
 *
 * $Log: eval.c,v $
 * Revision 1.4  2005/08/25 13:44:22  vierinen
 * windoze stuff
 *
 * Revision 1.3  2005/05/27 12:26:20  vierinen
 * changed header install location
 *
 * Revision 1.2  2005/05/26 12:34:53  vierinen
 * windows stuff
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:35  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"

VARIABLE *evaltree(root) TREE *root;
/*======================================================================
?  Evaluate the equation tree (real hard). The tree might be a list of 
|  trees. If there is more than one tree in the list, the return value
|  will be a single column vector, the result from each tree is 
|  concatenated to this vector. If there's only one tree, return
|  value is tree's evaluated result.
|
&  strlen(), fopen(), fclose(), var_temp_new(), var_temp_copy(), 
|  var_delete_temp(), com_check(), var_check(), fnc_check(), 
|  fnc_exec(), com_pointw(), com_source()
|
^=====================================================================*/
{

  VARIABLE *first,     /* used to hold start address for    */
                       /* temporary list of VARIABLES, the  */
                       /* evaluated trees results           */
            *subs,     /* indexes for tree elements         */
            *par,      /* parameters for tree elements      */
            *res,      /* result is returned in this var.   */ 
                       /* probably used as a temporary too  */ 
            *leftptr,  /* result from left leaf             */
            *rightptr, /* result from right leaf            */
            *tmp,      /* can't manage without these        */
            *tmp1;

  FUNCTION *fnc;       /* needed, if someone's using defined functions */

  COMMAND *com;        /* needed, maybe */

  MATRIX *opres;       /* don't know really */

  TREE *argptr;        /* temporary */

  int i,               /* there's always reason for these         */
      dim,             /* the final dimension of the return value */
      argcount;        /* parameter count                         */

  char *resbeg;        /* for memcpying */

  FILE *fp;            /* used to check if a file exists */


  if (root == NULL) return NULL;

  dim = 0;  first = NULL;

/*-------------------------------------------------------------------*/

  /* 
      while there's trees in the list.
  */
  while(root) 
  {

    subs = par = tmp = NULL;

    /*
        check if the result will be indexed
    */
    argptr = SUBS(root);
    if (argptr)
    {
      subs = tmp = evaltree(argptr);
      argptr = NEXT(argptr);
      while(argptr)
      {
        NEXT(tmp) = evaltree(argptr);
        argptr = NEXT(argptr); tmp = NEXT(tmp);
      }
    }

    switch(ETYPE(root))
    {

    /******************************************************
              some kind of existing identifier.
    *******************************************************/
    case ETYPE_NAME:
      /*
          is there parameters for this one ?
          and how many ?
      */
      argptr = ARGS(root);
      argcount = 0;
      if (argptr)
      {
        par = tmp = evaltree(argptr);
        argptr = NEXT(argptr);
        argcount++;
        while(argptr)
        {
          argcount++;
          NEXT(tmp) = evaltree(argptr);
          argptr = NEXT(argptr); tmp = NEXT(tmp);
        }
      }


      /* 
          one of the builtin commands (sin, ones, help, ...)
      */
      if ((com = com_check(SDATA(root))) != (COMMAND *)NULL)
      {

        if (argcount < com->minp || argcount > com->maxp)
        {
          if (com->minp == com->maxp)
          {
            error( "Builtin function [%s] requires %d argument(s).\n", SDATA(root), com->minp);
          }
          else
          {
             error("Builtin function [%s] takes from %d to %d argument(s).\n",
                       SDATA(root), com->minp, com->maxp);
          }
        }


        if (com->flags & CMDFLAG_PW)
        {
          tmp = com_pointw((double (*)())com->sub, par);
        }
        else
        {
          tmp = (*com->sub)(par);
        }
      }


      /*
          a variable name
      */
      else if ((tmp1 = var_check(SDATA(root))) != (VARIABLE *)NULL)
      {
        tmp = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
        tmp ->this = tmp1->this;
        REFCNT(tmp)++;
        if (par != NULL)
        {
          subs = par; par = NULL;
        }
      }

      /*
          user defined function
      */
      else if ((fnc = fnc_check(SDATA(root))) != (FUNCTION *)NULL)
      {
        tmp = fnc_exec(fnc, par); par = NULL;
      }


      /*
           maybe a file name ?
      */
      else if ((fp = fopen(SDATA(root),"r")) != (FILE *)NULL)
      {
        fclose(fp); 
        tmp = var_temp_new(TYPE_STRING, 1, strlen(SDATA(root)));
        for(i = 0; i < strlen(SDATA(root)); i++)
        {
          M(tmp, 0, i) = SDATA(root)[i];
        }
        (*com_source)(tmp);
        var_delete_temp(tmp);
        tmp = NULL;
      }

      /*
          troubles!
      */
      else
      {
        error("Undeclared identifier: [%s].\n", SDATA(root));
      }

      break;

    
    /******************************************************
                  single string constant
    *******************************************************/
    case ETYPE_STRING:
      tmp = var_temp_new(TYPE_STRING, 1, strlen(SDATA(root)));
      for(i = 0; i < strlen(SDATA(root)); i++)
      {
        M(tmp,0,i) = SDATA(root)[i];
      }
      break;

    /******************************************************
                  single numeric constant
    *******************************************************/
    case ETYPE_NUMBER:
      tmp = var_temp_new(TYPE_DOUBLE, 1, 1);
      M(tmp,0,0) = DDATA(root);
      break;

    /******************************************************
                    constant matrix
    *******************************************************/
    case ETYPE_CONST:
      tmp = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
      tmp->this = CDATA(root)->this; 
      REFCNT(tmp)++;
      break;

    /******************************************************
                           huh ?
    *******************************************************/
    case ETYPE_EQUAT:
      tmp = evaltree(LEFT(root));
      break;

    /******************************************************
         left oper [right]
         oper = divide, multiply, transpose, power,...
    *******************************************************/
    case ETYPE_OPER:
      /*
       *   evaluate both leafs.
       */
      leftptr = rightptr = NULL; opres = NULL;

      leftptr  = evaltree(LEFT(root));
      rightptr = evaltree(RIGHT(root));

      if (leftptr && rightptr) 
        opres = (*VDATA(root))(leftptr->this, rightptr->this);
      else if (leftptr)
        opres = (*VDATA(root))(leftptr->this, NULL);
      else if (rightptr)
        opres = (*VDATA(root))(rightptr->this, NULL);

      var_delete_temp(leftptr);
      var_delete_temp(rightptr);

      if (opres != NULL)
      {
        tmp = (VARIABLE *)ALLOCMEM(VARIABLESIZE);
        tmp->this = opres;
        REFCNT(tmp) = 1;
      }

      break;
    }

    /*
       if NULL result, don't really know what to do, so we try quitting.
    */
    /*
    if (tmp == (VARIABLE *)NULL)
    {
      if (subs) var_delete_temp(subs);
      if (par) var_delete_temp(par);
      if (first) var_delete_temp(first);
      return (VARIABLE *)NULL;
    }
    */

    /******************************************************
       deleteting temporaries, preparing for next loop
    *******************************************************/
    if (subs != NULL)
    {
      if (tmp != NULL)
      {
        tmp1 = tmp;
        NEXT(tmp1) = subs;
        tmp = com_el(tmp1);
        var_delete_temp(tmp1);
      }
      else
      {
        var_delete_temp(subs);
      }
      tmp1 = subs = NULL;
    }

    if (first == NULL)
    {
      first = res = tmp;
    }
    else if (tmp != NULL)
    {
      NEXT(res) = tmp; res = NEXT(res);
    }

    if (subs) var_delete_temp(subs);
    if (par) var_delete_temp(par); 

    if (tmp != NULL) dim += NROW(tmp) * NCOL(tmp);/* count size of the result */

    root = LINK(root);
  }
 /*-------------------------------------------------------------------*/

  /* 
      there was only one tree, so return  it's result.
  */
  if (tmp == first) return first;

  /* 
      it was a list of trees, concatenate the results
  */
  res = var_temp_new(TYPE(first), 1, dim);   

  resbeg = (char *)MATR(res);
  for(tmp = first; tmp; tmp = NEXT(tmp))
  {
    memcpy(resbeg, (char *)MATR(tmp), MATSIZE(tmp));
    resbeg += MATSIZE(tmp);
  }
  /* 
      and delete rest of the temporaries.
  */
  var_delete_temp(first);

  return res;
}

VARIABLE *evalclause(root) CLAUSE *root;
/*======================================================================
?  Evaluate the operations list. The list contains equations trees
|  (actually lists of trees :-).
|
&  ALLOCMEM, FREEMEM, STRCOPY, var_temp_new(), var_delete_temp(),
&  var_check(), com_check(), fnc_check(), evaltree(), 
^=====================================================================*/
{       
  VARIABLE *ptr = NULL,
           *res, *par, *tmp;      /* used for something magic */

  TREE *argptr;                   /* pointer to parameters    */

  double *d;

  int i;

/*-------------------------------------------------------------------*/

  while(root)
  {
    if (root->data == endsym) return ptr;
 
    switch(root->data)
    {
      /************************************************************
                       System call
      ************************************************************/
      case systemcall:
      {
#if defined(WIN32) || defined(MINGW32)
           FILE *fp = _popen( SDATA(root->this), "r" );
#else
           FILE *fp = popen( SDATA(root->this), "r" );
#endif
           static char s[101];

           if ( !fp ) error( "systemcall: open failure: [%s].\n", SDATA(root->this) );

           while( fgets( s, 120, fp ) ) PrintOut( s );
  
#if defined(WIN32) || defined(MINGW32)
           _pclose( fp );
#else
           pclose( fp );
#endif
      }
      break;
      /************************************************************
                       Function definition 
      ************************************************************/
      case funcsym:
      {

        FUNCTION *fnc;                    /* pointer to created function */
        TREE *tptr;                       /*  function parameters        */
        char *name = SDATA(root->this);   /* function name */
        int i, n, argcount;               /* well... */

        /*
            check for name conflicts
        */
        if (var_check(name) || com_check(name))
        {
          error( "Function not created [%s], identifier in use.\n",name);
        }
          
        /*
         *   allocate mem for FUNCTION structure and add it
         *   to the the FUNCTIONS list
         */
        if (fnc = fnc_check(name)) fnc_free_entry(fnc);
        fnc = (FUNCTION *)ALLOCMEM(sizeof(FUNCTION));
        NAME(fnc) = STRCOPY(name);
        lst_add(FUNCTIONS, (LIST *)fnc);

        /*
         *   copy parameter names to the structure.
         */
        argcount = 0;
        for(tptr = ARGS(root->this); tptr; tptr = NEXT(tptr)) argcount++; 
        if (argcount > 0)
        {
          fnc -> parnames = (char **)ALLOCMEM(argcount * sizeof(char *));
          for(i = 0, tptr = ARGS(root->this); tptr; tptr = NEXT(tptr))
            fnc -> parnames[i++] = STRCOPY(SDATA(tptr));
        }
        fnc -> parcount = argcount;

        /*
         *   copy help text if any to the structure.
         */
        argcount = n = 0;
        for( tptr = SUBS(root->this); tptr; tptr = NEXT(tptr))
        {
            if ( SDATA(tptr) )
            {
                argcount++;
                n += strlen( SDATA(tptr) );
            }
        }
        if ( argcount > 0 && n > 0)
        {
            fnc->help = (char *)ALLOCMEM( (n+argcount+1)*sizeof(char) );
            for( tptr = SUBS(root->this); tptr; tptr = NEXT(tptr) )
            {
                if ( SDATA(tptr) )
                {
                    strcat( fnc->help, SDATA(tptr) );
                    strcat( fnc->help, "\n" );
                }
            }
        }
  
        /*
         *   copy imported variable names to the structure.
         */
        argcount = 0;
        for(tptr = LEFT(root->this); tptr; tptr = NEXT(tptr)) argcount++; 
        if (argcount > 0)
        {
          fnc -> imports = (char **)ALLOCMEM((argcount+1) * sizeof(char *));
          for(i = 0, tptr = LEFT(root->this); tptr; tptr = NEXT(tptr))
            fnc -> imports[i++] = STRCOPY(SDATA(tptr));
          fnc -> imports[i] = NULL;
        }
        else
          fnc -> imports = NULL;

        /*
         *   copy exported variable names to the structure.
         */
        argcount = 0;
        for(tptr = RIGHT(root->this); tptr; tptr = NEXT(tptr)) argcount++; 
        if (argcount > 0)
        {
          fnc -> exports = (char **)ALLOCMEM((argcount+1) * sizeof(char *));
          for(i = 0, tptr = RIGHT(root->this); tptr; tptr = NEXT(tptr))
            fnc -> exports[i++] = STRCOPY(SDATA(tptr));
          fnc -> exports[i] = NULL;
        }
        else
          fnc -> exports = NULL;
 
        /*
            and finally the function body
        */
        fnc -> body = LINK(root); LINK(root) = NULL;
  
/*        PrintOut( "Function defined: [%s].\n", name); */
 
        return NULL;
      }
     

      /***************************************************************
                             statement
      ****************************************************************/
      case assignsym:
      {

        int  iflg = FALSE,          /* is the result indexed   */
             pflg = TRUE;           /* should it be printed to */
                                    /* output stream           */

        char *r = "ans";            /* resulted VARIABLE name  */

        par = NULL;

        /*
         *   there is an explicit assigment  (ie. x = x + 1, not x + 1)
         */
        if (root->this)
        {
          /*
           *  VARIABLE name
           */
          r = SDATA( root->this );

          /*
           *  check for name conflicts 
           */
          if ( fnc_check(r) || com_check(r) || lst_find(CONSTANTS, r) )
          {
              error( "VARIABLE not created [%s], identifier in use.\n", r );
          }
  
          /*
           *   is it indexed ?
           */
          pflg = FALSE; 
          argptr = ARGS(root->this);
          if (argptr)
          {
            iflg = TRUE;
            if ((par = tmp = evaltree(argptr)) != NULL)
            { 
              argptr = NEXT(argptr);
              while(argptr)
              {
                if ((NEXT(tmp) = evaltree(argptr)) == NULL) break;
                argptr = NEXT(argptr);
                tmp = NEXT(tmp);
              }
            }
          }
        }
  
        /*
         *  evaluate the right side of the statement
         *  and put the result where it belongs
         */
        ptr = evaltree( LINK(root)->this );
        ptr = put_result( ptr, r, par, iflg, pflg );
   
        if ( par ) var_delete_temp( par );
        root = LINK(root);
   
      break;
      }

      /***************************************************************
                             if statement
      ****************************************************************/
      case ifsym:

        if ((res = evaltree(root->this)) != NULL)
        {
          for(i = 0, d=MATR(res); i < NROW(res)*NCOL(res); i++)
            if (*d++ == 0) break;

          /*
           *   condition false
           */
          if (*--d == 0)
          {
            root = root->jmp;
            if (root->data == elsesym)
            {
              ptr = evalclause(LINK(root));
              root = root -> jmp;
            }
          }
  
          /*
           *   condition true
           */
          else
          {
            ptr = evalclause(LINK(root));
            root = root->jmp;
            if (root->data == elsesym)
            {
              root = root->jmp;
            }
          }
          var_delete_temp(res);
        }
        else
        {
          root = root -> jmp;
          if (root->data == elsesym)
          {
            root = root->jmp;
          }
        }
      break;
  
      /***************************************************************
                             while statement
      ****************************************************************/
      case whilesym:
  
        while(TRUE)
        {
          if ((res = evaltree(root->this)) == NULL) break;

          for(i=0, d=MATR(res); i < NROW(res)*NCOL(res); i++)
            if (*d++ == 0) break;
  
          /*
              condition true, go for another loop
          */
          if (*--d != 0)
          {
            ptr = evalclause(LINK(root));
            var_delete_temp(res);
          }
  
          /*
              condition false, done with while-loop
          */
          else
          {
            var_delete_temp(res);
            break;
          }
        }
        root = root->jmp;
      break;
  
      /***************************************************************
                             for statement
      ****************************************************************/
      case forsym:
      {
        VARIABLE *var;
        char *r;
  
        /*
         *  VARIABLE name
         */
        r = SDATA(root->this);
  
        /*
         *  check for name conflicts 
         */
        if (fnc_check(r) || com_check(r) || lst_find(CONSTANTS, r))
        {
          error( "VARIABLE not created [%s], identifier in use.\n ", r);
        }
  
        if ((res = evaltree(LINK(root->this))) != NULL) 
        {
          var_delete(r);
          var = var_new(r,TYPE(res),1,1);
  
          d = MATR(res);
          for(i = 0; i < NCOL(res)*NROW(res); i++) 
          {
            *MATR(var) = *d++;
            ptr = evalclause(LINK(root));
          }

          var_delete_temp(res);
        }
        root = root->jmp;
      break;
      }
    } 
    root = LINK(root);
  }
  return ptr;
}

VARIABLE *put_values(ptr, resname, par) VARIABLE *ptr, *par; char *resname;
/*======================================================================
?  extract values from ptr, indexes in par, and put them to res
^=====================================================================*/
{
  static double defind = 0.0;

  double *ind1,
         *ind2,
         *dtmp;

  int i, j, k,
      rows, cols,
      size1, size2,
      imax1, imax2, 
      ind;

  VARIABLE *res;
  
  res = (VARIABLE *)lst_find(VARIABLES, resname);

  if (NEXT(par) == NULL)
  {
    if (
      res != NULL && NROW(par) == NROW(res) && NCOL(par) == NCOL(res)
                  && !(NROW(res) == 1 && NCOL(res) == 1)
      )
    {
      int logical = TRUE,
          csize   = 0;

      dtmp = MATR(par);
      for(i = 0; i < NROW(par)*NCOL(par); i++)
        if (dtmp[i] != 0 && dtmp[i] != 1)
        {
          logical = FALSE;
          break;
        }

      if (logical)
      {
        imax1 = NROW(ptr) * NCOL(ptr);
        dtmp  = MATR(ptr);
        for(i = 0, k = 0; i < NROW(res); i++)
          for(j = 0,csize=0; j < NCOL(res); j++)
          {
            while(M(par,i,j)==1 && j+csize<NCOL(res) && k+csize<imax1) csize++;
            if (csize > 0)
            {
              memcpy(&M(res,i,j),&dtmp[k],csize*sizeof(double));
              j += csize-1;
              k += csize;
              csize = 0; 
              if (k >= imax1) k = 0;
            }
          }
         var_delete_temp(ptr);
         return res;
      }
      else
      {
        ind1  = &defind;
        size1 = 1;
        ind2  = MATR(par);
        size2 = NCOL(par); 
      }
    }
    else
    {
      ind1  = &defind;
      size1 = 1;
      ind2  = MATR(par);
      size2 = NCOL(par); 
    }
  }
  else
  {
    ind1  = MATR(par);
    size1 = NCOL(par); 
    ind2  = MATR(NEXT(par));
    size2 = NCOL(NEXT(par));
  }

  imax1 = (int)ind1[0];
  for(i = 1; i < size1; i++)
  {
    imax1 = max(imax1, (int)ind1[i]);
  }

  imax2 = (int)ind2[0];
  for(i = 1; i < size2; i++)
  {
    imax2 = max(imax2, (int)ind2[i]);
  }

  if (res == NULL)
  {
    res = var_new(resname, TYPE(ptr), imax1+1, imax2+1);
  }
  else if (NROW(res) <= imax1 || NCOL(res) <= imax2)
  {
    int ir = NROW(res), jc = NCOL(res);
    MATRIX *t;

    imax1 = max(ir, imax1 + 1);
    imax2 = max(jc, imax2 + 1);

    t = mat_new(TYPE(res), imax1, imax2);
    dtmp = t->data;

    for(i = 0; i < ir; i++)
    {
      memcpy(&dtmp[i*imax2],&M(res,i,0),jc*sizeof(double));
    }

    if (--REFCNT(res) == 0)
    {
      mat_free(res->this);
    }
    res->this = t;
    REFCNT(res) = 1;
  }
  else if (REFCNT(res) > 1)
  {
    --REFCNT(res);
    res->this = mat_copy(res->this);
  }

  imax1 = NROW(ptr) * NCOL(ptr);
  dtmp = MATR(ptr);
  for(i = 0, k = 0; i < size1; i++)
  {
    ind = (int)ind1[i];
    for(j = 0; j < size2; j++)
    {
      memcpy(&M(res,ind,(int)ind2[j]),&dtmp[k++],sizeof(double));
      if (k >= imax1) k = 0;
    }
  }

  var_delete_temp(ptr);

  return res;
}

VARIABLE *put_result(ptr, resname, par, indexflag, printflag)
/*======================================================================
?  copy VARIABLE from one place to another
|  and conditionally print the result  
^=====================================================================*/
   int indexflag, printflag;
   VARIABLE *ptr, *par;
   char *resname;
{                       
   VARIABLE *res = NULL;

  var_delete( "ans" );
  if (indexflag && par)
  {
    res = put_values(ptr, resname, par);
  }
  else
  {
    res = var_rename( ptr, resname );
  }

  if ( res ) res->changed = 1;
  if (printflag) var_print(res);

  return res;
}
