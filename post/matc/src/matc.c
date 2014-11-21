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
 *     MATC main module.
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
|  MATC - Last Edited 9. 8. 1988
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
 * $Id: matc.c,v 1.7 2007/06/08 08:12:17 jpr Exp $ 
 *
 * $Log: matc.c,v $
 * Revision 1.7  2007/06/08 08:12:17  jpr
 * *** empty log message ***
 *
 * Revision 1.6  2006/02/07 10:21:42  jpr
 * Changed visibility of some variables to local scope.
 *
 * Revision 1.5  2006/02/02 06:54:44  jpr
 * small formatting changes.
 *
 * Revision 1.3  2005/08/25 13:44:22  vierinen
 * windoze stuff
 *
 * Revision 1.2  2005/05/27 12:26:20  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:48  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#define MODULE_MATC
#include "elmer/matc.h"
#include "str.h"
#include "../config.h"

#ifdef DEBUG
      static FILE *fplog;
      static int tot;
#endif
/*======================================================================
? main program, initialize few constants and go for it.
^=====================================================================*/
void mtc_init( FILE *input_file, FILE *output_file, FILE *error_file )
{
   VARIABLE *ptr;

   char str[256];

   int i;                    /* i'm getting tired with all these i's */

   static char *evalHelp =
   {
       "eval( str )\n\n"
       "Evaluate content variable. Another form of this command is @str.\n"
   };

   static char *sourceHelp =
   {
       "source( name )\n\n"
       "Execute commands from file given name.\n"
   };

   static char *helpHelp =
   {
       "help or help(\"symbol\")\n\n"
       "First form of the command gives list of available commands.\n"
       "Second form gives help on specific routine.\n"
   };

#ifdef DEBUG
      fplog = fopen("matcdbg","w");
#endif
  ALLOC_HEAD = (LIST *)NULL;

  /*
   *   input & output & error streams
   */
  math_in  = input_file;
  math_err = error_file;
  math_out = output_file;

  mtr_com_init();          /* initialize matrix handling commands   */
  var_com_init();          /*     ""     VARIABLE  ""      ""       */
  fnc_com_init();          /*     ""     function handling commands */
  fil_com_init();          /*     ""     file handling commands     */
  gra_com_init();          /*     ""     graphics commands          */
  str_com_init();          /*     ""     string handling            */

  /*
   *    and few others.
   */
  com_init( "eval"   , FALSE, FALSE, com_apply,   1, 1, evalHelp   );
  com_init( "source" , FALSE, FALSE, com_source,  1, 1, sourceHelp );
  com_init( "help"   , FALSE, FALSE, com_help   , 0, 1, helpHelp   );
  com_init( "quit"   , FALSE, FALSE, com_quit   , 0, 0, "quit\n" );
  com_init( "exit"   , FALSE, FALSE, com_quit   , 0, 0, "exit\n" );
  
  /*
   *    these constants will always be there for you.
   */
  ptr = const_new("true", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 1.0;

  ptr = const_new("false", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 0.0;

  ptr = const_new("stdin", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 0;

  ptr = const_new("stdout", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 1;

  ptr = const_new("stderr", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 2;

  ptr = const_new("pi", TYPE_DOUBLE, 1, 1);
  M(ptr,0,0) = 2*acos(0.0);

#if 0
  /*
   *  trap INTERRUPT and Floating Point Exeption signals
   */ 
   signal(SIGFPE, sig_trap);

   sprintf( str, "%s/lib/mc.ini", getenv("ELMER_POST_HOME") );

   if ( (math_in = fopen( str, "r" ) ) != (FILE *)NULL)
   {
       doread();
       fclose( math_in );
   }

  /*
   *   and finally standard input.
   */
  math_in = stdin;

  doread();

  var_free();
  com_free();
  fnc_free();
  const_free();

  mem_free_all();

#ifdef DEBUG
      fclose(fplog);
#endif
#endif

  return;   /* done */
}

char * mtc_domath( char *str )
{
  VARIABLE *headsave;            /* this should not be here */

  jmp_buf jmp, *savejmp;         /* save program context */

  void (*sigfunc)() =  (void (*)())signal( SIGINT, sig_trap );

  if ( !str || !*str )
  {
      str = (char *)doread();
      signal( SIGINT, sigfunc );
      return math_out_str;
  }

  savejmp = jmpbuf;
  jmpbuf = &jmp;

#ifdef DEBUG 
  fprintf( stderr, "got [%s]\n", str );
#endif 
  if ( math_out_str ) math_out_str[0] = '\0';
  math_out_count  = 0;

  /*
   *   try it
   */
  if (*str != '\0')
  {
     ALLOC_HEAD = (LIST *)NULL;
     headsave = (VARIABLE *)VAR_HEAD;

    /*
     *  normal return takes branch 1,
     *  error() takes branch 2, 
     *  quit() takes branch 3.
     */
    switch (setjmp(*jmpbuf))
    {
      case 0:
        (void)doit( str );
        longjmp(*jmpbuf, 1);
      break;

      case 1:
      break;

      case 2:
        VAR_HEAD = (LIST *)headsave;
      break;

      case 3:
      break;
    }
  }

  jmpbuf = savejmp;

  signal( SIGINT, sigfunc );

  return math_out_str;
}

char *doread()
/*======================================================================
?  doread() is really the main loop of this program. Function reads
|  it's input as strings and gives them to function doit() for
|  execution. setjmp() function is used for error recovery.
|
|  Memory allocated during the lifetime of this function is 
|  collected to a list represented by the global VARIABLE
|  ALLOCLIST *alloc_list. If function error() is called, this
|  list is used to deallocte memory. Normally (well I certainly 
|  hope so) functions which allocate memory deallocate it themselves.
|
|  Program stays in this function until an end of file -condition
|  is reached or exit- or quit-commands are given.
|
@  jmp_buf *jmpbuf, ALLOC_LIST *alloc_list
&  ALLOCMEM, FREEMEM, setjmp(), longjmp()
~  doit(), quit(), error()
^=====================================================================*/
{
  VARIABLE *headsave;            /* this should not be here */

  jmp_buf jmp, *savejmp;         /* save program context */

  char *p, *q;                   /* buffer for input stream */

  savejmp = jmpbuf;
  jmpbuf = &jmp;

  if ( math_out_str ) math_out_str[0] = '\0';
  math_out_count  = 0;

  p  = q = ALLOCMEM(4096);
  /*
   *   try it
   */
  while(dogets(p, PMODE_MAIN))
  {
    if (*p != '\0')
    {
      ALLOC_HEAD = (LIST *)NULL;
      headsave = (VARIABLE *)VAR_HEAD;

      /*
       *  normal return takes branch 1,
       *  error() takes branch 2, 
       *  quit() takes branch 3.
       */
      switch (setjmp(*jmpbuf))
      {
        case 0:
          (void)doit(p);
          longjmp(*jmpbuf, 1);
	break;

        case 1:
	break;

        case 2:
          VAR_HEAD = (LIST *)headsave;
        break;

        case 3:
         goto ret;
        break;
      }
    }
  }

 ret:

  jmpbuf = savejmp;

  FREEMEM(q);

  return math_out_str;
}

VARIABLE *com_quit()
/*======================================================================
?  Quit current doread entry by longjumping back to it (nasty).
&  longjmp
~  doread
^=====================================================================*/
{
  longjmp(*jmpbuf, 3);

  return (VARIABLE *)NULL;   /* won't be executed (hopefully) */
}

int dogets(buff, prompt) char *buff; char *prompt;
/*======================================================================
?  Get line from input stream. If both input & output streams are 
|  connected to terminal, this function gives user one of three 
|  (default) prompts:
|
|     MATC> 
|        - normal prompt                      (PMODE_MAIN)
|     ....> 
|        - begin end- block is beign defined  (PMODE_BLOCK)
|     ####>                                   (PMODE_CONT)
|        - user has given a #-sign as a last character of 
|          previous line, this line will be catenated with it
|
|  If current comment character is found from input stream, the
|  line after this character is discarded. Likewise if current
|  system command character is found, the rest of the line is
|  passed to system()-call.
|
=  line got -> TRUE, EOF -> FALSE
!  There should be a way to get an echo when reading from file.
&  fprintf(), isatty(), fileno(), strlen(), fgets(), system()
^=====================================================================*/
{
   char *ptr = buff, *p;    /* Can't get rid of these. */

   if ( !math_in ) return FALSE;

   /* 
       Try figuring out if input & output streams are
       terminals, if they both are, give user a prompt.
   */
   if (isatty(fileno(math_in)) && isatty(fileno(math_out)))
     PrintOut( "%s", prompt );

   /*
      i'm not in the mood to explain this.
   */
   *ptr++ = ' ';

   /* 
        Go for it.
   */
   while((ptr = fgets(ptr, 256, math_in)) != NULL)
   {

     ptr[strlen(ptr)-1] = '\0';

     /*
      *    Check if the user wants to continue with this line.
      */
     while(ptr[strlen(ptr)-1] == '\\')
     {
       ptr += strlen(ptr) - 1;
       dogets(ptr, PMODE_CONT);
     }

     /* 
      *   if there is only spaces in this line, 
      *   don't bother returning it, instead
      *   let's read afresh, otherwise return.
      */
     p = ptr; while(isspace(*p)) p++;

     if (*p != '\0')   /* GOOD EXIT HERE */
     {
#if 0
       /* 
        *   Look for the system character, if found
        *   pass rest of the line to system()-call
        */
       for(p = buff; *p; p++)
       {
         switch(*p)
         {
           case SYSTEM:
             system(p + 1);
             PrintOut("\n");
             *p = '\0'; p--;
           break;
         }
       }
#endif
       if (*buff != '\0')   
         return TRUE;          /* OR IF WE ARE HONEST, IT'S HERE */
     }

     /* 
          if it's terminal give a prompt.
     */
     if (isatty(fileno(math_in)) && isatty(fileno(math_out)))
         PrintOut("%s", prompt);
   }

   return FALSE;
}


void com_init(word, flag_pw, flag_ce, sub, minp, maxp, help_text ) 
/*======================================================================
?  Adds commands to global command list. 
|
|  Parameters:
|      char *word 
|         - the keyword user gives for this command to be executed. 
|      int flag_pw 
|         - flag telling if the command can be executed element
|           by element using function *(*sub)().
|      int flag_ce
|         - flag telling if the command can be executed when
|           preprosessing if constant arguments
|      double *(*sub)()
|         - function to be executed when this command is given
|      int minp, maxp
|         - maximum and minimum number of parameters to command
|
|  The global list of available commands is updated (or created if
|  nonexistent).
|
&  lst_add()
~  *_com_init()
^=====================================================================*/
     char *word; 
     VARIABLE *(*sub)();
     int minp, maxp, flag_pw, flag_ce;
     char *help_text;
{
  COMMAND *ptr;        /* can't get rid of this */
  

  /* 
     Fill the structure... 
  */
  ptr = (COMMAND *)ALLOCMEM(COMSIZE);
  NAME(ptr) = STRCOPY(word);
  if (flag_pw) 
    ptr->flags |= CMDFLAG_PW;
  if (flag_ce) 
    ptr->flags |= CMDFLAG_CE;
  ptr->minp = minp;
  ptr->maxp = maxp;
  ptr->sub  = sub;
  ptr->help = help_text;

  /* 
     ...and update the list. 
  */
  lst_add(COMMANDS, (LIST *)ptr);

  return;
}

void com_free()  
/*======================================================================
?  Deletes the list of commands and frees associated memory.
|
&  lst_purge()
^=====================================================================*/
{
  /* 
       Give memory back to system 
  */
  lst_purge(COMMANDS);

  return;
}

COMMAND *com_check(str) char *str;
/*======================================================================
?  Look for command from COMMANDS list by name.
|
=  COMMAND *NULL if does not exist, pointer to command otherwise
&  lst_find()
^=====================================================================*/
{
  return (COMMAND *)lst_find(COMMANDS, str);
}

VARIABLE *com_help( VARIABLE *ptr )
/*======================================================================
?  Print list of commands and user defined functions from global lists.
|
!  The command to get here is "help" but it really is not very helpful.
|
&  lst_print()
^=====================================================================*/
{
   COMMAND  *cmd;
   FUNCTION *fnc;
   char *name;

   if ( !ptr )
   {

       lst_print(COMMANDS);
       lst_print(FUNCTIONS);

   } else {

      name = var_to_string( ptr );

      if ( (cmd = com_check( name ) ) !=  (COMMAND *)NULL )
      {

          if ( cmd->help )
              PrintOut( "\n%s\n", cmd->help );
          else
              PrintOut( "\nSorry: no help available on [%s].\n", name );

      } else if ( (fnc = fnc_check( name ) ) !=  (FUNCTION *)NULL )
      {

          if ( fnc->help )
              PrintOut( "\n%s", fnc->help );
          else
              PrintOut( "\nSorry: no help available on [%s].\n", name );

      } else {

          error( "help: symbol not found: [%s]\n", name );

      }

      FREEMEM( name );
   }

  return (VARIABLE *)NULL;
}

VARIABLE *com_pointw(sub, ptr)  double (*sub)(); VARIABLE *ptr;
/*======================================================================
?  This routine does a function call (*sub)(), for each element in
|  matrix given by ptr.  
|
=  a temporay VARIABLE for which M(res, i, j) = (*sub)(M(ptr, i, j)
&  var_temp_new(), *(sub)()
^=====================================================================*/
{
  VARIABLE *res;     /*  pointer to result structure */

  double *a, *b;     /* pointer to matrices */
  int n, m;          /* matrix dimensions   */

  int i;             /* loop index          */

  /* 
     Get space for result and ...
  */
  n = NROW(ptr); m = NCOL(ptr);
  res = var_temp_new(TYPE(ptr) ,n , m);

  n *= m;
  a = MATR(ptr); b = MATR(res);
  /* 
      ...to action.
  */
  for(i = 0; i < n; i++) *b++ = (*sub)(*a++);

  return res;
}

VARIABLE *com_el(ptr) VARIABLE *ptr;
/*======================================================================
?  Extracts specified elements from a matrix. Indexes are given by two
|  column vectors. The values of the elements of these vectors give
|  the required indexes. If there is only one index vector given
|  it is assumed to be column index and row index is set to scalar 0.
|
|  If matrix x is, for example, 
|
|        1 2
|        3 4 
|
|  you get the first row by 
|
|        x[0, 0 1]
|
|  or by 
| 
|        x(0 1)
|
=  A new temporary VARIABLE, whose size equals to
|  number of row indexes times number of column indexes.
|  
&  var_temp_new(), var_delete_temp()
^=====================================================================*/
{
  VARIABLE *res,               /* result ... */
           *par = NEXT(ptr);   /* pointer to list of VARIABLES */
                               /* containig indexes            */

  static double defind = 0.0;

  double *ind1 = &defind, *ind2;

  int      i, j, k,      /* loop indexes */
           rows, cols,   /* no. of rows and columns in the matrix */
                         /* to be indexed.                        */
           size1 = 1, size2,
           ind;

  rows = NROW(ptr); cols = NCOL(ptr);

  /*
   * check if scalar ....
   */
  if (rows == 1 && cols == 1)
  {
    if (*MATR(par) != 0) error("Index out of bounds.\n");
    if (NEXT(par) != NULL)
      if (*MATR(NEXT(par)) != 0) error("Index out of bounds.\n"); 
    res = var_temp_new(TYPE(ptr),1,1);
    *MATR(res) = *MATR(ptr);
    return res;
  }

  /* 
     The matrix will be indexed by two column vectors. 
     If there is just one assume it's column index and
     make rowindex 0.
  */
  if (NEXT(par) == NULL)
  {
    if (NROW(par) == rows && NCOL(par) == cols)
    {
      int logical = TRUE,
          onecount=0;

      double *dtmp;

      dtmp = MATR(par);
      for(i = 0; i < NROW(par)*NCOL(par); i++)
        if (dtmp[i] == 0)
        {
        }
        else if (dtmp[i] == 1)
        {
          onecount++;
        }
        else
        {
          logical = FALSE;
          break;
        }

      if (logical)
      {
        if (onecount == 0) return NULL;

        res = var_temp_new(TYPE(ptr),1,onecount); 
        for(i=0,k=0; i < rows; i++)
          for(j=0; j < cols; j++)
            if (M(par,i,j) == 1)
            {
              memcpy(&M(res,0,k++),&M(ptr,i,j),sizeof(double));
            }
        return res;
      }
    }

    ind2 = MATR(par); size2 = NCOL(par); 
    cols *= rows; rows = 1;
  }
  else
  {
    ind1 = MATR(par); size1 = NCOL(par);
    size2 = NCOL(NEXT(par));
    ind2  = MATR(NEXT(par));
  }

  /* 
      Space for result  
  */
  res = var_temp_new(TYPE(ptr), size1, size2);

  /* 
     Extract the values (try making sense out of that
     if you feel like it).
  */
  for(i = 0; i < size1; i++)
  {
    ind = (int)ind1[i];
    for(j = 0; j < size2; j++)
      if (ind < rows && (int)ind2[j] < cols)
        memcpy(&M(res,i,j),&M(ptr,ind,(int)ind2[j]),sizeof(double));
      else
        error("Index out of bounds.\n");
  }

  return res;
}

VARIABLE *com_source(ptr) VARIABLE *ptr;
/*======================================================================
?  Redirect input stream to a file, whose name is given.
|
@  FILE *math_in
&  ALLOCMEM, FREEMEM, fopen(), fclose(), error()
^=====================================================================*/
{
  char *name;               /* Hold converted string (file name) */

  FILE *save_in = math_in;  /* Save previous input stream until  */
                            /* we are done with the new file.    */

  /*
      convert the file name from ptr.
  */
  name = var_to_string(ptr);

  /*
      Execute the file.
  */
  if ((math_in = fopen(name,"r")) != NULL)
  {
/*   PrintOut("Executing commands from file, %s...\n", name); */
     doread();
     fclose(math_in);
  }
  else
  {
    PrintOut( "Source: Can't open file, %s.\n",name );
  }

  math_in = save_in;
  FREEMEM(name);

  return (VARIABLE *)NULL;
}


VARIABLE *com_apply(ptr) VARIABLE *ptr;
/*======================================================================
?  Executes given string.
|
&  ALLOCMEM, FREEMEM, doit()
^=====================================================================*/
{
  VARIABLE *res;        /* result pointer */

  char *p, *q;          /* holds the string to be executed, after */
                        /* conversion from structure VARIABLE *   */

  int i, j;             /* just loop indexes */


  /* 
      Allocate space for the string... 
  */
  p = q = ALLOCMEM(NROW(ptr) * NCOL(ptr) + 1);

  /*
      ... convert it ... 
  */
  for(i = 0; i < NROW(ptr); i++)
    for(j = 0; j < NCOL(ptr); j++) 
      *p++ = (char)M(ptr,i,j);

  *p = '\0';

  /*
      ... and try executing it. 
  */
  res = doit( q ); 
  
  FREEMEM(q);

  return res;
}

void mem_free(void *mem)
/*======================================================================
?  Free memory given by argument, and unlink it from alloction list.
|  Currently FREEMEM(ptr) is defined to be mem_free(ptr).
|
&  free()
~  mem_alloc(), mem_free_all()
^=====================================================================*/
{
  ALLOC_LIST *lst;

#ifdef DEBUG
  tot--; fprintf(fplog,"free addr: %d total: %d\n", ALLOC_LST(mem), tot);
  fflush( fplog );
#endif
  /*
      if the list is empty return
  */
  if ( (lst = (ALLOC_LIST *)ALLOC_HEAD) == (ALLOC_LIST *)NULL )
  {
#if 1
/* ????? */
     free( ALLOC_LST(mem) );
#else
fprintf( stderr, "SHOULD THIS HAPPEN ????\n" );
#endif
    return;
  }
 
  /*
   *  it's not the header, look if it's in list at all  
   */
  if (ALLOC_PTR(lst) != mem)
  {

    for(; NEXT(lst); lst = NEXT(lst))
    {
      if (ALLOC_PTR(NEXT(lst)) == mem) break;
    }

    /*
     *  item was not found from the list. free ptr and return.
     */
    if (NEXT(lst) == (ALLOC_LIST *)NULL)
    {
      free(ALLOC_LST(mem));
      return;
    }

    /*
     *   unlink
     */
    NEXT(lst) = NEXT(NEXT(lst));
  }

  /*
   *  item was the header, unlink it    
   */
  else
    ALLOC_HEAD = NEXT(ALLOC_HEAD);

  /*
   *  and at last return memory back to system
   */
  free(ALLOC_LST(mem));

  return;
}

void mem_free_all()
/*======================================================================
?  Free all memory allocated since last entry of doread.
|  (actually free all memory from list ALLOCATIONS).
|
~  mem_alloc(), mem_free(), doread(), error()
^=====================================================================*/
{
  ALLOC_LIST *lst, *lstn; 

  for(lst = (ALLOC_LIST *)ALLOC_HEAD; lst;)
  {
#ifdef DEBUG
  tot--; fprintf(fplog,"freeall addr: %d total: %d\n", lst, tot);
  fflush( fplog );
#endif
    lstn = NEXT(lst);
    free( (char *)lst );
    lst = lstn;
  }

  ALLOC_HEAD = (LIST *)NULL;  /* security */

  return;
}

void *mem_alloc(size) size_t size;
/*======================================================================
?  Allocate memory and link it to  memory allocation list.
|
~   calloc(), free(), error()
^=====================================================================*/
{
  ALLOC_LIST *lst;

  /* 
   *    try allocating memory 
   */
  if ((lst = (ALLOC_LIST *)calloc(size+sizeof(ALLOC_LIST), 1)) != NULL)
  {
    NEXT(lst) = (ALLOC_LIST *)ALLOC_HEAD; ALLOC_HEAD = (LIST *)lst;
  }
  else
    error("Can't alloc mem.\n");

#ifdef DEBUG
  tot++; fprintf(fplog,"alloc addr: %d size: %d total: %d\n", 
                          lst, size, tot);
  fflush( fplog );
#endif
  return ALLOC_PTR(lst);
}
