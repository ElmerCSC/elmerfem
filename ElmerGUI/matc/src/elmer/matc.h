/***********************************************************************
|
|  MATC.H - Last Edited 7. 8. 1988
|
************************************************************************/


/*
 * $Id: matc.h,v 1.2 2007/06/08 08:12:19 jpr Exp $ 
 *
 * $Log: matc.h,v $
 * Revision 1.2  2007/06/08 08:12:19  jpr
 * *** empty log message ***
 *
 * Revision 1.1  2005/05/27 12:26:22  vierinen
 * changed header install location
 *
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.3  2001/06/08 09:20:29  jpr
 * *** empty log message ***
 *
 * Revision 1.2  1998/08/01 12:34:49  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <string.h>
#include <signal.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>


#ifdef MODULE_MATC
#define EXT
#else
#define EXT extern
#endif

/*******************************************************************
                 LIST HANDLING DEFINITIONS
*******************************************************************/

typedef struct list {
  struct list *next;              /* pointer to next item in list */
  char *name;                     /* name of list item            */
} LIST;

/*
    pointers to start of global lists
*/
#ifdef MODULE_MATC

  EXT LIST listheaders[5] = {
    {
      NULL, "Allocations"                    /* memory allocations     */
    }, {
      NULL, "Constants"                     /* global CONSTANTS        */
    }, {
      NULL, "Currently defined VARIABLES"    /* global VARIABLES       */
    }, { 
      NULL, "Builtin Functions"              /* internal commands      */
    }, {
      NULL, "User Functions"                 /* user defined functions */
    }
  };
#else

  EXT LIST listheaders[];

#endif

#define ALLOCATIONS 0
#define CONSTANTS   1
#define VARIABLES   2
#define COMMANDS    3
#define FUNCTIONS   4

#define MAX_HEADERS 4

#define ALLOC_HEAD listheaders[ALLOCATIONS].next
#define CONST_HEAD listheaders[CONSTANTS].next
#define VAR_HEAD   listheaders[VARIABLES].next
#define COM_HEAD   listheaders[COMMANDS].next
#define FUNC_HEAD  listheaders[FUNCTIONS].next

#define NEXT(lst) (lst)->next
#define NAME(lst) (lst)->name

/*******************************************************************
                      MEMORY HANDLING
********************************************************************/

/*
    memory allocation and deallocation routines
*/
#define ALLOCMEM(size) mem_alloc(size)
#define FREEMEM(ptr) mem_free(ptr)

/*
    we use a lot of string copying.
*/
#define STRCOPY(str) strcpy((char *)ALLOCMEM(strlen(str)+1),(str))

typedef struct alloc_list {
    struct alloc_list *next;
    char *mem;
} ALLOC_LIST;

#define ALLOC_LST(mem) (ALLOC_LIST *)((char *)mem-sizeof(ALLOC_LIST))
#define ALLOC_PTR(lst) (char *)((char *)lst+sizeof(ALLOC_LIST))

/*******************************************************************
                           VARIABLES
*******************************************************************/

/*
 *    MATC matrix is internally represented by this structure.
 */
typedef struct MATRIX
{
   int type,                    /* TYPE_DOUBLE or TYPE_STRING       */
       refcount,                /* reference count                  */
       nrow, ncol;              /* number of rows and columns       */
   double *data;                /* pointer to double array          */
} MATRIX;


/*
 *   list of VARIABLES
 */

typedef struct variable
{
    struct variable *next;       /* pointer to next item in list     */
    char *name;                  /* name of the item                 */
    int changed;
    MATRIX *this;
} VARIABLE;

/*
     shortcuts for accsessing structure MATRIX
*/
#define MATR(ptr)    (ptr)->this->data
#define TYPE(ptr)    (ptr)->this->type
#define NROW(ptr)    (ptr)->this->nrow
#define NCOL(ptr)    (ptr)->this->ncol
#define REFCNT(ptr)  (ptr)->this->refcount
#define M(ptr,i,j)   (ptr)->this->data[(i) * NCOL(ptr) + (j)]

#define VARIABLESIZE sizeof(VARIABLE)
#define MATRIXSIZE   sizeof(MATRIX)
#define MATSIZE(ptr) NROW(ptr)*NCOL(ptr)*sizeof(double)

#define TYPE_DOUBLE  0
#define TYPE_COMPLEX 1       /* this is not */
#define TYPE_STRING  2

/*******************************************************************
               INTERNAL COMMANDS AND USER FUNCTIONS
*******************************************************************/

typedef struct command
{
  struct command *next;        /* pointer to next item in list    */
  char *name;                  /* name of the item                */
  int flags,                   /* CMDFLAG_PW & CMDFLAG_CE         */
      minp, maxp;              /* min. and max. no. of parameters */
   VARIABLE *(*sub)();         /* function to execute             */
  char *help;                  /* help string... */
} COMMAND;

#define COMSIZE sizeof(COMMAND)

#define CMDFLAG_PW 1           /* element by element operation    */
#define CMDFLAG_CE 2           /* command can be executed when 
                                  preprosessing if constant
                                  arguments.                      */

/*******************************************************************
               USER DEFINED FUNCTIONS
*******************************************************************/

typedef struct function
{
  struct function *next;     /* pointer to next function in list  */
  char *name,                /* name of the function              */
       **parnames,           /* function parameter names (if any) */
       **exports,            /* functions exported variables      */
       **imports,            /* functions imported variables      */
       *help;                /* functions help text               */
  int parcount;              /* defined number of parameters      */
  struct clause *body;       /* function body                     */
} FUNCTION;

#define FUNCSIZE sizeof(FUNCTION)

/*******************************************************************
               MISC DEFINITONS FOR PARSER
*******************************************************************/

typedef enum symbols {
  nullsym,  leftpar,  rightpar, indopen, indclose, power, times, ptimes, divide,
  plus, minus, reduction,  transpose, eq, neq, lt, gt, le, ge, and, or, not,
  assignsym, apply, resize, vector, statemend, argsep, name, number, string, 
  funcsym, import, export, ifsym, thensym, elsesym, whilesym, forsym,
  beginsym, endsym, breaksym, comment, systemcall
} SYMTYPE;

#ifdef MODULE_MATC

/*--------------------------------------------------------------------*/

SYMTYPE ssymbols[] = {
  leftpar, rightpar, indopen, indclose, beginsym, endsym, power, times, ptimes,
  divide, plus,  minus,  reduction,   transpose, lt,  gt, and, or, not,
  assignsym, apply, resize, vector, statemend, argsep, comment, systemcall
};

char csymbols[] = {
  '(', ')', '[', ']', '{', '}', '^', '*', '#', '/', '+', '-', '\?',
  '\'', '<', '>', '&', '|', '~', '=', '@', '%', ':', ';', ',', '!', '$'
};

/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/

char *reswords[] = {
  "function", "import", "export", "if", "then", "else", "while", "for", 
  "begin", "end", "break", NULL
};

SYMTYPE rsymbols[] = {
  funcsym, import, export, ifsym, thensym, elsesym, whilesym, forsym,
  beginsym, endsym, breaksym
};

/*--------------------------------------------------------------------*/

char *symchars = "`._";

#else

   EXT SYMTYPE ssymbols[];
   EXT char csymbols[];
   EXT char *reswords[];
   EXT SYMTYPE rsymbols[];
   EXT char *symchars;

#endif

/*
    dataentry for expression trees
*/
typedef struct treeentry
{
  struct tree *args,       /* parameters for functions                   */
              *subs;       /* indexes for VARIABLES, equation results    */
  int entrytype;           /* type of entrydata                          */

  union data_entry
  {
     char *s_data;          /* function or VARIABLE names or string constants */
     double d_data;         /* numeric constant                               */
     VARIABLE *c_data;      /* real constant, with no name references         */
     MATRIX *(*v_data)();   /* function address (for builtin operations)      */
  } entrydata;

} TREEENTRY;

#define ETYPE_NAME   0
#define ETYPE_NUMBER 1
#define ETYPE_STRING 2
#define ETYPE_OPER   3
#define ETYPE_CONST  4
#define ETYPE_EQUAT  5

/*
 *   four leaf tree, isn't that odd
 */
typedef struct tree {
  struct tree *next;
  struct tree *link;
  struct tree *left, *right;
  TREEENTRY tentry;
} TREE;

/*
    shortcuts for accsessing above structures
*/
#define SDATA(ptr) (ptr)->tentry.entrydata.s_data
#define DDATA(ptr) (ptr)->tentry.entrydata.d_data
#define CDATA(ptr) (ptr)->tentry.entrydata.c_data
#define VDATA(ptr) (ptr)->tentry.entrydata.v_data
#define ETYPE(ptr) (ptr)->tentry.entrytype
#define SUBS(ptr)  (ptr)->tentry.subs
#define ARGS(ptr)  (ptr)->tentry.args
#define LEFT(ptr)  (ptr)->left
#define RIGHT(ptr) (ptr)->right

/*
    this is an operations list. data can be one of
    the following:

        ifsym, elsesym, whilesym, forsym, assignsym, funcsym

     every input line is compiled to this type of list,
     and it is used to hold function bodies.
*/

typedef struct clause
{
   struct clause *link;
   struct clause *jmp;
   TREE *this;
   SYMTYPE data;
} CLAUSE;

#define LINK(ptr)  (ptr)->link

/*******************************************************************
                           THIS AND THAT
*******************************************************************/

#ifdef sign
#  undef sign
#endif
#ifdef max
#  undef max
#endif
#ifdef min
#  undef min
#endif
#ifdef abs
#  undef abs
#endif

#define sign(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) > (y) ? (y) : (x))
#define abs(x) ((x) > 0 ? (x) : -(x))

#define FOREVER for(;;)
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
    promptmode flags:

       PMODE_MAIN  ===>   MATC>
       PMODE_BLOCK ===>   ....>
*/
#define PMODE_MAIN  "MATC> "
#define PMODE_BLOCK "....> "
#define PMODE_CONT  "####> "

EXT FILE *math_in, *math_out, *math_err;

/* 
     see doread(), error() in matc.c
*/
EXT jmp_buf *jmpbuf;

EXT int term;

#ifdef VAX
struct desc
{ 
  int length; 
  char *addr;
} ;
#endif

#define COMMENT '!'               /* comment introducer        */
#define SYSTEM '$'                /* system()-call introducer  */

#define STRING_OUTPUT

#ifdef STRING_OUTPUT
#ifdef MODULE_MATC
static char *math_out_str = NULL;
static int math_out_count;
#endif
#endif

void mem_free_all(void);


#ifdef MODULE_MATC
static int math_out_allocated  = 0;
void error( char *format, ... )
{
    va_list args;

    va_start( args, format );
#ifdef STRING_OUTPUT
    if ( math_out_count+512 > math_out_allocated )
    { 
        math_out_allocated += 512;
        math_out_str = (char *)realloc( math_out_str, math_out_allocated );
    }
    math_out_count += sprintf( &math_out_str[math_out_count], "MATC ERROR: " );
    math_out_count += vsprintf( &math_out_str[math_out_count], format, args );
#else
    fprintf( math_err, "MATC ERROR: " );
    vfprintf( math_err, format, args );
#endif
    va_end( args );

    (void)mem_free_all();
    longjmp( *jmpbuf, 2 );
}

void PrintOut( char *format, ... )
{

    va_list args;

    va_start( args, format );
#ifdef STRING_OUTPUT
    if ( math_out_count+512 > math_out_allocated )
    { 
        math_out_allocated += 512;
        math_out_str = (char *)realloc( math_out_str, math_out_allocated );
    }
    math_out_count += vsprintf( &math_out_str[math_out_count], format, args );
#else
    vfprintf( math_out, format, args );
#endif
    va_end( args );
}
#else
extern void error( char *format, ... );
extern void PrintOut( char *format, ... );
#endif

/*******************************************************************
                       function prototypes 
*******************************************************************/ 
#include "fnames.h"

/*******************************************************************
                  graphics package defitions
*******************************************************************/ 
#include "gra.h"
