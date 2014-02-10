#ifndef MC_H
#define MC_H

extern "C" {
typedef struct list {
  struct list *next;              /* pointer to next item in list */
  char *name;                     /* name of list item            */
} LIST;

extern LIST listheaders[];

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
#define ALLOC_SIZE(size) (size+sizeof(ALLOC_LIST))
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
    MATRIX *me;
} VARIABLE;

/*
     shortcuts for accsessing structure MATRIX
*/
#define MATR(ptr)    (ptr)->me->data
#define TYPE(ptr)    (ptr)->me->type
#define NROW(ptr)    (ptr)->me->nrow
#define NCOL(ptr)    (ptr)->me->ncol
#define REFCNT(ptr)  (ptr)->me->refcount
#define M(ptr,i,j)   (ptr)->me->data[(i) * NCOL(ptr) + (j)]

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

}

#endif // MC_H
