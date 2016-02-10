/*****************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! ******************************************************************************
! *
! *  Utilities for dynamic loading of user functions, and other operating
! *  system interfaces.
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
/* #include <elmer/matc.h> maybe in the future */

/* eg. FC_CHAR_PTR and FC_FUNC is defined here */

#include "../config.h"

#if defined(WIN32) | defined(MINGW32)
#  include <direct.h>
#  include <windows.h>
#define ELMER_PATH_SEPARATOR ";"
#else
#include <strings.h>
#  include <dlfcn.h>
#  include <sys/stat.h>
#define ELMER_PATH_SEPARATOR ":"
#endif

#define MAX_PATH_LEN 512
#define ERROR_BUF_LEN 10*MAX_PATH_LEN

#ifndef USE_ISO_C_BINDINGS
#ifdef SGI64
void corename_()
{
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/prctl.h>

 prctl( PR_COREPID,0,0 );
}
#endif
#endif

/* pc needs more bits on 64bit arch  */
#ifdef ARCH_32_BITS
#define f_ptr int32_t *
#else 
#define f_ptr int64_t *
#endif

/*#if defined(MINGW32)*/
/*--------------------------------------------------------------------------
  work around mingw rxvt shell stdio/err buffering troubles
  -------------------------------------------------------------------------*/
/*void STDCALLBULL FC_FUNC_(set_stdio_bufs,SET_STDIO_BUFS) ()*/
/*[>void set_stdio_bufs_()<]*/
/*{*/
   /*setvbuf( stdout, NULL, _IOLBF, 2048 );*/
   /*setvbuf( stderr, NULL, _IONBF, 2048 );*/
/*}*/
/*#endif*/

/*--------------------------------------------------------------------------
  This routine will return the home directory of elmer solver.
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL getsolverhome( char *solverDir, int *len)
#else
void STDCALLBULL FC_FUNC(getsolverhome,GETSOLVERHOME) 
     ( char *solverDir, int *len)
#endif
{
  *len = 0;

  char *elmer_home = getenv("ELMER_HOME");

  if(elmer_home != NULL) {
    /* Return solver home relative to ELMER_HOME*/
#if defined(WIN32) || defined(MINGW32)
    _snprintf(solverDir, MAX_PATH_LEN, "%s\\share\\elmersolver", elmer_home);
#else
    snprintf(solverDir, MAX_PATH_LEN, "%s/share/elmersolver", elmer_home);
#endif
    *len = strlen(elmer_home) + 18;
    if(*len > MAX_PATH_LEN) *len = MAX_PATH_LEN;
    return;
  }

#if defined(WIN32) || defined(MINGW32)
  static char appPath[MAX_PATH_LEN] = "";
  static char appDir[MAX_PATH_LEN] = "";
  char *exeName = NULL;
  int n = 0;

  /* Get the full module file name  */
  GetModuleFileName(NULL, appPath, MAX_PATH_LEN);
  if(appPath == NULL) return;
  exeName = strrchr(appPath, '\\');
  if(exeName == NULL) return;
  n = (int)(exeName - appPath);
  if(n < 0) return; /* play safe */
  if(n > MAX_PATH_LEN) n = MAX_PATH_LEN;

  /* This is where the executable resides */
  strncpy(appDir, appPath, n);

  /* Return solver home relative to appDir */
  _snprintf(solverDir, MAX_PATH_LEN, "%s\\..\\share\\elmersolver", appDir);
  *len = n + 21;
  if(*len > MAX_PATH_LEN) *len = MAX_PATH_LEN;
#else

  /* Use the directory defined in config.h */
  snprintf(solverDir, MAX_PATH_LEN, "%s", ELMER_SOLVER_HOME);
  *len = strlen(ELMER_SOLVER_HOME);
  if(*len > MAX_PATH_LEN) *len = MAX_PATH_LEN;
#endif
}

/*--------------------------------------------------------------------------
  This routine will create a directory given name of the directory.
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL makedirectory(char *Name)
#else
void STDCALLBULL FC_FUNC(makedirectory,MAKEDIRECTORY)
     (char *Name)
#endif
{
#if defined(WIN32) || defined(MINGW32)
    if ( _mkdir( Name ) != 0 ) {
#else
    if ( mkdir( Name, 0700 ) != 0 ) {
      chmod( Name, 0700 );
#endif
    }
}

#ifndef USE_ISO_C_BINDINGS
/*--------------------------------------------------------------------------
  This routine execute a operating system command.
  -------------------------------------------------------------------------*/
void STDCALLBULL FC_FUNC(systemc,SYSTEMC) ( char *str )
{
   system( str );
}

/*--------------------------------------------------------------------------
  This routine will return value of a environment variable to a
  given string variable.
  -------------------------------------------------------------------------*/
void STDCALLBULL FC_FUNC(envir,ENVIR) (char *Name, char *Value, int *len)
{
    if ( getenv( Name ) ) {
      strncpy( Value,(char *)getenv(Name), MAX_PATH_LEN );
      *len = strlen( Value );
    } else {
      *len = 0;
      *Value = '\0';
    }
}
#endif

/*--------------------------------------------------------------------------
  Internal: convert function names into to fortran mangled form for dynamical
  loading
  ---------------------------------------------------------------------------*/
static void STDCALLBULL fortranMangle(char *orig, char *mangled)
{
  int uscore, i;
  
  strcpy( mangled, orig );

  if(ELMER_LINKTYP == 1 || ELMER_LINKTYP == 3 || ELMER_LINKTYP == 4)
  {
    for( i=0 ; i<strlen(mangled) ; i++ ) /* to lower case */
    {
      if ( mangled[i] >= 'A'  && mangled[i] <= 'Z' ) 
	mangled[i] += 'a' - 'A';
    }
  }
  if(ELMER_LINKTYP == 2)
  {
    for( i=0; i<strlen(mangled); i++ ) /* to upper case */
    {
      if ( mangled[i] >= 'a'  && mangled[i] <= 'z' ) 
	mangled[i] += 'A' - 'a';
    }
  }
  
  if(ELMER_LINKTYP == 1) /* underscore */
  {
      strcat( mangled, "_" );
  }
  else if(ELMER_LINKTYP == 4) /* 1-2 underscores  */
  {
    uscore = 0;
    for( i=0; i<strlen(mangled); i++ )
      if(mangled[i] == '_')
	uscore++;
    
    if(uscore == 0)
    {
      strcat( mangled, "_" );
    } 
    else 
    {
      strcat( mangled, "__" );
    }
  }

}

/*--------------------------------------------------------------------------
  INTERNAL: Appends two paths with slash checking
  Args: path1, path2 - string to join
 -------------------------------------------------------------------------*/
static void STDCALLBULL append_path(char *path1, char *path2)
{
    size_t len1;

    len1 = strnlen(path1, 2*MAX_PATH_LEN);
#if defined(WIN32) || defined(MINGW)
    if (path1[len1-1] != '\\') {
        strncat(path1, "\\", 2*MAX_PATH_LEN);
    }
#else
    if (path1[len1-1] != '/') {
        strncat(path1, "/", 2*MAX_PATH_LEN);
    }
#endif
    strncat(path1, path2, 2*MAX_PATH_LEN);
}

/*--------------------------------------------------------------------------
  INTERNAL: Tries to open library with dlopen, first without
            any extensions and then with SHL_EXTENSION.
  Args: Libname - name of the library file
        Handle - handle to the dl, NULL if fails
        error_buf - string buffer for error messages
 -------------------------------------------------------------------------*/
static void STDCALLBULL try_dlopen(char *LibName, void **Handle, char *errorBuf)
{
    static char dl_names[2][2*MAX_PATH_LEN];
    char error_tmp[MAX_PATH_LEN];
    int i;

    strncpy(dl_names[0], LibName, 2*MAX_PATH_LEN);
    strncpy(dl_names[1], LibName, 2*MAX_PATH_LEN);

    strncat(dl_names[1], SHL_EXTENSION, MAX_PATH_LEN);

    for (i = 0; i < 2; i++) {
#ifdef HAVE_DLOPEN_API
        if ((*Handle = dlopen(dl_names[i], RTLD_NOW)) == NULL) {
            strncat(errorBuf, dlerror(), MAX_PATH_LEN);
            strncat(errorBuf, "\n", MAX_PATH_LEN);
        } else {
            break;
        }
#elif defined(HAVE_LOADLIBRARY_API)
        if ((*Handle = LoadLibrary(dl_names[i])) == NULL) {
            sprintf(error_tmp, "Can not find %s.\n", dl_names[i]);
            strncat(errorBuf, error_tmp, ERROR_BUF_LEN);
        } else {
            break;
        }
#endif
    }
}

/*--------------------------------------------------------------------------
  INTERNAL: Parses the search path and tries to open a solver.
            First search is done without any path prefixes.
  Args: SearchPath - colon separated list of searhc paths
        Library - name of the library file to be opened
        Handle - handle to the dl file, NULL if fails
        error_buf - string buffer for error messages
 --------------------------------------------------------------------------*/
static void STDCALLBULL
try_open_solver(char *SearchPath, char *Library, void **Handle, char *errorBuf)
{
    static char CurrentLib[2*MAX_PATH_LEN];
    char *tok;

    /* Try to open first without any prefixes */
    try_dlopen(Library, Handle, errorBuf);

    /* and then using the provided paths */
    if (*Handle == NULL) {

        tok = strtok(SearchPath, ELMER_PATH_SEPARATOR);
        while (tok != NULL) {
            strncpy(CurrentLib, tok, 2*MAX_PATH_LEN);
            append_path(CurrentLib, Library);

            try_dlopen(CurrentLib, Handle, errorBuf);
            if (*Handle != NULL)
                break;
            tok = strtok(NULL, ELMER_PATH_SEPARATOR);
        }
    }
}

/*--------------------------------------------------------------------------
  This routine will return address of a function given path to a dynamically
  loaded library and name of the routine.
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void *STDCALLBULL loadfunction_c( int *Quiet, int *abort_not_found, char *Library, char *Name )
#else
void *STDCALLBULL FC_FUNC(loadfunction,LOADFUNCTION) ( int *Quiet, int *abort_not_found,
                                                       char *Library, char *Name )
#endif
{
/*--------------------------------------------------------------------------*/
   void (*Function)(),*Handle;
   char *cptr;
   static char ElmerLib[2*MAX_PATH_LEN], NewLibName[3*MAX_PATH_LEN],
               NewName[MAX_PATH_LEN], ErrorBuffer[ERROR_BUF_LEN];
/*--------------------------------------------------------------------------*/
   static char appPath[MAX_PATH_LEN] = "";
   char *exeName = NULL;
   int n = 0;
/*--------------------------------------------------------------------------*/
   memset(appPath, 0, MAX_PATH_LEN);
   memset(ElmerLib, 0, 2*MAX_PATH_LEN);
   memset(NewLibName, 0, 3*MAX_PATH_LEN);
   memset(NewName, 0, MAX_PATH_LEN);
   memset(ErrorBuffer, 0, ERROR_BUF_LEN);
/*--------------------------------------------------------------------------*/
   fortranMangle( Name, NewName );
   strncpy( NewLibName, Library, 3*MAX_PATH_LEN );

   if ( *Quiet==0 ) {
     fprintf(stdout,"Loading user function library: [%s]...[%s]\n", Library, Name );
     fflush(stdout);
   }

   /* First path is always current directory (.) */
   strncpy(ElmerLib, ".", 2*MAX_PATH_LEN);
   cptr = (char *)getenv( "ELMER_LIB" );
   if ( cptr != NULL ) {
      strncat( ElmerLib, ELMER_PATH_SEPARATOR, 2*MAX_PATH_LEN );
      strncat( ElmerLib, cptr, 2*MAX_PATH_LEN );
   } else {
      cptr = (char *)getenv("ELMER_HOME");
      if ( cptr != NULL  ) {
         strncat( ElmerLib, ELMER_PATH_SEPARATOR, 2*MAX_PATH_LEN);
         strncat( ElmerLib, cptr, 2*MAX_PATH_LEN );
         strncat( ElmerLib, "/share/elmersolver/lib", 2*MAX_PATH_LEN );
      } else {
#if defined(WIN32) || defined(MINGW32)
	/* Should not get here unless WIN32 implements DLOPEN_API */
	GetModuleFileName(NULL, appPath, MAX_PATH_LEN);
	exeName = strrchr(appPath, '\\');
	n = (int)(exeName - appPath);
	if(n < 0) n = 0;
	if(n > MAX_PATH_LEN) n = MAX_PATH_LEN;
        strncat(ElmerLib, ELMER_PATH_SEPARATOR, 2*MAX_PATH_LEN);
	strncat(ElmerLib, appPath, n);
	strncat(ElmerLib, "\\..\\share\\elmersolver\\lib", 2*MAX_PATH_LEN);
#else
        strncat( ElmerLib, ELMER_PATH_SEPARATOR, 2*MAX_PATH_LEN );
	strncat( ElmerLib, ELMER_SOLVER_HOME, 2*MAX_PATH_LEN );
	strncat( ElmerLib, "/lib", 2*MAX_PATH_LEN );
#endif
      }
   }

   cptr = (char *)getenv( "ELMER_MODULES_PATH" );
   if ( cptr != NULL ) {
      strncat( ElmerLib, ELMER_PATH_SEPARATOR, 2*MAX_PATH_LEN);
      strncat( ElmerLib, cptr, 2*MAX_PATH_LEN);
   }

   try_open_solver(ElmerLib, Library, &Handle, ErrorBuffer);
   if ( Handle == NULL ) {
      fprintf(stderr, "%s", ErrorBuffer);
      exit(0);
   }
   
#ifdef HAVE_DLOPEN_API

   if ( (Function = (void(*)())dlsym( Handle,NewName)) == NULL && *abort_not_found )
   {
      fprintf( stderr, "Load: FATAL: Can't find procedure [%s]\n", NewName );
      exit(0);
   }

#elif defined(HAVE_LOADLIBRARY_API)

   if ( (Function = (void *)GetProcAddress(Handle,NewName)) == NULL && *abort_not_found )
   {
     fprintf( stderr,"Load: FATAL: Can't find procedure [%s]\n", NewName );
     exit(0);
   }

#endif

   return (void *)Function;
}

/*--------------------------------------------------------------------------
  INTERNAL: Execute given function returning integer value
  -------------------------------------------------------------------------*/
static int IntExec( int (STDCALLBULL *Function)(),void *Model )
{
   return (*Function)( Model );
}

/*--------------------------------------------------------------------------
   Execute given function returning integer value
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
int STDCALLBULL execintfunction_c( f_ptr Function,void *Model )
#else
int STDCALLBULL FC_FUNC(execintfunction,EXECINTFUNCTION) ( f_ptr Function,void *Model )
#endif
{
  return IntExec( (int (STDCALLBULL *)())*Function,Model );
}

/*--------------------------------------------------------------------------
   INTERNAL: Execute given function returning double value
  -------------------------------------------------------------------------*/
static void DoubleArrayExec( double *(STDCALLBULL *Function)(), void *Model,
               int *Node, double *Value, double *Array )
{
   (*Function)( Model,Node,Value,Array );
}

/*--------------------------------------------------------------------------
   Execute given function returning double value
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL execrealarrayfunction_c( f_ptr Function, void *Model,
										int *Node, double *Value, double *Array )
#else
void STDCALLBULL FC_FUNC(execrealarrayfunction,EXECREALARRAYFUNCTION)
     ( f_ptr Function, void *Model,
       int *Node, double *Value, double *Array )
#endif
{
   DoubleArrayExec( (double*(STDCALLBULL *)())*Function,Model,Node,Value, Array );
}

/*--------------------------------------------------------------------------
   INTERNAL: Execute given function returning double value
  -------------------------------------------------------------------------*/
static double DoubleExec( double (STDCALLBULL *Function)(), void *Model,
               int *Node, double *Value )
{
   return (*Function)( Model,Node,Value );
}

/*--------------------------------------------------------------------------
   Execute given function returning double value
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
double STDCALLBULL execrealfunction_c( f_ptr Function, void *Model,
									 int *Node, double *Value )
#else
double STDCALLBULL FC_FUNC(execrealfunction,EXECREALFUNCTION)
     ( f_ptr Function, void *Model,
       int *Node, double *Value )
#endif
{
   return DoubleExec( (double (STDCALLBULL *)())*Function,Model,Node,Value );
}

/*--------------------------------------------------------------------------
   INTERNAL: Execute given function returning double value
  -------------------------------------------------------------------------*/
static double ConstDoubleExec( double (STDCALLBULL *Function)(), void *Model,
			       double *x, double *y, double *z )
{
   return (*Function)( Model, x,y,z );
}

/*--------------------------------------------------------------------------
   Execute given function returning double value
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
double STDCALLBULL execconstrealfunction_c( f_ptr Function, void *Model,
										  double *x, double *y, double *z )
#else
double STDCALLBULL FC_FUNC(execconstrealfunction,EXECCONSTREALFUNCTION)
     ( f_ptr Function, void *Model,
       double *x, double *y, double *z )
#endif
{
   return ConstDoubleExec( (double (STDCALLBULL *)())*Function,Model,x,y,z );
}


/*--------------------------------------------------------------------------
   Return argument (just to fool Fortran type checking)
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void *STDCALLBULL addrfunc_c( void *Function )
#else
void *STDCALLBULL FC_FUNC(addrfunc,ADDRFUNC) ( void *Function )
#endif
{
   return (void *)Function;
}

/*--------------------------------------------------------------------------
   INTERNAL: Call solver routines at given address
  -------------------------------------------------------------------------*/
static void DoExecSolver(
  void (STDCALLBULL *SolverProc)(), void *Model, void *Solver, void *dt, void *Transient)
{
  (*SolverProc)( Model,Solver,dt,Transient ); 
  return;
}

/*--------------------------------------------------------------------------
   Call solver routines at given address
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL execsolver_c( f_ptr *SolverProc, void *Model, void *Solver,
								void *dt, void *Transient )
#else
void STDCALLBULL FC_FUNC(execsolver,EXECSOLVER)
     ( f_ptr *SolverProc, void *Model, void *Solver, void *dt, void *Transient )
#endif
{
  DoExecSolver( (void (STDCALLBULL *)())*SolverProc,Model,Solver,dt,Transient );
}

/*--------------------------------------------------------------------------
   INTERNAL: Call lin. solve routines at given address
  -------------------------------------------------------------------------*/
static int DoLinSolveProcs(
  int (STDCALLBULL *SolverProc)(), void *Model, void *Solver, void *Matrix, void *b, 
                void *x, void *n, void *DOFs, void *Norm )
{
   return (*SolverProc)( Model,Solver,Matrix,b,x,n, DOFs,Norm );
}


/*--------------------------------------------------------------------------
   Call lin. solver routines at given address
   -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
int STDCALLBULL execlinsolveprocs_c( f_ptr *SolverProc, void *Model, void *Solver,
						void *Matrix, void *b, void *x, void *n, void *DOFs, void *Norm )
#else
int STDCALLBULL FC_FUNC(execlinsolveprocs,EXECLINSOLVEPROCS)
     ( f_ptr *SolverProc, void *Model, void *Solver, void *Matrix, void *b, void *x, void *n, void *DOFs, void *Norm )
#endif
{
   return DoLinSolveProcs( (int (STDCALLBULL *)())*SolverProc,Model,Solver,Matrix,b,x,n,DOFs,Norm );
}

char *mtc_domath(char *);
void mtc_init(FILE *,FILE *, FILE *);

/*--------------------------------------------------------------------------
  This routine will call matc and return matc variable array values
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL matc_get_array(char *name, double *values, int *nrows, int *ncols )
#else
void STDCALLBULL FC_FUNC_(matc_get_array,MATC_GET_ARRAY) (char *name, 
           double *values, int *nrows, int *ncols )
#endif
{
  var_copy_transpose(name,values,*nrows,*ncols);
}

/*--------------------------------------------------------------------------
  This routine will call matc and return matc result
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL matc( char *cmd, char *Value, int *len )
#else
void STDCALLBULL FC_FUNC(matc,MATC) ( char *cmd, char *Value, int *len )
#endif
{
#define MAXLEN 8192

  static int been_here = 0;
  char *ptr, c, cc[32];
  int slen, start;
#pragma omp threadprivate(been_here)

  /* MB: Critical section removed since Matc library
   * modified to be thread safe */

   slen = *len;
   if ( been_here==0 ) {
     mtc_init( NULL, stdout, stderr ); 
     strcpy( cc, "format( 12,\"rowform\")" );
     mtc_domath( cc );
     been_here = 1;
   }

  c = cmd[slen];
  cmd[slen] = '\0';

  start = 0;
  if (strncmp(cmd,"nc:",3)==0) start=3;

  ptr = (char *)mtc_domath(&cmd[start]);
  if ( ptr )
  {
    strcpy( Value, (char *)ptr );
    *len = strlen(Value)-1; /* ignore linefeed! */

    if ( strncmp(Value, "MATC ERROR:",11)==0 || strncmp(Value,"WARNING:",8)==0 ) {
      if (start==0) {
          fprintf( stderr, "Solver input file error: %s\n", Value );
          fprintf( stderr, "...offending input line: %s\n", cmd );
          exit(0);
      } else {
        Value[0]=' ';
        *len = 0;
        }
    }
  } else {
    *len = 0;
    *Value = ' ';
  }
  cmd[slen]=c;
  }

/*--------------------------------------------------------------------------
  INTERNAL: execute user material function
  -------------------------------------------------------------------------*/
static double DoViscFunction(double (STDCALLBULL *SolverProc)(), void *Model, void *Element, void *Nodes, void *n,
     void *Basis, void *GradBasis, void *Viscosity, void *Velo, void *GradV )
{
   double s;
   s = (*SolverProc)( Model,Element,Nodes,n,Basis,GradBasis,
                   Viscosity, Velo, GradV );
   return s;
}

/*--------------------------------------------------------------------------
  This routine will call user defined material def. function
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
double STDCALLBULL materialuserfunction_c( f_ptr Function, void *Model, void *Element,
		void *Nodes, void *n, void *nd, void *Basis, void *GradBasis, void *Viscosity, void *Velo, void *gradV )
#else
double STDCALLBULL FC_FUNC(materialuserfunction,MATERIALUSERFUNCTION)
  ( f_ptr Function, void *Model, void *Element, void *Nodes, void *n, void *nd, void *Basis, void *GradBasis, void *Viscosity, void *Velo, void *gradV )
#endif
{
   return DoViscFunction( (double (STDCALLBULL *)())*Function,Model,Element,Nodes,n,Basis,
                  GradBasis,Viscosity,Velo,gradV );
}

/*--------------------------------------------------------------------------
  INTERNAL: execute user material function
  -------------------------------------------------------------------------*/
static void DoSimulationProc( void (STDCALLBULL *SimulationProc)(), void *Model )
{ 
  (*SimulationProc)( Model ); 
}

/*--------------------------------------------------------------------------
  This routine will call user defined material def. function
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL execsimulationproc_c( f_ptr Function, void *Model )
#else
void STDCALLBULL FC_FUNC(execsimulationproc,EXECSIMULATIONPROC)
     ( f_ptr Function, void *Model )
#endif
{
   DoSimulationProc( (void (STDCALLBULL *)())*Function,Model );
}


/*--------------------------------------------------------------------------
  INTERNAL: execute (Krylov) iterator 
  -------------------------------------------------------------------------*/
static void DoIterCall( void (STDCALLBULL *iterProc)(),
       void *x,void *b,void *ipar,void *dpar,void *work,
       void (STDCALLBULL *mvProc)(),
       void (STDCALLBULL *pcondProc)(),
       void (STDCALLBULL *pcondrProc)(),
       void (STDCALLBULL *dotProc)(),
       void (STDCALLBULL *normProc)(),
       void (STDCALLBULL *STOPC)() )
{ 
  (*iterProc)( x,b,ipar,dpar,work,mvProc,pcondProc, 
       pcondrProc,dotProc,normProc,STOPC );
}

/*--------------------------------------------------------------------------
  This routine will call (Krylov) iterator
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL itercall_c( f_ptr iterProc, void *x, void *b, void *ipar, void *dpar, void *work,
       f_ptr mvProc, f_ptr pcondProc, f_ptr pcondrProc, f_ptr dotProc, f_ptr normProc, f_ptr STOPC )
#else
void STDCALLBULL FC_FUNC(itercall,ITERCALL)
     ( f_ptr iterProc, void *x, void *b, void *ipar, void *dpar, void *work, 
       f_ptr mvProc, f_ptr pcondProc, f_ptr pcondrProc, f_ptr dotProc, f_ptr normProc, f_ptr STOPC )
#endif
{
   DoIterCall( (void (STDCALLBULL *)())*iterProc,x,b,ipar,dpar,work,
       (void (STDCALLBULL *)())*mvProc, 
       (void (STDCALLBULL *)())*pcondProc,
       (void (STDCALLBULL *)())*pcondrProc,
       (void (STDCALLBULL *)())*dotProc,
       (void (STDCALLBULL *)())*normProc,
       (void (STDCALLBULL *)())*STOPC );
}

/*--------------------------------------------------------------------------
  INTERNAL: execute localmatrix call
  -------------------------------------------------------------------------*/
static void DoLocalCall( void (STDCALLBULL *localProc)(),
  void *Model,void *Solver,void *G, void *F, void *Element,void *n,void *nd )
{ 
  (*localProc)( Model, Solver, G, F, Element, n, nd );
}

/*--------------------------------------------------------------------------
  This routine will call local matrix add-on
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL execlocalproc_c( f_ptr localProc, void *Model,void *Solver,
								void *G, void *F, void *Element,void *n,void *nd )
#else
void STDCALLBULL FC_FUNC(execlocalproc, EXECLOCALPROC )
     ( f_ptr localProc, void *Model,void *Solver,void *G, void *F, void *Element,void *n,void *nd )
#endif
{
   DoLocalCall( (void (STDCALLBULL *)())*localProc,Model,Solver,G,F,Element,n,nd );
}



/*--------------------------------------------------------------------------
  INTERNAL: execute complete localmatrix call
  -------------------------------------------------------------------------*/
static void DoLocalAssembly( void (STDCALLBULL *LocalAssembly)(),
  void *Model,void *Solver,void *dt,void *transient,void *M, void *D, void *S,void *F, void *Element,void *n,void *nd )
{ 
  (*LocalAssembly)( Model, Solver, dt, transient, M, D, S, F, Element, n, nd );
}

/*--------------------------------------------------------------------------
  This routine will call complete local matrix add-on
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL execlocalassembly_c( f_ptr LocalAssembly, void *Model,
		         void *Solver,void *dt,void *transient,
		         void *M, void *D, void *S,void *F,
		         void *Element,void *n,void *nd )
#else
void STDCALLBULL FC_FUNC(execlocalassembly, EXECLOCALASSEMBLY )
     ( f_ptr LocalAssembly, void *Model,void *Solver,void *dt,void *transient,void *M, void *D, void *S,void *F,void *Element,void *n,void *nd )
#endif
{
   DoLocalAssembly( (void (STDCALLBULL *)())*LocalAssembly,Model,Solver,dt,transient,M,D,S,F,Element,n,nd );
}



/*--------------------------------------------------------------------------
  INTERNAL: execute complete localmatrix call
  -------------------------------------------------------------------------*/
static void DoMatVecSubr( void (STDCALLBULL *matvec)(),
  void **SpMV,void *n,void *rows,void *cols,void *vals,void *u, void *v, void *reinit )
{ 
  (*matvec)( SpMV,n,rows,cols,vals,u,v,reinit);
}

/*--------------------------------------------------------------------------
  This routine will call complete local matrix add-on
  -------------------------------------------------------------------------*/
#ifdef USE_ISO_C_BINDINGS
void STDCALLBULL matvecsubrext_c( f_ptr matvec, void **SpMV, void *n, void *rows,
		                          void *cols, void *vals, void *u, void *v,void *reinit )
#else
void STDCALLBULL FC_FUNC(matvecsubr, MMATVECSUBR)
     ( f_ptr matvec, void **SpMV, void *n, void *rows, void *cols, void *vals, void *u, void *v,void *reinit )
#endif
{
   DoMatVecSubr( (void (STDCALLBULL *)())*matvec,SpMV,n,rows,cols,vals,u,v,reinit);
}
