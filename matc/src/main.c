/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library (in file ../LGPL-2.1); if not, write 
 * to the Free Software Foundation, Inc., 51 Franklin Street, 
 * Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/


#include <stdio.h>
#include <signal.h>
#include "../config.h"

#ifdef USE_READLINE
# ifdef HAVE_READLINE_READLINE_H
#  include <readline/readline.h>
#  include <readline/history.h>
# else
#  ifdef HAVE_READLINE_H
#   include <readline.h>
#   include <history.h>
#  endif
# endif
#endif 

/* prototype */
char *mtc_domath(char *);

int main( int argc, char **argv )
{
  char strt[2000];
  char *str;

#ifdef _OPENMP
  /* Set number of threads to 1, computations are single threaded anyway */
  omp_set_num_threads(1);
#endif

  (void)mtc_init( stdin, stdout, stderr );
  str = mtc_domath( "source(\"mc.ini\")" );

  signal( SIGINT, SIG_IGN );

  while( 1 )
  {
#ifdef USE_READLINE
      str = readline ("MATC> ");
      /* add to history */
      if (str && *str)
	add_history (str);

#else
      fgets( strt,  2000 , stdin);
      str = strt;      
#endif
      
/* kludge to enable exit. */
#if defined(WIN32) || defined(MINGW32)
      if( stricmp(str,"exit") == 0  || stricmp(str,"quit") == 0 )
#else
      if( strcasecmp(str,"exit") == 0  || strcasecmp(str,"quit") == 0 )
#endif
      {
	return 0;
      }
      if ( *str ) fprintf( stdout, "%s\n", mtc_domath( str ) );
      
#ifdef USE_READLINE
      free(str);
#endif
    }  
    return 0;
}
