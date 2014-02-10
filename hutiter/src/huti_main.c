/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 *
 * huti_main.c - HUTIter libarary auxiliary routines
 *
 * $Id: huti_main.c,v 1.2 2005/05/04 20:18:42 vierinen Exp $
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "huti_defs.h"


#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Global HUTI variables */

int huti_num_of_procs;

/*
 *
 * HUTI_Init - Initialize HUTI environment
 *
 */
void HUTI_Init()
{
  char *evname;
  static int huti_init_done = FALSE;

  if (huti_init_done)
    return;

  /* Read environment variables */

  if ((evname = getenv(NUMBER_OF_PROCESSORS)) == NULL)
    huti_num_of_procs = 1;
  else
    if ((huti_num_of_procs = atoi(evname)) == 0) {
#if 0
      (void) HUTI_Errexit("Environment variable ",NUMBER_OF_PROCESSORS,
		   " has an illegal value ", evname, NULL);
#endif
  }

  huti_init_done = TRUE;
}

/*
 *
 *
 *
 */

void HUTI_Exit()
{
}

/*
 *
 * HUTI_Errexit - Print error messages and exit. There can be several
 *                messages as an argument.
 *
 */

#if 0
HUTI_Errexit(va_alist)
     va_dcl
{
  va_list ap;
  char *buffer[MAX_ERRMSGS];
  int msgno;

  msgno = 0;
  va_start( ap );
  while ( (buffer[msgno] = va_arg(ap, char *)) != NULL )
    fprintf(stderr, "%s", buffer[msgno++]);
  fprintf(stderr, "\n");
  fflush(stderr);
  va_end( ap );

  exit(1);
}
#endif

/*
 *
 *
 *
 */
