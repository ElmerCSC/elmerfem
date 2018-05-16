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

/*******************************************************************************
 *
 *     String handling user functions.
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

/*
 * $Id: str.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: str.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:55  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"
#include "str.h"

VARIABLE *str_sprintf(var) VARIABLE *var;
{
  char *fmt = var_to_string(var);
  VARIABLE *res;
  int i;
 
  if (NEXT(var) != NULL)
  {
    for(i = 0; i < NCOL(NEXT(var)); i++)
    {
      str_p[i] = M(NEXT(var),0,i);
    }
    sprintf(str_pstr, fmt, 
            str_p[0],  str_p[1],  str_p[2],  str_p[3],  str_p[4],  str_p[5],
            str_p[6],  str_p[7],  str_p[8],  str_p[9],  str_p[10], str_p[11],
            str_p[12], str_p[13], str_p[14], str_p[15], str_p[16], str_p[17],
            str_p[18], str_p[19], str_p[20], str_p[21], str_p[22], str_p[23],
            str_p[24], str_p[25], str_p[26], str_p[27], str_p[28], str_p[29]);
  }
  else 
  {
    sprintf(str_pstr, fmt);
  }

  FREEMEM(fmt);

  res = var_temp_new(TYPE_STRING,1,strlen(str_pstr));
  for(i = 0; i < NCOL(res); i++)
  {
    M(res,0,i) = str_pstr[i];
  }

  return res;
}

VARIABLE *str_sscanf(var) VARIABLE *var;
{
  char *fmt = var_to_string(NEXT(var));
  char *str = var_to_string(var);
  VARIABLE *res;
  int i, got;
 
  got = sscanf(str, fmt, 
      &str_p[0],  &str_p[1],  &str_p[2],  &str_p[3],  &str_p[4],  &str_p[5],
      &str_p[6],  &str_p[7],  &str_p[8],  &str_p[9],  &str_p[10], &str_p[11],
      &str_p[12], &str_p[13], &str_p[14], &str_p[15], &str_p[16], &str_p[17],
      &str_p[18], &str_p[19], &str_p[20], &str_p[21], &str_p[22], &str_p[23],
      &str_p[24], &str_p[25], &str_p[26], &str_p[27], &str_p[28], &str_p[29]);

  FREEMEM(str);
  FREEMEM(fmt);

  res = NULL;
  if (got > 0) {
    res = var_temp_new(TYPE_DOUBLE,1,got);
    for(i = 0; i < got; i++)
    {
      M(res,0,i) = str_p[i];
    }
  }

  return res;
}

VARIABLE *str_matcvt(var) VARIABLE *var;
{
  VARIABLE *res = NULL;

  char *type = var_to_string(NEXT(var));
  double *d = MATR(var);

  int i, rlen;

  if (strcmp(type, "float")==0)
  {
    float *f;

    rlen = (MATSIZE(var)/2+7)/8;
    res = var_temp_new(TYPE(var), 1, rlen);
    f = (float *)MATR(res);
  
    for(i = 0; i < NCOL(var)*NROW(var); i++)
    {
      *f++ = (float)*d++;
    }
  }
  else if (strcmp(type, "int")==0)
  {
    int *n;

    rlen = (MATSIZE(var)/2+7)/8;
    res = var_temp_new(TYPE(var), 1, rlen);
    n = (int *)MATR(res);
  
    for(i = 0; i < NCOL(var)*NROW(var); i++)
    {
      *n++ = (int)*d++;
    }
  }
  else if (strcmp(type, "char")==0)
  {
    char *c;

    rlen = (MATSIZE(var)/8+7)/8;
    res = var_temp_new(TYPE(var), 1, rlen);
    c = (char *)MATR(res);
  
    for(i = 0; i < NCOL(var)*NROW(var); i++)
    {
      *c++ = (char)*d++;
    }
  }
  else 
  {
    fprintf(math_err, "matcvt: unknown result type specified.\n");
  }

  FREEMEM(type);

  return res;
}

VARIABLE *str_cvtmat(var) VARIABLE *var;
{
  VARIABLE *res = NULL;
  double *d;

  char *type = var_to_string(NEXT(var));

  int i, rlen;

  if (strcmp(type, "float")==0)
  {
    float *f = (float *)MATR(var);

    rlen = MATSIZE(var)/4; 
    res = var_temp_new(TYPE(var), 1, rlen);
    d = MATR(res);
  
    for(i = 0; i < rlen; i++)
    {
      *d++ = (double)*f++;
    }
  }
  else if (strcmp(type, "int")==0)
  {
    int *n = (int *)MATR(var);

    rlen = MATSIZE(var)/4;
    res = var_temp_new(TYPE(var), 1, rlen);
    d = MATR(res);
  
    for(i = 0; i < rlen; i++)
    {
      *d++ = (double)*n++;
    }
  }
  else if (strcmp(type, "char")==0)
  {
    char *c = (char *)MATR(var);

    rlen = MATSIZE(var);
    res = var_temp_new(TYPE(var), 1, rlen);
    d = MATR(res);
  
    for(i = 0; i < rlen; i++)
    {
      *d++ = (double)*c++;
    }
  }
  else 
  {
    fprintf(math_err, "matcvt: unknown result type specified.\n");
  }

  FREEMEM(type);

  return res;
}



VARIABLE *str_env(var) VARIABLE *var;
{
  VARIABLE *res = NULL;
  int i;
  char *name = var_to_string(var), *str;

  str = getenv(name);

  if ( str ) {
    res = var_temp_new(TYPE_STRING,1,strlen(str));
    for(i = 0; i < NCOL(res); i++)
    {
      M(res,0,i) = str[i];
    }
  }

  return res;
}

void str_com_init()
{
  static char *sprintfHelp =
  {
     "str = sprintf( fmt[, vec] )\n"
     "Return a string formated using fmt and values from vec. A call to\n"
     "corresponding C-language function is made.\n\n"
  };

  static char *sscanfHelp =
  {
     "vec = sscanf( str,fmt )\n"
     "Return values from str using format fmt. A call to corresponding C-language\n"
     "function is made.\n\n"
  };

  static char *matcvtHelp =
  {
     "special = matcvt( matrix, type )\n"
     "Makes a type conversion from MATC matrix double precision array to given\n"
     "type, which can be one of the following: \"int\", \"char\" or \"float\"\n\n"
     "\n"
     "SEE ALSO: cvtmat, fwrite\n"
  };

  static char *cvtmatHelp =
  {
     "matrix = cvtmat( special, type )\n"
     "Makes a type conversion from given type to MATC matrix.\n"
     "Type can be one of the following: \"int\", \"char\" or \"float\".\n\n"
     "\n"
     "SEE ALSO: fread, matcvt.\n"
  };

  static char *envHelp =
  {
     "str = env(name)\n"
     "return environment variable value.\n"
  };

  com_init( "sprintf", FALSE, TRUE, str_sprintf, 1, 2, sprintfHelp );
  com_init( "sscanf",  FALSE, TRUE, str_sscanf,  2, 2, sscanfHelp  );
  com_init( "matcvt",  FALSE, TRUE, str_matcvt,  2, 2, matcvtHelp  );
  com_init( "cvtmat",  FALSE, TRUE, str_cvtmat,  2, 2, cvtmatHelp  );
  com_init( "env",     FALSE, TRUE, str_env,  1, 1, envHelp  );
}
