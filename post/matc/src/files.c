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
 *     IO handling for MATC.
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
 * $Id: files.c,v 1.1.1.1 2005/04/14 13:29:14 vierinen Exp $ 
 *
 * $Log: files.c,v $
 * Revision 1.1.1.1  2005/04/14 13:29:14  vierinen
 * initial matc automake package
 *
 * Revision 1.2  1998/08/01 12:34:36  jpr
 *
 * Added Id, started Log.
 * 
 *
 */

#include "elmer/matc.h"
#include "str.h"

/***********************************************************************
|
|  FILES.C - Last Edited 10. 8. 1988
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

#define FILE_BINARY 0
#define FILE_ASCII 1

#define MAXFILES 32
static FILE *fil_fps[MAXFILES];
static FILE *fil_fps_save[3];

VARIABLE *fil_fread(var) VARIABLE *var;
{
  VARIABLE *res;
  FILE *fp;

  int i, ind, len;

  ind = *MATR(var);
  if (ind < 0 || ind >= MAXFILES)
  {
    error("fread: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fread: File not open.\n");
  }
  fp = fil_fps[ind];

  if (feof(fp))
  {
    clearerr(fp);
    error("fread: end of file detected.\n");
  }

  len = *MATR(NEXT(var));
  if (len <= 0) 
  {
    error("fread: invalid length specified.\n");
  }
  res = var_temp_new(TYPE_DOUBLE, 1, (len+sizeof(double)-1)>>3);
  fread(MATR(res), 1, len, fp);

  if (feof(fp))
  {
    clearerr(fp);
    error("fread: end of file detected.\n");
  }

  if (ferror(fp))
  {
    clearerr(fp);
    error("fread: error reading file.\n");
  }

  return res;
}

VARIABLE *fil_fwrite(var) VARIABLE *var;
{
  int i, ind, len;
  FILE *fp;

  ind = *MATR(var);
  if (ind < 0 || ind >= MAXFILES)
  {
    error("fwrite: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fwrite: File not open.\n");
  }
  fp = fil_fps[ind];

  if (NEXT(NEXT(var)) != NULL)
  {
    len = *MATR(NEXT(NEXT(var)));
    if (len > MATSIZE(NEXT(var)))
    {
      error("fwrite: attempt to write more data than provided.\n");
    }
  }
  else
  {
    len = MATSIZE(NEXT(var));
  }
  fwrite(MATR(NEXT(var)), 1, len, fp);

  if (ferror(fp))
  {
    clearerr(fp);
    error("fwrite: error writing file.\n");
  }

  return (VARIABLE *)NULL;
}

VARIABLE *fil_fscanf(var) VARIABLE *var;
{ 
  VARIABLE *res;
  FILE *fp;

  char *fmt = var_to_string(NEXT(var));
  int i, ind, got;

  ind = *MATR(var);
  if (ind < 0 || ind >= MAXFILES)
  {
    error("fscanf: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fscanf: File not open.\n");
  }
  fp = fil_fps[ind];

  if (feof(fp))
  {
    clearerr(fp);
    error("fscanf: end of file detected.\n");
  }

  got = fscanf(fp, fmt, 
      &str_p[0],  &str_p[1],  &str_p[2],  &str_p[3],  &str_p[4],  &str_p[5],
      &str_p[6],  &str_p[7],  &str_p[8],  &str_p[9],  &str_p[10], &str_p[11],
      &str_p[12], &str_p[13], &str_p[14], &str_p[15], &str_p[16], &str_p[17],
      &str_p[18], &str_p[19], &str_p[20], &str_p[21], &str_p[22], &str_p[23],
      &str_p[24], &str_p[25], &str_p[26], &str_p[27], &str_p[28], &str_p[29]);

  res = NULL;
  if (got > 0) {
    res = var_temp_new(TYPE_DOUBLE,1,got);
    for(i = 0; i < got; i++)
    {
      M(res,0,i) = str_p[i];
    }
  }
  FREEMEM(fmt);

  if (feof(fp))
  {
    clearerr(fp);
    error("fscanf: end of file detected.\n");
  }

  if (ferror(fp))
  {
    clearerr(fp);
    error("fscanf: error reading file.\n");
  }

  return res;
}

VARIABLE *fil_fgets(var) VARIABLE *var;
{
  VARIABLE *res;
  FILE *fp;

  int i, ind;

  ind = *MATR(var);
  if (ind < 0 || ind >= MAXFILES)
  {
    error("fgets: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fgets: File not open.\n");
  }
  fp = fil_fps[ind];

  if (feof(fp))
  {
    clearerr(fp);
    error("fgets: end of file detected.\n");
  }

  fgets(str_pstr, STR_MAXLEN, fp);

  if (feof(fp))
  {
    clearerr(fp);
    error("fgets: end of file detected.\n");
  }

  if (ferror(fp))
  {
    clearerr(fp);
    error("fgets: error reading file.\n");
  }

  res = var_temp_new(TYPE_STRING, 1, strlen(str_pstr)-1);
  for(i = 0; i < strlen(str_pstr)-1; i++)
    M(res,0,i) = str_pstr[i];

  return res;
}

VARIABLE *fil_fprintf(var) VARIABLE *var;
{
  int i, ind;
  char *str;
  FILE *fp;

  ind = *MATR(var);
  if (ind < 0 || ind >= MAXFILES)
  {
    error("fprintf: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fprintf: File not open.\n");
  }
  fp = fil_fps[ind];

  var = str_sprintf(NEXT(var));
  str = var_to_string(var);  
  fprintf(fp, "%s",str);

  var_delete_temp(var);
  FREEMEM(str);

  if (ferror(fp))
  {
    clearerr(fp);
    error("fprintf: error writing file.\n");
  }

  return (VARIABLE *)NULL;
}

VARIABLE *fil_fputs(var) VARIABLE *var;
{
  char *str = var_to_string(NEXT(var)); 
  int ind = *MATR(var);
  FILE *fp;

  if (ind < 0 || ind >= MAXFILES)
  {
    error("fputs: Invalid file number.\n");
  }
  else if (fil_fps[ind] == NULL)
  {
    error("fputs: File not open.\n");
  }
  fp = fil_fps[ind];

  fprintf(fp, "%s", str);

  FREEMEM(str);

  if (ferror(fp))
  {
    clearerr(fp);
    error("fprintf: error writing file.\n");
  }

  return (VARIABLE *)NULL;
}

VARIABLE *fil_fopen(var) VARIABLE *var;
{
  VARIABLE *res;

  char *name, *mode;
  int file;

  mode = var_to_string(NEXT(var));
  name = var_to_string(var);

  for(file = 0; file < MAXFILES; file++)
  { 
    if (fil_fps[file] == NULL) break;
  }

  if (file >= MAXFILES)
  {
    error("fopen: maximum number of files already open.\n");
  }

  if ((fil_fps[file] = fopen(name, mode)) == (FILE *)NULL)
  {
    error("fopen: can't open file: %s.\n", name);
  }

  switch(file)
  {
    case 0:
      fil_fps_save[0] = math_in;
      math_in = fil_fps[0];
    break;

    case 1:
      fil_fps_save[1] = math_out;
      math_out = fil_fps[1];
    break;

    case 2:
      fil_fps_save[2] = math_err;
      math_err = fil_fps[2];
    break;
  }

  res = var_temp_new(TYPE_DOUBLE, 1, 1);
  M(res,0,0) = file;

  FREEMEM(name);
  FREEMEM(mode);

  return res;
}

VARIABLE *fil_fclose(var) VARIABLE *var;
{
  int file = *MATR(var);

  if (file < 0 || file >= MAXFILES)
  {
    error("fclose: Invalid file number.\n");
  }

  switch(file)
  {
    case 0:
      math_in = fil_fps_save[0];
      if (
           fil_fps[0] != math_out && fil_fps[0] != NULL
         ) fclose(fil_fps[0]);
      fil_fps[0] = math_in;
    break;

    case 1:
      math_out = fil_fps_save[1];
      if (
           fil_fps[1] != math_out && fil_fps[1] != NULL
         ) fclose(fil_fps[1]);
      fil_fps[1] = math_out;
    break;

    case 2:
      math_err = fil_fps_save[2];
      if (
           fil_fps[2] != math_err && fil_fps[2] != NULL
         ) fclose(fil_fps[2]);
      fil_fps[2] = math_err;
    break;

    default:
      if (fil_fps[file] != NULL)
        fclose(fil_fps[file]);
      fil_fps[file] = NULL;
    break;
  }

  return (VARIABLE *)NULL;
}

VARIABLE *fil_freopen(var) VARIABLE *var;
{
  int file = *MATR(var);

  fil_fclose(var);
  fil_fps[file] = NULL;
  var = fil_fopen(NEXT(var));
  var_delete_temp(var);

  return (VARIABLE *)NULL;
}

VARIABLE *fil_save(ptr) VARIABLE *ptr;
{
  VARIABLE *tmp;
  char *file;
  FILE *fp;

  int i, j, ascflg = FALSE;

  file = var_to_string(ptr);
   
  if ((fp = fopen(file, "w")) == (FILE *)NULL)
  {
    error( "save: can't open file: %s.\n", file );
  }
   
  tmp = NEXT(ptr);

  if (NEXT(NEXT(ptr)) != (VARIABLE *)NULL) 
    ascflg = M(NEXT(NEXT(ptr)), 0, 0);

  if (ascflg)
  {

    fprintf(fp, "%d %d %d %d\n", FILE_ASCII, TYPE(tmp), NROW(tmp), NCOL(tmp));
    if (ferror(fp))
    {
       fclose(fp); error("save: error writing file.\n");
    }

    for(i = 0; i < NROW(tmp); i++) 
      for(j = 0; j < NCOL(tmp); j++)
      {
        fprintf(fp, "%e\n", M(tmp, i, j));
        if (ferror(fp))
        {
           fclose(fp); error("save: error writing file.\n");
        }
      }
  }

  else
  {

    fprintf(fp, "%d %d %d %d\n", FILE_BINARY, TYPE(tmp), NROW(tmp), NCOL(tmp));
    if (ferror(fp))
    {
       fclose(fp); error("save: error writing file.\n");
    }

    fwrite(MATR(tmp), 1, MATSIZE(tmp), fp);
    if (ferror(fp))
    {
       fclose(fp); error("save: error writing file.\n");
    }
  }

  fclose(fp); FREEMEM(file);

  return NULL;
}

VARIABLE *fil_load(ptr) VARIABLE *ptr;
{
  int i, j, ftype, type, ncol, nrow;

  VARIABLE *res;

  char *file;
  FILE *fp;

  file = var_to_string(ptr);
  
  if ((fp = fopen(file, "r")) == (FILE *)NULL)
  {
    error( "load: can't open file: %s.\n", file );
  }

  fscanf(fp, "%d %d %d %d", &ftype, &type, &nrow, &ncol);

  if (ferror(fp)) {
    fclose(fp); error("load: error reading file.n");
  }
   
  res = var_temp_new(type, nrow, ncol);     

  if (ftype == FILE_ASCII)
  {
    for(i = 0; i < nrow; i++) 
      for(j = 0; j < ncol; j++)
      {
        fscanf(fp, "%lf", &M(res, i, j));
        if (ferror(fp))
        {
           fclose(fp); error("load: error reading file.\n");
         }
      }
  }
  else
  {
    fgetc(fp);
    fread(MATR(res), 1, MATSIZE(res), fp);
    if (ferror(fp))
    {
        fclose(fp); error("load: error reading file.\n");
     }
  }

  fclose(fp); FREEMEM(file);

  return res;
}
   
void fil_com_init()
{
  static char *freadHelp =
  {
      "str = fread( fp,len )\n\n"
      "Read len character  from file fp.  File pointer fp shoud have been\n"
      "obtained from a call to fopen or freopen, or be the standard input\n"
      "file stdin. Characters are returned as function value.\n"
      "\n"
      "SEE ALSO: fopen,freopen,fgets,fscanf,matcvt,cvtmat.\n"
  };

  static char *fscanfHelp =
  {
      "vec = fscanf( fp,format )\n\n"
      "Read file fp as given in format. Format is equal to C-language format\n"
      "File pointer fp shoud have been obtained from a call to fopen or freopen,\n"
      "or be the standard input.\n"
      "\n"
      "SEE ALSO: fopen,freopen,fgets,fread,matcvt,cvtmat.\n"
  };

  static char *fgetsHelp =
  {
      "str = fgets( fp )\n\n"
      "Read next line from fp. File pointer fp shoud have been obtained from a call\n"
      "to fopen or freopen or be the standard input.\n"
      "\n"
      "SEE ALSO: fopen,freopen,fread,fscanf,matcvt,cvtmat.\n"
  };

  static char *fwriteHelp =
  {
      "n = fwrite( fp, buf,len )\n\n"
      "Write len bytes form buf to file fp. File pointer fp shoud have been obtained\n"
      "from a call to fopen or freopen or be the standard output (stdout) or standard\n"
      "error (stderr). Return value is number of characters actually written.\n"
      "\n"
      "SEE ALSO: fopen,freopen,fputs,fprintf,matcvt,cvtmat.\n"
  };

  static char *fprintfHelp =
  {
      "n = fprintf( fp, format[, vec] )\n\n"
      "Write formatted string to file fp. File pointer fp shoud have been obtained\n"
      "from a call to fopen or freopen or be the standard output (stdout) or standard\n"
      "error (stderr). The format is equal to C-language format.\n"
      "\n"
      "SEE ALSO: fopen,freopen,fputs,fwrite,matcvt,cvtmat.\n"
  };

  static char *fputsHelp =
  {
      "fputs( fp, str )\n\n"
      "Write line to file fp. File pointer fp should have been obtained from a call\n"
      "to fopen or freopen or be the standard input (stdin).\n"
      "\n"
      "SEE ALSO: fopen,freopen,fwrite,matcvt,cvtmat.\n"
  };

  static char *fopenHelp =
  {
      "fp = fopen( name, mode )\n\n"
      "Open file given name and access mode. The most usual modes are \"r\" for reading\n"
      "and \"w\" for writing. Return value fp is used in functions reading and writing\n"
      "the file.\n"
      "\n"
      "SEE ALSO: freopen.\n"
  };

  static char *freopenHelp =
  {
      "fp = freopen( fp, name, mode )\n\n"
      "Reopen file given previous file pointer, name and access mode. The most usual modes\n"
      "are \"r\" for reading and \"w\" for writing. Return value fp is used in functions  \n"
      "reading and writing the file.\n"
      "\n"
      "SEE ALSO: fopen.\n"
  };

  static char *fcloseHelp =
  {
      "fclose( fp )\n\n"
      "Close file previously opened with fopen or freopen.\n"
      "\n"
      "SEE ALSO: fopen, freopen.\n"
  };

  static char *saveHelp =
  {
      "save( name, matrix[, ascii_flag] )\n\n"
      "Save matrix in file with name given as first parameter. If ascii_flag is\n"
      "given and is not zero the file will be in ascii format, otherwise matrix\n"
      "is saved in double precsision binary format. In either case the first line\n"
      "of the file contains four digits (in ascii):\n\n"
      "ascii_flag 0 NROW(matrix) NCOL(matrix).\n"
      "\n"
      "SEE ALSO: load.\n"
  };

  static char *loadHelp =
  {
      "matrix = load( name )\n\n"
      "Load matrix from a file given name and in format used by save-command.\n"
      "\n"
      "SEE ALSO: save.\n"
  };

  com_init( "fread",   FALSE, FALSE, fil_fread,   2, 2, freadHelp   );
  com_init( "fscanf",  FALSE, FALSE, fil_fscanf,  2, 2, fscanfHelp  );
  com_init( "fgets",   FALSE, FALSE, fil_fgets,   1, 1, fgetsHelp   );
  com_init( "fwrite",  FALSE, FALSE, fil_fwrite,  2, 3, fwriteHelp  );
  com_init( "fprintf", FALSE, FALSE, fil_fprintf, 2, 3, fprintfHelp );
  com_init( "fputs",   FALSE, FALSE, fil_fputs,   2, 2, fputsHelp   );
  com_init( "fopen",   FALSE, FALSE, fil_fopen,   2, 2, fopenHelp   );
  com_init( "freopen", FALSE, FALSE, fil_freopen, 3, 3, freopenHelp );
  com_init( "fclose",  FALSE, FALSE, fil_fclose,  1, 1, fcloseHelp  );
  com_init( "save",    FALSE, FALSE, fil_save,    2, 3, saveHelp    );
  com_init( "load",    FALSE, FALSE, fil_load,    1, 1, loadHelp    );

  fil_fps[0] = fil_fps_save[0] = stdin;
  fil_fps[1] = fil_fps_save[1] = stdout; 
  fil_fps[2] = fil_fps_save[2] = stderr;
}
