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

/***********************************************************************
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#include <string.h>
#include <ctype.h>

#include "EIOConstraintAgent.h"

extern void make_filename(char *buf, const char *model, const char *suffix);

static char *extension[] = {
  "constraints.header",
  "constraints.table"
};

enum { HEADER = 0, TABLE };

EIOConstraintAgent::EIOConstraintAgent(EIOModelManager *mm)
{
  manager = mm;
}

EIOConstraintAgent::~EIOConstraintAgent()
{
}

int EIOConstraintAgent::
openConstraints()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < constraintFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(constraintFileStream[i], extension[i], std::ios::in);
    }

  fstream& str = constraintFileStream[HEADER];
  str >> constraintCount;
  return 0;
}

int EIOConstraintAgent::
closeConstraints()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < constraintFiles; ++i)
    {
      manager->closeStream(constraintFileStream[i]);
    } 
  return 0;
}

int EIOConstraintAgent::
descriptor(int& cCount)
{
  cCount = constraintCount;
  return 0;
}

static int cstep = 0;
int EIOConstraintAgent::
nextConstraint(int& tag, int& field, 
	       int *constraintType, double *constraintValue)
{
  int i;
  fstream& str = constraintFileStream[TABLE];
  if(cstep == constraintCount)
    {
      streampos pos = 0;
      filebuf *fbuf = str.rdbuf();
      fbuf->pubseekpos(pos, std::ios::in);
      cstep = 0;
      return -1;
    }
  str >> tag >> field;
  for(i = 0; i < field; ++i)
    {
      char name[10];
      str >> name;
      int j, len;
      len = strlen(name);
      for(j = 0; j < len; ++j) name [j] = toupper(name[j]);
      if(!strcmp(name, "U"))
	{
	  constraintType[i] = 0;
	}
      else if(!strcmp(name, "V"))
	{
	  constraintType[i] = 1;
	}
      else if(!strcmp(name, "T"))
	{
	  constraintType[i] = 100;
	}
      else
	{
	  constraintType[i] = 999;
	}
      str >> constraintValue[i];
    }
  ++cstep;
  return 0;
}

