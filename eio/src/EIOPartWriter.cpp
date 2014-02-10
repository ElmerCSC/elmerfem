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

#include "EIOPartWriter.h"
#include <stdio.h>
#include <string.h>

extern void make_filename(char *buf, const char *model, const char *suffix);

static char *extension[] = {
  "%s/part.%d.header",
  "%s/part.%d.nodes",
  "%s/part.%d.shared",
  "%s/part.%d.elements",
  "%s/part.%d.border"
};

enum { HEADER = 0, NODES, SHARED, ELEMENTS, BORDER};

EIOPartWriter::
EIOPartWriter(int& partCount, EIOModelManager *mm)
{
  parts = partCount;
  me = -1;
  manager = mm;
}

EIOPartWriter::~EIOPartWriter()
{
}

int EIOPartWriter::
createPartitioning(const char *dir)
{
  sprintf(newdir, "%s/partitioning.%d", dir, parts);
  return manager->makeDirectory(newdir);  
}

int EIOPartWriter::
activatePart(int part)
{
  me = part;
  openStreams();
  return 0;
}

int EIOPartWriter::
deactivatePart()
{
  closeStreams();
  me = -1;
  return 0;
}

int EIOPartWriter::
closePartitioning()
{
  if(me != -1)
    {
      closeStreams();
    }
  return 0;
}

int EIOPartWriter::
write_descriptor(int& nodeC,                 /* nodes */
		 int& sharedC,               /* shared nodes */
		 int& elementC,              /* elements (inner) */
		 int& borderC,               /* elements bordering the part */
		 int& boundaryC,
		 int& usedElementTypes, 
		 int* elementTypeTagsH,
		 int* elementCountByType)
{
  int i;
  fstream& str = meshFileStream[HEADER];
  str << nodeC << ' ' << elementC << ' ' << boundaryC << '\n';
  str << usedElementTypes << '\n';
  for(i = 0; i < usedElementTypes; ++i)
    str << elementTypeTagsH[i] << ' ' << elementCountByType[i] << '\n';

  str << sharedC << ' ' << borderC << '\n';
  return 0;
}

int EIOPartWriter::
write_node(int& tag, int& type, double *coord, int& partC, int *partsH)
{
  fstream& str = meshFileStream[NODES];
  fstream& str2 = meshFileStream[SHARED];

  str << tag << ' ' << type << ' ';

  str.setf(std::ios::scientific);
  str.precision(16);

  str << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
  if(partC > 1)
    {
      int i;
      str2 << tag << ' ' << partC << ' ';
      for(i = 0; i < partC; ++i)
	{
	  str2 << partsH[i] << ' ';
	}
      str2 << std::endl;
    }
  return 0;
}
int EIOPartWriter::
write_element(int& tag, int& body, int& type, int *nodes, int& border)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  fstream& str2 = meshFileStream[BORDER];

  str << tag << ' ' << body << ' ' << type << ' ';
  if(type == 303)
    {
      for(i = 0; i < 3; ++i)
	{
	  str << nodes[i] << ' ';
	}
    }
  str << std::endl;

  if(border)
    {
      str2 << tag << std::endl;
    }
  return 0;
}

void EIOPartWriter::
openStreams()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < partWriterFiles; ++i)
    {
      sprintf(filename, extension[i], newdir, me);
      manager->openStream(meshFileStream[i], filename, std::ios::out);
    }
}

void EIOPartWriter::
closeStreams()
{
  int i;
  for(i = 0; i < partWriterFiles; ++i)
    {
      manager->closeStream(meshFileStream[i]);
    }
}
