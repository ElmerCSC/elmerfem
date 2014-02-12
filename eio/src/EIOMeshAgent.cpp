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

#include "EIOMeshAgent.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <iostream>


//extern "C" void *
//void*
//bsearch(const void *key, const void *base, size_t nmemb,
//        size_t size, int (*compar)(const void *, const void *));

static int step = 0;

int CDECL nodecomp(const void *a, const void *b)
{
  cache_node *aptr = (cache_node*)a;
  cache_node *bptr = (cache_node*)b;

  if(aptr->tag <  bptr->tag) return -1;
  else if(aptr->tag >  bptr->tag) return 1;
  return 0;
}

void rewind_stream(fstream& str)
{
  streampos pos = 0;
  filebuf *fbuf = str.rdbuf();
  fbuf->pubseekpos(pos, std::ios::in);
}

int elementNodes(const int type)
{
  int cnt;

  cnt = type - 100*(type/100);
  return cnt;
}

void make_filename(char *buf, const char *model, const char *suffix)
{
  buf[0] = '\0';
  strcat(buf, model);
  strcat(buf, suffix);
}

static char *sequential_extensions[] = {
  "/mesh.header",
  "/mesh.nodes",
  "/mesh.elements",
  "/mesh.boundary"
};
static char *parallel_extensions[] = {
  "%s/part.%d.header",
  "%s/part.%d.nodes",
  "%s/part.%d.elements",
  "%s/part.%d.boundary",              //"%s/%s.part.%d.boundary",
  "%s/part.%d.shared"
};

static char **extension = (char **)0;

enum { HEADER = 0, NODES, ELEMENTS, BOUNDARY, SHARED };

EIOMeshAgent::EIOMeshAgent(EIOModelManager *mm, int split, int part)
{
  manager = mm;

  parts = split;
  me = part;
  if(me > 0) parallel = 1;
  else parallel = 0;

  dim = 3;
  clist = (cache_node *) NULL;

  elementTypeTags = (int*) 0;
  elementTypeCount = (int*) 0;

  if(parallel)
    {
      meshFiles = 5;
      extension = parallel_extensions;
    }
  else
    {
      meshFiles = 4;
      extension = sequential_extensions;      
    }
  meshFileStream = new fstream[meshFiles];
}

EIOMeshAgent::~EIOMeshAgent()
{
}

int EIOMeshAgent::
createMesh(const char *dir)
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < meshFiles; ++i)
    {
      make_filename(filename, dir, extension[i]);
      manager->openStream(meshFileStream[i], filename, std::ios::out);
    }

  return 0;
}

int EIOMeshAgent::
openMesh(const char *dir)
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < meshFiles; ++i)
    {
      if(parallel) { //  && (i != BOUNDARY)) {
	sprintf(newdir, "%s/partitioning.%d", dir, parts);
	sprintf(filename, extension[i], newdir, me);
      }
      else
	make_filename(filename, dir, extension[i]);
      if ( !manager->openStream(meshFileStream[i],
        filename, std::ios::in ) )
          return -1;
    }

  // Read header
  fstream& str = meshFileStream[HEADER];
  str >> nodeCount;
  str >> elementCount;
  str >> boundaryElementCount;
  str >> elementTypes;
  
  elementTypeTags = new int[elementTypes];
  elementTypeCount = new int[elementTypes];

  for(i = 0; i < elementTypes; ++i)
    {
      int etype, ecount;
      str >> etype >> ecount;
      elementTypeTags[i] = etype;
      elementTypeCount[i] = ecount;
    }

  if(parallel)
    {
      str >> sharedNodeCount >> borderElementCount;
    }


  step = 0;
  clist = (cache_node *)NULL;

  return 0;
}

int EIOMeshAgent::
closeMesh()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < meshFiles; ++i)
    {
      manager->closeStream(meshFileStream[i]);
    } 

  if (clist) delete []clist;
  clist = (cache_node *)NULL;

  delete []elementTypeTags;
  delete []elementTypeCount;

  return 0;
}

int EIOMeshAgent::
read_descriptor(int &nodeC, int& elementC, int& boundaryElementC, 
		int& usedElementTypes, int* usedElementTypeTags,	
		int *usedElementTypeCount)
{
  int i;
  nodeC = nodeCount;
  elementC = elementCount;
  boundaryElementC = boundaryElementCount;
  usedElementTypes = elementTypes;

  for(i = 0; i < elementTypes; ++i)
    {
      usedElementTypeTags[i] = elementTypeTags[i];
      usedElementTypeCount[i] = elementTypeCount[i];
    }
  return 0;
}

int EIOMeshAgent::
read_partDescriptor(int& shared)
{
  shared = sharedNodeCount;
  return 0;
}

int EIOMeshAgent::
read_nextElementConnections(int& tag, int& part, int& body, int& type, int* pdofs, int* nodes )
{
  int i, gotnodal;
  char typestr[32], tagstr[32];
  fstream& str = meshFileStream[ELEMENTS];
  if(step == elementCount)
    {
      rewind_stream(str);
      step = 0;
      return -1;
    }
  for( i=0; i<6; i++ ) pdofs[i] = 0;
  str >> tagstr >> body >> typestr ;
  part = 0;
  sscanf( tagstr, "%d/%d", &tag, &part );

  gotnodal = 0;
  for( i=0; i<strlen(typestr); i++ )
  {
     switch( typestr[i] ) {
     case('n'):
        sscanf( &typestr[i+1], "%d",  &pdofs[0] );
        gotnodal = 1;
        break;
     case('e'):
        sscanf( &typestr[i+1], "%d",  &pdofs[1] );
        break;
     case('f'):
        sscanf( &typestr[i+1], "%d",  &pdofs[2] );
        break;
     case('d'):
        sscanf( &typestr[i+1], "%d",  &pdofs[3] );
        break;
     case('b'):
        sscanf( &typestr[i+1], "%d",  &pdofs[4] );
        break;
     case('p'):
        sscanf( &typestr[i+1], "%d",  &pdofs[5] );
        break;
     }
  }
  typestr[3] = '\0';
  sscanf( typestr, "%d",  &type );
                                                                                                             
  int elNodes = elementNodes(type);
  for(i = 0; i < elNodes; ++i)
   {
     str >> nodes[i];
   }
   if ( !gotnodal ) pdofs[0] = 1;

  ++step;
  return 0;
}

int EIOMeshAgent::
read_nextElementCoordinates(int& tag, int& body, int& type, int* nodes, double *coord)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  if(step == elementCount)
    {
      rewind_stream(str);
      step = 0;
      return -1;
    }
  else if(step == 0)
    {
      cache_nodes();
    }

  str >> tag >> body >> type;
  int elNodes = elementNodes(type);
  for(i = 0; i < elNodes; ++i)
    {
      str >> nodes[i];
    }
  for(i = 0; i < elNodes; ++i)
    {
      if(!copy_coords(coord+i*3, nodes[i]))
	{
	  std::cout << tag << " exiting" << std::endl;
	  exit(14);
	}
    }
  ++step;
  return 0;
}


int EIOMeshAgent::
read_nextBoundaryElement(int& tag, int& part, int& boundary, int& leftElement,
                         int& rightElement, int& type, int* nodes, double* coord)
{
  int i;
  char tagstr[32];
  fstream& str = meshFileStream[BOUNDARY];

  if(step == boundaryElementCount)
    {
      rewind_stream(str);
      step = 0;
      return -1;
    }
  else if(step == 0)
    {
      cache_nodes();
    }

  str >> tagstr >> boundary >> leftElement >> rightElement;
  part = 0;
  sscanf( tagstr, "%d/%d", &tag, &part );

  str >> type;
  int elNodes = elementNodes(type);
  for(i = 0; i < elNodes; ++i)
    {
      str >> nodes[i];
    }

  if(parallel)
    {
      int ok = 1;
      for(i = 0; i < elNodes; ++i)
	if(search_node(nodes[i]) == NULL)
	  {
	    ok = 0;
	    break;
	  }
      if(!ok)
	{
	  ++step;
	  return read_nextBoundaryElement(tag, part, boundary, 
					  leftElement, rightElement, 
					  type, nodes, coord);
	}
    }

  for(i = 0; i < elNodes; ++i)
    {
      if(!copy_coords(coord+i*3, nodes[i]))
	{
	  exit(14);
	}
    }

  ++step;
  return 0;
}

int EIOMeshAgent::
write_descriptor(int& nodeC, int& elementC, int& boundaryElementC, 
		 int& usedElementTypes, int* elementTypeTags,
		 int* elementCountByType)
{
  int i;
  fstream& str = meshFileStream[HEADER];
  str << nodeC << ' ' << elementC << ' ' << boundaryElementC << '\n';
  str << usedElementTypes << '\n';
  for(i = 0; i < usedElementTypes; ++i)
    str << elementTypeTags[i] << ' ' << elementCountByType[i] << '\n';
  return 0;
}

int EIOMeshAgent::
write_node(int& tag, int& type, double* coord)
{
  int i;
  fstream& str = meshFileStream[NODES];
  str << tag << ' ' << type << ' ';
  str.setf(std::ios::scientific);
  str.precision(16);

  for(i = 0; i < dim; ++i)
    {
      str << coord[i] << ' ';
    }

  str << std::endl;
  return 0;
}

int EIOMeshAgent::
write_elementConnections(int& tag, int& body, int& type, int* nodes)
{
  int i;
  fstream& str = meshFileStream[ELEMENTS];
  str << tag << ' ' << body << ' ' << type << ' ';  
  int elNodes = elementNodes(type);
  for(i = 0; i < elNodes; ++i)
    {
      str << nodes[i] << ' ';
    }
  str << std::endl;
  return 0;
}
 

int EIOMeshAgent::
write_boundaryElement(int& tag, int& boundary, int& leftElement, int& rightElement,
		                  int& type, int* nodes)
{
  int i;
  fstream& str = meshFileStream[BOUNDARY];
  int elNodes = elementNodes(type);

  str << tag << ' ' << boundary << ' ';
  str << leftElement << ' ';
  str << rightElement << ' ';
  str << type << ' ';
  for(i = 0; i < elNodes; ++i)
    {
      str << nodes[i] << ' ';
    }
  str << std::endl;
  return 0;
}

int EIOMeshAgent::
read_allNodes(int *tags,double* coord)
{
  int i = 0;
  int pt = 0;

  cache_nodes();
  for(i = 0; i < nodeCount; ++i)
    {
      tags[i] = clist[i].tag;
      coord[pt]   = clist[i].x;
      coord[pt+1] = clist[i].y;
      coord[pt+2] = clist[i].z;
      pt += 3;
    }
  return 0;
}

int EIOMeshAgent::
read_sharedNode(int& tag, 
		int& constraint,      
		double *coord, 
		int& partcount, 
		int *parts)
{
  int i;
  fstream& str = meshFileStream[SHARED];

  if(step == sharedNodeCount)
    {
      rewind_stream(str);
      step = 0;
      return -1;
    }
  else if(step == 0)
    {
      cache_nodes();
    }

  str >> tag >> partcount;
  for(i = 0; i < partcount; ++i) str >> parts[i];
 
  cache_node *retval = search_node(tag);
  if(retval == NULL) 
    {
      std::cout << "Partition error: PANIC PANIC!!! "<< tag << std::endl;
      exit(23);
    }
  else
    {
      constraint = retval->constraint;
      coord[0] = retval->x;
      coord[1] = retval->y;
      coord[2] = retval->z;
    }
  ++step;
  return 0;
}

void EIOMeshAgent::
cache_nodes()
{
  if(!clist)
    {
      clist = new cache_node[nodeCount];
      fstream& str = meshFileStream[NODES];
      for(int i = 0; i < nodeCount; ++i)
	{
	  if(parallel) // assume that everything is sorted by splitter
	    {
	      str >> clist[i].tag >> clist[i].constraint
		  >> clist[i].x >> clist[i].y >> clist[i].z;
	    }
	  else
	    {
	      int tag;
	      str >> tag;
	      clist[tag-1].tag = tag;
	      str >> clist[tag-1].constraint
                  >> clist[tag-1].x >> clist[tag-1].y >> clist[tag-1].z;
	    }
	}
      rewind_stream(str);
    }
}

int EIOMeshAgent::
copy_coords(double *target, const int address)
{
  int found = 1;
  if(parallel)
    {
      cache_node *retval = search_node(address);
      if(retval == NULL) 
	{
	  std::cout << "Partition error: PANIC PANIC!!! "<< address << std::endl;
	  found = 0;
	}
      else
	{
	  target[0] = retval->x;
	  target[1] = retval->y;
	  target[2] = retval->z;
	}
    }
  else
    {
      int offset = address - 1;
      target[0] = clist[offset].x;
      target[1] = clist[offset].y;
      target[2] = clist[offset].z;      
    }
  return found;
}


cache_node* EIOMeshAgent::
search_node(const int address)
{
  cache_node entry;
  cache_node *retval;

  entry.tag = address;
  retval = (cache_node *) bsearch((void *)&entry, (void *)clist, 
				  nodeCount, 
				  sizeof(cache_node),
				  nodecomp);
  return retval;
}

