#ifndef ELMERGRID_API_H
#define ELMERGRID_API_H

#include "src/meshtype.h"

class ElmergridAPI
{
 public:
  ElmergridAPI();
  ~ElmergridAPI();
  
  int loadElmerMeshStructure(const char*);
  int createElmerMeshStructure(mesh_t *mesh,const char *options);
};

#endif // #ifndef ELMERGRID_API_H
 
