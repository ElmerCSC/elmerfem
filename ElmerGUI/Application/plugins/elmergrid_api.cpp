#include <iostream>
#include "elmergrid_api.h"
#include "egmain.h"


using namespace std;

ElmergridAPI::ElmergridAPI()
{
  cout << "Constructing ElmergridAPI... done" << endl;
  cout.flush();
}


ElmergridAPI::~ElmergridAPI()
{
  cout << "Destructing ElmergridAPI... done" << endl;
  cout.flush();
}
 
int ElmergridAPI::loadElmerMeshStructure(const char *filename)
{
  int retval = eg_loadmesh(filename); 
  return retval;
}


int ElmergridAPI::createElmerMeshStructure(mesh_t *mesh,const char *options)
{
#if 1
  int retval = eg_transfermesh(mesh,options);
  return retval;
#else
    // nodes:
  mesh->nodes = 4;
  mesh->node = new node_t[4];
 node_t *n;

  n = &mesh->node[0];
  n->x[0] = 0.0;
  n->x[1] = 0.0;
  n->x[2] = 0.0;
  n->index = -1;
  
  n = &mesh->node[1];
  n->x[0] = 1.0;
  n->x[1] = 0.0;
  n->x[2] = 0.0;
  n->index = -1;
  
  n = &mesh->node[2];
  n->x[0] = 0.0;
  n->x[1] = 1.0;
  n->x[2] = 0.0;
  n->index = -1;
  
  n = &mesh->node[3];
  n->x[0] = 0.0;
  n->x[1] = 0.0;
  n->x[2] = 1.0;
  n->index = -1;

  // elements:
  mesh->elements = 1;
  mesh->element = new element_t[1];
  
  element_t *e;

  e = &mesh->element[0];
  e->code = 504;
  e->nodes = 4;
  e->node = new int[4];
  e->node[0] = 0;
  e->node[1] = 1;
  e->node[2] = 2;
  e->node[3] = 3;
  e->index = 1;
  
  // boundaryelements:
  mesh->boundaryelements = 4;
  mesh->boundaryelement = new boundaryelement_t[4];

  boundaryelement_t *b;

  b = &mesh->boundaryelement[0];
  b->code = 303;
  b->nodes = 3;
  b->node = new int[3];
  b->node[0] = 0;
  b->node[1] = 1;
  b->node[2] = 2;
  b->edges = 3;
  b->edge = new int[3];
  b->edge[0] = -1;
  b->edge[1] = -1;
  b->edge[2] = -1;
  b->elements = 0;
  b->element = new int[2]; // reserve space for 2 (will be corrected)
  b->index = 1;

  b = &mesh->boundaryelement[1];
  b->code = 303;
  b->nodes = 3;
  b->node = new int[3];
  b->node[0] = 0;
  b->node[1] = 1;
  b->node[2] = 3;
  b->edges = 3;
  b->edge = new int[3];
  b->edge[0] = -1;
  b->edge[1] = -1;
  b->edge[2] = -1;
  b->elements = 0;
  b->element = new int[2]; // reserve space for 2 (will be corrected)
  b->index = 1;

  b = &mesh->boundaryelement[2];
  b->code = 303;
  b->nodes = 3;
  b->node = new int[3];
  b->node[0] = 0;
  b->node[1] = 2;
  b->node[2] = 3;
  b->edges = 3;
  b->edge = new int[3];
  b->edge[0] = -1;
  b->edge[1] = -1;
  b->edge[2] = -1; 
  b->elements = 0;
  b->element = new int[2]; // reserve space for 2 (will be corrected)
 b->index = 1;

  b = &mesh->boundaryelement[3];
  b->code = 303;
  b->nodes = 3;
  b->node = new int[3];
  b->node[0] = 1;
  b->node[1] = 2;
  b->node[2] = 3;
  b->edges = 3;
  b->edge = new int[3];
  b->edge[0] = -1;
  b->edge[1] = -1;
  b->edge[2] = -1;
  b->elements = 0;
  b->element = new int[2]; // reserve space for 2 (will be corrected)
  b->index = 1;

  cout << "ok, palautan nyt uuden verkon" << endl;
  cout.flush();
#endif
}
