/*  
   ElmerFront - A graphical user interface of Elmer software
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Front
Module:     ecif_bodyLayer.cpp
Language:   C++
Date:       01.01.03
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_bodyLayer.h"
#include "ecif_body.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_model.h"
#include "ecif_parameter.h"


//Initialize static class variables.
int BodyLayer::last_tag = 0;


// BodyLayer class
// ===================

BodyLayer::BodyLayer()
{
  tag = ++last_tag;
  init();
  initName();
}


BodyLayer::BodyLayer(ecif_BodyLayer_X& tx)
{
  enum bodyLayerType tp;

  if ( tx.tag == NO_INDEX ) {
    tag = ++last_tag;
    tp = IMPLICIT_LAYER;
  } else {
    tag = tx.tag;
    tp = EXPLICIT_LAYER;
  }

  if (last_tag < tag) {
    last_tag = tag;
  }

  init();

  // NOTE: This must be set after init
  type = tp;

  int i;

  if ( tx.is_open ) {
    tplgType = OPEN_LAYER;
  }

  if ( tx.has_color ) {
    colorIndex = ef_nodefault;
    for (i = 0; i < 4; i++)
      color[i] = tx.color[i];

  } else {
    colorIndex = DEFAULT_COLOR_INDEX;
  }

  bodyId = tx.body_id;
  bodyTag = tx.body_tag;

  update_dyna_string(name, tx.name);

  int nof_ids;
  int* gids;
  int* mids;
  int* xids;

  // Possible grid parameter ids and mesh indices for layer-0
  nof_ids = tx.nof_grid_param_ids;
  gids = tx.grid_param_ids;
  mids = tx.grid_param_mesh_indices;
  setGridParameterData(nof_ids, gids, mids);

  // Possible excluded mesh indices for layer-0
  nof_ids = tx.nof_excluded_meshes;
  xids = tx.excluded_mesh_indices;
  setExcludedMeshData(nof_ids, xids);
}


BodyLayer::~BodyLayer()
{
  delete[] gridParameterIds;
  delete[] gridParameterMeshIndices;
  delete[] excludedMeshIndices;
}

// Correct this !!!###!!!
bool
BodyLayer::acceptsStructuredMesh()
{
  return false;
}


const Body*
BodyLayer::getBody()
{
  return model->getBodyById(bodyId);
}


void
BodyLayer::getColor(Color4& clr)
{
  for (int i = 0; i < 4; i++)
    clr[i] = color[i];
}


int
BodyLayer::getGridParameterId(int mesh_index)
{
  for (int i = 0; i < nofGridParameterIds; i++) {
    if ( gridParameterMeshIndices[i] == mesh_index ) {
      return gridParameterIds[i];
    }
  }

  return NO_INDEX;
}


bool
BodyLayer::getMeshDensityValue(int mesh_index, char& dtype, double& dvalue)
{
  if ( mesh_index < 0 || mesh_index >= model->getNofMeshes() )
    return false;

  dtype = ' ';
  dvalue = 0;

  // Check if parameter is defined for the mesh
  for (int i = 0; i < nofGridParameterIds; i++) {

    // Mesh found
    if ( mesh_index == gridParameterMeshIndices[i] ) {

      Parameter* param = model->getParameterById(ECIF_GRID_PARAMETER, gridParameterIds[i]);

      if ( param == NULL )
        return false;

      // If value type given
      char type[2] = " ";
      if ( param->getFieldValueBySifName("Mesh Density Type", 1, (char*)type) ) {

        double value;

        // MeshH
        if ( type[0] == 'H' && param->getFieldValueBySifName("Mesh H", value) ) {
          dvalue = value;
          dtype = 'H';
          return true;

        // Mesh R
        } else if ( type[0] == 'R' && param->getFieldValueBySifName("Mesh R", value) ) {
          dvalue = value;
          dtype = 'R';
          return true;
        }

      }

      break;
    }

  } // for all parameters

  return false;
}


// Correct this!!!###!!!
int
BodyLayer::getMeshQuadGridN(int mesh_index, int element_id)
{
  if ( !acceptsStructuredMesh() ) return 0;

  if ( mesh_index < 0 || mesh_index >= model->getNofMeshes() )
    return 0;

  Parameter* param = NULL;

  for (int i = 0; i < nofGridParameterIds; i++) {

    if ( gridParameterMeshIndices[i] == mesh_index ) {
      param = model->getParameterById(ECIF_GRID_PARAMETER, gridParameterIds[i]);
      break;
    }
  }

  if ( param == NULL ) {
    return 0;
  }

  char value_buffer[1 + 128];

  // Check if Quadrilateral mesh defined for the grid layer
  if ( !( param->getFieldValueBySifName("Mesh Element Type", 128, value_buffer) &&
          LibFront::ncEqual(value_buffer, "Quad")
        )
     ) {
    return 0;
  }


  Body* body = model->getBodyById(bodyId);

  if (body == NULL) return 0;

  // Check element index
  int index = NO_INDEX;

  int layer = 0; // Correct this !!!###!!!

  int be_index = 0;
  while (true) {
    BodyElement* be = body->getElement(layer, be_index++);
    if (be==NULL) break;
    if ( be->Id() == element_id ) {
      index++;
      break;
    }
  }

  if ( index == NO_INDEX ) {
    return 0;
  }

  int n;
  bool found = false;

  if ( index == 0 || index == 2 ) {
    found = param->getFieldValueBySifName("Mesh Quadgrid N1", n);
  } else if ( index == 1 || index == 3 ) {
    found = param->getFieldValueBySifName("Mesh Quadgrid N2", n);
  } else {
    return 0;
  }

  if (found) {
    return n;
  }

  return 0;

}


// Total number of mif-file layers which must be created
// from one body layer
//
int
BodyLayer::getNofMifLayers(const IdList* elem_loop_ids)
{
  int count = 0;

  // Nof geometries in the first element in the first
  // loop defines number of mif layers
  //
  IdList::iterator itr = ((IdList*)elem_loop_ids)->begin();

  int did = *itr;
  int direction = (did < 0 )?-1:1;
  BodyElementLoop* bel = model->getBodyElementLoopById(direction * did);

  if ( bel != NULL ) {
    BodyElement* be = bel->getElement(0);
    if ( be != NULL ) {
      count += be->getNofMifGeometries();
    }
  }

  return count;
}


// Total number of mif-file layer loops which must be created
// for one mif-layer
// NOTE: This is for 'normal' cases when we have a single layer
//
int
BodyLayer::getNofMifLayerLoops(const IdList* elem_loop_ids)
{
  int count = 0;

  // The sum of Nof geometries in the first elements in each bodyelement loop
  // defines the number of mif-layers loops
  //

  int nof_bels = elem_loop_ids->size();
  IdList::iterator itr = ((IdList*)elem_loop_ids)->begin();

  for (int i = 0; i < nof_bels; i++, itr++) {

    int did = *itr;
    int direction = (did < 0 )?-1:1;
    BodyElementLoop* bel = model->getBodyElementLoopById(direction * did);

    if ( bel != NULL ) {
      BodyElement* be = bel->getElement(0);
      if ( be != NULL ) {
        count += be->getNofMifGeometries();
      }
    }
  }

  return count;
}


// Total number of mif-file layer loops which must be created
// for one mif-layer
// NOTE: This is for situations when we have a multigeometry
// element (as the first element in the first layer) creating
// multiple mif-layers.
// index: the multi-geometry (mif-layer index) index
//
int
BodyLayer::getNofMifLayerLoops(int gmtr_index, const IdList* elem_loop_ids)
{
  int count = 0;

  // Nof geometries in the first element in the first
  // loop defines number of mif layers
  //
  int nof_bels = elem_loop_ids->size();
  IdList::iterator itr = ((IdList*)elem_loop_ids)->begin();

  for (int i = 0; i < nof_bels; i++, itr++) {

    int did = *itr;
    int direction = (did < 0 )?-1:1;
    BodyElementLoop* bel = model->getBodyElementLoopById(direction * did);

    if ( bel != NULL ) {
      BodyElement* be = bel->getElement(0);
      if ( be != NULL ) {
        if ( NO_INDEX != be->getMifGeometryTag(gmtr_index) ) {
          count += 1;
        }
      }
    }
  }

  return count;
}


bool
BodyLayer::hasBody(int bd_id)
{
  if ( bodyId == bd_id ) {
   return true;
  }

  return false;
}


void
BodyLayer::initClass(Model* mdl)
{
  BodyLayer::last_tag = 0;
}


void
BodyLayer::init()
{
  int i;

  model->addModelObject(this, OT_BODY_LAYER);

  type = NONE_LAYER;

  tplgType = CLOSED_LAYER;

  name = NULL;

  bodyId = NO_INDEX;
  bodyTag = NO_INDEX;

  colorIndex = DEFAULT_COLOR_INDEX;
  for (i = 0; i< 4; i++) {
    color[i] = colorValues[DEFAULT_COLOR_INDEX][i];
  }

  nofExcludedMeshes = 0;
  excludedMeshIndices = NULL;

  nofGridParameterIds = 0;
  gridParameterMeshIndices = NULL;
  gridParameterIds = NULL;
}


void
BodyLayer::initName()
{
  if ( name == NULL || name[0] == '\0' ) {
    strstream strm;
    strm << "Layer" << tag << ends;

    update_dyna_string(name, strm.str());
  }
}


// If body layer excluded from the currently drawn mesh
bool
BodyLayer::isExcludedFromMesh(int mesh_index)
{
  if ( mesh_index == NO_INDEX ) {
    return false;
  }

  for (int i = 0; i < nofExcludedMeshes; i++) {
    if ( excludedMeshIndices[i] == mesh_index ) {
      return true;
      break;
    }
  }

  return false;
}


// Layer Front model file output
//
ostream&
BodyLayer::output_emf(ostream& out, short indent_size, short indent_level)
{
  char* QM = "\"";
  short is = indent_size;
  short il = indent_level;

  int nof_ids;
  int* gids;
  int* mids;

  //--Excluded mesh indices
  nof_ids = nofExcludedMeshes;
  mids = excludedMeshIndices;

  if ( nof_ids > 0 ) {
    LibFront::output_vector(out, is, il, "Excluded Mesh Indices", NULL, nof_ids, mids, false);
  }

  //--GridParameter ids
  nof_ids = nofGridParameterIds;
  gids = gridParameterIds;
  mids = gridParameterMeshIndices;

  if ( nof_ids > 0 ) {
    LibFront::output_vector(out, is, il, "Grid Parameter Ids", NULL, nof_ids, gids, false);
    LibFront::output_vector(out, is, il, "Grid Parameter Mesh Indices", NULL, nof_ids, mids, false);
  }

  return out;
}


// Layer Mesh2D file output
//
ostream&
BodyLayer::output_mif(ostream& out, int& next_mif_id, const IdList* elem_loop_ids)
{
  char* QM = "\"";

  int i;
  char value_buffer[1 + 128];

  int mesh_index = model->getCurrentMeshIndex();
  int pid = getGridParameterId(mesh_index);
  Parameter* grid_param = model->getParameterById(ECIF_GRID_PARAMETER, pid);

  char dtype;
  double dvalue;
  bool dgiven = getMeshDensityValue(mesh_index, dtype, dvalue);

  bool is_open;
  bool is_quadGrid;

  if ( tplgType == OPEN_LAYER ) {
    is_open = true;
  } else {
    is_open = false;
  }

  // Check if Quadrilateral elements
  if ( grid_param != NULL &&
       grid_param->getFieldValueBySifName("Mesh Element Type", 128, value_buffer) &&
       LibFront::ncEqual(value_buffer, "Quad")
     ) {
    is_quadGrid = true;
  } else {
    is_quadGrid = false;
  }

  int nof_mif_layers = getNofMifLayers(elem_loop_ids);

  for (int layer = 0; layer < nof_mif_layers; layer++) {

    // Layer id
    // --------
    indent(out, 4) << "LayerId: " << next_mif_id++ << endl;

    // Layer mesh density value
    // ------------------------
    if ( dgiven && dvalue > 0 )
      indent(out, 4) << dtype << ": " << dvalue << endl;

    // Layer type (meshing method)
    // ---------------------------
    indent(out, 6) << "LayerType: ";

    //-Structured grid
    if ( is_quadGrid ) {
      out << "QuadGrid";

      indent(out, 6) << "GridSize: ";
      grid_param->getFieldValueBySifName("Mesh Quadgrid N1", 128, value_buffer);
      out << value_buffer << " ";
      grid_param->getFieldValueBySifName("Mesh Quadgrid N2", 128, value_buffer);
      out << value_buffer << endl;

    //-Triangles
    } else if ( grid_param != NULL &&
                grid_param->getFieldValueBySifName("Mesh Layer Type", 128, value_buffer)
              ) {
      out << value_buffer;

    //-Open layer (1D body!)
    } else if ( is_open ) {
      out << "BoundaryMesh";

    //-Default
    } else {
      out << "MovingFront";
    }

    out << endl;

    // Background mesh stuff
    // ---------------------
    if ( !(is_open || is_quadGrid) ) {

      //-Fixed nodes
      indent(out, 6) << "FixedNodes: " << 0 << endl;  // NOTE: Currently always 0!, MVe 16.02.00

      //-Background mesh type
      indent(out, 6) << "BGMesh: ";

      //--Background mesh method specified
      if ( grid_param != NULL &&
           grid_param->getFieldValueBySifName("Mesh Bg Mesh", 128, value_buffer)
         ) {

        //-External bg-mesh file
        if ( LibFront::in(value_buffer, "External") ) {

          //-File name given
          if ( grid_param->getFieldValueBySifName("Mesh Bg Mesh File", 128, value_buffer) &&
               value_buffer != NULL && value_buffer[0] != '\0'
             ) {
            out << "External " << QM << value_buffer << QM;

          //-File name missing, use default method!
          } else {
            out << "Delaunay";
          }

        //-Some of the predefined methods selected by the user
        } else {
          out << value_buffer;
        }

      //--Use default method (no bg mesh)
      } else {
        out << "Delaunay";
      }

      out << endl;

      // Mesh seed (optional)
      // --------------------
      if ( grid_param != NULL &&
           grid_param->getFieldValueBySifName("Mesh Seed Type", 128, value_buffer)
         ) {

        // Only Implicit seed currently supported, MVe 16.02.99
        if ( LibFront::ncEqual(value_buffer, "Implicit") &&
             grid_param->getFieldValueBySifName("Mesh Seed Edge", 128, value_buffer)
           ) {
          indent(out,6) << "Seed: Implicit" << endl;
          indent(out,8) << "Edge: " << value_buffer << endl;
        }
      }
    }

    // Layer loops
    // -----------
    int nof_mif_loops = 0;

    if ( nof_mif_layers == 1 ) {
      nof_mif_loops = getNofMifLayerLoops(elem_loop_ids);
    } else {
      nof_mif_loops = getNofMifLayerLoops(layer, elem_loop_ids);
    }

    indent(out, 6) << "Loops: " << nof_mif_loops << endl;

    int next_mif_loop = 1;

    int nof_loops = ((IdList*)elem_loop_ids)->size();
    IdList::iterator itr = ((IdList*)elem_loop_ids)->begin();

    for (int loop = 0; loop < nof_loops; loop++, itr++) {

      int did = *itr;
      int direction = (did < 0 )?-1:1;

      BodyElementLoop* bel = model->getBodyElementLoopById(direction * did);

      if ( bel == NULL ) continue;

      int nof_bel_mif_loops = 0;

      if ( nof_mif_layers == 1 ) {
        nof_bel_mif_loops = bel->getNofMifLoops(NO_INDEX);
      } else {
        nof_bel_mif_loops = bel->getNofMifLoops(layer);
      }

      for (int idx = 0; idx < nof_bel_mif_loops; idx++) {

        // Print positions before and after printing fixed data
        // (for getting the indent size for ids)
        int pos1, pos2;

        indent(out, 8) << "LoopId: " << next_mif_loop++ << endl;
        indent(out, 10) << "Direction: " << direction << endl;

        pos1 = out.tellp();

#if 0
        indent(out, 10) << "Edges: " << bel->getNofElementMifTags(idx) << "  ";
        pos2 = out.tellp();
        bel->outputDirectedElementMifTags(out, pos2 - pos1, idx);
#endif

#if 1
        // Open layer are output in both directions:
        // Edges: 6  1 2 3 -3 -2 -1
        // Owing to Mesh2D disability to handle open bodies!!!
        //
        int factor = isOpen()?2:1;

        // Output 'Edges: nof-edges'
        if ( nof_mif_layers == 1 ) {
          indent(out, 10) << "Edges: " << factor * bel->getNofElementMifTags(idx) << "  ";
        } else {
          indent(out, 10) << "Edges: " << factor * bel->getNofElementMifTags(layer) << "  ";
        }

        // Pick current output position
        pos2 = out.tellp();

        // Output edge tags
        if ( nof_mif_layers == 1 ) {
          bel->outputDirectedElementMifTags(out, isOpen(), pos2 - pos1, idx);
        } else {
          bel->outputDirectedElementMifTags(out, isOpen(), pos2 - pos1, layer);
        }
#endif

        out << endl;
      }
    }

  }

  return out;
}


void
BodyLayer::setBodyId(int body_id)
{
  bodyId = body_id;
}


void
BodyLayer::setBodyTag(int body_tag)
{
  bodyTag = body_tag;
}


void
BodyLayer::setColorIndex(colorIndices color_index)
{
  // Set color index
  colorIndex = color_index;
  // Set actual color
  for (int i = 0; i< 4; i++) {
    color[i] = colorValues[color_index][i];
  }
}


void
BodyLayer::setExcludedMeshData(int nof_indices, int* excluded_mesh_indices)
{
  nofExcludedMeshes = nof_indices;

  delete[] excludedMeshIndices;

  excludedMeshIndices = NULL;

  if ( nof_indices == 0 )
    return;

  excludedMeshIndices = new int[nof_indices];

  for (int i = 0; i < nof_indices; i++) {
    excludedMeshIndices[i] = excluded_mesh_indices[i];
  }
}


void
BodyLayer::setGridParameterData(int nof_ids, int* gids, int* mesh_indices)
{
  nofGridParameterIds = nof_ids;

  delete[] gridParameterIds;
  delete[] gridParameterMeshIndices;

  gridParameterIds = NULL;
  gridParameterMeshIndices = NULL;

  if ( nof_ids == 0 )
    return;

  gridParameterIds = new int[nof_ids];
  gridParameterMeshIndices = new int[nof_ids];

  for (int i = 0; i < nof_ids; i++) {
    gridParameterIds[i] = gids[i];
    gridParameterMeshIndices[i] = mesh_indices[i];
  }
}


void
BodyLayer::setGridParameterIds(int nof_ids, int* gids)
{
  nofGridParameterIds = nof_ids;

  delete[] gridParameterIds;

  gridParameterIds = NULL;

  if ( nof_ids == 0 )
    return;

  gridParameterIds = new int[nof_ids];

  for (int i = 0; i < nof_ids; i++) {
    gridParameterIds[i] = gids[i];
  }
}


void
BodyLayer::setGridParameterMeshIndices(int nof_ids, int* mesh_indices)
{
  delete[] gridParameterMeshIndices;

  gridParameterMeshIndices = NULL;

  if ( nof_ids == 0 )
    return;

  gridParameterMeshIndices = new int[nof_ids];

  for (int i = 0; i < nof_ids; i++) {
    gridParameterMeshIndices[i] = mesh_indices[i];
  }
}
