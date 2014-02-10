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
Module:     ecif_body2D.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation

************************************************************************/

#include "ecif_body2D.h"
#include "ecif_bodyElement.h"
#include "ecif_bodyElement2D.h"
#include "ecif_bodyElementLoop.h"
#include "ecif_model.h"
#include "ecif_renderer.h"
#include "ecif_userinterface.h"


// Constructors.

Body2D::Body2D()
  : Body()
{
}

//--Body-elements etc are not known
Body2D::Body2D(bodyGmtrType body_type, int ext_tag, char* name, colorIndices color)
  : Body(body_type, ext_tag, name, color)
{
}

//--Component ids are known, but not more.
Body2D::Body2D(ecif_Body_X& trx_body, bool add_default_layer)
  : Body(trx_body, add_default_layer)
{
}


Body2D::Body2D(bodyGmtrType body_type, int int_tag, int ext_tag,
               int nof_mesh_elements, int* mesh_element_ids)
               :Body(body_type, int_tag, ext_tag, NULL)
{
  maxNofMeshElements = nof_mesh_elements;
  nofMeshElements = nof_mesh_elements;
  meshElementIds = mesh_element_ids;
}


bool
Body2D::acceptsStructuredMesh(int layer)
{
  if (!checkLayerIndex(layer)) return false;

  // If more than one loops, no grid!
  if ( nofElementLoops[layer] == 0 || nofElementLoops[layer] > 1 ) {
    return false;
  }

  // If more or less than four elements, n grid!
  if ( nofElements[layer] != 4 ) {
    return false;
  }

  // Ok, gridable
  return true;
}


// Add all pending vertices as elements
// Returns nof elements added
int
Body2D::addAllPendingVertices(int layer)
{
  if (!checkLayerIndex(layer)) return 0;

  int nof_vertex_tags = pendingVertexTags[layer]->size();

  if ( nof_vertex_tags == 0 )
    return 0;

  // Vertex ids are stored in this table
  // Allocate space also for the possible closing
  // vetex
  int* vertex_groups = new int[1 + nof_vertex_tags];
  int* vertex_ids = new int[1 + nof_vertex_tags];

  int nof_vertex_ids = 0;

  int i;

  for (i = 0; i < nof_vertex_tags; i++) {

    // Read and remove last element
    int v_group = (*pendingVertexGroups[layer])[i];
    int v_tag = (*pendingVertexTags[layer])[i];

    // Find vertex by tag and get id
    BodyElement* v = model->getVertexByTag(v_tag);

    if (v != NULL) {
      vertex_groups[nof_vertex_ids] = v_group;
      vertex_ids[nof_vertex_ids] = v->Id();
      nof_vertex_ids++;
    }
  }

  int nof_elements = 0;

  int groupd_id = vertex_groups[0];

  for (i = 1; i < nof_vertex_ids; i++) {

    int v1 = vertex_ids[i-1];
    int v2 = vertex_ids[i];

    // Inside group, create new element
    if ( groupd_id == vertex_groups[i] ) {

      BodyElement* be = new BodyElement2D(v1, v2);
      addElement(layer, be);
      model->addBodyElement(be, false);
      nof_elements++;
    }

    groupd_id = vertex_groups[i];
  }

  pendingVertexGroups[layer]->clear();
  pendingVertexTags[layer]->clear();

  return nof_elements;
}



#if 0
// Method checks that edges create a closed loop and it stores the
// edge-loop as table of consecutive edge-ids. A negative id-number
// in the  table means that edge is in the reversed order.
// Return true if body is ok, otherwise false.
bool
Body2D::check()
{
  checked = true;
  bool result = false;

  nofElements = belements->size();
  if (nofElements == 0)
    return false;
  // A table where row is: (edge-id, vertex1-id, vertex2-id)
  int (*vertex_table)[3] = new int[nofElements][3];
  BodyElement* be;
  BodyElement* vertex;
  int e_id, v1_id, v2_id;
  BodyElementTable::iterator pos = belements->begin();
  // Now the (edge,vertex1,vertex2)-table is filled.
  int row = 0;
  while (pos != belements->end()) {
    be = (*pos++).second;
    e_id = be->ID();
    vertex = be->getSubElement(1);
    v1_id = vertex->ID();
    vertex = be->getSubElement(2);
    v2_id = vertex->ID();
    vertex_table[row][0] = e_id;
    vertex_table[row][1] = v1_id;
    vertex_table[row][2] = v2_id;
    row++;
  }
  // Next we create a edgeLoop-table which is simply a list of edge-id
  // numbers. A negative id-number in the table means that edge's
  // vertices are in opposite order compared to whole loop order.
  IdList* edgeLoop = new IdList;

  // Starting edge is selcted.
  int table_row = 0;//we take the first edge as a staring edge.
  int vrtx_nbr = 1;
  // These are needed to check at the end that loop is closed.
  int start_vertex = vertex_table[0][1];
  int end_vertex = vertex_table[0][2];
  // OUTER LOOP: once for each edge in the resulting edge-loop table.
  for (int i = 0; i < nofElements; i++) {
    // order of vertices compared to whole loop order
    int sign = (vrtx_nbr == 1)? 1: -1;
    edgeLoop->push_back(sign * vertex_table[table_row][0]);
    // id for the other vertex in the edge's vertex pair.
    //  Next edge is
    int v_id = vertex_table[table_row][(vrtx_nbr==1)?2:1];
    end_vertex = v_id;
    // we don't want to search for this table-row any more, lets mark it off!
    vertex_table[table_row][0] *= -1;
    // now next row in *vertex_table* is searhced.
    // Note: updating of *vrtx_nbr* is based on short-circuit property
    // of the logical operators in C++!!
    // INNER LOOP: we search for the 'next' edge in the (edge,vrtx1,vrtx2)-table
    for (int j = 0; j < nofElements; j++) {
      vrtx_nbr = 1;
      if ( (vertex_table[j][0] > 0)
           &&
          // here short-circuit is used in logical-or (||) when
          // vtrx_nbr is selcted based on first matching column.
          (vertex_table[j][vrtx_nbr] == v_id ||
           vertex_table[j][++vrtx_nbr] == v_id )) {
        table_row = j;
        break;
      }
    }
  }

  if (start_vertex == end_vertex)
    result = true;
  else
    return false;

  // If we have a closed loop of edges, we can create the final
  // element-loop to be stored in the Model.
  // We also ccheck loops ccw/cw orientation
  int loopDirection;
  BodyElementLoop* bel;
  if (result == true) {
    bel = new BodyElementLoop(edgeLoop);
    model->addBodyElementLoop(bel);
    elementLoopIds->push_back(bel->ID());
    loopDirection = calcDirection();
  }
  // If for some reason edge-loop orientation couldn't
  // be calulated, status is not ok.
  if (loopDirection == 0)
    result = false;
  // If loop direction was clock-wise, we have to
  // change the loop-ids sign to indicate this
  else if (loopDirection < 0) {
    elementLoopIds->pop_back();
    elementLoopIds->push_back(-1 * bel->ID());
  }

  return result;
}
#endif


// New version checking all element-loops defined for the body
// Method checks that edges create a closed loop and it stores the
// edge-loop as table of consecutive edge-ids. A negative id-number
// in the  table means that edge is in the reversed order.
// Return true if body is ok, otherwise false.
bool
Body2D::check()
{
  if (nofLayers == 0) return false;

  initName();

  // For a mesh, open or virtual body
  // ================================
  if ( isMeshBody() ||
       isBemBody()  ||
       isOpen()     ||
       isVirtual()
     ) {
    return Body::check();
  }

  // For a closed cad body
  // =====================
  UserInterface* gui = (UserInterface*)model->getGui();


  // Check all layers
  // ----------------
  int layer = -1;
  while (true) {

    if (!selectLayer(++layer)) break;

    int nof_pending_vertices = addAllPendingVertices(layer);

    int nof_pending_elements = addAllPendingElements(layer);

    checked = true;

    nofElements[layer] = belements[layer]->size();

    if ( nofElements[layer] == 0 ) {

      // We have elementLoops constructed already,
      // Store element ids from the loop to the belements-array
      if ( nofElementLoops[layer] > 0 ) {

        if ( !addAllLoopElements(layer) ) {
          status = false;
          return false;
        }

        nofElements[layer] = belements[layer]->size();

        // All input loops were checked (emf-file!), the body
        // shoud be ok
        if ( nof_pending_elements == 0 && status == true ) {
          continue;
        }

      // We have neither elements nor elementLoops,
      // something is really wrong!
      } else {

        strstream strm1, strm2;
        strm1 << "***ERROR in geometry for Body " << tag << ":" << ends;
        strm2 << "---No Edges or Edge Loops defined!" << ends;
        gui->showMsg(strm1.str());
        gui->showMsg(strm2.str(), 1);

        status = false;
        return false;
      }
    }

    // Construct loops from "scratch"
    // ==============================

    //--Remove existing loops
    removeElementLoops(layer);
    status = false;

    //--Create a table where rows are: (edge-id, vertex1-id, vertex2-id)
    int (*vertex_table)[3] = new int[nofElements[layer]][3];

    BodyElement* edge;
    BodyElement* vertex;
    int e_id, v1_id, v2_id;

    BodyElementTable::iterator pos = belements[layer]->begin();
    int row = 0;

    //--Fill the table.
    while ( pos != belements[layer]->end() ) {

      edge = (*pos++).second;

      e_id = edge->Id();

      // If no vertices in the element or we have a closed element like
      // a  circle
      if ( edge->isClosedU() || edge->getNofSubElements() == 0 ) {
        v1_id = 0;
        v2_id = 0;

      } else {
        v1_id = edge->getFirstSubElement()->Id();
        v2_id = edge->getLastSubElement()->Id();
      }

      vertex_table[row][0] = e_id;
      vertex_table[row][1] = v1_id;
      vertex_table[row][2] = v2_id;
      row++;
    }

    // Next we create a edgeLoop-table which is simply a list of edge-id
    // numbers. A negative id-number in the table means that edge's
    // vertices are in opposite order compared to whole loop order.

    // Starting edge for the loop is selected.
    int nof_free_elements = nofElements[layer];

    while (nof_free_elements > 0) {

      // OUTERMOST For: once for each edge in the resulting edge-loop table.
      for (int k = 0; k < nofElements[layer]; k++) {

        if (vertex_table[k][0] < 0)
          continue;

        IdList edge_loop;

        int table_row = k;//we take the first free edge as the starting edge.
        int vrtx_nbr = 1;

        // These are needed to finally check that the element loop is closed.
        int start_vertex, end_vertex;
        start_vertex = vertex_table[k][1];
        end_vertex = -1;

        while (nof_free_elements > 0 && start_vertex != end_vertex) {

          // OUTER For: once for each edge in the resulting edge-loop table.
          for (int i = 0; i < nofElements[layer]; i++) {

            if ( vertex_table[i][0] < 0 ) continue;

            // Order of vertices compared to whole loop order
            int sign = (vrtx_nbr == 1)? 1: -1;

            // Add new edge
            edge_loop.push_back(sign * vertex_table[table_row][0]);
            nof_free_elements--;

            // We don't want to search for this table-row any more, lets mark it off!
            vertex_table[table_row][0] *= -1;

            // Id for the other vertex in the edge's vertex pair.
            int v_id = vertex_table[table_row][(vrtx_nbr==1)?2:1];
            end_vertex = v_id;

            if (start_vertex == end_vertex) break;

            // Now next row in *vertex_table* is searhced.
            // Note: updating of *vrtx_nbr* is based on short-circuit property
            // of the logical operators in C!!

            // INNER For: we search for the 'next' edge in the (edge,vrtx1,vrtx2)-table
            for (int j = 0; j < nofElements[layer] ; j++) {

              if (vertex_table[j][0] < 0) continue;

              vrtx_nbr = 1;

              // NOTE: Here short-circuit is used in logical-or (||) when
              // vtrx_nbr is selcted based on first matching column.
              if (vertex_table[j][vrtx_nbr] == v_id ||
                  vertex_table[j][++vrtx_nbr] == v_id ) {
                table_row = j;
                break;
              }

            } // End INNER For

            if (nof_free_elements <= 0 ) break;

          } // End OUTER For

        } // End OUTER While

        // We could not close the loop, error!
        if (start_vertex != end_vertex) {

          strstream strm;
          strm << "***ERROR in geometry for Body " << tag << " in layer "
               << layer + 1 << ends;
          gui->showMsg(strm.str());
          gui->showMsg("---Body not closed!", 1);

          return false;
        }

        // All but the first loops are cw-directed
        int loop_direction = ( nofElementLoops[layer] > 0 )?-1:1;

        BodyElementLoop* bel = NULL;
        int bel_id;
        int direction; // Direction relative to the possible existing loop

        //-Existing loop, just pick id and relative direction
        //-NOTE All loops are stored in ccw-direction in the model, so
        // we do not use the 'direction' argument here
        //
        if ( model->getBodyElementLoopId(&edge_loop, bel_id, direction) ) {
          bel = model->getBodyElementLoopById(bel_id);
          if ( bel == NULL ) return false;

        //-Otherwise create new loop and add to the model
        //
        } else {

          bel = new BodyElementLoop(&edge_loop, false, OT_BOUNDARY);
          bel->check();
          model->addBodyElementLoop(bel);
        }

        // Add loop id to the body
        elementLoopIds[layer]->push_back(loop_direction * bel->Id());
        elementLoopTags[layer]->push_back(bel->Tag());
        nofElementLoops[layer]++;

        // No more free edges
        if (nof_free_elements <= 0) break;

      } // End OUTERMOST For

    } // End OUTERMOST While

  } // Each body grid layer

  status = true;
  return true;
}


