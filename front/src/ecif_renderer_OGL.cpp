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
Module:     ecif_renderer_ogl.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation
  NOTE: Platform specific (WIN32, UNIX) parts are implemented
  in WIN32/Unix specific *.hpp files which are included later
  in this file depending on the macro value WIN32.

************************************************************************/

#include "ecif_body.h"
#include "ecif_bodyElement.h"
#include "ecif_control.h"
#include "ecif_geometry.h"
#include "ecif_mesh.h"
#include "ecif_model.h"
#include "ecif_modelObject.h"
#include "ecif_renderer_OGL.h"
#include "ecif_userinterface.h"
#include "ecif_timer.h"

static char oglWinClass [] = "ModelWindow";

GLParam::GLParam()
{
  lineWidthCoordAxis = 1.0f;
  lineWidthThin = 0.0f;
  lineWidthNormal = 1.0f;
  lineWidthSelected = 3.0f;
};

void Callback
Renderer_OGL::errorCallback(GLenum errorCode)
{
  const GLubyte *estring;
  estring = gluErrorString(errorCode);
  strstream strm;
  strm << "OpenGL Error: \n" << estring << ends;
  theControlCenter->getUI()->showMsg(strm.str());
}


// Constructor
Renderer_OGL::Renderer_OGL(Hinst app_instance, ecif_modelDimension dim)
{
  appInstance = app_instance;

  if (dim == ECIF_2D) {
    is2D = true;
    is2DSimulation = true;

  } else {
    is2D = false;
    is2DSimulation = false;
  }

  useOrthoProjection = true;

  visible = false;

  createData();

  // Init environment
  attachRenderer();
  setRendererInfo();
  init();

  // Displaylist for printing digits
  makeRasterFont();
}


Renderer_OGL::~Renderer_OGL()
{
  destroyWindow(rendererDisplay, rendererWindow);
  deleteData();
}


// Method attaches the renderer object with a window.
void
Renderer_OGL::attachRenderer()
{
  if (status == HAS_NO_WINDOW){

    // Create GL-window.
    WindowInfo winfo;

    createGLWindow(appInstance, "Model", 0, 0, winX, winY, winfo);

    rendererDisplay = winfo.display;
    rendererWindow = winfo.window;

    status = HAS_WINDOW;
  }
}


void
Renderer_OGL::bodySelectionHits()
{
  if ( selectionHits == 0 ) {
    return;
  }

  int picked_body_id = NO_INDEX;
  int picked_lr_id = NO_INDEX;

  GLuint tot_min_z = 0xffffffff;
  GLuint tot_max_z = 0;

  int record_start = 0;

  for (int i = 0; i < selectionHits; i++) {

    int nof_names = selectionBuffer[record_start];

    if (nof_names == 0) {
      break;
    }

    GLuint min_z   = selectionBuffer[record_start + 1];
    GLuint max_z   = selectionBuffer[record_start + 2];
    int body_id = selectionBuffer[record_start + 3];
    int lr_id = selectionBuffer[record_start + 4];

    // NOTE: 3 covers: nofNames, minZ, maxZ  and
    // nof_names covers all ids found for pick
    //
    record_start += 3 + nof_names;

    if ( min_z < tot_min_z ) {
      tot_min_z = min_z;
      picked_body_id = body_id;
      picked_lr_id = lr_id;
    }
  }

  if ( picked_body_id != NO_INDEX ) {
    model->bodySelected(this, picked_body_id, picked_lr_id);
  }

  selectionHits = 0;
}


void
Renderer_OGL::boundarySelectionHits()
{
  if ( selectionHits == 0) {
    return;
  }

  int max_nof_names = -1;

  int picked_bd1_id = NO_INDEX;
  int picked_bd2_id = NO_INDEX;
  int picked_lr1_id = NO_INDEX;
  int picked_lr2_id = NO_INDEX;
  int picked_bndr_id  = NO_INDEX;
  bool vertex_picked = false;

  GLuint tot_min_z = 0xffffffff;
  GLuint tot_max_z = 0;
  int record_start = 0;

  for (int i = 0; i < selectionHits; i++) {

    int nof_names = selectionBuffer[record_start];

    if (nof_names == 0) {
      break;
    }

    bool higer_found = false;

    if ( is2D && nof_names == 4 ) {
      vertex_picked = true;

    } else if ( nof_names == 5 ) {
      vertex_picked = true;
    }

    GLuint min_z = selectionBuffer[record_start + 1];
    GLuint max_z = selectionBuffer[record_start + 2];

    if ( max_z > tot_max_z ) {
      tot_max_z = max_z;
      higer_found = true;
    }

    // First parent
    //
    if ( picked_bd1_id == NO_INDEX || higer_found ) {

      picked_bd1_id = selectionBuffer[record_start + 3];
      picked_lr1_id = selectionBuffer[record_start + 4];

      if ( nof_names > 2 ) {
        picked_bndr_id = selectionBuffer[record_start + 2 + nof_names];
      }

    }

    // Second parent
    //
    if ( nof_names > 2 &&
         picked_bd1_id != selectionBuffer[record_start + 3] &&
         picked_bndr_id == selectionBuffer[record_start + 2 + nof_names] &&
        picked_bd2_id == NO_INDEX
      ) {
      picked_bd2_id = selectionBuffer[record_start + 3];
      picked_lr2_id = selectionBuffer[record_start + 4];
    }

    // If a deeper level name encountered
    if ( nof_names > max_nof_names ) {
      max_nof_names = nof_names;
      picked_bndr_id = selectionBuffer[record_start + 2 + nof_names];
    }

    // NOTE: 3 covers: nofNames, minZ, maxZ  and
    // nof_names covers all ids found for pick
    //
    record_start += 3 + nof_names;

  }

  // If picked found
  if ( picked_bndr_id != NO_INDEX ) {

    if ( vertex_picked ) {
      //model->setCurrentAnchorPoint(picked_bndr_id);
    }

    model->boundarySelected(this, picked_bndr_id,
                            picked_bd1_id, picked_lr1_id,
                            picked_bd2_id, picked_lr2_id);
  }

  selectionHits = 0;
}


void
Renderer_OGL::boundaryEdgeSelectionHits()
{
  // Nothing so far!
  return;
}


void
Renderer_OGL::boundaryVertexSelectionHits()
{
  // Nothing so far!
  return;
}


// Clean buffers
void
Renderer_OGL::clean()
{
  // Background color to white. Clear buffers.
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


// Clear screen
void
Renderer_OGL::clear()
{
  clean();
  draw();
}


// Constructor
void
Renderer_OGL::createData()
{

  displayListIds = new IdList;

  tesselator = gluNewTess();
  tesselatorPoints = new Point3List;

  gluTessCallback(tesselator, (GLenum)GLU_BEGIN, (GLvoid (Callback*)()) &glBegin);
  gluTessCallback(tesselator, (GLenum)GLU_VERTEX, (GLvoid (Callback*)()) &glVertex3dv);
  gluTessCallback(tesselator, (GLenum)GLU_END, &glEnd);
  gluTessCallback(tesselator, (GLenum)GLU_ERROR, (GLvoid (Callback*) ()) &errorCallback);

  pointBuffer = new double*[MAX_NOF_NODES];

  selectionBufferSize = 1024;
  selectionBuffer = new GLuint[selectionBufferSize];
  selectionHits = 0;

  Color4 current_color = {0.0f, 0.0f, 0.0f, 1.0f};  // black
  Color4 mesh_color    = {0.0f, 0.0f, 0.0f, 1.0f};  // black
  Color4 select_color  = {1.0f, 1.0f, 1.0f, 1.0f};  // white

  // Set initial values colors
  for (int i = 0; i < 4; i++) {
    currentColor[i] = current_color[i];
    meshColor[i]    = mesh_color[i];
    selectColor[i]  = select_color[i];
    //selectColor[i] = float(colorValues[0][i]) / MAX_NOF_COLOR_LEVELS;
  }

  inBoxDrawMode = false;
  inLineDrawMode = false;
  inMeshEditingMode = false;
  inMeshPickingMode = false;
  inVectorDrawMode = true;

  mkState = 0;

  // Initial window size
  winX   = 450;
  winY   = 400;
  winZ   = 2; // (-1,1);

  // Initially no rotation priorities, so rotation
  // by the mouse movements in the screen is for X,Y axis
  rotateAxisX = false;
  rotateAxisY = false;
  rotateAxisZ = false;

  // Timers
  doubleClickTimer = new Timer;
  mouseMoveTimer = new Timer;

  doubleClickTimer->start();
  mouseMoveTimer->start();
}


void
Renderer_OGL::deleteData()
{
  delete[] pointBuffer;

  gluDeleteTess(tesselator);

  deleteDisplayLists();
  delete displayListIds; // container itself

  deleteTesselatorPoints();
  delete tesselatorPoints; // container itself

  delete doubleClickTimer;
  delete mouseMoveTimer;
}


void
Renderer_OGL::deleteDisplayList(int list_id)
{
  glDeleteLists(list_id, 1);
}


void
Renderer_OGL::deleteDisplayLists()
{
  while ( !displayListIds->empty() ) {
    deleteDisplayList(displayListIds->back());
    displayListIds->pop_back();
  }
}


void
Renderer_OGL::deleteTesselatorPoints()
{
  while ( !tesselatorPoints->empty() ) {
    Point3* point = tesselatorPoints->back();
    delete[] point;
    tesselatorPoints->pop_back();
  }
}


// *** This function does the actual painting of the OpenGL-window.
void Callback
Renderer_OGL::display(Renderer* renderer)
{
#if 0
  static int counter = 0;
  strstream strm;
  strm << "Now doing display " << ++counter << endl << ends;
  theUI->showMsg(strm.str());
#endif
  renderer->clean();
  renderer->transform_scene();
  renderer->drawAllBodies();
  renderer->drawCoordinateAxis(0, 0);
  renderer->draw();
}


void
Renderer_OGL::displayRenderer()
{
  attachRenderer();
  status = SHOW;
  paintRenderer();
}


// **************************
// *** Model drawing loop ***
// **************************

// Draw bodies as filled surfaces or
// filled boundary mesh elements
void
Renderer_OGL::drawAllBodies()
{
  bool is_drawing_cad = false;
  bool is_drawing_mesh = false;

  flagName geometries[2];
  short nof_geometries = 0;

  if (model->getFlagValue(DRAW_SOURCE_CAD)) {
    geometries[nof_geometries] = DRAW_SOURCE_CAD;
    nof_geometries++;
    is_drawing_cad = true;
  }

  if (model->getFlagValue(DRAW_SOURCE_MESH)) {
    geometries[nof_geometries] = DRAW_SOURCE_MESH;
    nof_geometries++;
    is_drawing_mesh = true;
  }

  if ( nof_geometries == 0) return;

  // Possible selection hits
  if ( selectionHits > 0 ) {

    if ( model->getFlagValue(DRAW_TARGET_BODIES) ) {
        bodySelectionHits();

    } else if ( model->getFlagValue(DRAW_TARGET_BOUNDARIES) ) {
        boundarySelectionHits();

    } else if ( model->getFlagValue(DRAW_TARGET_BOUNDARY_EDGES) ) {
        boundaryEdgeSelectionHits();

    } else if ( model->getFlagValue(DRAW_TARGET_BOUNDARY_VERTICES) ) {
        boundaryVertexSelectionHits();
    }

    selectionHits = 0;
  }

  setDrawMode();

  bool draw_done = false;

  // Init possible mesh elements for renderering
  if (is_drawing_mesh) {
    model->resetMeshRendered();
    model->resetMeshEdgesSelected();
  }

  //---Draw all selected geometries
  short geometry = 0;
  while ( geometry < nof_geometries ) {

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //--Start body loop
    objectDrawingMode dmode;
    objectDrawingState dstate;

    int index = 0;
    while (true) {

      Body* body = model->getBody(index++);

      if (body==NULL) break;

      dmode = body->getDrawMode();
      dstate = body->getDrawState();

      if ( dmode == DM_HIDDEN) continue;

      // If we draw both Cad and Mesh, we draw Cad
      // as wireframe bodies to make mesh visible!
      if ( (nof_geometries == 2 && geometries[geometry] == DRAW_SOURCE_CAD) ||
            model->getFlagValue(DRAW_TARGET_BOUNDARIES)
         ) {
        dmode = DM_WIREFRAME;
      }

      //-Set body color etc.
      // Normal drawing mode
      if ( dstate == DM_NORMAL || model->getFlagValue(DRAW_TARGET_BOUNDARIES) ) {
        body->getColor(currentColor);

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, currentColor);
        //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, currentColor);
        //glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

      // Selected drawing mode
      } else {
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, selectColor);
      }

      //-Draw body
      //glEnable(GL_BLEND);
      body->draw(this, geometries[geometry], dmode);
      //glDisable(GL_BLEND);

    } //End body loop

    geometry++;
  } // End geometries loop


  //---Draw element labels
  // NOTE This must done last (not to be overwritten!)
  bool only_active = true;
  bool draw_sub_labels = true;

  int index = 0;
  while (true)  {

    BodyElement* be = model->getBoundary(index++, only_active);

    if (be==NULL) break;

    Body* bd1 = model->getBodyById(be->getParentId(1));
    Body* bd2 = model->getBodyById(be->getParentId(2));

    bool bd1_draw = true; // Parent body1 draw flag
    bool bd2_draw = true; // Parent body2 draw flag
    bool be_draw = true;  // Boundary draw flag

    if ( bd1 == NULL || DM_HIDDEN == bd1->getDrawMode() ) {
      bd1_draw = false;
    }

    if ( bd2 == NULL || DM_HIDDEN == bd2->getDrawMode() ) {
      bd2_draw = false;
    }

    // A selected element is not 'hidden'
    if ( DM_HIDDEN == be->getDrawMode() &&
         DS_SELECTED != be->getDrawState()
       ) {
      be_draw = false;
    }

    // Draw label (if element not hidden)
    if ( be_draw && (bd1_draw || bd2_draw) ) {
      be->drawLabel(this, draw_sub_labels);
    }
  }

  // Just a test!!!
  //model->drawCurrentPickVector();
}


// Draw bodies as bulk mesh elements
void
Renderer_OGL::drawAllMeshBodies()
{
  if ( !model->getFlagValue(GEOMETRY_TYPE_MESH) ) {
    return;
  }

  // Non-filled
  //glPolygonMode(GL_FRONT, GL_FILL);
  //glPolygonMode(GL_BACK, GL_FILL);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LINE_SMOOTH);
  glLineWidth(glParam.lineWidthNormal);

  //-Use black color
  Color4f meshColor = {0.0f, 0.0f, 0.0f, 1.0f};
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);


#if 0
  if ( glIsList(500) ) {
    glCallList(500);
    return;
  }

  glNewList(500, GL_COMPILE );
#endif

#if 0
  model->drawMesh();
  return;
#endif

  model->resetMeshRendered();

  //---Start body loop
  int index = 0;
  while (true) {

    Body* body = model->getBody(index++);

    if (body==NULL) break;

    // Set body color.
    //Color4f color;
    //body->getColor(color);
    //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color);
    //glColor4fv(color);
    if (DM_HIDDEN != body->getDrawMode()) {
      body->draw(this, DRAW_SOURCE_MESH, DM_NORMAL);
    }
  }

#if 0
  glEndList();
#endif
  //End body loop
  //glEnable(GL_LINE_SMOOTH);
}



void
Renderer_OGL::mouse_clck_action(int mk_state, int x1, int y1)
{
  if ( mk_state & MK_CONTROL )
    theControlCenter->getRenderer()->processSelection(MOUSE_CTRL_CLCK, mk_state, x1, y1);
  //else if ( mk_state & MK_SHIFT)
  //  theControlCenter->getRenderer()->processSelection(MOUSE_SHFT_CLCK, mk_state, x1, y1);
}


void
Renderer_OGL::mouse_dblclck_action(int mk_state, int x1, int y1)
{
  if ( mk_state & MK_CONTROL )
    theControlCenter->getRenderer()->processSelection(MOUSE_CTRL_DBL_CLCK, mk_state, x1, y1);
  else
    theControlCenter->getRenderer()->processSelection(MOUSE_DBL_CLCK, mk_state, x1, y1);
}


void
Renderer_OGL::mouse_move_action(int mk_state, int x1, int x2, int y1, int y2)
{
  static int curr_axis = -1;
  static int prev_pos_x = 0;
  static int prev_pos_y = 0;

  int amount;
  int axis;
  int code;
  int del_x = x2 - x1;
  int del_y = y2 - y1;
  int pos_x, pos_y;
  char* time_buffer;
  double time;
  int decimal, sign;

  switch (mk_state) {

  case MK_SHIFT:
#if 0
    // Time in milliseconds
    time = mouseMoveTimer->getLapTime();
    //time_buffer = _fcvt(time, 10, &decimal, &sign);
    //theUI->showMsg(time_buffer);
    if (time > 0.05) {
      mouseMoveTimer->stop();
      theControlCenter->getRenderer()->processSelection(MOUSE_SHFT_MOVE, mk_state, x1, y1);
      mouseMoveTimer->start();
    }
#endif
#if 1
    break; // Pure shift+move does not select any more, it is now shift+button+move
    theControlCenter->getRenderer()->processSelection(MOUSE_SHFT_MOVE, mk_state,  x1, y1);
#endif
    break;

  // Selection (shift+left down) or
  // Move (translate) (only left down)
  case MK_LBUTTON | MK_SHIFT:
  case MK_RBUTTON | MK_SHIFT:
  case MK_LBUTTON :

    if ( mk_state & MK_SHIFT ) {
      theControlCenter->getRenderer()->processSelection(MOUSE_SHFT_MOVE, mk_state, x1, y1);
    }
    else {
      theControlCenter->getRenderer()->translate(del_x, del_y, 0);
    }
    break;

  // Rotation (right button down)
  case MK_RBUTTON:

    // If some priority is set, we use it
    axis = -1;
    theControlCenter->getRenderer()->getRotateAxis(axis);

    if ( axis != -1 ) {
      curr_axis = axis;
    }

    pos_x = (del_x < 0)?-del_x:del_x;
    pos_y = (del_y < 0)?-del_y:del_y;

    // No axis set
    if ( axis == -1 ) {

      // Change axis if mouse direction is changed relevantly!
      if ( abs(pos_y - pos_x) > 1.05 * abs(prev_pos_y - prev_pos_x) ) {
        curr_axis = (pos_x >= pos_y)?1:0;
      }

      amount = (curr_axis == 1)?del_x:-del_y;

    // Axis set
    } else {
      amount = (curr_axis == 1)?del_x:-del_y;
    }

    theControlCenter->getRenderer()->rotate(curr_axis, amount);
    break;

  // Scale (both buttons down)
  case (MK_LBUTTON | MK_RBUTTON):
  case MK_MBUTTON:
    theControlCenter->getRenderer()->scale(del_y);
    break;

  default:
    break;
  }

  prev_pos_x = pos_x;
  prev_pos_y = pos_y;
}


void
Renderer_OGL::startDisplayList(int list_id)
{
  //glNewList(list_id, GL_COMPILE );
  glNewList(list_id, GL_COMPILE_AND_EXECUTE);
  storeDisplayListId(list_id);
}


void
Renderer_OGL::startDrawingCadBody()
{
  if (!is2D) {
    return;
  }

  gluBeginPolygon(tesselator);
}


void
Renderer_OGL::startDrawingCadBodyElementLoop(bool is_first_loop)
{
  if (!is2D) {
    return;
  }

  if (is_first_loop)
    gluNextContour(tesselator, (GLenum)GLU_EXTERIOR);
  else
    gluNextContour(tesselator, (GLenum)GLU_INTERIOR);
}


void
Renderer_OGL::startDrawingCadSurface()
{
  if (!is2D) {
    glDepthFunc(GL_LEQUAL);
  }
}



void
Renderer_OGL::startDrawingMeshSurface()
{
  if (!is2D) {
    glDepthFunc(GL_LEQUAL);
  }

  //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);
}


void
Renderer_OGL::startDrawingMeshSurfaceEdges()
{
  if (!is2D) {
    glDepthFunc(GL_LEQUAL);
  }

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);
}


void
Renderer_OGL::stopDisplayList(int list_id)
{
  glEndList();
}



void
Renderer_OGL::stopDrawingCadBody()
{
  if (!is2D) {
    return;
  }

  glFlush();

  gluEndPolygon(tesselator);
}


void
Renderer_OGL::stopDrawingCadBodyElementLoop()
{
  if (!is2D) {
    return;
  }
}



void
Renderer_OGL::stopDrawingCadSurface()
{
}


void
Renderer_OGL::stopDrawingMeshSurface()
{
  //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
}


void
Renderer_OGL::stopDrawingMeshSurfaceEdges()
{
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
}


void
Renderer_OGL::storeDisplayListId(int id)
{
  displayListIds->push_back(id);
}


void
Renderer_OGL::storeTesselatorPoint(Point3* point)
{
  tesselatorPoints->push_back(point);
}



void
Renderer_OGL::drawCoordinateAxis(int x_pos, int y_pos)
{
  if ( winX == 0 || winY == 0 )
    return;

  static double len = 0.15;       // length of the axis
  static double fac = 1.05;       // for the position of the labels
  static double pos = -1.0 + len; // Position from the SW-corner

  static double origo[3] = {0.0, 0.0, 0.0};
  static double axis_x[3] = {1.0, 0.0, 0.0};
  static double axis_y[3] = {0.0, 1.0, 0.0};
  static double axis_z[3] = {0.0, 0.0, 1.0};

  static float color_x[4] = {1.0f, 0.0f, 0.0f, 1.0f}; // Red
  static float color_y[4] = {0.0f, 1.0f, 0.0f, 1.0f}; // Green
  static float color_z[4] = {0.0f, 0.0f, 1.0f, 1.0f}; // Blue
  static float color_lbl[4] = {0.0f, 0.0f, 0.0f, 1.0f}; // Black

  static char label_x[21];
  static char label_y[21];
  static char label_z[21];

  model->getCoordinateLabels(20, label_x, label_y, label_z);

  // Correct y length by the window aspct ratio
  double len_y = len * (double(winX) / winY);

  axis_x[0] = len;
  axis_y[1] = len_y;
  axis_z[2] = len;

   // Save current transformation and projection matrix
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslated(pos, pos, 0.0);

  // Use model rotation matrix for coordinates
  rotate(0);

  glLineWidth(glParam.lineWidthCoordAxis);

  glDisable(GL_DEPTH_TEST);

  glBegin(GL_LINES);
    // x-axis
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color_x);
    glVertex3dv(origo);
    glVertex3dv(axis_x);
    // y-axis
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color_y);
    glVertex3dv(origo);
    glVertex3dv(axis_y);
    // z-axis
    if (!is2D) {
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color_z);
      glVertex3dv(origo);
      glVertex3dv(axis_z);
    }
  glEnd();

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color_lbl);

  glRasterPos3d(axis_x[0] * fac, 0.0, 0.0);
  printString(label_x);

  glRasterPos3d(0.0, axis_y[1] * fac, 0.0);
  printString(label_y);

  if (!is2DSimulation) {
    glRasterPos3d(0.0, 0.0, axis_z[2] * fac);
    printString(label_z);
  }

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
}


void
Renderer_OGL::drawElementLabel(char* lbl, Point3& lbl_p)
{
#if 0
  //-Test if label is invisible
  static GLfloat fb_buffer[4];

  glFeedbackBuffer(4, GL_3D, fb_buffer);

  glRenderMode(GL_FEEDBACK);

  glBegin(GL_POINTS);
    glVertex3d(lbl_p[0], lbl_p[1], lbl_p[2]);
  glEnd();

  int buffer_hits = glRenderMode(GL_RENDER);

  if ( buffer_hits == 0) {
    return;
  }
#endif

  //-Store state
  glPushAttrib(GL_LIGHTING_BIT);


  //-Use black color for the label
  static Color4f label_color = {0.0f, 0.0f, 0.0f, 1.0f};
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, label_color);

  //-Label position.
  glRasterPos3d(lbl_p[0], lbl_p[1], lbl_p[2]);

  //-Render label
  //glDisable(GL_DEPTH_TEST);
  printString(lbl);
  //glEnable(GL_DEPTH_TEST);

  //-Restore state
  glPopAttrib();
}


void
Renderer_OGL::drawLine(objectDrawingMode dmode, objectDrawingState dstate, short direction,
                       const Point3* start, const Point3* end, int elem_id)
{
  // We must take into account bodyelement's orientation within body!.
  // Here direction == 1 means ccw-order.
  // Calling functions must take care of this ***!!!***

  const Point3* p1;
  const Point3* p2;

  if (direction == 1) {
    p1 = start;
    p2 = end;

  } else {
    p2 = start;
    p1 = end;
  }

  // Just a line, no trimming
  if ( !is2D              ||
       dstate == DS_SELECTED  ||
       dmode == DM_WIREFRAME ||
       model->getFlagValue(DRAW_TARGET_EDGES) ||
       model->getFlagValue(DRAW_SOURCE_MESH)
     ) {

    if (dstate == DS_SELECTED)
      glLineWidth(glParam.lineWidthSelected);
    else
      glLineWidth(glParam.lineWidthNormal);

    if ( dmode == DM_INTRA_LAYER ) {
      glColor4d(1.0, 1.0, 1.0, 1.0);
    }

    glBegin(GL_LINES);
      glVertex3dv(*p1);
      glVertex3dv(*p2);
    glEnd();

  // Tesselating a polygon
  } else if ( is2D && dmode == DM_NORMAL ) {
    // NOTE: Tesselator object does not make a copy of the data!
    double* point = new double[3];
    point[0] = (*p1)[0];
    point[1] = (*p1)[1];
    point[2] = (*p1)[2];

    gluTessVertex(tesselator, point, point);

  // Drawing as a  polygon (only first vertex!)
  } else {

    if ( dmode == DM_INTRA_LAYER )
      glLineWidth(glParam.lineWidthThin);

    glVertex3dv(*p1);
  }
}


// Draw one mesh boundary element
//
// NOTE 1: Line element is the only element type where
// possible middle node (ie. for 203) is drawn. This means that
// middle nodes are seen in wireframe modes, because bulk elements
// are finally drawn usin there edges (line elements). When drawing shaded
// geometry (surfaces), we do not use the middle nodes.
//
// NOTE 2: Possible center node is never drawn!
// In wireframe mode this would need special handling, because we draw in
// practice only line elements and they know nothing about the center node, which
// is strictly a bulk element level concept! To draw the center node, bulk element
// should consequently handle it directly!
void
Renderer_OGL::drawMeshElement(int elem_type,
                              const int* node_ids,
                              const Point3* nodeNormals,
                              const Point3* nodes,
                              short direction,
                              bool selected)
{
  // Vertex
  // ======
  if (elem_type == 101 ) {
    drawMeshVertexElement(elem_type, node_ids, nodes, selected);

  // Lines
  // =====
  } else if (elem_type >= 202 && elem_type < 300) {

    if (nodeNormals != NULL) {
      double* n = (double*)nodeNormals;
      glNormal3d(n[0], n[1], n[2]);
    }
    drawMeshLineElement(elem_type, node_ids, nodes, selected);

  // Triangles
  // =========
  } else if (elem_type >= 303 && elem_type < 400) {

    drawMeshTriangleElement(elem_type, node_ids, nodes, nodeNormals, direction, selected);

  // Quadrilaterals
  // ==============
  } else if (elem_type >= 404 && elem_type < 500) {

    drawMeshQuadriElement(elem_type, node_ids, nodes, nodeNormals, direction, selected);
  }

}


// Draw one mesh boundary element
void
Renderer_OGL::drawMeshBoundaryElement(int elem_type, const int* elem_nodes,
                              objectDrawingState dstate,
                              const Point3* nodeNormals,
                              const Point3* node_data,
                              short direction)
{
  static double color[4];

  // Set highlite color mode
  if (dstate == DS_SELECTED) {
    glGetDoublev(GL_CURRENT_COLOR, color);
    color[3] = 0.5;
    glColor4dv(color);
  }

  drawMeshElement(elem_type, elem_nodes, nodeNormals, node_data, direction);

  // Reset color mode
  if (dstate == DS_SELECTED) {
    color[3] = 1.0;
    glColor4dv(color);
  }

}


// Draw one mesh bulk element
void
Renderer_OGL::drawMeshBulkElement(int elem_type, const int* elem_nodes,
                              objectDrawingState dstate,
                              const Point3* node_data)
{
  if (dstate == DS_SELECTED) {
    glPolygonMode(GL_FRONT, GL_FILL);
  }

  drawMeshElement(elem_type, elem_nodes, NULL, node_data);

  if (dstate == DS_SELECTED) {
    glPolygonMode(GL_FRONT, GL_LINE);
  }

}


void
Renderer_OGL::drawMeshLineElement(int elem_type, const int* node_ids, const Point3*  nodes,
                                  bool selected)
{

  if (selected) {
    glLineWidth(glParam.lineWidthSelected);
  }

  glBegin(GL_LINE_STRIP);

  // Three nodes
  if ( elem_type == 203 ) {
    glVertex3dv(nodes[node_ids[0]]);
    glVertex3dv(nodes[node_ids[2]]);
    glVertex3dv(nodes[node_ids[1]]);

  // Two nodes
  } else  {
    glVertex3dv(nodes[node_ids[0]]);
    glVertex3dv(nodes[node_ids[1]]);
  }

  glEnd();

  // Draw possible middle point
  if ( elem_type == 203 ) {

    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);
    glPointSize(1.5f);

    glBegin(GL_POINTS);
      glVertex3dv(nodes[node_ids[2]]);
    glEnd();

    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
    glPointSize(1.0f);
  }

  if (selected) {
    glLineWidth(glParam.lineWidthNormal);
  }

}


void
Renderer_OGL::drawMeshTriangleElement(int elem_type, const int* node_ids,
                                      const Point3*  nodes, const Point3*  nodeNormals,
                                      short direction, bool selected)
{
  // Set select color alpha
  if (selected) {
    //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, selectColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, selectColor);
  }

  //glBegin(GL_TRIANGLES);
  glBegin(GL_POLYGON);

    for (int i = 0; i < 3; i++) {
      if ( nodeNormals != NULL ) {
        glNormal3dv(nodeNormals[i]);
      }
      glVertex3dv(nodes[node_ids[i]]);
    }

  glEnd();

  // Reset select color alpha
  if (selected) {
    //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, currentColor);
  }

// For debugging, draw normal vectors!
#if 0
  if ( nodeNormals != NULL ) {
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);
    drawVector((double*)nodes[node_ids[0]], (double*)nodeNormals[0], 0.05 * modelLength[0]);
    drawVector((double*)nodes[node_ids[1]], (double*)nodeNormals[1], 0.05 * modelLength[0]);
    drawVector((double*)nodes[node_ids[2]], (double*)nodeNormals[2], 0.05 * modelLength[0]);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
  }
#endif

}


void
Renderer_OGL::drawMeshQuadriElement(int elem_type, const int* node_ids,
                                    const Point3*  nodes, const Point3*  nodeNormals,
                                    short direction, bool selected)
{
  // Set select color alpha
  if (selected) {
    //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, selectColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, selectColor);
  }

  glBegin(GL_POLYGON);

    for (int i = 0; i < 4; i++) {
      if ( nodeNormals != NULL ) {
        glNormal3dv(nodeNormals[i]);
      }
      glVertex3dv(nodes[node_ids[i]]);
    }

  glEnd();

  // Reset select color alpha
  if (selected) {
    //glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, currentColor);
  }

// For debugging, draw normal vectors!
#if 0
  if ( nodeNormals != NULL ) {
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);
    drawVector((double*)nodes[node_ids[0]], (double*)nodeNormals[0], 0.05 * modelLength[0]);
    drawVector((double*)nodes[node_ids[1]], (double*)nodeNormals[1], 0.05 * modelLength[0]);
    drawVector((double*)nodes[node_ids[2]], (double*)nodeNormals[2], 0.05 * modelLength[0]);
    drawVector((double*)nodes[node_ids[3]], (double*)nodeNormals[3], 0.05 * modelLength[0]);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
  }
#endif

}


void
Renderer_OGL::drawMeshVertexElement(int elem_type, const int* node_ids, const Point3*  nodes,
                                    bool selected)
{
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, meshColor);

  if (selected) {
    glPointSize(5.0f);
  } else {
    glPointSize(1.5f);
  }

  glBegin(GL_POINTS);
    glVertex3dv(nodes[node_ids[0]]);
  glEnd();

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currentColor);
  glPointSize(1.0f);
}


void
Renderer_OGL::drawNurbsCrv(objectDrawingMode dmode, objectDrawingState dstate, short direction,
                              ecif_NurbsCurve& data, int elem_id)
{
  static GLfloat knPt[ecif_MAX_NOF_KNOTS];
  static GLfloat ctPt[3 * ecif_MAX_NOF_CPOINTS];

  // Tables for knot- and control-points.
  // NOTE: Memory for these arries must be allocated
  //    in continues blocks for OpenGL. So no two-phase
  //    dynamic allocation doesn't work !!!###!!!
  // Knot-points are already in the [0,1] domain. A new dynamically allocated
  // structure would be needed only if we have to change their order. However they
  // are in double-format, but we must feed them as GLfloat into OPenGL, so
  // we will copy them inot a new structure.
//GLfloat* knPt = new GLfloat[data.nof_knots];
  // Control points must be transformed into [0,1] domain for trimming!
//GLfloat (*ctPt)[3] = new GLfloat[data.nof_cpoints][3];

  // Trimming curve must be counter-clockwise creature in OPenGL.
  // We must take into account bodyelement's orientation within body!.
  // Here direction == 1 means ccw-order.
  // Calling functions must take care of this ***!!!***
  //
  int i;
  // *** Copy knot-points and adjust the order (if necessary)
  for (i = 0; i < data.nofKnots; i++) {
    if (direction == 1)
      // original values
      knPt[i] = GLfloat(data.knots[i]);
    else
      // order is reversed (read backwards and subtract 1.
      knPt[i] = GLfloat(1 - data.knots[(data.nofKnots - 1) - i]);
  }

  // *** Adjust control-points into (0,1) domain (and revere order if necessary).
  int pos = (direction == 1) ? 0 : data.nofCpoints - 1;
  int incr = (direction == 1) ? 1 : -1;

  // Set control points according to the direction
  for (i = 0; i < data.nofCpoints; i++, pos += incr) {
    int t_cpos = 3 * i;
    ctPt[t_cpos + 0] = data.cpoints[pos][0];
    ctPt[t_cpos + 1] = data.cpoints[pos][1];
    ctPt[t_cpos + 2] = 0.0f;
    //ctPt[i][2] = data.cpoints[pos][2]);
  }


  GLUnurbsObj* nrb = gluNewNurbsRenderer();
  gluBeginCurve(nrb);
    gluNurbsCurve(nrb, data.nofKnots, knPt, 3, ctPt,
                  data.degree + 1, GL_MAP1_VERTEX_3);
  gluEndCurve(nrb);
  gluDeleteNurbsRenderer(nrb);
}


// Draw geometry point (vertex)
void
Renderer_OGL::drawPoint(objectDrawingMode dmode, objectDrawingState dstate, const Point3* point)
{
  //-Use black color
  static Color4f vertex_color = {0.0f, 0.0f, 0.0f, 1.0f};

  //---Rendering
  glPushAttrib(GL_LIGHTING_BIT);

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vertex_color);

  if (dstate == DS_SELECTED) {
    glPointSize(5.0f);
  } else {
    glPointSize(2.5f);
  }

  glBegin(GL_POINTS);
    glVertex3dv(*point);
  glEnd();

  glPopAttrib();
  glPointSize(1.0f);
}


void
Renderer_OGL::drawTwoSidedPolygons(bool turn_on)
{
  if ( turn_on ) {
    glDisable(GL_CULL_FACE);

  } else {
    glEnable(GL_CULL_FACE);
  }
}


void
Renderer_OGL::drawVector(double start[3], double dir[3], double scale, bool draw_point, bool draw_arrow)
{
  static double p1[3];
  static double p2[3];
  int i;

  for (i = 0; i < 3; i++) {
    p1[i] = start[i];
    p2[i] = p1[i] + scale * dir[i];
  }

  glBegin(GL_LINES);
   glVertex3dv(p1);
   glVertex3dv(p2);
  glEnd();

  if ( draw_point ) {
   glPointSize(5.0f);
   glBegin(GL_POINTS);
    glVertex3dv(start);
   glEnd();
  }

  if ( draw_arrow ) {
   glPointSize(2.0f);
   glBegin(GL_POINTS);
    glVertex3dv(start);
   glEnd();
  }
}


bool
Renderer_OGL::hasDisplayList(int list_id)
{
  if ( GL_TRUE == glIsList(list_id) )
    return true;
  else
    return false;
}


void
Renderer_OGL::hideRenderer()
{
  destroyWindow(rendererDisplay, rendererWindow);
}


// Method initialises the GL-environment
void
Renderer_OGL::init()
{
  int i;

  // Reset model transformation
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Reset projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  for (i = 0; i < 3; i++) {
    scaleVector[i]      = 1.0f;
    translateVector[i]  = 0.0f;
  }

  rotateAxis = 0;

  for (i = 0; i < 16; i++) {
    rotateMatrix[i] = 0.0;
  }

  rotateMatrix[0] = 1.0;
  rotateMatrix[5] = 1.0;
  rotateMatrix[10] = 1.0;
  rotateMatrix[15] = 1.0;

  //glEnable(GL_NORMALIZE);


  // Add lights to the scene.

  // Light-0, Undirected, strong ambient light
  GLfloat position0[4] = { 0.0, 0.0, 1.0, 0.0 };
  GLfloat ambient0[4]  = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat diffuse0[4]  = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1.0 };
  //glLightfv(GL_LIGHT0, GL_POSITION, position0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);

  // Light-1, Undirected, weak ambient light
  GLfloat position1[4] = { 0.0, 0.0, 1.0, 0.0 };
  GLfloat ambient1[4]  = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat diffuse1[4]  = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat specular1[4] = { 0.0, 0.0, 0.0, 1.0 };
  //glLightfv(GL_LIGHT1, GL_POSITION, position1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse1);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular1);

  // Light-2, Diffuse light
  GLfloat position2[4] = { 0.0, 0.0, 1.0, 0.0 };
  GLfloat ambient2[4]  = { 0.35, 0.35, 0.35, 1.0 };
  GLfloat diffuse2[4]  = { 0.65, 0.65, 0.65, 1.0 };
  GLfloat specular2[4] = { 0.0, 0.0, 0.0, 1.0 };
  //glLightfv(GL_LIGHT2, GL_POSITION, position2);
  glLightfv(GL_LIGHT2, GL_AMBIENT, ambient2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse2);
  glLightfv(GL_LIGHT2, GL_SPECULAR, specular2);

  // Light-3, Directed, left diffuse light
  GLfloat position3[4] = { -0.5, 0.0, 1.0, 0.0 };
  GLfloat ambient3[4]  = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat diffuse3[4]  = { 0.8, 0.8, 0.8, 1.0 };
  GLfloat specular3[4] = { 0.3, 0.3, 0.3, 1.0 };
  glLightfv(GL_LIGHT3, GL_POSITION, position3);
  glLightfv(GL_LIGHT3, GL_AMBIENT, ambient3);
  glLightfv(GL_LIGHT3, GL_DIFFUSE, diffuse3);
  glLightfv(GL_LIGHT3, GL_SPECULAR, specular3);

  // Right-3, Directed, right diffuse light
  GLfloat position4[4] = { 0.5, 0.0, 1.0, 0.0 };
  GLfloat ambient4[4]  = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat diffuse4[4]  = { 0.8, 0.8, 0.8, 1.0 };
  GLfloat specular4[4] = { 1.0, 1.0, 1.0, 1.0 };
  glLightfv(GL_LIGHT4, GL_POSITION, position4);
  glLightfv(GL_LIGHT4, GL_AMBIENT, ambient4);
  glLightfv(GL_LIGHT4, GL_DIFFUSE, diffuse4);
  glLightfv(GL_LIGHT4, GL_SPECULAR, specular4);


  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  //glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

  // Select active light(s)
  //
  setLightDirected();

  // Turn on the lights
  glEnable(GL_LIGHTING);

  // Shading mode
  //glShadeModel(GL_FLAT);
  glShadeModel(GL_SMOOTH);

  // Set default normal vector for 2D geometry
  glNormal3d(0.0, 0.0, 1.0);

  setParameters();

  reshape();
}


void
Renderer_OGL::reset()
{
  // If not visible
  if ( winX == 0 || winY == 0 ) {
    return;
  }

  init();
  refresh();
}


void
Renderer_OGL::resetData()
{
  deleteData();
  createData();
  init();
}


void
Renderer_OGL::makeRasterFont()
{
  GLuint i;
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  fontOffset = glGenLists (128);

  for (i = 32; i < 127; i++) {
    glNewList(i+fontOffset, GL_COMPILE);
    glBitmap(8, 13, 0.0f, 2.0f, 10.0f, 0.0f, rasterFont[i-32]);
    glEndList();
  }
}


void
Renderer_OGL::name_delete(int name_id)
{
  glPopName();
}


void
Renderer_OGL::name_replace(int name_id)
{
  glLoadName((GLuint)name_id);
}


void
Renderer_OGL::name_save(int name_id)
{
  GLint stack[1];
  stack[0] = 0;
  glPushName((GLuint)name_id);
  glGetIntegerv(GL_NAME_STACK_DEPTH, stack);
}


// Scale double-point to [0,1]
void
Renderer_OGL::normalize1_point(Point3& p)
{
  for (int i = 0; i < 3; i++) {
    if ( modelLength[i] > 0 )
      p[i] = (p[i] - modelStart[i]) / modelLength[i];
  }
}


// Scale float-point to [0,1]
void
Renderer_OGL::normalize1_point(Point3f& p)
{
  for (int i = 0; i < 3; i++) {
    if ( modelLength[i] > 0 )
      p[i] = (p[i] - modelStart[i]) / modelLength[i];
  }
}


// Scale double-point to [-1,1]
void
Renderer_OGL::normalize2_point(Point3& p)
{
  for (int i = 0; i < 3; i++) {
    p[i] = 2 * (p[i] - modelCenter[i]) / modelLength[i];
  }
}


// Scale float-point to [-1,1]
void
Renderer_OGL::normalize2_point(Point3f& p)
{
  for (int i = 0; i < 3; i++) {
    p[i] = 2 * (p[i] - modelCenter[i]) / modelLength[i];
  }
}


void
Renderer_OGL::printString(char* s)
{
#if defined(WIN32)
  glPushAttrib(GL_LIST_BIT);
  glListBase(fontOffset);
  glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte*) s);
  glPopAttrib();

#else
  // NOTE: Push/Pop attributes crashes (after scaling)
  // in Linux, so we skip them in Unix!!!
  // Consequences...?
  glListBase(fontOffset);
  glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte*) s);
#endif
}


// Process mouse selections when editing mesh boundaries, selecting boundaries etc.
//
void
Renderer_OGL::processSelection(mouseAction maction,int mk_state, int x_pos, int y_pos)
{
  if ( winX == 0 || winY == 0) {
    return;
  }

  // Update mouse-keyboard state attribute
  mkState = mk_state;

  bool in_mesh_draw_mode = model->getFlagValue(DRAW_SOURCE_MESH);
  bool in_draw_mode = inBoxDrawMode || inLineDrawMode || inVectorDrawMode;
  bool do_picking = false;

  // Select by mouse action
  switch (maction) {

  case MOUSE_DBL_CLCK:
    // Turn on picking mode, turn off extended object selection
    do_picking = true;
    model->setFlagValue(SELECT_OBJECTS, SELECT_OBJECTS_EXTEND, false);
    break;

  case MOUSE_CTRL_DBL_CLCK:
    // Turn on picking mode, turn on extend selection
    do_picking = true;
    model->setFlagValue(SELECT_OBJECTS, SELECT_OBJECTS_EXTEND, true);
    break;

  case MOUSE_CTRL_CLCK:
    // Effective only in editing mode
    if ( !inMeshEditingMode )
      return;
    break;

  case MOUSE_SHFT_MOVE:
    // Effective only in editing mode
    if ( !(inMeshEditingMode || in_draw_mode) )
      return;
    break;

  default:
    return;
    break;
  }

  //----Ray-cast picking for the mesh boundary/bulk elements

  // Editing bodies/boundaries
  if ( inMeshEditingMode || in_mesh_draw_mode || in_draw_mode
     ) {

    // Find model coordinates for the window coordinates
    Point3 p1, p2;
    int res1 = gluUnProject( x_pos, y_pos, 0.0, modelMatrix, projectionMatrix, viewport,
                             &p1[0], &p1[1], &p1[2]);
    int res2 = gluUnProject( x_pos, y_pos, 1.0, modelMatrix, projectionMatrix, viewport,
                             &p2[0], &p2[1], &p2[2]);

    // Line direction for the casted ray in the form:
    // start-point + t*dir-vec
    Point3 ldir;
    diff3(p2, p1, ldir);

    model->setCurrentPickInfo(p1,ldir);

    if ( inMeshEditingMode || in_mesh_draw_mode ) {

      int bndr_id;
      int bd1_id, lr1_id;
      int bd2_id, lr2_id;

      bool use_cur_bndr = inMeshEditingMode;

      int fem_id = model->findSelectedMeshBoundaryElement(this, p1, ldir, use_cur_bndr,
                                                          bndr_id,
                                                          bd1_id, lr1_id,
                                                          bd2_id, lr2_id);

      if ( fem_id == NO_INDEX ) {
        return;
      }

      if ( inMeshEditingMode ) {
        model->meshBoundaryElementSelected(this, fem_id);

      } else {
        model->boundarySelected(this, bndr_id, bd1_id, lr1_id, bd2_id, lr2_id, true);
      }

      display(this);
      return;
    }
  }

  if (!do_picking) return;

  //---Save current projection
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  //---Build picking "cube"

  //-Window center
  // Convert mouse's location in the window to the location in the model world
  double x_ratio = 1.0;
  double y_ratio = 1.0;

  if ( viewport[2] > 0 ) {
    x_ratio = (double)x_pos / (viewport[2] - viewport[0]);
  }

  if ( viewport[3] > 0 ) {
    y_ratio = (double)(y_pos) / (viewport[3] - viewport[1]);
  }

  double x_center = projection[0] + x_ratio * (projection[1] - projection[0]);
  double y_center = projection[2] + y_ratio * (projection[3] - projection[2]);

  //-Picking window size

  double pixels;
  // NOTE:If this is small the gluEndSurface(background)
  // in 2D nurbs trimming slows down dramatically!!!
  if ( model->getDimension() == ECIF_2D         &&
       model->getFlagValue(DRAW_SOURCE_CAD)     &&
       model->getFlagValue(DRAW_TARGET_BODIES)

     ) {
    pixels = 10.0;
  } else {
    pixels = 1.0;
  }

  double x_del = pixels *  (projection[1] - projection[0]) / winX;
  double y_del = pixels *  (projection[3] - projection[2]) / winY;

  x_del *= scaleVector[0];
  y_del *= scaleVector[1];

  //----OpenGL "selection picking" for boundaries/bodies
  // "Double click selection mode"
  //---Set selection buffer
  glSelectBuffer(selectionBufferSize, selectionBuffer);

  //---Goto selection mode
  // Comment this call if you want to see the screen after
  // selection (for debugging!)
  glRenderMode(GL_SELECT);
  renderMode = RENDER_SELECTION;

  //---Init selection name buffer
  glInitNames();

  glLoadIdentity();

  //-Clipping viewing volume
  glOrtho(x_center - x_del, x_center + x_del,
          y_center - y_del, y_center + y_del,
          projection[4], projection[5]);


  //---Get selections
  //inMeshPickingMode = 1;
  selectionHits = 0;
  display(this);
  selectionHits = glRenderMode(GL_RENDER);

  inMeshPickingMode = 0;
  renderMode = RENDER_NORMAL;

  //---Restore original projection
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW); // Default matrix mode from now on

  display(this);
}


void
Renderer_OGL::meshBoundaryElementSelectionHits()
{
  if ( selectionHits == 0)
    return;

  int picked_fem_id = NO_INDEX;
  GLuint tot_min_z = 0xffffffff;
  GLuint tot_max_z = 0;
  int record_start = 0;

  for (int i = 0; i < selectionHits; i++) {

    int nof_names = selectionBuffer[record_start];
    GLuint min_z = selectionBuffer[record_start + 1];
    GLuint max_z = selectionBuffer[record_start + 2];
    int body_id = selectionBuffer[record_start + 3];
    int bndr_id = selectionBuffer[record_start + 4];

    int fem_id = NO_INDEX;

    if (nof_names > 2) {
      fem_id = selectionBuffer[record_start + 5];
    }

    record_start += 3 + nof_names;

    if ( fem_id != NO_INDEX && min_z < tot_min_z ) {
      tot_min_z = min_z;
      picked_fem_id = fem_id;
    }
  }

  if ( picked_fem_id != NO_INDEX ) {
    model->meshBoundaryElementSelectionHit(this, picked_fem_id);
  }

  selectionHits = 0;
}


void
Renderer_OGL::meshBulkElementSelectionHits()
{
  if ( selectionHits == 0)
    return;

  int picked_fem_id = NO_INDEX;
  GLuint tot_min_z = 0xffffffff;
  GLuint tot_max_z = 0;
  int record_start = 0;
  for (int i = 0; i < selectionHits; i++) {
    int nof_names = selectionBuffer[record_start];
    GLuint min_z = selectionBuffer[record_start + 1];
    GLuint max_z = selectionBuffer[record_start + 2];
    int body_id = selectionBuffer[record_start + 3];

    int fem_id = NO_INDEX;
    if (nof_names > 1)
      fem_id = selectionBuffer[record_start + 4];

    record_start += 3 + nof_names;

    if ( fem_id != NO_INDEX &&
         max_z > tot_max_z
       ) {
      tot_max_z = max_z;
      picked_fem_id = fem_id;
    }
  }

  if ( picked_fem_id != NO_INDEX )
    model->meshBulkElementSelectionHit(this, picked_fem_id);

  selectionHits = 0;
}


void
Renderer_OGL::refresh()
{
  if (!visible) {
    displayRenderer();

  } else {
    display(this);
  }
}


void
Renderer_OGL::removeDisplayLists()
{
  IdList::iterator itr = displayListIds->begin();

  while ( itr != displayListIds->end() ) {
    int id = (*itr++);
    glDeleteLists(id, 1);
  }

  displayListIds->clear();
}


// Setup renderer
void
Renderer_OGL::reshape()
{
  static GLdouble ratio[3];
  static GLdouble marginal[3];
  static GLdouble aspect;

  // Window (client) size
  findRendererWindowSize(winX, winY);

  // If window is minimized
  if (winX == 0 || winY == 0 ) {
    return;
  }

  // Set viewport to full window (client) size
  glViewport(0, 0, winX, winY);

  // Correct for window aspect ratio
  GLdouble xw_ratio = (GLfloat)winY / winX;
  GLdouble yw_ratio = 1;

  // Correct model geometry aspect ratio from the -1, 1 normalized world
  // NOTE: This is in "opposite" direction compared to window ratio correction!
  GLdouble xg_ratio = 1;

  if ( modelLength[1] > 0 ) {
    xg_ratio = modelLength[0] / modelLength[1];
  }

  GLdouble yg_ratio = 1;

  aspect = yw_ratio / xw_ratio;

  // Calculate the final aspect-ratios for an ortho projection
  // Here the larger aspect is set to value 1.0 and the smaller
  // is normed > 1.0, (this way smaller is made smaller in ortho projection!)
  double x_ratio = xw_ratio * xg_ratio;
  double y_ratio = yw_ratio * yg_ratio;
  double z_ratio = 1.0f;

  if (x_ratio > y_ratio && y_ratio > 0 ) {
    y_ratio = x_ratio / y_ratio;
    x_ratio = 1.0f;

  } else if ( x_ratio > 0 ) {
    x_ratio = y_ratio / x_ratio;
    y_ratio = 1.0f;
  }

  // Set data into arries (to make filling projection[] easier)
  ratio[0] = x_ratio;
  ratio[1] = y_ratio;
  ratio[2] = z_ratio;

  double max_length = modelLength[0];

  if (modelLength[1] > max_length) {
    max_length = modelLength[1];
  }

  if (modelLength[2] > max_length) {
    max_length = modelLength[2];
  }

  marginal[0] = 0.10 * modelLength[0];  // 10% x marginal
  marginal[1] = 0.10 * modelLength[1];  // 10% y marginal
  marginal[2] = 10.0 * scaleVector[2] * max_length; // That's a marginal!
  //marginal[2] = 10.0 * max_length; // That's also a marginal!

  // Clipping values for orthogonal projection
  for (int i = 0; i < 3; i++) {

    // correct with the ratio
    double half = 0.5 * modelLength[i] * ratio[i];
    double marg = marginal[i] * ratio[i];

    // add marginal
    projection[0 + 2 * i] = modelCenter[i] - half - marg;
    projection[1 + 2 * i] = modelCenter[i] + half + marg;
  }

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();


  // Ortho projection
  if ( useOrthoProjection ) {

    glOrtho(projection[0], projection[1],
            projection[2], projection[3],
            projection[4], projection[5]
           );

    //glOrtho(-aspect * max_length, aspect * max_length, -1.0 * max_length, 1.0 * max_length,
    //        //projection[4], projection[5]
    //        scaleVector[2] * -max_length, scaleVector[2] * 3 * max_length
    //       );

    //translateVector[2] = - scaleVector[2] * max_length;

  // Perspective projection
  // NOTE: This is not correct, do NOT use this!
  } else {

    //glFrustum(projection[0], projection[1],
    //          projection[2], projection[3],
              //2.0, 2.0 + 5 * max_length
    //          0.0, 1.0 + 3 * max_length
    //       );
    gluPerspective(15.0, aspect, 1.0, 1.0 + 3 * scaleVector[2] * max_length);

    translateVector[2] = -1.0 - scaleVector[2] * max_length;

  }

  glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);
}


void
Renderer_OGL::rotate(double angle)
{

  // Note: Reaction to rotation angle is slower when model is
  // scaled to a large one!
  //
  double factor = 2.0 / (scaleVector[0] + scaleVector[1]);

  if ( factor < 1.0 ) {
    angle *= sqrt(factor);
  } else {
    angle *= factor;
  }

  glMatrixMode(GL_MODELVIEW);

  glPushMatrix();

  glLoadIdentity();

  GLdouble* rm = rotateMatrix;

  // 2D model
  if ( is2D ) {
    glRotated(angle, 0, 0, 1);

  // 3D model
  } else if (rotateAxisX) {
    glRotated(angle, rm[0], rm[1], rm[2]);

  } else if (rotateAxisY) {
    glRotated(angle, rm[4], rm[5], rm[6]);

  } else if (rotateAxisZ) {
    glRotated(angle, rm[8], rm[9], rm[10]);

  } else if ( rotateAxis == 0 ) {
    glRotated(angle, 1.0f, 0.0f, 0.0f);

  } else if ( rotateAxis == 1 ) {
    glRotated(angle, 0.0f, 1.0f, 0.0f);

  } else if ( rotateAxis == 2 ) {
    glRotated(angle, 0.0f, 0.0f, 1.0f);
  }

  glMultMatrixd(rotateMatrix);

  glGetDoublev(GL_MODELVIEW_MATRIX, rotateMatrix);

  glPopMatrix();

  glMultMatrixd(rotateMatrix);

}


// Rotates a fixed amount to given direction around given axis
void
Renderer_OGL::rotate(int axis, short direction)
{
  if ( winX == 0 || winY == 0 )
    return;

  int base_rt[] = {winX, winY, winZ};

  if (axis == 2)
    direction *= -1;

  // Rotate 1.5 pixel equivalent
  double rt = direction * 1.5;
  rotate(axis, rt);
}


// Rotate a given "screen amount" around given axis
// This is used mainly when using mouse movements
// for rotation
void
Renderer_OGL::rotate(int axis, int amount)
{
  // In 2D we rotate only around z-axis, and
  // only if z-rotation is turned on!
  if (is2D && axis != 2) {
    return;
  }

  //int base[3] = {winX, winY, winZ};
  int base[3] = {winX, winY, (winX + winY) /2};

  GLfloat degrees = 0.0;
  if ( base[axis] != 0 )
    degrees = 90 * (GLfloat)amount / base[axis];

  if (axis == 2)
    degrees *= -1;

  rotate(axis, degrees);
}


void
Renderer_OGL::rotate(int axis, double degrees)
{
  if ( is2D ) {
    rotate(degrees);
  }

  else {
    rotate(degrees);
    rotateAxis  = axis;
  }

  display(this);
}


void
Renderer_OGL::scale()
{
  glScaled(scaleVector[0], scaleVector[1], scaleVector[2]);
}


void
Renderer_OGL::scale(short direction)
{
  if ( winX == 0 || winY == 0 )
    return;

  // Scale 5% of the avg. window size
  int sc = direction * 0.05 * 0.5 * (winX + winY);
  scale(sc);
}


void
Renderer_OGL::setDrawMode()
{
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_POLYGON_OFFSET_FILL);

  // In 3D
  if (!is2D) {

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    //glDepthFunc(GL_LESS);
    //glDisable(GL_DEPTH_TEST);

    //glPolygonOffset(1.0f, 1.1f);
    glPolygonOffset(2.1f, 2.1f);
    glEnable(GL_POLYGON_OFFSET_FILL);
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glPolygonMode(GL_FRONT, GL_FILL);
  //glPolygonMode(GL_BACK, GL_LINE);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  //glDisable(GL_CULL_FACE);
}


void
Renderer_OGL::setLightDirected()
{
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LIGHT2);
  glDisable(GL_LIGHT3);
  glDisable(GL_LIGHT4);

  if (is2D) {
    glEnable(GL_LIGHT0);   // Undirected strong ambient

  } else {
    glEnable(GL_LIGHT1); // Undirected weak ambient
    //glEnable(GL_LIGHT2); // Undirected diffuse
    glEnable(GL_LIGHT3); // Directed ,left diffuse
    //glEnable(GL_LIGHT4); // Directed ,right diffuse
  }
 }


void
Renderer_OGL::setLightUndirected()
{
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LIGHT2);
  glDisable(GL_LIGHT3);
  glDisable(GL_LIGHT4);

  if (is2D) {
    glEnable(GL_LIGHT0);   // Undirected strong ambient

  } else {
    glEnable(GL_LIGHT0);   // Undirected strong ambient
    glEnable(GL_LIGHT1); // Undirected weak ambient
    glEnable(GL_LIGHT2); // Undirected diffuse
  }
}


// Set general drawing paramteters.
void
Renderer_OGL::setParameters()
{
  int i;

  // Model's scale for normalization parameters (class attributes)
  // Transformation (x + shift) / norm transforms data into [-1, 1]
  RangeVector rv;
  model->getBoundingBox(rv);

  // If 2D model, set z-dimensions to [-1,1] for proper
  // clipping in orto projection
  if ( is2D ) {
    rv[4] = -1.0;
    rv[5] =  1.0;
  }

  for (i = 0; i < 3; i++) {
    int i1 = 2 * i;
    int i2 = i1 + 1;
    modelStart[i]  = rv[i1];
    modelEnd[i]    = rv[i2];
    modelCenter[i] = (rv[i2] + rv[i1]) / 2;
    modelLength[i] = rv[i2] - rv[i1];
  }

  if ( isZero(modelLength[0]) ) {
    modelLength[0] = max(modelLength[1], modelLength[2]);
  }

  if ( isZero(modelLength[1]) ) {
    modelLength[1] = max(modelLength[0], modelLength[2]);
  }

  if ( !is2D && isZero(modelLength[2]) ) {
    modelLength[2] = max(modelLength[0], modelLength[1]);
  }

}


void
Renderer_OGL::setRendererInfo()
{
  WindowInfo winfo;

  if (status == HAS_NO_WINDOW){
    createGLWindow(NULL, "Test", 0, 0, 1, 1, winfo);
  }

  GLfloat f1[1];
  GLfloat f2[2];

  glGetFloatv(GL_LINE_WIDTH_GRANULARITY, f1);
  rendererInfo.LINE_WIDTH_GRANULARITY = f1[0];

  glGetFloatv(GL_LINE_WIDTH_RANGE, f2);
  rendererInfo.LINE_WIDTH_RANGE[0] = f2[0];
  rendererInfo.LINE_WIDTH_RANGE[1] = f2[1];

  theControlCenter->getUI()->updateRendererInfo(rendererInfo);

  if (status == HAS_NO_WINDOW){
    destroyWindow(winfo.display, winfo.window);
  }

}


void
Renderer_OGL::scale(int amount)
{
  GLfloat sc_amount = 0.0;

  if ( winY > 0 ) {
    sc_amount = 1.0f + (GLfloat)amount / winY;
  }

  //glMatrixMode(GL_MODELVIEW);
  //glScaled(sc_amount, sc_amount, sc_amount);
  scaleVector[0] *= sc_amount;
  scaleVector[1] *= sc_amount;
  scaleVector[2] *= sc_amount;

  reshape();

  display(this);
}


// Generic test-procedure for testing renderer related stuff
void
Renderer_OGL::test()
{
}


void
Renderer_OGL::transform_scene()
{
  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity();


  translate();

  // Re-center from origo for translation
  glTranslated(modelCenter[0], modelCenter[1], modelCenter[2]);

  scale();
  rotate();

  // Center to origo before rotate and scale
  glTranslated(-modelCenter[0], -modelCenter[1], -modelCenter[2]);

  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
}


void
Renderer_OGL::translate()
{
  glTranslated(translateVector[0]*scaleVector[0],
               translateVector[1]*scaleVector[1],
               translateVector[2]*scaleVector[2]
              );
}


void
Renderer_OGL::translate(int coordinate, short direction)
{
  if ( winX == 0 || winY == 0 )
    return;

  int base_tv[] = {winX, winY, winZ};
  int tv[] = {0, 0, 0};

  // Translate 2.5%
  tv[coordinate] = (int) direction * 0.025 * base_tv[coordinate];
  translate(tv[0], tv[1], tv[2]);
}


void
Renderer_OGL::translate(int x_delta, int y_delta, int z_delta)
{
  if ( winX == 0 || winY == 0 )
    return;

  GLdouble avg = 0.5 * (modelLength[0] + modelLength[1]);

  GLdouble x_amount = avg * (GLdouble)x_delta / winX;
  GLdouble y_amount = avg * (GLdouble)y_delta / winY;
  GLdouble z_amount = modelLength[2] * (GLdouble)z_delta / winZ;

  // The Macig Speedup Factor of Translations
  // (from The Handbook of Contants)
  double factor = 1.5;

  translateVector[0] += factor * x_amount / scaleVector[0];
  translateVector[1] += factor * y_amount / scaleVector[1];
  translateVector[2] += z_amount / scaleVector[2];

  display(this);
}


void
Renderer_OGL::useDisplayList(int list_id)
{
   glCallList(list_id);
}


// *** System dependent section ***
#ifdef WIN32
#include "ecif_renderer_OGL_WIN32.hpp"
#else
#include "ecif_renderer_OGL_UNIX.hpp"
#endif
// *** End system dependent ***


//************************
//*** FONTS DEFINITION ***
//************************
// Set value for the static attribute:
GLubyte Renderer_OGL::rasterFont[][13] = {
//GLubyte rasterFont[][13] = {
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x36, 0x36, 0x36},
  {0x00, 0x00, 0x00, 0x66, 0x66, 0xff, 0x66, 0x66, 0xff, 0x66, 0x66, 0x00, 0x00},
  {0x00, 0x00, 0x18, 0x7e, 0xff, 0x1b, 0x1f, 0x7e, 0xf8, 0xd8, 0xff, 0x7e, 0x18},
  {0x00, 0x00, 0x0e, 0x1b, 0xdb, 0x6e, 0x30, 0x18, 0x0c, 0x76, 0xdb, 0xd8, 0x70},
  {0x00, 0x00, 0x7f, 0xc6, 0xcf, 0xd8, 0x70, 0x70, 0xd8, 0xcc, 0xcc, 0x6c, 0x38},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x1c, 0x0c, 0x0e},
  {0x00, 0x00, 0x0c, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c},
  {0x00, 0x00, 0x30, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x30},
  {0x00, 0x00, 0x00, 0x00, 0x99, 0x5a, 0x3c, 0xff, 0x3c, 0x5a, 0x99, 0x00, 0x00},
  {0x00, 0x00, 0x00, 0x18, 0x18, 0x18, 0xff, 0xff, 0x18, 0x18, 0x18, 0x00, 0x00},
  {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06, 0x03, 0x03},
  {0x00, 0x00, 0x3c, 0x66, 0xc3, 0xe3, 0xf3, 0xdb, 0xcf, 0xc7, 0xc3, 0x66, 0x3c},
  {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0xe7, 0x7e},
  {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0x07, 0x03, 0x03, 0xe7, 0x7e},
  {0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0xff, 0xcc, 0x6c, 0x3c, 0x1c, 0x0c},
  {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x03, 0x03, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
  {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x03, 0x7f, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e},
  {0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x1c, 0x1c, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x06, 0x0c, 0x18, 0x30, 0x60, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06},
  {0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x60},
  {0x00, 0x00, 0x18, 0x00, 0x00, 0x18, 0x18, 0x0c, 0x06, 0x03, 0xc3, 0xc3, 0x7e},
  {0x00, 0x00, 0x3f, 0x60, 0xcf, 0xdb, 0xd3, 0xdd, 0xc3, 0x7e, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18},
  {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
  {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e},
  {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06},
  {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3},
  {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e},
  {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c},
  {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
  {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff},
  {0x00, 0x00, 0x3c, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x3c},
  {0x00, 0x03, 0x03, 0x06, 0x06, 0x0c, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x60, 0x60},
  {0x00, 0x00, 0x3c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x3c},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18},
  {0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x38, 0x30, 0x70},
  {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0x7f, 0x03, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
  {0x00, 0x00, 0x7e, 0xc3, 0xc0, 0xc0, 0xc0, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03, 0x03},
  {0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x33, 0x1e},
  {0x7e, 0xc3, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0},
  {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x00, 0x18, 0x00},
  {0x38, 0x6c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x00, 0x00, 0x0c, 0x00},
  {0x00, 0x00, 0xc6, 0xcc, 0xf8, 0xf0, 0xd8, 0xcc, 0xc6, 0xc0, 0xc0, 0xc0, 0xc0},
  {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78},
  {0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xfc, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x7c, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x7c, 0x00, 0x00, 0x00, 0x00},
  {0xc0, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00},
  {0x03, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xfe, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x1c, 0x36, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x00},
  {0x00, 0x00, 0x7e, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xdb, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18, 0x3c, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
  {0xc0, 0x60, 0x60, 0x30, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0xff, 0x60, 0x30, 0x18, 0x0c, 0x06, 0xff, 0x00, 0x00, 0x00, 0x00},
  {0x00, 0x00, 0x0f, 0x18, 0x18, 0x18, 0x38, 0xf0, 0x38, 0x18, 0x18, 0x18, 0x0f},
  {0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18},
  {0x00, 0x00, 0xf0, 0x18, 0x18, 0x18, 0x1c, 0x0f, 0x1c, 0x18, 0x18, 0x18, 0xf0},
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x8f, 0xf1, 0x60, 0x00, 0x00, 0x00}
  };
