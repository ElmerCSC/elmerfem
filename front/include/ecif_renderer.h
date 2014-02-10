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

/***********************************************************************
Program:    ELMER Front
Module:     ecif_renderer.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   An abstract base class for geometry renderers.

************************************************************************/


#ifndef _ECIF_RENDERER_
#define _ECIF_RENDERER_

#include "ecif_def.h"
#include "ecif_geometry.h"


// A structure for storing renderer environmment info
struct RendererInfo {
  double LINE_WIDTH_GRANULARITY;
  double LINE_WIDTH_RANGE[2];
};


class Renderer {
friend class Control;
friend class Model;
public:
  virtual ~Renderer() {};
  virtual void clean() {};  // Clean buffers etc.
  virtual void clear() {};  // Clear screen
  virtual void deleteDisplayList(int list_id) {}
  virtual void displayRenderer() {};
  virtual void draw() = 0;
  virtual void drawAllBodies() = 0;
  virtual void drawAllMeshBodies() = 0;
  virtual void drawCoordinateAxis(int x_pos, int y_pos) {}
  virtual void drawElementLabel(char* lbl, Point3& lbl_p) {}
  virtual void drawLine(objectDrawingMode dmode, objectDrawingState dstate, short direction, const Point3* start, const Point3* end, int elem_id = 0) {}
  virtual void drawMeshBoundaryElement(int elem_type, const int* elem_nodes,
                                       objectDrawingState dstate,
                                       const Point3* elem_normal,
                                       const Point3* node_data,
                                       short direction) {}
  virtual void drawMeshBulkElement(int elem_type, const int* elem_nodes,
                                   objectDrawingState dstate,
                                   const Point3* node_data) {}
  virtual void drawMeshElement(int elem_type,
                               const int* elem_nodes,
                               const Point3* elem_normal,
                               const Point3* node_data,
                               short direction = 1,
                               bool selected = false) {}
  virtual void drawMeshElement2(int elem_type,
                               const int* elem_nodes,
                               const Point3* elem_normal,
                               const Point3* node_data,
                               short direction = 1,
                               bool selected = false) {}
  virtual void drawNurbsCrv(objectDrawingMode dmode, objectDrawingState dstate, short direction, ecif_NurbsCurve& data, int elem_id = 0) {}
  virtual void drawNurbsSrf(objectDrawingMode dmode, objectDrawingState dstate, short direction, ecif_NurbsSurface& data, int elem_id = 0) {}
  virtual void drawPoint(objectDrawingMode dmode, objectDrawingState dstate, const Point3* point) {}
  virtual void drawTwoSidedPolygons(bool value) {}
  virtual void drawVector(double start[3], double dir[3], double scale, bool draw_point = true, bool draw_arrow = true) {}
  virtual void dummyWindowProc() {}
  void findRendererWindowSize(int& w, int& h);
  int getMouseKeyboardState() { return mkState; }
  void getRotateAxis(int& default_axis);
  renderingMode getRenderingMode() {return renderMode;}
  void getRendererWindowSize(int& w, int& h) { w = winX; h = winY;}
  virtual bool hasDisplayList(int list_id) {return false;}
  virtual void hideRenderer() {};
  static void initClass(Model* model);
  virtual bool isEditingMesh() { return inMeshEditingMode;}
  virtual bool isPickingMesh() { return inMeshPickingMode;}
  virtual bool isVisible() { return visible;}
  virtual void name_delete(int name_id) {}
  virtual void name_replace(int name_id) {}
  virtual void name_save(int name_id) {}
  virtual void printString(char* str) {}
  virtual void processSelection(mouseAction maction, int mk_state, int x_pos, int y_pos) {};
  virtual void refresh() {}
  virtual void removeDisplayLists() {}
  virtual void reset() {}
  virtual void resetData() {}
  virtual void reshape() {}
  virtual void rotate(int axis, int amount) {}
  virtual void rotate(int axis, short direction) {}
  virtual void rotate(int axis, double degrees) {}
  virtual void scale(int amount) {}
  virtual void scale(short direction) {}
  virtual void selectBoundary(int elem_id,
                              int body1_id, int layer1_id,
                              int body2_id, int layer2_id,
                              bool accept_body_change = false, bool update_gui = true);
  virtual void selectBoundaries(int nof_elements, int* elem_ids,
                                int* body1_ids, int* layer1_ids,
                                int* body2_ids, int* layer2_ids,
                                bool accept_body_change = false, bool update_gui = true);
  virtual void setParameters() {}
  void setDimension(ecif_modelDimension dimension);
  void setDrawBox(bool in_draw_mode);
  void setDrawVector(bool in_draw_mode);
  void setEditBoundaries(bool in_edit_mode);
  virtual void setLightDirected() {}
  virtual void setLightUndirected() {}
  void setRenderingMode(renderingMode render_mode) {renderMode = render_mode;}
  void setRotatePriorities(bool rot_x, bool rot_y, bool rot_z);
  void setSimulationDimension(ecif_modelDimension simulation_dimension);
  virtual void setWindowTitle(char* title) {}
  virtual void startDisplayList(int list_id) {}
  virtual void startDrawingCadBody() = 0;
  virtual void startDrawingCadBodyElementLoop(bool is_first_loop = true) = 0;
  virtual void startDrawingCadSurface() = 0;
  virtual void startDrawingMeshSurface() = 0;
  virtual void startDrawingMeshSurfaceEdges() = 0;
  virtual void stopDisplayList(int list_id) {}
  virtual void stopDrawingCadBody() = 0;
  virtual void stopDrawingCadBodyElementLoop() = 0;
  virtual void stopDrawingCadSurface() = 0;
  virtual void stopDrawingMeshSurface() = 0;
  virtual void stopDrawingMeshSurfaceEdges() = 0;
  virtual void test() {}
  virtual void transform_scene() {}
  virtual void translate(int x_amount, int y_amount, int z_amount) {}
  virtual void translate(int coordinate, short direction) {}
  virtual void useDisplayList(int list_id) {}
protected:
  static Timer* doubleClickTimer;
  static Model* model;
  static Timer* mouseMoveTimer;
  static RendererInfo rendererInfo;
  static rendererStatus status;
  static Control* theControlCenter;
  static bool visible;
  Hinst appInstance;
  Color4f currentColor;
  bool inBoxDrawMode;
  bool inLineDrawMode;
  bool inMeshEditingMode;
  bool inMeshPickingMode;
  bool inVectorDrawMode;
  bool is2D;
  bool is2DSimulation;
  Color4f meshColor;
  int mkState; // mouse-keyborad state
  Hdisplay* rendererDisplay;
  Hwindow rendererWindow;
  bool rotateAxisX;
  bool rotateAxisY;
  bool rotateAxisZ;
  Color4f selectColor;
  renderingMode renderMode;
  bool useOrthoProjection;
  int winX;
  int winY;
  int winZ;
};

#endif
