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
Module:     ecif_renderer_ogl.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   A base class for geometry renderers. 
  Currently only OpenGL renderer implemented.

************************************************************************/


#ifndef _ECIF_RENDERER_OGL_
#define _ECIF_RENDERER_OGL_

#include <GL/gl.h>
#include <GL/glu.h>

#include "ecif_renderer.h"


struct WindowInfo {
  Hdisplay* display;          
  Hwindow window;
};

// A structure for transmitting GL-specific parameters
struct GLParam {
  GLParam();
  GLfloat lineWidthCoordAxis;
  GLfloat lineWidthThin;
  GLfloat lineWidthNormal;
  GLfloat lineWidthSelected;
};


class Body;
class BodyElement;
class Model;
class GcPoint;
struct MeshInfo;


class Renderer_OGL : public Renderer {
friend class Model;
public:
  Renderer_OGL(Hinst appInstance, enum ecif_modelDimension dim);
  ~Renderer_OGL();
  static void setRendererInfo();
  void clear(); // Clear screen
  void deleteDisplayList(int list_id);
  void deleteDisplayLists();
  void deleteTesselatorPoints();
  void displayRenderer();
  void draw() { swapBuffers(); }
  void drawAllBodies();
  void drawAllMeshBodies();
  //void drawBody(Body* body);
  virtual void drawCoordinateAxis(int x_pos, int y_pos);
  void drawElementLabel(char* lbl, Point3& lbl_p);
  void drawLine(objectDrawingMode dmode, objectDrawingState dstate, short direction, const Point3* start, const Point3* end, int id = 0);
  void drawMeshBoundaryElement(int elem_type, const int* elem_nodes,
                               objectDrawingState dstate,
                               const Point3* elem_normal,
                               const Point3* node_data,
                               short direction);
  void drawMeshBulkElement(int elem_type, const int* elem_nodes,
                           objectDrawingState dstate,
                           const Point3* node_data);
  void drawMeshElement(int elem_type,
                       const int* elem_nodes,
                       const Point3* elem_normal,
                       const Point3* node_data,
                       short direction = 1,
                       bool selected = false);
  void drawNurbsCrv(objectDrawingMode dmode, objectDrawingState dstate, short direction, ecif_NurbsCurve& data, int id = 0);
  void drawPoint(objectDrawingMode dmode, objectDrawingState dstate, const Point3* point);
  void drawTwoSidedPolygons(bool turn_on);
  void drawVector(double start[3], double dir[3], double scale, bool draw_point = true, bool draw_arrow = true);
  void dummyWindowProc();
  bool hasDisplayList(int list_id);
  void hideRenderer();
  void name_delete(int name_id);
  void name_replace(int name_id);
  void name_save(int name_id);
  void printString(char *s);
  void processSelection(mouseAction maction, int mk_state, int x_pos, int y_pos);
  void refresh();
  void removeDisplayLists();
  void reset();
  void resetData();
  void reshape();
  void rotate(double angle = 0.0);
  void rotate(int axis, int amount);
  void rotate(int axis, double degrees);
  void rotate(int axis, short direction);
  void scale();
  void scale(int amount);
  void scale(short direction);
  void setLightDirected();
  void setLightUndirected();
  void setParameters();
  void setWindowTitle(char* title);
  void startDisplayList(int list_id);
  void startDrawingCadBody();
  void startDrawingCadBodyElementLoop(bool is_first_loop = true);
  void startDrawingCadSurface();
  void startDrawingMeshSurface();
  void startDrawingMeshSurfaceEdges();
  void stopDisplayList(int list_id);
  void stopDrawingCadBody();
  void stopDrawingCadBodyElementLoop();
  void stopDrawingCadSurface();
  void stopDrawingMeshSurface();
  void stopDrawingMeshSurfaceEdges();
  void test();
  void translate();
  void translate(int x_amount, int y_amount, int z_amount);
  void translate(int coordinate, short direction);
  void transform_scene();
  void useDisplayList(int list_id);
protected:
  //
  void attachRenderer();
  void bodySelectionHits();
  void boundarySelectionHits();
  void boundaryEdgeSelectionHits();
  void boundaryVertexSelectionHits();
  void clean();  // Clear buffers
  void createData();
  static void createGLWindow(Hinst appInstance, const char* pName,
                             int xPos, int yPos, int pWidth, int pHeight,
                             WindowInfo& winfo);
  void deleteData();
  static void destroyWindow(Hdisplay* display, Hwindow window);
  static void Callback display(Renderer* renderer);
  void drawMeshLineElement(int elem_type, const int* node_ids, const Point3*  points, bool selected);
  void drawMeshQuadriElement(int elem_type, const int* node_ids, const Point3*  points,  const Point3*  normals, short direction, bool selected);
  void drawMeshTriangleElement(int elem_type, const int* node_ids, const Point3*  points, const Point3*  normals, short direction, bool selected);
  void drawMeshVertexElement(int elem_type, const int* node_ids, const Point3*  points, bool selected);
  static void Callback errorCallback(GLenum errorCode);
  void init();
  void makeRasterFont();
  void meshBoundaryElementSelectionHits();
  void meshBulkElementSelectionHits();
  void normalize1_point(Point3& p);
  void normalize1_point(Point3f& p);
  void normalize2_point(Point3& p);
  void normalize2_point(Point3f& p);
  static void mouse_clck_action(int mk_state, int x1, int y1);
  static void mouse_dblclck_action(int mk_state, int x1, int y1);
  static void mouse_move_action(int mk_state, int x1, int x2, int y1, int y2);
  void paintRenderer();
  void setDrawMode();
  void storeDisplayListId(int id);
  void storeTesselatorPoint(Point3* point);
  void swapBuffers();
  static LresCallback windowProcedure(Hwindow hWnd, Event wMsg,
                    Wparam wParam = NULL, Lparam lParam = NULL);

  IdList* displayListIds;
  GLuint fontOffset;
  GLParam glParam;
  static GLubyte rasterFont[][13];
  GLdouble modelStart[3];
  GLdouble modelEnd[3];
  GLdouble modelCenter[3];
  GLdouble modelLength[3];
  GLdouble modelMatrix[16];
  double** pointBuffer;
  GLdouble projectionMatrix[16];
  GLdouble projection[6];
  GLdouble rotateMatrix[16];
  int rotateAxis;
  GLdouble scaleVector[3];
  GLuint* selectionBuffer;
  GLuint selectionBufferSize;
  GLuint selectionHits;
  //GLUtesselator* tesselator;
  GLUtriangulatorObj* tesselator;
  Point3List* tesselatorPoints;
  GLdouble translateVector[3];
  GLint viewport[4];

};


#endif
