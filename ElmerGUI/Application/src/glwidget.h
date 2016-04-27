/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland    *
 *                                                                           *
 *  This program is free software; you can redistribute it and/or            *
 *  modify it under the terms of the GNU General Public License              *
 *  as published by the Free Software Foundation; either version 2           *
 *  of the License, or (at your option) any later version.                   *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program (in file fem/GPL-2); if not, write to the        *
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,         *
 *  Boston, MA 02110-1301, USA.                                              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 *  ElmerGUI glwidget                                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                   *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                 *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 15 Mar 2008                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef GLWIDGET_H
#define GLWIDGET_H

enum ListTypes {
  POINTLIST,
  EDGELIST,
  SURFACELIST,
  SURFACEMESHLIST,
  SHARPEDGELIST,
  VOLUMEMESHLIST,
  UNKNOWNLIST
};

#ifndef WIN32
#ifndef __APPLE__
#include <GL/glu.h>
#else
#include <OpenGL/glu.h>
#endif
#endif

#ifdef __MINGW32__
#include <GL/glu.h>
#endif

#include <QGLWidget>
#include <QHash>
#include <QVector>
#include "helpers.h"
#include "meshutils.h"

#define DUMMY_NAME 0xffffffff

class list_t {
 public:
  list_t();
  ~list_t();

  void setNature(int);
  int getNature() const;
  void setType(int);
  int getType() const;
  void setIndex(int);
  int getIndex() const;
  void setObject(GLuint);
  GLuint getObject() const;
  void setChild(int);
  int getChild() const;
  void setParent(int);
  int getParent() const;
  void setSelected(bool);
  bool isSelected() const;
  void setVisible(bool);
  bool isVisible() const;

 private:
  int nature;        // PDE_UNKNOWN, PDE_BOUNDARY, PDE_BULK, ...
  int type;          // POINTLIST, EDGELIST, SURFACELIST, ...
  int index;         // Boundary condition as defined in input file
  GLuint object;     // GL list index as returned by glGenLists()
  int child;         // Index to the child list (-1 = no child)
  int parent;        // Index to the parent list (-1 = no parent)
  bool selected;     // Currently selected?
  bool visible;      // Currently visible?
};

class GLWidget : public QGLWidget
{
  Q_OBJECT
    
public:
  GLWidget(QWidget *parent = 0);
  ~GLWidget();
  
  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  void setMesh(mesh_t*);
  mesh_t* getMesh() const;
  void newMesh();
  void deleteMesh();
  bool hasMesh() const;

  list_t* getList(int) const;
  int getLists() const;

  void rebuildLists();
  void rebuildSurfaceLists();
  void rebuildEdgeLists();
  void changeProjection();

  bool toggleCoordinates();

  static void indexColors(double *, int);
  static void indexColors(int *, int);
  
  // public state variables:
  bool stateOrtho;
  bool stateFlatShade;
  bool stateDrawSurfaceMesh;
  bool stateDrawVolumeMesh;
  bool stateDrawSharpEdges;
  bool stateDrawSurfaceElements;
  bool stateDrawEdgeElements;
  bool stateDrawCoordinates;
  bool stateDrawSurfaceNumbers;
  bool stateDrawEdgeNumbers;
  bool stateDrawNodeNumbers;
  bool stateDrawBoundaryIndex;
  bool stateDrawBodyIndex;
  bool stateBcColors;
  bool stateBodyColors;
  bool ctrlPressed;
  bool shiftPressed;
  bool altPressed;
  bool bodyEditActive;
  bool stateUseBgImage;
  bool stateStretchBgImage;
  bool stateAlignRightBgImage;
  QString bgImageFileName;
  int currentlySelectedBody;
  QColor backgroundColor;
  QColor surfaceColor;
  QColor edgeColor;
  QColor surfaceMeshColor;
  QColor sharpEdgeColor;
  
  // public hash tables:
  QHash<int, int> boundaryMap;
  QHash<int, int> bodyMap;

public slots:

signals:
  void signalBoundarySelected(list_t*);
  void escPressed();

protected:
  void initializeGL();
  void paintGL();
  void resizeGL(int, int);
  
  void focusInEvent(QFocusEvent*);
  void mouseDoubleClickEvent(QMouseEvent*);
  void mousePressEvent(QMouseEvent*);
  void mouseMoveEvent(QMouseEvent*);
  void wheelEvent(QWheelEvent*);
  void keyPressEvent(QKeyEvent*);
  void keyReleaseEvent(QKeyEvent*);
  
private:
  QVector<list_t*> list;

  mesh_t *mesh;
  
  Helpers *helpers;
  Meshutils *meshutils;

  GLuint makeLists();
  
  qreal matrix[16];
  qreal invmatrix[16];
  void getMatrix();
  
  QPoint lastPos;
  
  GLuint generateSurfaceList(int, QColor);
  GLuint generateSurfaceMeshList(int, QColor);
  GLuint generateVolumeMeshList(QColor);
  GLuint generateEdgeList(int, QColor);
  GLuint generateSharpEdgeList(QColor);
  
  GLUquadricObj *quadric_axis;
  void drawCoordinates();

  double drawTranslate[3];
  double drawScale;

  int bgSizeX;
  int bgSizeY;
  GLuint bgTexture;
  void drawBgImage();

  void changeNormalDirection(double*, double*);
};

#endif
