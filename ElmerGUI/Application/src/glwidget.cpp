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

#include <QtGui>
#include <QtOpenGL>
#include <QWheelEvent>
#include <QKeyEvent>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "glwidget.h"
#include "mainwindow.h"

using namespace std;

#define MY_PI 3.14159265
#define ZSHIFT -5.0

// Get qreal regardless of whether it's float or double
static inline void glGetQrealv(GLenum e, GLfloat* data) { glGetFloatv(e,data); }
static inline void glGetQrealv(GLenum e, GLdouble* data) { glGetDoublev(e,data); }
static inline void glMultMatrixq( const GLdouble *m ) { glMultMatrixd(m); }
static inline void glMultMatrixq( const GLfloat *m ) { glMultMatrixf(m); }

list_t::list_t()
{
  nature = PDE_UNKNOWN;
  type = UNKNOWNLIST;
  index = -1;
  object = 0;
  child = -1;
  parent = -1;
  selected = false;
  visible = false;
}

list_t::~list_t()
{
}

void list_t::setNature(int n)
{
  this->nature = n;
}

int list_t::getNature(void) const
{
  return this->nature;
}

void list_t::setType(int n)
{
  this->type = n;
}

int list_t::getType(void) const
{
  return this->type;
}

void list_t::setIndex(int n)
{
  this->index = n;
}

int list_t::getIndex(void) const
{
  return this->index;
}

void list_t::setObject(GLuint n)
{
  this->object = n;
}

GLuint list_t::getObject(void) const
{
  return this->object;
}

void list_t::setChild(int n)
{
  this->child = n;
}

int list_t::getChild(void) const
{
  return this->child;
}

void list_t::setParent(int n)
{
  this->parent = n;
}

int list_t::getParent(void) const
{
  return this->parent;
}

void list_t::setSelected(bool b)
{
  this->selected = b;
}

bool list_t::isSelected(void) const
{
  return this->selected;
}

void list_t::setVisible(bool b)
{
  this->visible = b;
}

bool list_t::isVisible(void) const
{
  return this->visible;
}

// Construct glWidget...
//-----------------------------------------------------------------------------
#if WITH_QT6
GLWidget::GLWidget(QWidget *parent)
  : QOpenGLWidget(parent)
#else
GLWidget::GLWidget(QWidget *parent)
  : QGLWidget(parent)
#endif
{
  backgroundColor = Qt::white;
  surfaceColor = Qt::cyan;
  edgeColor = Qt::green;
  surfaceMeshColor = Qt::black;
  sharpEdgeColor = Qt::black;
  selectionColor = Qt::red;

  stateOrtho = false;
  stateFlatShade = true;
  stateDrawSurfaceMesh = true;
  stateDrawVolumeMesh = false;
  stateDrawSharpEdges = true;
  stateDrawSurfaceElements = true;
  stateDrawEdgeElements = true;
  stateDrawCoordinates = false;
  stateDrawSurfaceNumbers = false;
  stateDrawEdgeNumbers = false;
  stateDrawNodeNumbers = false;
  stateDrawBoundaryIndex = false;
  stateDrawBodyIndex = false;
  stateBcColors = false;
  stateBodyColors = false;

  currentlySelectedBody = -1;

  drawScale = 1.0;
  drawTranslate[0] = 0.0;
  drawTranslate[1] = 0.0;
  drawTranslate[2] = 0.0;

  mesh = NULL;

  helpers = new Helpers;
  meshutils = new Meshutils;

  // Coordinate axis:
  quadric_axis = gluNewQuadric();

  // Background image:
  stateUseBgImage = false;
  stateStretchBgImage = false;
  stateAlignRightBgImage = false;
  bgImageFileName = "";
  bgTexture = 0;
  bgSizeX = 0;
  bgSizeY = 0;
}


// dtor...
//-----------------------------------------------------------------------------
GLWidget::~GLWidget()
{
}


// Min size hint...
//-----------------------------------------------------------------------------
QSize GLWidget::minimumSizeHint() const
{
  return QSize(64, 64);
}


// Default size...
//-----------------------------------------------------------------------------
QSize GLWidget::sizeHint() const
{
  return QSize(720, 576);
}

void GLWidget::setMesh(mesh_t *m)
{
  this->mesh = m;
}

mesh_t* GLWidget::getMesh(void) const
{
  return this->mesh;
}

void GLWidget::newMesh(void)
{
  this->mesh = new mesh_t;
}

void GLWidget::deleteMesh(void)
{
  delete this->mesh;
}

bool GLWidget::hasMesh(void) const
{
  if(this->mesh)
    return true;

  return false;
}

// Init GL...
//-----------------------------------------------------------------------------
void GLWidget::initializeGL()
{
  cout << "Initialize GL" << endl;
  
#if WITH_QT6
  initializeOpenGLFunctions();
#endif

  cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
  cout << "Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Version: " << glGetString(GL_VERSION) << endl;
  cout.flush();

  static GLfloat light_ambient[]  = {0.2, 0.2, 0.2, 1.0};
  static GLfloat light_diffuse[]  = {0.6, 0.6, 0.6, 1.0};
  static GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] = {0.0, 0.0, 5.0, 0.0};

  static GLfloat mat_ambient[]    = {0.2, 0.2, 0.2, 1.0};
  static GLfloat mat_diffuse[]    = {1.0, 1.0, 1.0, 1.0};
  static GLfloat mat_specular[]   = {0.9, 0.9, 0.9, 1.0};
  static GLfloat high_shininess[] = {20.0};

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0);
  glEnable(GL_LIGHTING);

  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER,1.0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);
  
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  // glDepthRange(-10.0, 10.0);

  glShadeModel(GL_SMOOTH);
  // glEnable(GL_LINE_SMOOTH);

  glEnable(GL_NORMALIZE);

#if WITH_QT6
  glClearColor( backgroundColor.redF(), backgroundColor.greenF(), backgroundColor.blueF(), backgroundColor.alphaF());
#else
  qglClearColor(backgroundColor);
#endif

  glEnable(GL_TEXTURE_2D);
}



// Paint event...
//-----------------------------------------------------------------------------
void GLWidget::paintGL()
{
#if WITH_QT6
  glClearColor( backgroundColor.redF(), backgroundColor.greenF(), backgroundColor.blueF(), backgroundColor.alphaF());
#endif
  float xabs[3], xrel[3];

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Background image:
  if(stateUseBgImage)
    drawBgImage();

  // FE objects:
  if(getLists() > 0) {

    for(int i = 0; i < getLists(); i++) {
      list_t *l = getList(i);

      if(l->isVisible()) {
	glPushName(i);

	if((l->getType() == SURFACEMESHLIST) && stateDrawSurfaceMesh) {	  
	  // translate slightly towards viewer
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, 0.01);
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject()); 
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);
	  
	} else if((l->getType() == VOLUMEMESHLIST) && stateDrawVolumeMesh) {
	  // translate slightly towards viewer
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, 0.01);
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject()); 
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);
	  
	} else if ((l->getType() == SHARPEDGELIST) && stateDrawSharpEdges) {
	  // translate slightly towards viewer
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, 0.01);
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject()); 
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);

	} else if((l->getType() == EDGELIST) && stateDrawEdgeElements ) {	  
	  // translate slightly towards viewer
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, 0.02);
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject()); 
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);

	} else if((l->getType() == SURFACELIST) && stateDrawSurfaceElements ) {
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject());
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);

	} else {
	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  glTranslated(0, 0, ZSHIFT);
	  glCallList(l->getObject()); 
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);
	}

	glPopName();
      }
    }
  }

  if(stateDrawCoordinates) {
    // push a dummy name
    glPushName(DUMMY_NAME);
    drawCoordinates();
    glPopName();
  }

  if(mesh) {
    if(stateDrawSurfaceNumbers) {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glTranslated(0, 0, ZSHIFT);
      glTranslated(0, 0, 0.1);
      glColor3d(0.5, 0, 0);
      
      for(int i=0; i < mesh->getSurfaces(); i++) {
	surface_t *surface = mesh->getSurface(i);
	int nodes = surface->getCode() / 100;
	
	xabs[0] = xabs[1] = xabs[2] = 0.0;
	
	for(int j = 0; j < nodes; j++) {
	  int ind = surface->getNodeIndex(j);
	  xabs[0] = xabs[0] + mesh->getNode(ind)->getX(0);
	  xabs[1] = xabs[1] + mesh->getNode(ind)->getX(1);
	  xabs[2] = xabs[2] + mesh->getNode(ind)->getX(2);
	}
	
	xrel[0] = (xabs[0]/nodes - drawTranslate[0]) / drawScale;
	xrel[1] = (xabs[1]/nodes - drawTranslate[1]) / drawScale;
	xrel[2] = (xabs[2]/nodes - drawTranslate[2]) / drawScale;
	
	renderText(xrel[0], xrel[1], xrel[2], QString::number(i+1) );
      }

      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
    }       
    
    if(stateDrawEdgeNumbers) {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glTranslated(0, 0, ZSHIFT);
      glTranslated(0, 0, 0.1);
      glColor3d(0.0, 0.5, 0);
      
      for(int i=0; i < mesh->getEdges(); i++) {
	edge_t *edge = mesh->getEdge(i);
	int nodes = edge->getCode() / 100;
	
	xabs[0] = xabs[1] = xabs[2] = 0.0;
	
	for(int j = 0; j < nodes; j++) {
	  int ind = edge->getNodeIndex(j);
	  xabs[0] = xabs[0] + mesh->getNode(ind)->getX(0);
	  xabs[1] = xabs[1] + mesh->getNode(ind)->getX(1);
	  xabs[2] = xabs[2] + mesh->getNode(ind)->getX(2);
	}
	xrel[0] = (xabs[0]/nodes - drawTranslate[0]) / drawScale;
	xrel[1] = (xabs[1]/nodes - drawTranslate[1]) / drawScale;
	xrel[2] = (xabs[2]/nodes - drawTranslate[2]) / drawScale;
	
	renderText(xrel[0], xrel[1], xrel[2], QString::number(i+1) );
      }

      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
    }       
    
    if(stateDrawNodeNumbers) {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glTranslated(0, 0, ZSHIFT);
      glTranslated(0, 0, 0.1);
      glColor3d(0, 0, 0.5);

      for(int i = 0; i < mesh->getNodes(); i++) {
	xabs[0] = mesh->getNode(i)->getX(0);
	xabs[1] = mesh->getNode(i)->getX(1);
	xabs[2] = mesh->getNode(i)->getX(2);
	
	xrel[0] = (xabs[0] - drawTranslate[0]) / drawScale;
	xrel[1] = (xabs[1] - drawTranslate[1]) / drawScale;
	xrel[2] = (xabs[2] - drawTranslate[2]) / drawScale;
	
	renderText(xrel[0], xrel[1], xrel[2], QString::number(i+1) );
      }

      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
    }

    if(stateDrawBoundaryIndex || stateDrawBodyIndex) {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glTranslated(0, 0, ZSHIFT);
      glTranslated(0, 0, 0.1);
      glColor3d(0.5, 0, 0);

      for(int i = 0; i < mesh->getEdges(); i++) {
	edge_t *edge = mesh->getEdge(i);
	int nodes = edge->getCode() / 100;
	
	xabs[0] = xabs[1] = xabs[2] = 0.0;
	
	for(int j = 0; j < nodes;j++) {
	  int ind = edge->getNodeIndex(j);
	  xabs[0] = xabs[0] + mesh->getNode(ind)->getX(0);
	  xabs[1] = xabs[1] + mesh->getNode(ind)->getX(1);
	  xabs[2] = xabs[2] + mesh->getNode(ind)->getX(2);
	}

	xrel[0] = (xabs[0]/nodes - drawTranslate[0]) / drawScale;
	xrel[1] = (xabs[1]/nodes - drawTranslate[1]) / drawScale;
	xrel[2] = (xabs[2]/nodes - drawTranslate[2]) / drawScale;
	
	if(stateDrawBoundaryIndex && (edge->getNature() == PDE_BOUNDARY))
	  renderText(xrel[0], xrel[1], xrel[2], QString::number(edge->getIndex()));

	if(stateDrawBodyIndex && (edge->getNature() == PDE_BULK))
	  renderText(xrel[0], xrel[1], xrel[2], QString::number(edge->getIndex()));
      }
      
      for(int i = 0; i < mesh->getSurfaces(); i++) {
	surface_t *surface = mesh->getSurface(i);
	int nodes = surface->getCode() / 100;

	xabs[0] = xabs[1] = xabs[2] = 0.0;
	
	for(int j = 0; j < nodes; j++) {
	  int ind = surface->getNodeIndex(j);
	  xabs[0] = xabs[0] + mesh->getNode(ind)->getX(0);
	  xabs[1] = xabs[1] + mesh->getNode(ind)->getX(1);
	  xabs[2] = xabs[2] + mesh->getNode(ind)->getX(2);
	}

	xrel[0] = (xabs[0]/nodes - drawTranslate[0]) / drawScale;
	xrel[1] = (xabs[1]/nodes - drawTranslate[1]) / drawScale;
	xrel[2] = (xabs[2]/nodes - drawTranslate[2]) / drawScale;
	
	if(stateDrawBoundaryIndex && (surface->getNature() == PDE_BOUNDARY))
	  renderText(xrel[0], xrel[1], xrel[2], QString::number(surface->getIndex()));

	if(stateDrawBodyIndex && (surface->getNature() == PDE_BULK))
	  renderText(xrel[0], xrel[1], xrel[2], QString::number(surface->getIndex()));

	// case 3d:
	if(stateDrawBodyIndex && (surface->getNature() == PDE_BOUNDARY)) {
	  for(int i = 0; i < surface->getElements(); i++) {
	    int j = surface->getElementIndex(i);
	    if(j >= 0) {
	      element_t *element = mesh->getElement(j);
	      renderText(xrel[0], xrel[1], xrel[2], QString::number(element->getIndex()));
	    }
	  }
	}

      }

      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
    }
  }
}


// Change projection...
//-----------------------------------------------------------------------------
void GLWidget::changeProjection()
{
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  int width = viewport[2];
  int height = viewport[3];
  double top = 1.0;
  double bottom = -1.0;
  double left = -(double)width / (double)height;
  double right = (double)width / (double)height;
  double _near = -10.0;
  double _far = 10.0;

  if(stateOrtho) {
    glViewport(0, 0, (GLint)width, (GLint)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, _near, _far);
    glMatrixMode(GL_MODELVIEW);
  } else {
    glViewport(0, 0, (GLint)width, (GLint)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(54.0, (float)width/(float)height, 0.1, 10.0);
    glMatrixMode(GL_MODELVIEW);
  }
}

// Resize window...
//-----------------------------------------------------------------------------
void GLWidget::resizeGL(int width, int height)
{
  double top = 1.0;
  double bottom = -1.0;
  double left = -(double)width / (double)height;
  double right = (double)width / (double)height;
  double _near = -10.0;
  double _far = 10.0;

  if(stateOrtho) {
    glViewport(0, 0, (GLint)width, (GLint)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, _near, _far);
    glMatrixMode(GL_MODELVIEW);
  } else {
    glViewport(0, 0, (GLint)width, (GLint)height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (float)width/(float)height, 0.1, 10.0);
    glMatrixMode(GL_MODELVIEW);  
  }
}


// Focus in event...
//-----------------------------------------------------------------------------
void GLWidget::focusInEvent(QFocusEvent *event)
{
  Q_UNUSED(event)

  // Should we check the key pressed status here?
}


// Key pressed...
//-----------------------------------------------------------------------------
void GLWidget::keyPressEvent(QKeyEvent *event)
{
  if(event->key() == Qt::Key_Escape)
    emit(escPressed());
}


// Key released...
//-----------------------------------------------------------------------------
void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
}



// Mouse button pressed...
//-----------------------------------------------------------------------------
void GLWidget::mousePressEvent(QMouseEvent *event)
{
  lastPos = event->pos();
  lastPressPos = event->pos();
  setFocus();  // for tracing keyboard events
}

// Mouse button released...
//-----------------------------------------------------------------------------
void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
  if(event->button() == Qt::RightButton & event->pos() == lastPressPos){
#if WITH_QT6
    ((MainWindow*)parent())->showContextMenu(event->globalPosition().toPoint());
#else
    ((MainWindow*)parent())->showContextMenu(event->globalPos());
#endif	
  }
}

// Mouse wheel rotates...
//-----------------------------------------------------------------------------
void GLWidget::wheelEvent(QWheelEvent *event)
{
#if WITH_QT6
  makeCurrent();
  double s = exp((double)(event->angleDelta().y())*0.001);
#else
  double s = exp((double)(event->delta())*0.001);
#endif
  glScaled(s, s, s);
  updateGL();
#if WITH_QT6  
  lastPos = event->position().toPoint();
#else
  lastPos = event->pos();
#endif
  getMatrix();
}



// Mouse moves...
//-----------------------------------------------------------------------------
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
#if WITH_QT6
  makeCurrent();
#endif

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

#if WITH_QT6
  int dx = event->position().x() - lastPos.x();
  int dy = event->position().y() - lastPos.y();
#else
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();
#endif

  dy = -dy;
  
  if (
#if WITH_QT6
  ((event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::MiddleButton))
#else
  ((event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::MidButton))
#endif
        ||
    event->buttons() == Qt::RightButton  // added for easy scaling
       ) {

    // Scale:
    double s = exp(dy*0.01);
    glScaled(s, s, s);
    updateGL();

  } else if (event->buttons() == Qt::LeftButton) {
    
    // Rotation:
    double ax = -(double)dy;
    double ay =  (double)dx;
    double az = 0.0;

    double s = 180.0*sqrt(ax*ax+ay*ay+az*az)/(double)(viewport[3]+1);
    double bx = invmatrix[0]*ax + invmatrix[4]*ay + invmatrix[8]*az;
    double by = invmatrix[1]*ax + invmatrix[5]*ay + invmatrix[9]*az;
    double bz = invmatrix[2]*ax + invmatrix[6]*ay + invmatrix[10]*az;
    glRotated(s, bx, by, bz);
    updateGL();

  } else if (
#if WITH_QT6
  (event->buttons() == Qt::MiddleButton)
#else
  (event->buttons() == Qt::MidButton)
#endif
  
  ||
  (event->buttons() ==  (Qt::LeftButton | Qt::RightButton)) // added for 2 button mouse
    ){

    // Translation:
    double s = 2.0/(double)(viewport[3]+1);
    double ax = s*dx;
    double ay = s*dy;
    double az = 0.0;
    glLoadIdentity();
    glTranslated(ax, ay, az);
    glMultMatrixq(matrix);
    updateGL();
  }

  lastPos = event->pos();
  getMatrix();
}



// Mouse button double clicked...
//-----------------------------------------------------------------------------
void GLWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
#if WITH_QT6
  makeCurrent();
#endif
  if(getLists() == 0) 
    return;
  
  /*
  To avoid segmentation fault in MSYS2 environment, compass, numbers and indexes are hidden.
  These will be restored at the end of this function. Do not return before restoring these.
  */
  bool prevStateDrawCoordinates = stateDrawCoordinates;
  bool prevStateDrawSurfaceNumbers = stateDrawSurfaceNumbers;
  bool prevStateDrawEdgeNumbers = stateDrawEdgeNumbers;
  bool prevStateDrawNodeNumbers = stateDrawNodeNumbers;
  bool prevStateDrawBoundaryIndex = stateDrawBoundaryIndex;
  bool prevStateDrawBodyIndex = stateDrawBodyIndex;  
  stateDrawCoordinates = false;
  stateDrawSurfaceNumbers = false;
  stateDrawEdgeNumbers = false;
  stateDrawNodeNumbers = false;
  stateDrawBoundaryIndex = false;
  stateDrawBodyIndex = false;
  
  
  static list_t dummylist;
  static GLuint buffer[1024];
  const int bufferSize = sizeof(buffer)/sizeof(GLuint);
  
  GLint viewport[4];
  GLdouble projection[16];

  GLint hits;
  GLint i, j;

  updateGL();
#if WITH_QT6
  makeCurrent();
#endif
  
  glSelectBuffer(bufferSize, buffer);
  glRenderMode(GL_SELECT);
  glInitNames();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glLoadIdentity();

#if WITH_QT6
  GLdouble x = event->position().x();
  GLdouble y = (double)viewport[3]-event->position().y()-1;
#else
  GLdouble x = event->x();
  GLdouble y = (double)viewport[3]-event->y()-1;
#endif  
  
  GLdouble deltaX = 3.0;
  GLdouble deltaY = 3.0;
  
  gluPickMatrix(x, y, deltaX, deltaY, viewport);
  glMultMatrixd(projection); 
  
  glMatrixMode(GL_MODELVIEW);


  /*This is to avoid segmentation fault in Linux with old hardware*/
  setMeshVisibility(false, false, false);
  
#if WITH_QT6
  paintGL();
#else
  updateGL();
#endif

  /*Again, this is to avoid segmentation fault in Linux with old hardware and to suppress blinking*/
  setMeshVisibility(stateDrawSurfaceMesh, stateDrawVolumeMesh, stateDrawSharpEdges); 

  hits = glRenderMode(GL_RENDER);
  
  bool badDriver = true;
  GLuint smallestz = DUMMY_NAME;
  GLuint nearest = DUMMY_NAME;

  if(hits != 0) {
    for (i=0, j=0; i<hits; i++) {
      GLuint minz = buffer[j+1];
      GLuint resultz = buffer[j+3];

      badDriver = (badDriver && (minz == 0x80000000));

      if(minz < smallestz) {
	nearest = resultz;
	smallestz = minz;
      }

      j += 3 + buffer[j];
    }
  }
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  if(badDriver) {
    cerr << "Detected bad GL-context or broken graphics driver" << endl;
    cerr.flush();
    cout << "glRenderMode(GL_RENDER) produces bad z-values" << endl;
    cout << "Unable to reliably select objects" << endl;
    cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
    cout << "Renderer: " << glGetString(GL_RENDERER) << endl;
    cout << "Version: " << glGetString(GL_VERSION) << endl;
    cout.flush();
  }
  
  // highlight the selected boundary:
  if(nearest != DUMMY_NAME) {
    list_t *l = getList(nearest);

    // skip sharp edge lists
   if(l->getType() == SHARPEDGELIST) {

	  /* 
	  Restoration of view settings. These were adjusted at the beginning of this function.
	  */
	  stateDrawCoordinates = prevStateDrawCoordinates;
	  stateDrawSurfaceNumbers = prevStateDrawSurfaceNumbers;
	  stateDrawEdgeNumbers = prevStateDrawEdgeNumbers;
	  stateDrawNodeNumbers = prevStateDrawNodeNumbers;
	  stateDrawBoundaryIndex = prevStateDrawBoundaryIndex;
	  stateDrawBodyIndex = prevStateDrawBodyIndex;  
	  

      updateGL();
	  return;
	}
    
    // substitute surfacemeshlist with the parent surfacelist:
    if(l->getType() == SURFACEMESHLIST)
      l = getList(l->getParent());

    // if not ctrl pressed, rebuild all selected lists except this one:
    if(!(event->modifiers() & Qt::ControlModifier)) {
      for(i = 0; i < getLists(); i++) {
	list_t *l2 = getList(i);
	if(l2->isSelected() && (l2->getIndex() != l->getIndex())) {
	  glDeleteLists(l2->getObject(), 1);
	  l2->setSelected(false);
	  if(l2->getType() == SURFACELIST) {
            for( int j = 0; j < mesh->getSurfaces(); j++ ) {
              surface_t *surf = mesh->getSurface(j);
              if ( surf->getIndex() == l2->getIndex() )
                surf->setSelected(l2->isSelected());
            }
	    l2->setObject(generateSurfaceList(l2->getIndex(), surfaceColor)); // cyan
	  } else if(l2->getType() == EDGELIST) {
            for( int j=0; j < mesh->getEdges(); j++ ) {
              edge_t *edge = mesh->getEdge(j);
              if ( edge->getIndex() == l2->getIndex() )
                edge->setSelected(l2->isSelected());
            }
	    l2->setObject(generateEdgeList(l2->getIndex(), edgeColor)); // green
	  }
	}
      }
    }

    // Toggle selection:
    l->setSelected(!l->isSelected());
    
    glDeleteLists(l->getObject(), 1);

    // Highlight current selection:
    if(l->getType() == SURFACELIST) {
      if(l->isSelected()) {
	l->setObject(generateSurfaceList(l->getIndex(), selectionColor)); // red
      } else {
	l->setObject(generateSurfaceList(l->getIndex(), surfaceColor)); // cyan
      }

      for( int i=0; i<mesh->getSurfaces(); i++ ) {
        surface_t *surf = mesh->getSurface(i);
        if ( surf->getIndex() == l->getIndex() ) surf->setSelected(l->isSelected());
      }

    } else if(l->getType() == EDGELIST) {
      if(l->isSelected()) {
	l->setObject(generateEdgeList(l->getIndex(), selectionColor)); // red
      } else {
	l->setObject(generateEdgeList(l->getIndex(), edgeColor)); // green
      }
      for( int i=0; i < mesh->getEdges(); i++ ) {
        edge_t *edge = mesh->getEdge(i);
        if ( edge->getIndex() == l->getIndex() ) edge->setSelected(l->isSelected());
      }
    }

    // body selection:
    //----------------
    currentlySelectedBody = -1;
    if( (event->modifiers() & Qt::ShiftModifier) || bodyEditActive) {

      // determine the max bulk index
      int MAX_BULK_INDEX = -1;

      for(int i = 0; i < mesh->getElements(); i++) {
	element_t *elem = mesh->getElement(i);
	if(elem->getNature() != PDE_BULK)
	  break;
	if(elem->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = elem->getIndex();
      }

      for(int i = 0; i < mesh->getSurfaces(); i++) {
	surface_t *surf = mesh->getSurface(i);
	if(surf->getNature() != PDE_BULK)
	  break;
	if(surf->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = surf->getIndex();
      }

      for(int i = 0; i < mesh->getEdges(); i++) {
	edge_t *edge = mesh->getEdge(i);
	if(edge->getNature() != PDE_BULK)
	  break;
	if(edge->getIndex() > MAX_BULK_INDEX)
	  MAX_BULK_INDEX = edge->getIndex();
      }
      
      MAX_BULK_INDEX++;
      if(MAX_BULK_INDEX == 0) {
	cout << "Error in body selection: "
	  "There are no legal body indices from which to choose" << endl;
	cout.flush();
	goto body_selection_finished;
      }

      // allocate temp arrays:
      bool *tmp1 = new bool[MAX_BULK_INDEX];
      bool *tmp2 = new bool[MAX_BULK_INDEX];

      for(int i = 0; i < MAX_BULK_INDEX; i++) {
	tmp1[i] = true;
	tmp2[i] = false;
      }
      
      // check if the selected lists uniquely determine a bulk body:
      for(int i = 0; i < getLists(); i++) {
	list_t *l2 = getList(i);

	if(l2->isSelected() && (l2->getNature() == PDE_BULK)) {
	  for(int j = 0; j < MAX_BULK_INDEX; j++) {
	    if(j != l2->getIndex())
	      tmp1[j] = false;
	  }
	}
	
	if(l2->isSelected() && (l2->getNature() == PDE_BOUNDARY) && 
	   (l2->getType() == SURFACELIST)) {	  
	  for(int j = 0; j < mesh->getSurfaces(); j++) {
	    surface_t *surf = mesh->getSurface(j);
	    if(surf->getIndex() == l2->getIndex()) {
	      for(int k = 0; k < surf->getElements(); k++) {
		int l = surf->getElementIndex(k);
		if(l < 0) 
		  break;
		element_t *elem = mesh->getElement(l);
		if((elem->getIndex() < 0) || (elem->getIndex() >= MAX_BULK_INDEX))
		  break;
		tmp2[elem->getIndex()] = true;
	      }
	      for(int k = 0; k < MAX_BULK_INDEX; k++) {
		tmp1[k] &= tmp2[k];
		tmp2[k] = false;
	      }
	    }
	  }
	}
      }

      // array "tmp1" should contain only one entry with value "true"
      int count = 0;
      int found = -1;
      for(int i = 0; i < MAX_BULK_INDEX; i++) {
	if( tmp1[i] ) {
	  count++;
	  found = i;
	}
//	cout << i << ": tmp1,2= " << tmp1[i] << "," << tmp2[i] << endl; 
      }
      
      if((count == 1) && (found >= 0))
	currentlySelectedBody = found;
    else if((count > 1)){
	  int m = mostVisibleBody(MAX_BULK_INDEX, tmp1);
	  if(m >=0) currentlySelectedBody = m;
	}
     
      delete [] tmp1;
      delete [] tmp2;
    }
  body_selection_finished:
    
    // Emit result to mainwindow:
    emit(signalBoundarySelected(l, event->modifiers()));

  } else {

    // Emit "nothing selected":
    dummylist.setNature(-1);
    dummylist.setType(-1);
    dummylist.setIndex(-1);
    emit(signalBoundarySelected(&dummylist, event->modifiers()));

  }

  /* 
  Restoration of view settings. These were adjusted at the beginning of this function.
  */
  stateDrawCoordinates = prevStateDrawCoordinates;
  stateDrawSurfaceNumbers = prevStateDrawSurfaceNumbers;
  stateDrawEdgeNumbers = prevStateDrawEdgeNumbers;
  stateDrawNodeNumbers = prevStateDrawNodeNumbers;
  stateDrawBoundaryIndex = prevStateDrawBoundaryIndex;
  stateDrawBodyIndex = prevStateDrawBodyIndex;  
  
  
  updateGL();
}



// Get current matrix and its inverse...
//-----------------------------------------------------------------------------
void GLWidget::getMatrix()
{
  glGetQrealv(GL_MODELVIEW_MATRIX, matrix);
  helpers->invertMatrix(matrix, invmatrix);
}



// Rebuild lists...
//-----------------------------------------------------------------------------
void GLWidget::rebuildLists()
{
#if WITH_QT6
  makeCurrent();
#endif
  double *bb = mesh->boundingBox();
  
  drawTranslate[0] = bb[6]; // x-center
  drawTranslate[1] = bb[7]; // y-center
  drawTranslate[2] = bb[8]; // z-center
  drawScale = bb[9];        // scaling

  delete [] bb;
 
  if(getLists() > 0) {
    for(int i=0; i < getLists(); i++) {
      list_t *l = getList(i);
      glDeleteLists(l->getObject(), 1);
    }
  }

  makeLists();

  updateGL();
}

// Compose GL surface lists...
//-----------------------------------------------------------------------------
void GLWidget::rebuildSurfaceLists()
{
  for( int i = 0; i < getLists(); i++ )
  {
    list_t *l = getList(i);
     if( l->getType() == SURFACELIST )
     {
       glDeleteLists( l->getObject(), 1 );
       if(l->isSelected()) {
 	 l->setObject(generateSurfaceList(l->getIndex(), selectionColor)); // red
       } else {
 	 l->setObject(generateSurfaceList(l->getIndex(), surfaceColor)); // cyan
       }
     }
  }
}

// Compose GL edge lists...
//-----------------------------------------------------------------------------
void GLWidget::rebuildEdgeLists()
{
  for( int i = 0; i < getLists(); i++ )
  {
    list_t *l = getList(i);
     if ( l->getType() == EDGELIST )
     {
       glDeleteLists( l->getObject(), 1 );
       if(l->isSelected()) {
 	 l->setObject(generateEdgeList(l->getIndex(), selectionColor)); // red
       } else {
 	 l->setObject(generateEdgeList(l->getIndex(), edgeColor)); // green
       }
     }
  }
}



// Compose GL object lists...
//-----------------------------------------------------------------------------
GLuint GLWidget::makeLists()
{
#if WITH_QT6
  makeCurrent();
#endif

  int i;
  list_t *l;

  if((mesh == NULL) || mesh->isUndefined())
    return 0;

  // The rule for composing lists to display is the following:
  //---------------------------------------------------------------------------
  // - All surface elements with index >= 0 will be drawn - one list/index
  //   (list->type = SURFACELIST)
  // - For each surface element list, one auxiliary list will be drawn
  //   (list->type = SURFACEMESHLIST)
  // - All edge elements with index >= 0 will be drawn - one list/index
  //   (list->type = EDGELIST)
  // - All point elements with index >= 0 will be drawn - one list/index
  //   (list->type = POINTLIST)
  // - A list of sharp edges will always be drawn (even if it is empty)
  //---------------------------------------------------------------------------
  
  // Simultaneously, populate hash for mapping body & boundary indices:
  //--------------------------------------------------------------------
  boundaryMap.clear();
  bodyMap.clear();
  int boundaryCount = 0;
  int bodyCount = 0;

  // Scan volume elements to determine the number of material indices:
  //-------------------------------------------------------------------
  QMap<int, int> bodyNatures;

  for(i = 0; i < mesh->getElements(); i++) {
    element_t *element = mesh->getElement(i);
    int index = element->getIndex();

    if(index >= 0) {
      int nature = element->getNature();

      if(!bodyNatures.contains(index))
	bodyNatures.insert(index, nature);
    }
  }    

  QMapIterator<int, int> itrBody(bodyNatures);
  while (itrBody.hasNext()) {
    itrBody.next();
    int index = itrBody.key();
    int nature = itrBody.value();

    if(nature == PDE_BULK) 
      bodyMap.insert(index, bodyCount++);
  }

  // Scan surface elements to determine the number of bcs. / mat. indices:
  //-----------------------------------------------------------------------
  int surface_bcs = 0;

  QMap<int, int> surfaceNatures;

  for(i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);
    int index = surface->getIndex();

    if(index > 0) {
      int nature = surface->getNature();

      if(!surfaceNatures.contains(index))
	surfaceNatures.insert(index, nature);
    }
  }    

  QMapIterator<int, int> itrSurface(surfaceNatures);
  while (itrSurface.hasNext()) {
    itrSurface.next();
    int index = itrSurface.key();
    int nature = itrSurface.value();

    if(nature > 0) {
      surface_bcs++;

      if(nature == PDE_BULK)
	bodyMap.insert(index, bodyCount++);
      
      if(nature == PDE_BOUNDARY)
	boundaryMap.insert(index, boundaryCount++);
    }
  }

  cout << "Bcs/materials on surface elements: " << surface_bcs << endl;
  cout.flush();

  // Scan edge elements to determine the number of bcs. / mat. indices:
  //--------------------------------------------------------------------
  int edge_bcs = 0;
  
  QMap<int, int> edgeNatures;
  
  for(i = 0; i < mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);
    int index = edge->getIndex();
    
    if(index > 0) {
      int nature = edge->getNature();

      if(!edgeNatures.contains(index))
	edgeNatures.insert(index, nature);
    }
  }    

  QMapIterator<int, int> itrEdge(edgeNatures);
  while (itrEdge.hasNext()) {
    itrEdge.next();
    int index = itrEdge.key();
    int nature = itrEdge.value();

    if(nature > 0) {
      edge_bcs++;

      if(nature == PDE_BULK)
	bodyMap.insert(index, bodyCount++);
      
      if(nature == PDE_BOUNDARY)
	boundaryMap.insert(index, boundaryCount++);
    }
  }

  cout << "Bcs/materials on edge elements: " << edge_bcs << endl;  
  cout.flush();

  // Scan point elements to determine the number of bcs. / mat. indices:
  //---------------------------------------------------------------------
  int point_bcs = 0;

  // TODO

  cout << "Bcs/materials on point elements: " << point_bcs << endl;  
  cout.flush();

  // Generate lists:
  //---------------------------------------------------------------------
  for(i = 0; i < getLists(); i++)
    delete list.at(i);

  list.clear();

  cout << "Generating  lists to display" << endl;
  cout.flush();

  // Surface lists:
  //----------------
  for(i = 0; i < mesh->getSurfaces(); i++)
    mesh->getSurface(i)->setSelected(false);

  itrSurface.toFront();
  while (itrSurface.hasNext()) {
    itrSurface.next();
    int index = itrSurface.key();
    int nature = itrSurface.value();

    if(nature > 0) {
      l = new list_t;
      list.push_back(l);

      l->setNature(nature);
      l->setType(SURFACELIST);
      l->setIndex(index);
      l->setObject(generateSurfaceList(l->getIndex(), surfaceColor)); // cyan
      l->setChild(getLists());
      l->setParent(-1);
      l->setSelected(false);
      l->setVisible(stateDrawSurfaceElements);

      // edges of surface elements (just for visual):
      l = new list_t;
      list.push_back(l);

      l->setNature(PDE_UNKNOWN);
      l->setType(SURFACEMESHLIST);
      l->setIndex(index);
      l->setObject(generateSurfaceMeshList(l->getIndex(), surfaceMeshColor)); // black
      l->setChild(-1);
      l->setParent(getLists() - 2);
      l->setSelected(false);
      l->setVisible(stateDrawSurfaceMesh);
    }
  }
  
  // Edge lists (only PDE_BOUNDARY):
  //---------------------------------
  for(i = 0; i < mesh->getEdges(); i++)
    mesh->getEdge(i)->setSelected(false);
  
  itrEdge.toFront();
  while (itrEdge.hasNext()) {
    itrEdge.next();
    int index = itrEdge.key();
    int nature = itrEdge.value();
    
    if(nature > 0) {
      l = new list_t;
      list.push_back(l);

      l->setNature(nature); 
      l->setType(EDGELIST);
      l->setIndex(index);
      l->setObject(generateEdgeList(l->getIndex(), edgeColor)); // green
      l->setChild(-1);
      l->setParent(-1);
      l->setSelected(false);
      l->setVisible(stateDrawEdgeElements);
    }
  }

  // Point lists: TODO

  // Sharp edges (just for visual):
  //--------------------------------
  l = new list_t;
  list.push_back(l);

  l->setNature(PDE_UNKNOWN);
  l->setType(SHARPEDGELIST);
  l->setIndex(-1);
  l->setObject(generateSharpEdgeList(sharpEdgeColor)); // black
  l->setChild(-1);
  l->setParent(-1);
  l->setSelected(false);
  l->setVisible(stateDrawSharpEdges);

  // Volume mesh (visual only):
  //----------------------------
  l = new list_t;
  list.push_back(l);

  l->setNature(PDE_UNKNOWN);
  l->setType(VOLUMEMESHLIST);
  l->setIndex(-1);
  l->setObject(generateVolumeMeshList(Qt::black)); // black
  l->setChild(-1);
  l->setParent(-1);
  l->setSelected(false);
  l->setVisible(stateDrawVolumeMesh);

  // Clean up:
  //-----------
  edgeNatures.clear();
  surfaceNatures.clear();
  bodyNatures.clear();  

  updateGL();
  getMatrix();

  cout << "Generated " << getLists() << " lists" << endl;
  cout.flush();

  return getLists();
}


// Generate volume mesh list...
//-----------------------------------------------------------------------------
GLuint GLWidget::generateVolumeMeshList(QColor qColor)
{
  static int tetmap[6][2] = {{0, 1}, {0, 2}, {0, 3}, 
			     {1, 2}, {1, 3}, {2, 3}};

  static int wedgemap[9][2] = {{0, 1}, {1, 2}, {2, 0},
			       {0, 3}, {1, 4}, {2, 5},
			       {3, 4}, {4, 5}, {5, 3}};


  static int hexmap[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
			      {0, 4}, {1, 5}, {2, 6}, {3, 7},
			      {4, 5}, {5, 6}, {6, 7}, {7, 4}};

  double R = qColor.red() / 255.0;
  double G = qColor.green() / 255.0;
  double B = qColor.blue() / 255.0;

  GLuint current = glGenLists(1);
  glNewList(current, GL_COMPILE);

  glBegin(GL_LINES);

  for(int i = 0; i < mesh->getElements(); i++) {
    element_t *element = mesh->getElement(i);

    glColor3d(R, G, B);

    int nofEdges = 0;
    int *edgeMap = 0;

    switch((int)(element->getCode() / 100)) {
    case 5:
      nofEdges = 6;
      edgeMap = &tetmap[0][0];
      break;
    case 7:
      nofEdges = 9;
      edgeMap = &wedgemap[0][0];
      break;
    case 8:
      nofEdges = 12;
      edgeMap = &hexmap[0][0];
      break;
    }
    
    // draw edges:
    for(int j = 0; j < nofEdges; j++) {
      int p0 = *edgeMap++;
      int p1 = *edgeMap++;
      
      int q0 = element->getNodeIndex(p0);
      int q1 = element->getNodeIndex(p1);
      
      node_t *n0 = mesh->getNode(q0);
      node_t *n1 = mesh->getNode(q1);
      
      double x0 = ( n0->getX(0) - drawTranslate[0] ) / drawScale;
      double y0 = ( n0->getX(1) - drawTranslate[1] ) / drawScale;
      double z0 = ( n0->getX(2) - drawTranslate[2] ) / drawScale;
      
      double x1 = ( n1->getX(0) - drawTranslate[0] ) / drawScale;
      double y1 = ( n1->getX(1) - drawTranslate[1] ) / drawScale;
      double z1 = ( n1->getX(2) - drawTranslate[2] ) / drawScale;
      
      glVertex3d(x0, y0, z0);
      glVertex3d(x1, y1, z1);
    }
  }

  glEnd();

  glEndList();

  return current;
}



// Generate surface list...
//-----------------------------------------------------------------------------
GLuint GLWidget::generateSurfaceList(int index, QColor qColor)
{
  double x0[3], x1[3], x2[3], x3[3], u[3];

  double R = qColor.red() / 255.0;
  double G = qColor.green() / 255.0;
  double B = qColor.blue() / 255.0;

  GLuint current = glGenLists(1);
  glNewList(current, GL_COMPILE);

  // Draw triangles:
  //-----------------
  glBegin(GL_TRIANGLES);

  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);

    if((surface->getIndex() == index) && ((int)(surface->getCode() / 100) == 3)) {

      glColor3d(R, G, B);

      if(stateBcColors && (surface->getNature() == PDE_BOUNDARY)) {
        double c[3];
        indexColors(c, index);
        glColor3d(c[0], c[1], c[2]);
      } 
	
      if(stateBodyColors) {
	int bodyIndex = surface->getIndex();
	if(surface->getNature() == PDE_BOUNDARY) {
	  int parentIndex = surface->getElementIndex(0);
	  element_t *parent = mesh->getElement(parentIndex);
	  bodyIndex = parent->getIndex();
	}
        double c[3];
        indexColors(c, bodyIndex);
        glColor3d(c[0], c[1], c[2]);
      } 
      
      // change normal direction:
      changeNormalDirection(u, surface->getNormalVec());
      glNormal3dv(u); 
      
      int n0 = surface->getNodeIndex(0);
      int n1 = surface->getNodeIndex(1);
      int n2 = surface->getNodeIndex(2);
      
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
      
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
      
      x2[0] = (mesh->getNode(n2)->getX(0) - drawTranslate[0]) / drawScale;
      x2[1] = (mesh->getNode(n2)->getX(1) - drawTranslate[1]) / drawScale;
      x2[2] = (mesh->getNode(n2)->getX(2) - drawTranslate[2]) / drawScale;

      changeNormalDirection(u, surface->getVertexNormalVec(0));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x0);

      changeNormalDirection(u, surface->getVertexNormalVec(1));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x1);

      changeNormalDirection(u, surface->getVertexNormalVec(2));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x2);
    }
  }

  glEnd();

  // Draw quads:
  //------------
  glBegin(GL_QUADS);
  
  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);

    if((surface->getIndex() == index) && ((int)(surface->getCode() / 100) == 4)) {

      glColor3d(R, G, B);

      if(stateBcColors && (surface->getNature() == PDE_BOUNDARY)) {
        double c[3];
        indexColors(c, index);
        glColor3d(c[0], c[1], c[2]);
      }

      if(stateBodyColors) {
	int bodyIndex = surface->getIndex();
	if(surface->getNature() == PDE_BOUNDARY) {
	  int parentIndex = surface->getElementIndex(0);
	  element_t *parent = mesh->getElement(parentIndex);
	  bodyIndex = parent->getIndex();
	}
        double c[3];
        indexColors(c, bodyIndex);
        glColor3d(c[0], c[1], c[2]);
      } 

      // change normal direction:
      changeNormalDirection(u, surface->getNormalVec());
      glNormal3dv(u); 
      
      int n0 = surface->getNodeIndex(0);
      int n1 = surface->getNodeIndex(1);
      int n2 = surface->getNodeIndex(2);
      int n3 = surface->getNodeIndex(3);
      
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
      
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
      
      x2[0] = (mesh->getNode(n2)->getX(0) - drawTranslate[0]) / drawScale;
      x2[1] = (mesh->getNode(n2)->getX(1) - drawTranslate[1]) / drawScale;
      x2[2] = (mesh->getNode(n2)->getX(2) - drawTranslate[2]) / drawScale;
      
      x3[0] = (mesh->getNode(n3)->getX(0) - drawTranslate[0]) / drawScale;
      x3[1] = (mesh->getNode(n3)->getX(1) - drawTranslate[1]) / drawScale;
      x3[2] = (mesh->getNode(n3)->getX(2) - drawTranslate[2]) / drawScale;
      
      changeNormalDirection(u, surface->getVertexNormalVec(0));
      if ( !stateFlatShade ) glNormal3dv(u); 
      glVertex3dv(x0);

      changeNormalDirection(u, surface->getVertexNormalVec(1));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x1);

      changeNormalDirection(u, surface->getVertexNormalVec(2));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x2);

      changeNormalDirection(u, surface->getVertexNormalVec(3));
      if ( !stateFlatShade ) glNormal3dv(u);
      glVertex3dv(x3);      
    }
  }

  glEnd();
  glEndList();

  return current;
}



// Generate surface edge list...
//-----------------------------------------------------------------------------
GLuint GLWidget::generateSurfaceMeshList(int index, QColor qColor)
{
  double x0[3], x1[3], x2[3], x3[3];

  double R = qColor.red() / 255.0;
  double G = qColor.green() / 255.0;
  double B = qColor.blue() / 255.0;

  GLuint current = glGenLists(1);
  glNewList(current, GL_COMPILE);
 
  // Draw lines:
  //------------
  glLineWidth(1.0);
  glDisable(GL_LIGHTING);
  glColor3d(R, G, B);
  glBegin(GL_LINES);
  
  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);
    
    if((surface->getIndex() == index) && ((int)(surface->getCode() / 100) == 3)) {
      int n0 = surface->getNodeIndex(0);
      int n1 = surface->getNodeIndex(1);
      int n2 = surface->getNodeIndex(2);
      
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
      
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
      
      x2[0] = (mesh->getNode(n2)->getX(0) - drawTranslate[0]) / drawScale;
      x2[1] = (mesh->getNode(n2)->getX(1) - drawTranslate[1]) / drawScale;
      x2[2] = (mesh->getNode(n2)->getX(2) - drawTranslate[2]) / drawScale;
      
      glVertex3dv(x0);
      glVertex3dv(x1);

      glVertex3dv(x1);
      glVertex3dv(x2);

      glVertex3dv(x2);
      glVertex3dv(x0);
    }
  }
  
  for(int i=0; i < mesh->getSurfaces(); i++) {
    surface_t *surface = mesh->getSurface(i);

    if((surface->getIndex() == index) && ((int)(surface->getCode() / 100) == 4)) {
      int n0 = surface->getNodeIndex(0);
      int n1 = surface->getNodeIndex(1);
      int n2 = surface->getNodeIndex(2);
      int n3 = surface->getNodeIndex(3);
      
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
      
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
      
      x2[0] = (mesh->getNode(n2)->getX(0) - drawTranslate[0]) / drawScale;
      x2[1] = (mesh->getNode(n2)->getX(1) - drawTranslate[1]) / drawScale;
      x2[2] = (mesh->getNode(n2)->getX(2) - drawTranslate[2]) / drawScale;
      
      x3[0] = (mesh->getNode(n3)->getX(0) - drawTranslate[0]) / drawScale;
      x3[1] = (mesh->getNode(n3)->getX(1) - drawTranslate[1]) / drawScale;
      x3[2] = (mesh->getNode(n3)->getX(2) - drawTranslate[2]) / drawScale;
      
      glVertex3dv(x0);
      glVertex3dv(x1);

      glVertex3dv(x1);
      glVertex3dv(x2);

      glVertex3dv(x2);
      glVertex3dv(x3);

      glVertex3dv(x3);
      glVertex3dv(x0);      
    }
  }
  glEnd();

  glEnable(GL_LIGHTING);
  glEndList();
  
  return current;
}


// Generate edge list...
//-----------------------------------------------------------------------------
GLuint GLWidget::generateEdgeList(int index, QColor qColor)
{
  double x0[3], x1[3];

  double R = qColor.red() / 255.0;
  double G = qColor.green() / 255.0;
  double B = qColor.blue() / 255.0;

  GLuint current = glGenLists(1);
  glNewList(current, GL_COMPILE);
  glColor3d(R, G, B);  
  glLineWidth(4.0);
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);

  for(int i=0; i < mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);

    if(edge->getIndex() == index) {
      int n0 = edge->getNodeIndex(0);
      int n1 = edge->getNodeIndex(1);
	
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
	
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
	
      glVertex3dv(x0);
      glVertex3dv(x1);
    }
  }
  
  glEnd();

  glEnable(GL_LIGHTING);  
  glEndList();
  
  return current;
}



// Generate sharp edge list...
//-----------------------------------------------------------------------------
GLuint GLWidget::generateSharpEdgeList(QColor qColor)
{
  double x0[3], x1[3];

  double R = qColor.red() / 255.0;
  double G = qColor.green() / 255.0;
  double B = qColor.blue() / 255.0;

  GLuint current = glGenLists(1);
  glNewList(current, GL_COMPILE);

  glColor3d(R, G, B);  
  glLineWidth(1.0);
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  
  for(int i=0; i < mesh->getEdges(); i++) {
    edge_t *edge = mesh->getEdge(i);

    if(edge->isSharp()) {
      int n0 = edge->getNodeIndex(0);
      int n1 = edge->getNodeIndex(1);
	
      x0[0] = (mesh->getNode(n0)->getX(0) - drawTranslate[0]) / drawScale;
      x0[1] = (mesh->getNode(n0)->getX(1) - drawTranslate[1]) / drawScale;
      x0[2] = (mesh->getNode(n0)->getX(2) - drawTranslate[2]) / drawScale;
	
      x1[0] = (mesh->getNode(n1)->getX(0) - drawTranslate[0]) / drawScale;
      x1[1] = (mesh->getNode(n1)->getX(1) - drawTranslate[1]) / drawScale;
      x1[2] = (mesh->getNode(n1)->getX(2) - drawTranslate[2]) / drawScale;
	
      glVertex3dv(x0);
      glVertex3dv(x1);
    }
  }
  
  glEnd();
  glEnable(GL_LIGHTING);  
  glEndList();
  
  return current;
}


// Draw coordinates:
//-----------------------------------------------------------------------------
void GLWidget::drawCoordinates()
{
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  // glTranslated(-0.8, -0.8, 5.0);
  glTranslated(-0.8, -0.8, ZSHIFT);

  glMatrixMode(GL_MODELVIEW);

  // z-axis
  glColor3d(0, 0, 1);
  gluCylinder(quadric_axis, 0.02, 0.0, 0.2, 8, 8);  
  renderText(0.0, 0.0, 0.25, "Z");

  // x-axis
  glColor3d(1, 0, 0);
  glRotated(90, 0, 1, 0);
  gluCylinder(quadric_axis, 0.02, 0.0, 0.2, 8, 8);  
  renderText(0.0, 0.0, 0.25, "X");
  glRotated(-90, 0, 1, 0);

  // y-axis
  glColor3d(0, 1, 0);
  glRotated(-90, 1, 0, 0);
  gluCylinder(quadric_axis, 0.02, 0.0, 0.2, 8, 8);  
  renderText(0.0, 0.0, 0.25, "Y");
  glRotated(90, 1, 0, 0);
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);

  return;
}

bool GLWidget::toggleCoordinates()
{
#if WITH_QT6
  makeCurrent();
#endif
  stateDrawCoordinates = !stateDrawCoordinates;
  updateGL();
  return stateDrawCoordinates;
}

// Draw background image...
//-----------------------------------------------------------------------------
void GLWidget::drawBgImage()
{ 
  GLint viewport[4];

  if(!bgTexture) {
#if WITH_QT5 || WITH_QT6
    cout << "Bind texture " << string(bgImageFileName.toLatin1()) << "... ";
	QOpenGLTexture texture(QImage(bgImageFileName).mirrored());
	bgSizeX = texture.width();
    bgSizeY = texture.height();
    bgTexture = texture.textureId();
    cout << "done" << endl;
#else
    cout << "Bind texture " << string(bgImageFileName.toAscii()) << "... ";
    QPixmap pixmap(bgImageFileName);
    bgSizeX = pixmap.width();
    bgSizeY = pixmap.height();
    bgTexture = bindTexture(pixmap, GL_TEXTURE_2D);
    cout << "done" << endl;
#endif
  }
  
  if(!bgTexture) {
    cout << "Failed to bind texture" << endl;
    stateUseBgImage = false;
    return;
  }

  glGetIntegerv(GL_VIEWPORT, viewport);
  
  double relativeSizeX = (double)viewport[2] / (double)viewport[3];
  double relativeSizeY = 1.0;

  double bgRelativeSizeX = (double)bgSizeX / (double)viewport[3];
  double bgRelativeSizeY = (double)bgSizeY / (double)viewport[3];

  double width = 1.0;
  double height = 1.0;
  double depth = 9.9;
  double xshift = 0.0;
  double yshift = 0.0;

  if(stateAlignRightBgImage) {
    width = bgRelativeSizeX;
    height = bgRelativeSizeY;
    xshift = relativeSizeX - bgRelativeSizeX;
    yshift = bgRelativeSizeY - relativeSizeY;
  }

  if(stateStretchBgImage) {
    width = (double)viewport[2] / (double)viewport[3];
    height = 1.0;
    xshift = 0.0;
    yshift = 0.0;
  }
  
  glDisable(GL_DEPTH_TEST);
  glPushMatrix();
  glLoadIdentity();  
  glColor3d(1, 1, 1);
  glDisable(GL_LIGHTING);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, bgTexture);
  
  glBegin(GL_QUADS);    
  glTexCoord2d(0, 0);
  glVertex3d(-width+xshift, -height+yshift, -depth);
  glTexCoord2d(1, 0);
  glVertex3d(+width+xshift, -height+yshift, -depth);
  glTexCoord2d(1, 1);
  glVertex3d(+width+xshift, +height+yshift, -depth);
  glTexCoord2d(0, 1);
  glVertex3d(-width+xshift, +height+yshift, -depth);
  glEnd();  
  
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_LIGHTING);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
}


// Auxiliary function for changing the direction of a vector...
//---------------------------------------------------------------------------
void GLWidget::changeNormalDirection(double *u, double *v)
{
  u[0] = -v[0];
  u[1] = -v[1];
  u[2] = -v[2];
}

list_t* GLWidget::getList(int i) const
{
  return list.at(i);
}

int GLWidget::getLists() const
{
  return list.count();
}

// Set 'c' to an RGB color corresponding to index 'i'.
// 'c' should be pre-allocated to a length of at least 3.
void GLWidget::indexColors(double *c, int i)
{
  c[0] = 0.5 + 0.5 * sin(1.0 * i);
  c[1] = 0.5 + 0.5 * cos(2.0 * i);
  c[2] = 0.5 + 0.5 * cos(3.0 * i); 
}

void GLWidget::indexColors(int *c, int i)
{
  double tmp[3];

  indexColors(tmp, i);
  for (int j = 0; j < 3; j++) c[j] = int(tmp[j]*255 + 0.5);
}


void GLWidget::setMeshVisibility(bool stateDrawSurfaceMesh, bool stateDrawVolumeMesh, bool stateDrawSharpEdges){
/*
  This function is used in mouseDoubleClickEvent(mouseEvent event) to avoid segmentation fault observed Linux
  environment with old hardware.
*/

  mesh_t *mesh = getMesh();
  int lists = getLists();

  if (mesh == NULL) {
    return;
  }


  for (int i = 0; i < lists; i++) {
    list_t *l = getList(i);
	int type = l->getType();
    if (type == SURFACEMESHLIST) {
      l->setVisible(stateDrawSurfaceMesh);

      // do not set visible if the parent surface list is hidden
      int p = l->getParent();
      if (p >= 0) {
        list_t *lp = getList(p);
        if (!lp->isVisible())
          l->setVisible(false);
      }
    }else if (type == VOLUMEMESHLIST) {
      l->setVisible(stateDrawVolumeMesh);
	}else if (type == SHARPEDGELIST) {
      l->setVisible(stateDrawSharpEdges);
	}
  }
}

int GLWidget::mostVisibleBody(int n, bool* tmp1){
/*
This function is called in GLWidget::mouseDoubleClickEvent(QMouseEvent *event) to
identify the most visible body when the double-clicked surface is shared by multiple bodies.
This is a public function to be called from ObjectBrowser. 
*/

  long *nElement = new long[n];
  long *nVisibleElement = new long[n];
  for(int i = 0; i < n; i++){
    nElement[i] = 0;
    nVisibleElement[i] = 0;
  }

  for(int i = 0; i < getLists(); i++) {
    list_t *l2 = getList(i);
	if(l2->getNature() == PDE_BOUNDARY && l2->getType() == SURFACELIST) {
      for(int j = 0; j < mesh->getSurfaces(); j++) {
	    surface_t *surf = mesh->getSurface(j);
	    if(surf->getIndex() == l2->getIndex()) { 
	      for(int k = 0; k < surf->getElements(); k++) {
		    int l = surf->getElementIndex(k);
 		    if(l < 0) 
		      break;
		    element_t *elem = mesh->getElement(l);
		    if((elem->getIndex() < 0) || (elem->getIndex() >= n))
		      break;
		  	nElement[elem->getIndex()]++;
		    if(l2->isVisible())nVisibleElement[elem->getIndex()]++;
          }
		}
      }
    }
  }
  
  double max = -1.0;
  int selected = -1;
  double visibility = -2.0;
  for(int i = 0; i < n; i++){
    if(tmp1[i] && nElement[i] > 0){
    visibility = ((double) nVisibleElement[i]) / nElement[i];
//cout << i << " visibility=" << visibility << " (" << nVisibleElement[i] << "/" << nElement[i] << ")" << endl; 
      if(visibility > max){
		max = visibility;
		selected = i;
	  }
    }
  }
  delete[] nElement;
  delete[] nVisibleElement;
//cout << "selected: " << selected << endl;  
  return selected;
}

#if WITH_QT6
void GLWidget::updateGL()
{
	update();
}



/*
 Reference http://stackoverflow.com/questions/28216001/how-to-render-text-with-qopenglwidget/33674071#33674071
*/
void GLWidget::renderText(double x, double y, double z, const QString & str, const QFont & font/* = QFont()*/, int listBase /*= 2000*/)
{
	int width = this->width();
    int height = this->height();

    GLdouble model[4][4], proj[4][4];
    GLint view[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, &model[0][0]);
    glGetDoublev(GL_PROJECTION_MATRIX, &proj[0][0]);
    glGetIntegerv(GL_VIEWPORT, &view[0]);
    GLdouble textPosX = 0, textPosY = 0, textPosZ = 0;

    project(x, y, z, 
                &model[0][0], &proj[0][0], &view[0],
                &textPosX, &textPosY, &textPosZ);

    textPosY = height - textPosY; // y is inverted

    QPainter painter(this);
    painter.setPen(Qt::black);
    painter.setFont(font);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    painter.drawText(textPosX, textPosY, str); // z = pointT4.z + distOverOp / 4
    painter.end();
}

inline GLint GLWidget::project(GLdouble objx, GLdouble objy, GLdouble objz,
    const GLdouble model[16], const GLdouble proj[16],
    const GLint viewport[4],
    GLdouble * winx, GLdouble * winy, GLdouble * winz)
{
    GLdouble in[4], out[4];

    in[0] = objx;
    in[1] = objy;
    in[2] = objz;
    in[3] = 1.0;
    transformPoint(out, model, in);
    transformPoint(in, proj, out);

    if (in[3] == 0.0)
        return GL_FALSE;

    in[0] /= in[3];
    in[1] /= in[3];
    in[2] /= in[3];

    *winx = viewport[0] + (1 + in[0]) * viewport[2] / 2;
    *winy = viewport[1] + (1 + in[1]) * viewport[3] / 2;

    *winz = (1 + in[2]) / 2;
    return GL_TRUE;
}

inline void GLWidget::transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4])
{
#define M(row,col)  m[col*4+row]
    out[0] =
        M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
    out[1] =
        M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
    out[2] =
        M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
    out[3] =
        M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
#undef M
}

#endif
