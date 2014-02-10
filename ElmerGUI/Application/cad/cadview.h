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
 *  ElmerGUI cadview                                                         *
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

#ifndef CADVIEW_H
#define CADVIEW_H

#include <QMainWindow>
#include <QHash>

#include "cadpreferences.h"

namespace nglib {
#include "nglib.h"
}

#include <TopoDS_Shape.hxx>

class QMenu;
class QAction;
class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkAppendPolyData;

class pt
{
 public:
  int n;
  double x;
  double y;
};

class seg
{
 public:
  int p0;
  int p1;
  int bc;
};

class CadView : public QMainWindow
{
  Q_OBJECT

public:
  CadView(QWidget *parent = 0);
  ~CadView();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  QVTKWidget* GetQVTKWidget();

  bool readFile(QString);
  void generateSTL();

  void setMesh(nglib::Ng_Mesh*);
  void setGeom(nglib::Ng_STL_Geometry*);
  void setMp(nglib::Ng_Meshing_Parameters*);
  void setDeflection(double);
  double lengthOf(double*);
  void differenceOf(double*, double*, double*);
  int getFaceNumber(vtkActor*);
  int getDim();
  void generateIn2dFile();

 private slots:
  void closeSlot();
  void generateSTLSlot();
  void cadPreferencesSlot();
  void reloadSlot();

 private:
  void createActions();
  void createMenus();
  void clearScreen();
  TopoDS_Shape readBrep(QString);
  TopoDS_Shape readStep(QString);
  TopoDS_Shape readIges(QString);
  void restrictMeshSizeLocal(nglib::Ng_Mesh*, vtkPolyData*, double, double);

  QMenu* fileMenu;
  QMenu* modelMenu;

  QAction* exitAct;
  QAction* reloadAct;
  QAction* cadPreferencesAct;

  QVTKWidget* qVTKWidget;
  vtkRenderer* renderer;

  vtkAppendPolyData* stlSurfaceData;
  vtkAppendPolyData* stlEdgeData;

  int numberOfFaces;
  double modelLength;

  nglib::Ng_Mesh* mesh;
  nglib::Ng_STL_Geometry* geom;
  nglib::Ng_Meshing_Parameters* mp;

  CadPreferences* cadPreferences;

  QString fileName;

  QHash<vtkActor*, int> actorToFace;

  int modelDim;

  TopoDS_Shape shape;
};

#endif // CADVIEW_H
