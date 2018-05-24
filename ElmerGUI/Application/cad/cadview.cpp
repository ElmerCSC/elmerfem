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

#ifdef WITH_QT5
#include <QtWidgets>
#else
#include <QtGui>
#endif

#include <iostream>

#include "cadview.h"

#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkProperty.h>
#include <vtkPropPicker.h>
#include <vtkCallbackCommand.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCleanPolyData.h>

#include <BRep_Builder.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <BRepMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <GeomLProp_SLProps.hxx>
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

using namespace std;

static void pickEventHandler(vtkObject* caller, unsigned long eid, 
			     void* clientdata, void* calldata)
{
  CadView* cadView = reinterpret_cast<CadView*>(clientdata);
  QVTKWidget* qvtkWidget = cadView->GetQVTKWidget();
  vtkAbstractPicker* picker = qvtkWidget->GetInteractor()->GetPicker();
  vtkPropPicker* propPicker = vtkPropPicker::SafeDownCast(picker);
  vtkActor* actor = propPicker->GetActor();
  int faceNumber = cadView->getFaceNumber(actor);

  if(faceNumber > 0) {
    vtkProperty* p = actor->GetProperty();

    double color[3];
    p->GetColor(color);

    // Toggle color:
    //--------------
    if(color[0] < 0.5) {
      cout << "Secected face: ";
      p->SetColor(1, 0, 0);
    } else {
      cout << "Unselected face: ";
      p->SetColor(0, 1, 1);
    }
    cout << faceNumber << endl;
  }
}

CadView::CadView(QWidget *parent)
  : QMainWindow(parent)
{
  setWindowTitle("ElmerGUI geometry viewer");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));

  createActions();
  createMenus();

  qVTKWidget = new QVTKWidget(this);
  setCentralWidget(qVTKWidget);

  renderer = vtkRenderer::New();
  renderer->SetBackground(1, 1, 1);
  qVTKWidget->GetRenderWindow()->AddRenderer(renderer);
  renderer->GetRenderWindow()->Render();

  vtkPropPicker* propPicker = vtkPropPicker::New();
  vtkCallbackCommand* cbcPick = vtkCallbackCommand::New();
  qVTKWidget->GetInteractor()->SetPicker(propPicker);
  cbcPick->SetClientData(this);
  cbcPick->SetCallback(pickEventHandler);
  qVTKWidget->GetInteractor()->GetPicker()->AddObserver(vtkCommand::PickEvent, cbcPick);
  propPicker->Delete();
  cbcPick->Delete();

  stlSurfaceData = vtkAppendPolyData::New();
  stlEdgeData = vtkAppendPolyData::New();

  cadPreferences = new CadPreferences(this);

  modelLength = 0.0;
  numberOfFaces = 0;
  fileName = "";
  modelDim = -1;
}

CadView::~CadView()
{
}

QSize CadView::minimumSizeHint() const
{
  return QSize(64, 64);
}

QSize CadView::sizeHint() const
{
  return QSize(720, 576);
}

void CadView::createActions()
{
  exitAct = new QAction(QIcon(""), tr("&Quit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(closeSlot()));

  cadPreferencesAct = new QAction(QIcon(""), tr("Preferences..."), this);
  cadPreferencesAct->setShortcut(tr("Ctrl+P"));
  connect(cadPreferencesAct, SIGNAL(triggered()), this, SLOT(cadPreferencesSlot()));

  reloadAct = new QAction(QIcon(""), tr("Reload geometry"), this);
  reloadAct->setShortcut(tr("Ctrl+R"));
  connect(reloadAct, SIGNAL(triggered()), this, SLOT(reloadSlot()));  
}

void CadView::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(reloadAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  modelMenu = menuBar()->addMenu(tr("&Model"));
  modelMenu->addAction(cadPreferencesAct);
}

void CadView::closeSlot()
{
  close();
}

void CadView::cadPreferencesSlot()
{
  cadPreferences->show();
}

void CadView::reloadSlot()
{
  if(this->fileName.isEmpty()) return;
  readFile(this->fileName);
}

bool CadView::readFile(QString fileName)
{
  double deflection = cadPreferences->ui.deflection->text().toDouble();
  double featureAngle = cadPreferences->ui.featureAngle->text().toDouble();
  bool mergePoints = cadPreferences->ui.mergePoints->isChecked();

  if(deflection < 0.0) {
    deflection = 0.0005;
    cout << "Bad value for deflection. Using: " << deflection << endl;
  }

  if(featureAngle < 0.0) {
    featureAngle = 30.0;
    cout << "Bad value for feature angle. Using: " << featureAngle << endl;
  }

  if(stlSurfaceData->GetOutput()->GetNumberOfPoints() > 0)
    stlSurfaceData->Delete();

  if(stlEdgeData->GetOutput()->GetNumberOfPoints() > 0)
    stlEdgeData->Delete();

  stlSurfaceData = vtkAppendPolyData::New();
  stlEdgeData = vtkAppendPolyData::New();

  if(fileName.isEmpty()) {
    cout << "File name is empty. Aborting." << endl;
    return false;
  }

  QFileInfo fileInfo(fileName);
  QString fileSuffix = fileInfo.suffix().toLower();
  
  // TopoDS_Shape shape;

  if(fileSuffix == "brep")
    shape = readBrep(fileName);

  if((fileSuffix == "step") || (fileSuffix == "stp"))
    shape = readStep(fileName);
  
  if((fileSuffix == "iges") || (fileSuffix == "igs"))
    shape = readIges(fileName);
  
  if(shape.IsNull()) {
    cout << "Cad import: No shapes. Aborting" << endl;
    return false;
  }

  BRepTools::Clean(shape);

  // Check 3D properties:
  //----------------------
  GProp_GProps System;
  BRepGProp::VolumeProperties(shape, System);
  double mass = System.Mass();

  if(mass < 1.0e-12) {
    QMessageBox message;
    message.setIcon(QMessageBox::Warning);
    message.setText("Non 3D-shape detected");
    message.setInformativeText("The cad import features of ElmerGUI are currently limited to 3D models. Please consider using external software or other formats for meshing 1D and 2D geometries.");
    message.exec();
  }

  // Go:
  //-----
  this->fileName = fileName;

  actorToFace.clear();

  clearScreen();

  // Compute bounding box:
  //----------------------
  Bnd_Box boundingBox;
  double min[3], max[3];
  
  BRepBndLib::Add(shape, boundingBox);
  boundingBox.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

  cout << "Bounding box: "
       << "[ " << min[0] << ", " << min[1] << ", " << min[2] << "] x "
       << "[ " << max[0] << ", " << max[1] << ", " << max[2] << "]" << endl;

  double length = sqrt((max[2]-min[2])*(max[2]-min[2])
		       +(max[1]-min[1])*(max[1]-min[1]) 
		       +(max[0]-min[0])*(max[0]-min[0]));

  deflection *= length; // use relative units

  double t0 = sqrt((max[0] - min[0])*(max[0] - min[0]));
  double t1 = sqrt((max[1] - min[1])*(max[1] - min[1]));
  double t2 = sqrt((max[2] - min[2])*(max[2] - min[2]));

  modelDim = 3;

  double tol = 1.0e-6 * length;

  if((t0 < tol) || (t1 < tol) || (t2 < tol)) {
    modelDim = 2;
    // cout << "Cad import: Shape seems to be 2D. Unable to proceed. Aborting." << endl;
    // return false;
  }
  
  // Construct model data and draw surfaces:
  //-----------------------------------------
#if OCC_VERSION_HEX >= 0x060800
  BRepMesh_IncrementalMesh(shape,deflection);
#else
  BRepMesh::Mesh(shape, deflection);
#endif

  numberOfFaces = 0;
  TopExp_Explorer expFace;
  for(expFace.Init(shape, TopAbs_FACE); expFace.More(); expFace.Next()) {
    TopoDS_Face Face = TopoDS::Face(expFace.Current());

    TopLoc_Location Location;
    Handle(Poly_Triangulation) Triangulation = BRep_Tool::Triangulation(Face, Location);

    if(Triangulation.IsNull()) {
      cout << "Encountered empty triangulation after face: " << numberOfFaces+1 << endl;
      continue;
    }

    const gp_Trsf& Transformation = Location.Transformation();

    const Poly_Array1OfTriangle& Triangles = Triangulation->Triangles();
    const TColgp_Array1OfPnt& Nodes = Triangulation->Nodes();

    int nofTriangles = Triangulation->NbTriangles();
    int nofNodes = Triangulation->NbNodes();

    if(nofTriangles < 1) {
      cout << "No triangles for mesh on face: " << numberOfFaces+1 << endl;
      continue;
    }

    if(nofNodes < 1) {
      cout << "No nodes for mesh on face: " << numberOfFaces+1 << endl;
      continue;
    }
    
    numberOfFaces++;

    int n0, n1, n2;
    vtkPolyData* partGrid = vtkPolyData::New();
    vtkTriangle* triangle = vtkTriangle::New();
    partGrid->Allocate(nofTriangles, nofTriangles);

    for(int i = Triangles.Lower(); i <= Triangles.Upper(); i++) {
      Triangles(i).Get(n0, n1, n2);

      if(Face.Orientation() != TopAbs_FORWARD) {
	int tmp = n2; n2 = n1; n1 = tmp;
      }

      triangle->GetPointIds()->SetId(0, n0 - Nodes.Lower());
      triangle->GetPointIds()->SetId(1, n1 - Nodes.Lower());
      triangle->GetPointIds()->SetId(2, n2 - Nodes.Lower());

      partGrid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());
    }

    double x[3];
    vtkPoints* partPoints = vtkPoints::New();
    for(int i = Nodes.Lower(); i <= Nodes.Upper(); i++) {
      gp_XYZ XYZ = Nodes(i).Coord();
      Transformation.Transforms(XYZ);
      x[0] = XYZ.X(); x[1] = XYZ.Y(); x[2] = XYZ.Z();
      partPoints->InsertPoint(i - Nodes.Lower(), x);
    }

    partGrid->SetPoints(partPoints);

    // Draw part:
    //-----------
    vtkCleanPolyData* partCleaner = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION <= 5
    partCleaner->SetInput(partGrid);
#else
    partCleaner->SetInputData(partGrid);
    partCleaner->Update();
#endif
    if(mergePoints) {
      partCleaner->PointMergingOn();
    } else {
      partCleaner->PointMergingOff();
    }

    vtkPolyDataNormals* partNormals = vtkPolyDataNormals::New();
    partNormals->SetInputConnection(partCleaner->GetOutputPort());
    partNormals->SetFeatureAngle(featureAngle);
    
    vtkDataSetMapper* partMapper = vtkDataSetMapper::New();
    partMapper->SetInputConnection(partNormals->GetOutputPort());
    partMapper->ScalarVisibilityOff();
    
    vtkActor* partActor = vtkActor::New();
    partActor->SetPickable(1);
    partActor->GetProperty()->SetColor(0, 1, 1);
    partActor->SetMapper(partMapper);
    renderer->AddActor(partActor);
    actorToFace.insert(partActor, numberOfFaces);

    vtkFeatureEdges* partFeature = vtkFeatureEdges::New();
    partFeature->SetInputConnection(partCleaner->GetOutputPort());
    partFeature->SetFeatureAngle(featureAngle);
    partFeature->FeatureEdgesOff();
    partFeature->BoundaryEdgesOn();
    partFeature->NonManifoldEdgesOn();
    partFeature->ManifoldEdgesOff();

    vtkDataSetMapper *partFeatureMapper = vtkDataSetMapper::New();
    partFeatureMapper->SetInputConnection(partFeature->GetOutputPort());
    partFeatureMapper->SetResolveCoincidentTopologyToPolygonOffset();
    partFeatureMapper->ScalarVisibilityOff();
    
    vtkActor* partFeatureActor = vtkActor::New();
    partFeatureActor->SetPickable(0);
    partFeatureActor->GetProperty()->SetColor(0, 0, 0);
    partFeatureActor->SetMapper(partFeatureMapper);
    renderer->AddActor(partFeatureActor);

    // Add triangles and edges to STL structures:
    //--------------------------------------------
#if VTK_MAJOR_VERSION <= 5
    stlSurfaceData->AddInput(partCleaner->GetOutput());
    stlEdgeData->AddInput(partFeature->GetOutput());
#else
    stlSurfaceData->AddInputData(partCleaner->GetOutput());
    stlEdgeData->AddInputData(partFeature->GetOutput());
#endif

    // Clean up:
    //----------
    partFeatureActor->Delete();
    partFeatureMapper->Delete();
    partFeature->Delete();
    partActor->Delete();
    partNormals->Delete();
    partMapper->Delete();
    partCleaner->Delete();
    partGrid->Delete();
    partPoints->Delete();
    triangle->Delete();
  }

  if(numberOfFaces < 1) {
    cout << "Cad import: error: no surface triangulation was generated. Aborting." << endl;
    return false;
  }


  stlSurfaceData->Update();
  stlEdgeData->Update();
  modelLength = stlSurfaceData->GetOutput()->GetLength();
  cout << "StlSurfaceData: points: " << stlSurfaceData->GetOutput()->GetNumberOfPoints() << endl;
  cout << "StlSurfaceData: cells: " << stlSurfaceData->GetOutput()->GetNumberOfCells() << endl;
  cout << "StlEdgeData: lines: " << stlEdgeData->GetOutput()->GetNumberOfLines() << endl;
  cout << "StlModelData: length: " << modelLength << endl;

  // Draw:
  //------
  qVTKWidget->GetRenderWindow()->Render();
  renderer->ResetCamera();

  QCoreApplication::processEvents();

  return true;
}

void CadView::generateSTLSlot()
{
  double meshMinSize = modelLength * cadPreferences->ui.minh->text().toDouble();
  double meshMaxSize = modelLength * cadPreferences->ui.maxh->text().toDouble();
  bool restrictBySTL = cadPreferences->ui.restrictBySTL->isChecked();

  // Check also the MainWindow meshing preferences:
  if(mp->maxh < meshMaxSize) meshMaxSize = mp->maxh;
  double meshFineness = mp->fineness;
  
  if(meshMaxSize > 0.1 * modelLength)
    meshMaxSize = 0.1 * modelLength;

  if(meshMinSize > meshMaxSize)
    meshMinSize = meshMaxSize;

  if(meshMinSize < 0)
    meshMinSize = modelLength * 0.005;

  if(meshMaxSize < 0)
    meshMaxSize = modelLength * 0.1;

  cout << "Cad import: max mesh size: " << meshMaxSize << endl;
  cout << "Cad import: mesh fineness: " << meshFineness << endl;

  // Add STL triangles to geometry:
  //--------------------------------
  vtkCleanPolyData* stlSurface = vtkCleanPolyData::New();
  stlSurface->PointMergingOn();
#if VTK_MAJOR_VERSION <= 5
  stlSurface->SetInput(stlSurfaceData->GetOutput());
#else
  stlSurface->SetInputConnection(stlSurfaceData->GetOutputPort());
#endif

  stlSurface->Update();

  if(stlSurface->GetOutput()->GetNumberOfCells() < 1) {
    cout << "Cad import: error: geometry undefined - no STL available" << endl;
    return;
  }

  double p0[3], p1[3], p2[3];
  for(int i = 0; i < stlSurface->GetOutput()->GetNumberOfCells(); i++) {
    vtkCell* cell = stlSurface->GetOutput()->GetCell(i);
    int nofCellPoints = cell->GetNumberOfPoints();
    vtkPoints* cellPoints = cell->GetPoints();

    if(nofCellPoints == 3) {
      cellPoints->GetPoint(0, p0);
      cellPoints->GetPoint(1, p1);
      cellPoints->GetPoint(2, p2);
      nglib::Ng_STL_AddTriangle(geom, p0, p1, p2, NULL);
    }
  }

  // Add STL edges to geometry:
  //----------------------------
  vtkCleanPolyData* stlEdge = vtkCleanPolyData::New();
  stlEdge->PointMergingOn();
#if VTK_MAJOR_VERSION <= 5
  stlEdge->SetInput(stlEdgeData->GetOutput());
#else
  stlEdge->SetInputConnection(stlEdgeData->GetOutputPort());
#endif

  stlEdge->Update();

  for(int i = 0; i < stlEdge->GetOutput()->GetNumberOfCells(); i++) {
    vtkCell* cell = stlEdge->GetOutput()->GetCell(i);
    int nofCellPoints = cell->GetNumberOfPoints();
    vtkPoints* cellPoints = cell->GetPoints();

    if(nofCellPoints == 2) {
      cellPoints->GetPoint(0, p0);
      cellPoints->GetPoint(1, p1);
      nglib::Ng_STL_AddEdge(geom, p0, p1);
    }
  }

  // Init STL geometry:
  //--------------------
  nglib::Ng_STL_InitSTLGeometry(geom);

  // Generate edges:
  //-----------------
  nglib::Ng_STL_MakeEdges(geom, mesh, mp);
  
  // Global mesh size restrictions:
  //--------------------------------
  nglib::Ng_RestrictMeshSizeGlobal(mesh, meshMaxSize);
  
  // Local mesh size restrictions:
  //-------------------------------
  if(restrictBySTL)
    restrictMeshSizeLocal(mesh, stlSurface->GetOutput(),
			  meshMaxSize, meshMinSize);
}

QVTKWidget* CadView::GetQVTKWidget()
{
  return this->qVTKWidget;
}

void CadView::clearScreen()
{
  cout << "Clear screen" << endl;
  vtkActorCollection* actors = renderer->GetActors();
  vtkActor* lastActor = actors->GetLastActor();
  while(lastActor != NULL) {
    renderer->RemoveActor(lastActor);
    lastActor = actors->GetLastActor();
  }
}

TopoDS_Shape CadView::readBrep(QString fileName)
{
  TopoDS_Shape shape;
  BRep_Builder builder;
  Standard_Boolean result;

  result = BRepTools::Read(shape, fileName.toLatin1().data(), builder);    

  if(!result)
    cout << "Read brep failed" << endl;
  
  return shape;
}

TopoDS_Shape CadView::readStep(QString fileName)
{
  TopoDS_Shape shape;
  Handle_TopTools_HSequenceOfShape shapes;
  STEPControl_Reader stepReader;
  IFSelect_ReturnStatus status;

  status = stepReader.ReadFile(fileName.toLatin1().data());
  
  if(status == IFSelect_RetDone) {	  
    bool failsonly = false;
    stepReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);
    
    int nbr = stepReader.NbRootsForTransfer();
    stepReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
    
    for(Standard_Integer n = 1; n <= nbr; n++) {
      bool ok = stepReader.TransferRoot(n);
      int nbs = stepReader.NbShapes();

      // Display warning if nbs > 1
      //----------------------------
      if(nbs > 1) {
	QMessageBox message;
	message.setIcon(QMessageBox::Warning);
	message.setText("Loading multiple shapes");
	message.setInformativeText("The mesh generators of ElmerGUI are currently unable to handle cad files with multiple shapes. Please consider using external software for mesh generation in this case.");
	message.exec();
      }
      
      if(nbs > 0) {
	shapes = new TopTools_HSequenceOfShape();
	for(int i = 1; i <= nbs; i++) {
	  cout << "Added shape: " << i << endl;
	  shape = stepReader.Shape(i);
	  shapes->Append(shape);
	}
      }
    }
  }

  return shape;
}

TopoDS_Shape CadView::readIges(QString fileName)
{
  TopoDS_Shape shape;
  IGESControl_Reader igesReader;
  IFSelect_ReturnStatus status;

  status = igesReader.ReadFile(fileName.toLatin1().data());
  
  if(status == IFSelect_RetDone) {
    igesReader.TransferRoots();
    shape = igesReader.OneShape();
  }

  return shape;
}

void CadView::restrictMeshSizeLocal(nglib::Ng_Mesh* mesh, vtkPolyData* stlData,
				    double meshMaxSize, double meshMinSize)
{
  int n0, n1, n2;
  double h, h0, h1, h2;
  double t[3], p0[3], p1[3], p2[3];
  vtkFloatArray* mshSize = vtkFloatArray::New();
  mshSize->SetNumberOfComponents(1);
  mshSize->SetNumberOfTuples(stlData->GetNumberOfPoints());
  
  for(int i = 0; i < stlData->GetNumberOfPoints(); i++) 
    mshSize->SetComponent(i, 0, meshMaxSize);

  for(int i = 0; i < stlData->GetNumberOfCells(); i++) {
    vtkCell* cell = stlData->GetCell(i);
    int nofCellPoints = cell->GetNumberOfPoints();
    vtkPoints* cellPoints = cell->GetPoints();
    
    if(nofCellPoints == 3) {
      n0 = cell->GetPointId(0);
      n1 = cell->GetPointId(1);
      n2 = cell->GetPointId(2);
      
      h0 = mshSize->GetComponent(n0, 0);
      h1 = mshSize->GetComponent(n1, 0);
      h2 = mshSize->GetComponent(n2, 0);
      
      cellPoints->GetPoint(0, p0);
      cellPoints->GetPoint(1, p1);
      cellPoints->GetPoint(2, p2);
      
      differenceOf(t, p0, p1);
      h = lengthOf(t);
      if(h < meshMinSize) h = meshMinSize;
      if(h < h1) mshSize->SetComponent(n1, 0, h);
      if(h < h0) mshSize->SetComponent(n0, 0, h);
      
      differenceOf(t, p2, p0);
      h = lengthOf(t);
      if(h < meshMinSize) h = meshMinSize;
      if(h < h2) mshSize->SetComponent(n2, 0, h);
      if(h < h0) mshSize->SetComponent(n0, 0, h);
      
      differenceOf(t, p2, p1);
      h = lengthOf(t);
      if(h < meshMinSize) h = meshMinSize;
      if(h < h2) mshSize->SetComponent(n2, 0, h);
      if(h < h1) mshSize->SetComponent(n1, 0, h);
    }
  }  

  for(int i = 0; i < stlData->GetNumberOfPoints(); i++) {
    h = mshSize->GetComponent(i, 0);
    if(h < meshMinSize) h = meshMinSize;
    if(h > meshMaxSize) h = meshMaxSize;
    stlData->GetPoint(i, p0);
    nglib::Ng_RestrictMeshSizePoint(mesh, p0, h);
  }
  
  mshSize->Delete();
}

void CadView::generateSTL()
{
  this->generateSTLSlot();
}

void CadView::setMesh(nglib::Ng_Mesh* mesh)
{
  this->mesh = mesh;
}

void CadView::setGeom(nglib::Ng_STL_Geometry* geom)
{
  this->geom = geom;
}

void CadView::setMp(nglib::Ng_Meshing_Parameters* mp)
{
  this->mp = mp;
}

void CadView::setDeflection(double deflection)
{
  if(deflection < 0) return;

  this->cadPreferences->ui.deflection->setText(QString::number(deflection));
}

double CadView::lengthOf(double* v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void CadView::differenceOf(double* u, double* v, double* w)
{
  u[0] = v[0] - w[0];
  u[1] = v[1] - w[1];
  u[2] = v[2] - w[2];
}

int CadView::getFaceNumber(vtkActor* actor)
{
  if(actor == NULL) return -1;
  return actorToFace.value(actor);
}

int CadView::getDim()
{
  return modelDim;
}

void CadView::generateIn2dFile()
{
  cout << "Generating In2D file from 2D Iges geometry"<< endl;

  QVector<pt> pts;
  QVector<seg> segs;

  bool firstPoint = true;
  int numberOfFaces = 0;
  int numberOfPts = 0;
  int numberOfSegs = 0;
  gp_Pnt previousPnt;

  // Loop over faces:
  //------------------
  TopExp_Explorer expFace(shape, TopAbs_FACE);
  for(expFace; expFace.More(); expFace.Next()) {
    TopoDS_Face Face = TopoDS::Face(expFace.Current());
    cout << "Face: " << ++numberOfFaces << endl;

    // Loop over edges:
    //------------------
    int numberOfEdges = 0;
    TopExp_Explorer expEdge(Face, TopAbs_EDGE);
    for(expEdge; expEdge.More(); expEdge.Next()) {
      TopoDS_Edge Edge = TopoDS::Edge(expEdge.Current());
      cout << " Edge: " << ++numberOfEdges << endl;

      // Divide edge into segments:
      //----------------------------
      double AngularDeflection = 0.1;
      double CurvatureDeflection = 0.1;
      int MinPoints = 2;
      double Tolerance = 1.0e-9;

      BRepAdaptor_Curve2d Curve(Edge, Face);
      GCPnts_TangentialDeflection TD(Curve, AngularDeflection, 
				     CurvatureDeflection, MinPoints,
				     Tolerance);

      int nofPoints = TD.NbPoints();
      cout << "  Points: " << nofPoints << endl;

      // Loop over points:
      //-------------------
      for(int i = 2; i <= nofPoints; i++) {
	gp_Pnt value;

	if(firstPoint) {
	  value = TD.Value(1);
	  firstPoint = false;
	  previousPnt = value;
	  pt p;
	  p.n = ++numberOfPts;
	  p.x = value.X();
	  p.y = value.Y();
	  pts.push_back(p);
	}

	double p0 = TD.Parameter(i-1);
	double p1 = TD.Parameter(i);

	double dist = sqrt((p1-p0)*(p1-p0));

	if(dist < Tolerance) {
	  cout << "   Skipped one (based on parameter)" << endl;
	  continue;
	}

	value = TD.Value(i);

	double dx = value.X() - previousPnt.X();
	double dy = value.Y() - previousPnt.Y();

	dist = sqrt(dx*dx + dy*dy);

	if(dist < Tolerance) {
	  cout << "   Skipped one (based on value)" << endl;
	  continue;
	}
	
	pt p;
	p.n = ++numberOfPts;
	p.x = value.X();
	p.y = value.Y();
	pts.push_back(p);

	seg s;
	s.p0 = numberOfPts - 1;
	s.p1 = numberOfPts;
	s.bc = numberOfEdges;
	segs.push_back(s);
	numberOfSegs++;

	previousPnt = value;
      }
    }
  }

  // Write the in2d file:
  //----------------------
  fstream file("iges2ng.in2d", ios::out);

  file << "splinecurves2dv2" << endl;
  file << "1" << endl << endl;

  file << "points" << endl;
  for(int i = 0; i < pts.size() - 1; i++) {
    pt p = pts[i];
    file << p.n << " " << p.x << " " << p.y << endl;
  }

  seg s;
  file << endl << "segments" << endl;
  for(int i = 0; i < segs.size() - 1; i++) {
    s = segs[i];
    file << "1 0 2 " << s.p0 << " " << s.p1 << " -bc=" << s.bc << endl;
  }

  file << "1 0 2 " << s.p1  << " 1 -bc=" << s.bc << endl;

  file << endl << "materials" << endl;
  file << "1 mat1 -maxh=100000" << endl;
				       
  file.close();
}
