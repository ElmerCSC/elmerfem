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
 *  ElmerGUI vtkpost                                                         *
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
#include <QScriptEngine>
#include <iostream>

#include "epmesh.h"
#include "vtkpost.h"
#include "surface.h"
#include "isocontour.h"
#include "isosurface.h"
#include "colorbar.h"
#include "preferences.h"
#include "vector.h"
#include "readepfile.h"
#include "streamline.h"
#include "timestep.h"
#include "axes.h"
#include "text.h"
#include "featureedge.h"
#include "meshpoint.h"
#include "meshedge.h"
#include "ecmaconsole.h"

#ifdef EG_MATC
#include "matc.h"
#endif

#include <QVTKWidget.h>
#include <vtkLookupTable.h>
#include <vtkActor.h>
#include <vtkTextActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkLine.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCellDerivatives.h>
#include <vtkCellDataToPointData.h>
#include <vtkPlane.h>
#include <vtkPropPicker.h>
#include <vtkCallbackCommand.h>
#include <vtkAbstractPicker.h>
#include <vtkObject.h>
#include <vtkCommand.h>
#include <vtkFollower.h>
#include <vtkImplicitPlaneWidget.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIndent.h>

using namespace std;

// Custom print for QtScript:
//----------------------------
#if QT_VERSION >= 0x040403
QScriptValue printFun(QScriptContext* context, QScriptEngine* engine)
{
  QString result;
  for (int i = 0; i < context->argumentCount(); ++i) {
    if (i > 0)
      result.append(" ");
    result.append(context->argument(i).toString());
  }
  
  QScriptValue calleeData = context->callee().data();
  EcmaConsole* ecmaConsole = qobject_cast<EcmaConsole*>(calleeData.toQObject());
  ecmaConsole->append(result);
  
  return engine->undefinedValue();
}
#endif

// Interaction event handler (press 'i' to interact):
//-------------------------------------------------------------------
static void iEventHandler(vtkObject* caller, unsigned long eid, 
			  void* clientdata, void* calldata)
{
  VtkPost* vtkPost = reinterpret_cast<VtkPost*>(clientdata);
  vtkImplicitPlaneWidget* planeWidget = vtkPost->GetPlaneWidget();

  vtkPost->SetClipPlaneOrigin(planeWidget->GetOrigin());
  vtkPost->SetClipPlaneNormal(planeWidget->GetNormal());
}

// Pick event handler (press 'p' to pick):
//-------------------------------------------------------------------
static void pEventHandler(vtkObject* caller, unsigned long eid, 
			  void* clientdata, void* calldata)
{
  VtkPost* vtkPost = reinterpret_cast<VtkPost*>(clientdata);
  vtkRenderer* renderer = vtkPost->GetRenderer();
  vtkActor* pickedPointActor = vtkPost->GetPickedPointActor();
  QVTKWidget* qvtkWidget = vtkPost->GetQVTKWidget();
  vtkAbstractPicker* picker = qvtkWidget->GetInteractor()->GetPicker();
  vtkPropPicker* propPicker = vtkPropPicker::SafeDownCast(picker);

  vtkActor* actor = propPicker->GetActor();
  double* pickPos = vtkPost->GetCurrentPickPosition();
  propPicker->GetPickPosition(pickPos);
  vtkPost->SetCurrentPickPosition(pickPos);

  if(!actor) {
    renderer->RemoveActor(pickedPointActor);

  } else {
    vtkDataSetMapper* mapper = vtkDataSetMapper::New();
    vtkUnstructuredGrid* cross = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkLine* line = vtkLine::New();
    double l = vtkPost->GetLength() / 15.0;

    points->SetNumberOfPoints(6);
    points->InsertPoint(0, +l,  0,  0);
    points->InsertPoint(1,  0, +l,  0);
    points->InsertPoint(2,  0,  0, +l);
    points->InsertPoint(3, -l,  0,  0);    
    points->InsertPoint(4,  0, -l,  0);    
    points->InsertPoint(5,  0,  0, -l);
    cross->SetPoints(points);

    line->GetPointIds()->SetId(0, 0);
    line->GetPointIds()->SetId(1, 3);
    cross->InsertNextCell(line->GetCellType(), line->GetPointIds());

    line->GetPointIds()->SetId(0, 1);
    line->GetPointIds()->SetId(1, 4);
    cross->InsertNextCell(line->GetCellType(), line->GetPointIds());

    line->GetPointIds()->SetId(0, 2);
    line->GetPointIds()->SetId(1, 5);
    cross->InsertNextCell(line->GetCellType(), line->GetPointIds());

    mapper->SetInput(cross);

    pickedPointActor->SetMapper(mapper);
    pickedPointActor->SetPosition(pickPos);
    pickedPointActor->GetProperty()->SetColor(1, 0, 0);

    renderer->AddActor(pickedPointActor);

    mapper->Delete();
    cross->Delete();
    points->Delete();
    line->Delete();
  }
}

// Class VtkPost:
//-----------------------------------------------------------------
VtkPost::VtkPost(QWidget *parent)
  : QMainWindow(parent)
{
  // Initialize:
  //------------
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
  setWindowTitle("ElmerGUI postprocessor");
  resize(800, 600);

  createActions();
  createMenus();
  createToolbars();
  createStatusBar();

  // VTK:
  //-----
  volumeGrid = vtkUnstructuredGrid::New();
  surfaceGrid = vtkUnstructuredGrid::New();
  lineGrid = vtkUnstructuredGrid::New();
  isoContourActor = vtkActor::New();
  isoSurfaceActor = vtkActor::New();
  surfaceActor = vtkActor::New();
  meshEdgeActor = vtkActor::New();
  meshPointActor = vtkActor::New();
  colorBarActor = vtkScalarBarActor::New();
  featureEdgeActor = vtkActor::New();
  vectorActor = vtkActor::New();
  streamLineActor = vtkActor::New();
  axesActor = vtkActor::New();
  axesXTextActor = vtkFollower::New();
  axesYTextActor = vtkFollower::New();
  axesZTextActor = vtkFollower::New();
  pickedPointActor = vtkActor::New();
  textActor = vtkTextActor::New();

  // Default color map (from blue to red):
  //--------------------------------------
  currentLut = vtkLookupTable::New();
  currentLut->SetHueRange(0.6667, 0);
  currentLut->SetNumberOfColors(128);
  currentLut->Build();

  // User interfaces, widgets, and draw routines:
  //---------------------------------------------
  surface = new Surface(this);
  connect(surface, SIGNAL(drawSurfaceSignal()), this, SLOT(drawSurfaceSlot()));
  connect(surface, SIGNAL(hideSurfaceSignal()), this, SLOT(hideSurfaceSlot()));

  vector = new Vector(this);
  connect(vector, SIGNAL(drawVectorSignal()), this, SLOT(drawVectorSlot()));
  connect(vector, SIGNAL(hideVectorSignal()), this, SLOT(hideVectorSlot()));
  
  isoContour = new IsoContour(this);
  connect(isoContour, SIGNAL(drawIsoContourSignal()), this, SLOT(drawIsoContourSlot()));
  connect(isoContour, SIGNAL(hideIsoContourSignal()), this, SLOT(hideIsoContourSlot()));

  isoSurface = new IsoSurface(this);
  connect(isoSurface, SIGNAL(drawIsoSurfaceSignal()), this, SLOT(drawIsoSurfaceSlot()));
  connect(isoSurface, SIGNAL(hideIsoSurfaceSignal()), this, SLOT(hideIsoSurfaceSlot()));

  colorBar = new ColorBar(this);
  connect(colorBar, SIGNAL(drawColorBarSignal()), this, SLOT(drawColorBarSlot()));
  connect(colorBar, SIGNAL(hideColorBarSignal()), this, SLOT(hideColorBarSlot()));

  streamLine = new StreamLine(this);
  connect(streamLine, SIGNAL(drawStreamLineSignal()), this, SLOT(drawStreamLineSlot()));
  connect(streamLine, SIGNAL(hideStreamLineSignal()), this, SLOT(hideStreamLineSlot()));

  preferences = new Preferences(this);
  connect(preferences, SIGNAL(redrawSignal()), this, SLOT(redrawSlot()));

  timeStep = new TimeStep(this);
  connect(timeStep, SIGNAL(timeStepChangedSignal()), this, SLOT(timeStepChangedSlot()));
  connect(this, SIGNAL(canProceedWithNextSignal(vtkRenderWindow*)), timeStep, SLOT(canProceedWithNextSlot(vtkRenderWindow*)));

  readEpFile = new ReadEpFile(this);
  connect(readEpFile, SIGNAL(readPostFileSignal(QString)), this, SLOT(ReadPostFile(QString)));

  axes = new Axes(this);

  text = new Text(this);
  connect(text, SIGNAL(drawTextSignal()), this, SLOT(drawTextSlot()));
  connect(text, SIGNAL(hideTextSignal()), this, SLOT(hideTextSlot()));

  featureEdge = new FeatureEdge(this);
  meshPoint = new MeshPoint(this);
  meshEdge = new MeshEdge(this);

#ifdef EG_MATC
  matc = new Matc(this);
  connect(matc->ui.mcEdit, SIGNAL(returnPressed()), 
        this, SLOT(domatcSlot()));
  connect(matc->ui.mcHistory, SIGNAL(selectionChanged()), 
        this, SLOT(matcCutPasteSlot()));
#endif

  // Ep-data:
  //----------
  epMesh = new EpMesh;
  postFileName = "";
  postFileRead = false;
  scalarFields = 0;
  scalarField = new ScalarField[MAX_SCALARS];

  // Central widget:
  //----------------
  qvtkWidget = new QVTKWidget(this);
  setCentralWidget(qvtkWidget);

  // VTK interaction:
  //------------------
  renderer = vtkRenderer::New();
  renderer->SetBackground(1, 1, 1);
  qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
  renderer->GetRenderWindow()->Render();

  // Create a cell picker and set the callback & observer:
  //------------------------------------------------------
  vtkPropPicker* propPicker = vtkPropPicker::New();
  qvtkWidget->GetInteractor()->SetPicker(propPicker);
  propPicker->Delete();

  vtkCallbackCommand* cbcPick = vtkCallbackCommand::New();
  cbcPick->SetClientData(this);
  cbcPick->SetCallback(pEventHandler);

  vtkAbstractPicker* picker = qvtkWidget->GetInteractor()->GetPicker();
  picker->AddObserver(vtkCommand::EndPickEvent, cbcPick);
  cbcPick->Delete();

  // Create the clip plane & implicit plane widget:
  //------------------------------------------------
  clipPlane = vtkPlane::New();

  vtkCallbackCommand* cbcPlane = vtkCallbackCommand::New();
  cbcPlane->SetClientData(this);
  cbcPlane->SetCallback(iEventHandler);

  planeWidget = vtkImplicitPlaneWidget::New();
  planeWidget->SetInteractor(qvtkWidget->GetInteractor());
  planeWidget->AddObserver(vtkCommand::InteractionEvent, cbcPlane);
  cbcPlane->Delete();

  SetClipPlaneOrigin(planeWidget->GetOrigin());
  SetClipPlaneNormal(planeWidget->GetNormal());

  // Python bindings:
  //-----------------
#ifdef EG_PYTHONQT
  PythonQt::init(PythonQt::IgnoreSiteModule | PythonQt::RedirectStdOut);
  mainContext = PythonQt::self()->getMainModule();
  mainContext.addObject("egp", this);
#ifdef EG_MATC
  mainContext.addObject("matc", matc);
#endif
  mainContext.addObject("surfaces", surface);
  mainContext.addObject("vectors", vector);
  mainContext.addObject("isoContours", isoContour);
  mainContext.addObject("isoSurfaces", isoSurface);
  mainContext.addObject("streamLines", streamLine);
  mainContext.addObject("preferences", preferences);
  mainContext.addObject("timeStep", timeStep);
  mainContext.addObject("colorBar", colorBar);

  console = new PythonQtScriptingConsole(NULL, mainContext);
  console->setWindowIcon(QIcon(":/icons/Mesh3D.png"));
  console->setWindowTitle("ElmerGUI PythonQt");
  console->resize(400, 300);
#endif

  // ECMAScript
  //-----------
  ecmaConsole = new EcmaConsole(NULL);
  ecmaConsole->setWindowIcon(QIcon(":/icons/Mesh3D.png"));
  ecmaConsole->setWindowTitle("ElmerGUI ECMAScript");
  ecmaConsole->resize(400, 300);

  connect(ecmaConsole, SIGNAL(cmd(QString)), this, SLOT(evaluateECMAScriptSlot(QString)));

  engine = new QScriptEngine(this);

#if QT_VERSION >= 0x040403
  QScriptValue fun = engine->newFunction(printFun);
  fun.setData(engine->newQObject(ecmaConsole));
  engine->globalObject().setProperty("print", fun);
#endif

  QScriptValue egpValue = engine->newQObject(this);
#ifdef EG_MATC
  QScriptValue matcValue = engine->newQObject(matc);
#endif
  QScriptValue surfacesValue = engine->newQObject(surface);
  QScriptValue vectorsValue = engine->newQObject(vector);
  QScriptValue isoContoursValue = engine->newQObject(isoContour);
  QScriptValue isoSurfacesValue = engine->newQObject(isoSurface);
  QScriptValue streamLinesValue = engine->newQObject(streamLine);
  QScriptValue preferencesValue = engine->newQObject(preferences);
  QScriptValue timeStepValue = engine->newQObject(timeStep);
  QScriptValue colorBarValue = engine->newQObject(colorBar);
  QScriptValue textValue = engine->newQObject(text);

  engine->globalObject().setProperty("egp", egpValue);
#ifdef EG_MATC
  engine->globalObject().setProperty("matc", matcValue);
#endif
  engine->globalObject().setProperty("surfaces", surfacesValue);
  engine->globalObject().setProperty("vectors", vectorsValue);
  engine->globalObject().setProperty("isoContours", isoContoursValue);
  engine->globalObject().setProperty("isoSurfaces", isoSurfacesValue);
  engine->globalObject().setProperty("streamLines", streamLinesValue);
  engine->globalObject().setProperty("preferences", preferencesValue);
  engine->globalObject().setProperty("timeStep", timeStepValue);
  engine->globalObject().setProperty("colorBar", colorBarValue);
  engine->globalObject().setProperty("text", textValue);

  ecmaConsole->addNames("egp", this->metaObject());
#ifdef EG_MATC
  ecmaConsole->addNames("matc", matc->metaObject());
#endif
  ecmaConsole->addNames("surfaces", surface->metaObject());
  ecmaConsole->addNames("vectors", vector->metaObject());
  ecmaConsole->addNames("isoContours", isoContour->metaObject());
  ecmaConsole->addNames("isoSurfaces", isoSurface->metaObject());
  ecmaConsole->addNames("streamLines", streamLine->metaObject());
  ecmaConsole->addNames("preferences", preferences->metaObject());
  ecmaConsole->addNames("timeStep", timeStep->metaObject());
  ecmaConsole->addNames("colorBar", colorBar->metaObject());
  ecmaConsole->addNames("text", text->metaObject());

  ecmaConsole->initCompleter();
}

VtkPost::~VtkPost()
{
}

void VtkPost::createActions()
{
  // File menu:
  //-----------
  exitAct = new QAction(QIcon(":/icons/application-exit.png"), tr("&Quit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  exitAct->setStatusTip("Quit VTK widget");
  connect(exitAct, SIGNAL(triggered()), this, SLOT(exitSlot()));

  savePictureAct = new QAction(QIcon(""), tr("Save picture as..."), this);
  savePictureAct->setStatusTip("Save picture in file");
  connect(savePictureAct, SIGNAL(triggered()), this, SLOT(savePictureSlot()));

  savePovrayAct = new QAction(QIcon(""), tr("Save povray data..."), this);
  savePovrayAct->setStatusTip("Save model data in povray-format");
  connect(savePovrayAct, SIGNAL(triggered()), this, SLOT(savePovraySlot()));

  reloadPostAct = new QAction(QIcon(""), tr("Reload"), this);
  reloadPostAct->setStatusTip("Reloads input file");
  connect(reloadPostAct, SIGNAL(triggered()), this, SLOT(reloadPostSlot()));

  readEpFileAct = new QAction(QIcon(":/icons/document-open.png"), tr("Open..."), this);
  readEpFileAct->setShortcut(tr("Ctrl+O"));
  readEpFileAct->setStatusTip("Read input file");
  connect(readEpFileAct, SIGNAL(triggered()), this, SLOT(readEpFileSlot()));

  // View menu:
  //------------
  drawMeshPointAct = new QAction(QIcon(""), tr("Mesh points"), this);
  drawMeshPointAct->setStatusTip("Draw mesh points");
  drawMeshPointAct->setCheckable(true);
  drawMeshPointAct->setChecked(false);
  connect(drawMeshPointAct, SIGNAL(triggered()), this, SLOT(drawMeshPointSlot()));
  connect(drawMeshPointAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawMeshEdgeAct = new QAction(QIcon(""), tr("Mesh edges"), this);
  drawMeshEdgeAct->setStatusTip("Draw mesh edges");
  drawMeshEdgeAct->setCheckable(true);
  drawMeshEdgeAct->setChecked(false);
  connect(drawMeshEdgeAct, SIGNAL(triggered()), this, SLOT(drawMeshEdgeSlot()));
  connect(drawMeshEdgeAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawFeatureEdgesAct = new QAction(QIcon(""), tr("Feature edges"), this);
  drawFeatureEdgesAct->setStatusTip("Draw feature edges");
  drawFeatureEdgesAct->setCheckable(true);
  drawFeatureEdgesAct->setChecked(true);
  connect(drawFeatureEdgesAct, SIGNAL(triggered()), this, SLOT(drawFeatureEdgesSlot()));
  connect(drawFeatureEdgesAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawAxesAct = new QAction(QIcon(""), tr("Coordinate axes"), this);
  drawAxesAct->setStatusTip("Draw cordinate axes");
  drawAxesAct->setCheckable(true);
  drawAxesAct->setChecked(false);
  connect(drawAxesAct, SIGNAL(triggered()), this, SLOT(drawAxesSlot()));

  drawTextAct = new QAction(QIcon(""), tr("Text..."), this);
  drawTextAct->setStatusTip("Annotate text");
  drawTextAct->setCheckable(true);
  drawTextAct->setChecked(false);
  connect(drawTextAct, SIGNAL(triggered()), this, SLOT(showTextDialogSlot()));
  connect(drawTextAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawColorBarAct = new QAction(QIcon(""), tr("Colorbar"), this);
  drawColorBarAct->setStatusTip("Draw color bar");
  drawColorBarAct->setCheckable(true);
  drawColorBarAct->setChecked(false);
  connect(drawColorBarAct, SIGNAL(triggered()), this, SLOT(showColorBarDialogSlot()));
  connect(drawColorBarAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawSurfaceAct = new QAction(QIcon(""), tr("Surfaces"), this);
  drawSurfaceAct->setStatusTip("Draw scalar fields on surfaces");
  drawSurfaceAct->setCheckable(true);
  drawSurfaceAct->setChecked(false);
  connect(drawSurfaceAct, SIGNAL(triggered()), this, SLOT(showSurfaceDialogSlot()));
  connect(drawSurfaceAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawVectorAct = new QAction(QIcon(""), tr("Vectors"), this);
  drawVectorAct->setStatusTip("Visualize vector fields by arrows");
  drawVectorAct->setCheckable(true);
  drawVectorAct->setChecked(false);
  connect(drawVectorAct, SIGNAL(triggered()), this, SLOT(showVectorDialogSlot()));
  connect(drawVectorAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawIsoContourAct = new QAction(QIcon(""), tr("Isocontours"), this);
  drawIsoContourAct->setStatusTip("Draw isocontours (2d)");
  drawIsoContourAct->setCheckable(true);
  drawIsoContourAct->setChecked(false);
  connect(drawIsoContourAct, SIGNAL(triggered()), this, SLOT(showIsoContourDialogSlot()));
  connect(drawIsoContourAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawIsoSurfaceAct = new QAction(QIcon(""), tr("Isosurfaces"), this);
  drawIsoSurfaceAct->setStatusTip("Draw isosurfaces (3d)");
  drawIsoSurfaceAct->setCheckable(true);
  drawIsoSurfaceAct->setChecked(false);
  connect(drawIsoSurfaceAct, SIGNAL(triggered()), this, SLOT(showIsoSurfaceDialogSlot()));
  connect(drawIsoSurfaceAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  drawStreamLineAct = new QAction(QIcon(""), tr("Streamlines"), this);
  drawStreamLineAct->setStatusTip("Draw stream lines");
  drawStreamLineAct->setCheckable(true);
  drawStreamLineAct->setChecked(false);
  connect(drawStreamLineAct, SIGNAL(triggered()), this, SLOT(showStreamLineDialogSlot()));
  connect(drawStreamLineAct, SIGNAL(toggled(bool)), this, SLOT(maybeRedrawSlot(bool)));

  redrawAct = new QAction(QIcon(""), tr("Redraw"), this);
  redrawAct->setShortcut(tr("Ctrl+R"));
  redrawAct->setStatusTip("Redraw");
  connect(redrawAct, SIGNAL(triggered()), this, SLOT(redrawSlot()));

  fitToWindowAct = new QAction(QIcon(""), tr("Fit to window"), this);
  fitToWindowAct->setStatusTip("Fit model to window");
  connect(fitToWindowAct, SIGNAL(triggered()), this, SLOT(fitToWindowSlot()));

  clipAllAct = new QAction(QIcon(""), tr("Clip all"), this);
  clipAllAct->setShortcut(tr("Ctrl+C"));
  clipAllAct->setStatusTip("Apply clip plane to all actors");
  clipAllAct->setCheckable(true);
  clipAllAct->setChecked(false);
  connect(clipAllAct, SIGNAL(toggled(bool)), this, SLOT(clipAllToggledSlot(bool)));

  resetModelViewAct = new QAction(QIcon(""), tr("Reset model view"), this);
  resetModelViewAct->setStatusTip("Reset model view");
  connect(resetModelViewAct, SIGNAL(triggered()), this, SLOT(resetModelViewSlot()));

  preferencesAct = new QAction(QIcon(""), tr("Preferences"), this);
  preferencesAct->setStatusTip("Show preferences");
  connect(preferencesAct, SIGNAL(triggered()), this, SLOT(showPreferencesDialogSlot()));

  // Edit menu:
  //------------
#ifdef EG_MATC
  matcAct = new QAction(QIcon(""), tr("Matc..."), this);
  matcAct->setStatusTip("Matc window");
  connect(matcAct, SIGNAL(triggered()), this, SLOT(matcOpenSlot()));
#endif

  regenerateGridsAct = new QAction(QIcon(""), tr("Regenerate all..."), this);
  regenerateGridsAct->setStatusTip("Regerate all meshes");
  connect(regenerateGridsAct, SIGNAL(triggered()), this, SLOT(regenerateGridsSlot()));

  timeStepAct = new QAction(QIcon(""), tr("Time step control"), this);
  timeStepAct->setStatusTip("Time step control");
  connect(timeStepAct, SIGNAL(triggered()), this, SLOT(showTimeStepDialogSlot()));

#ifdef EG_PYTHONQT
  showPythonQtConsoleAct = new QAction(QIcon(""), tr("PythonQt console..."), this);
  showPythonQtConsoleAct->setStatusTip("Show/hide PythonQt console");
  connect(showPythonQtConsoleAct, SIGNAL(triggered()), this, SLOT(showPythonQtConsoleSlot()));
#endif

  showECMAScriptConsoleAct = new QAction(QIcon(""), tr("ECMAScript console..."), this);
  showECMAScriptConsoleAct->setStatusTip("Show/hide ECMAScript console");
  connect(showECMAScriptConsoleAct, SIGNAL(triggered()), this, SLOT(showECMAScriptConsoleSlot()));

  // Help menu:
  //-----------
  showHelpAct = new QAction(QIcon(":/icons/help-about.png"), tr("Help..."), this);
  showHelpAct->setStatusTip("Show help dialog");
  connect(showHelpAct, SIGNAL(triggered()), this, SLOT(showHelpSlot()));
}

void VtkPost::createMenus()
{
  // File menu:
  //-----------
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(readEpFileAct);
  fileMenu->addAction(reloadPostAct);
  fileMenu->addSeparator();
  fileMenu->addAction(savePictureAct);
  fileMenu->addSeparator();
  fileMenu->addAction(savePovrayAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  // Edit menu:
  //-----------
  editMenu = menuBar()->addMenu(tr("&Edit"));
  editGroupsMenu = new QMenu(tr("Groups"));
  editMenu->addMenu(editGroupsMenu);
  editMenu->addSeparator();
  editMenu->addAction(timeStepAct);
#ifdef EG_MATC
  editMenu->addSeparator();
  editMenu->addAction( matcAct );
#endif

#ifdef EG_PYTHONQT
  editMenu->addSeparator();
  editMenu->addAction(showPythonQtConsoleAct);
#endif

  editMenu->addSeparator();
  editMenu->addAction(showECMAScriptConsoleAct);

  // View menu:
  //-----------
  viewMenu = menuBar()->addMenu(tr("&View"));
  viewMenu->addAction(drawMeshPointAct);
  viewMenu->addAction(drawMeshEdgeAct);
  viewMenu->addAction(drawFeatureEdgesAct);
  viewMenu->addAction(drawAxesAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawTextAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawSurfaceAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawIsoContourAct);
  viewMenu->addAction(drawIsoSurfaceAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawVectorAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawColorBarAct);
  viewMenu->addSeparator();
  viewMenu->addAction(drawStreamLineAct);
  viewMenu->addSeparator();
  viewMenu->addAction(preferencesAct);
  viewMenu->addSeparator();
  viewMenu->addAction(clipAllAct);
  viewMenu->addSeparator();
  viewMenu->addAction(fitToWindowAct);
  viewMenu->addAction(resetModelViewAct);
  viewMenu->addAction(redrawAct);

  // Help menu:
  //-----------
  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(showHelpAct);
}

void VtkPost::createToolbars()
{
  viewToolBar = addToolBar(tr("View"));
  viewToolBar->addAction(drawSurfaceAct);
  viewToolBar->addAction(drawVectorAct);
  viewToolBar->addAction(drawIsoContourAct);
  viewToolBar->addAction(drawIsoSurfaceAct);
  viewToolBar->addAction(drawStreamLineAct);
  viewToolBar->addSeparator();
  viewToolBar->addAction(drawColorBarAct);
  viewToolBar->addSeparator();
  viewToolBar->addAction(drawTextAct);
  viewToolBar->addSeparator();  
  viewToolBar->addAction(preferencesAct);
  viewToolBar->addSeparator();
  viewToolBar->addAction(redrawAct);
}

void VtkPost::createStatusBar()
{
}

#ifdef EG_PYTHONQT
void VtkPost::showPythonQtConsoleSlot()
{
  console->show();
}
#endif

void VtkPost::showECMAScriptConsoleSlot()
{
  ecmaConsole->clearHistory();
  ecmaConsole->show();
}

void VtkPost::evaluateECMAScriptSlot(QString cmd)
{
  if(cmd.isEmpty()) return;
  if(cmd.trimmed().toLower() == "quit") ecmaConsole->hide();
  QScriptValue val = engine->evaluate(cmd);

  QString msg = val.toString();

  if(engine->hasUncaughtException())
    msg = "Uncaught exception:" + msg;

  // do not show "undefined"
  if(msg.trimmed() != "undefined")
    ecmaConsole->append(msg);
}

#ifdef EG_MATC
QString VtkPost::MatcCmd(QString cmd)
{
   matc->ui.mcEdit->setText(cmd);
   return domatcSlot();
}

void VtkPost::matcOpenSlot()
{
  matc->show();
}

void VtkPost::matcCutPasteSlot()
{
  matc->ui.mcHistory->copy();
  matc->ui.mcEdit->clear();
  matc->ui.mcEdit->paste();
}

QString VtkPost::domatcSlot()
{
  QString res=matc->domatc(this);
  populateWidgetsSlot();
  return res;
}
#endif

void VtkPost::minMax(ScalarField *sf)
{
   sf->minVal =  9e99;
   sf->maxVal = -9e99;
   for( int i=0; i<sf->values; i++ ) {
     if ( sf->minVal>sf->value[i] ) sf->minVal=sf->value[i];
     if ( sf->maxVal<sf->value[i] ) sf->maxVal=sf->value[i];
   }

   if(sf->name == "nodes_x") {
     boundingBoxMinX = sf->minVal;
     boundingBoxMaxX = sf->maxVal;
   }

   if(sf->name == "nodes_y") {
     boundingBoxMinY = sf->minVal;
     boundingBoxMaxY = sf->maxVal;
   }

   if(sf->name == "nodes_z") {
     boundingBoxMinZ = sf->minVal;
     boundingBoxMaxZ = sf->maxVal;
   }
}

// Populate widgets in user interface dialogs:
//----------------------------------------------------------------------
void VtkPost::populateWidgetsSlot()
{
  surface->populateWidgets(this);
  vector->populateWidgets(this);
  isoContour->populateWidgets(this);
  isoSurface->populateWidgets(this);
  streamLine->populateWidgets(this);
  colorBar->populateWidgets(this);
}

// Save picture:
//----------------------------------------------------------------------
void VtkPost::savePictureSlot()
{
  QString fileName = QFileDialog::getSaveFileName(this,
	 tr("Save picture"), "", tr("Picture files (*.png)"));
  
  if(fileName.isEmpty()) {
    cout << "File name is empty" << endl;
    return;
  }

  this->SavePngFile(fileName);
}

// Save povray:
//----------------------------------------------------------------------
void VtkPost::savePovraySlot()
{
  QFile file("case.pov");

  if(!file.open(QIODevice::WriteOnly))
    return;

  QTextStream text(&file);

  // Convert to vtkPolyData with normals:
  //======================================
  vtkGeometryFilter *geometry = vtkGeometryFilter::New();
  geometry->SetInput(surfaceGrid);
  
  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
  normals->SetInputConnection(geometry->GetOutputPort());
  
  vtkPolyData *polyData = normals->GetOutput();

  polyData->Update();

  polyData->ComputeBounds();
  
  double *bounds = polyData->GetBounds();

  double x0 = (bounds[1] + bounds[0])/2.0;
  double y0 = (bounds[3] + bounds[2])/2.0;
  double z0 = (bounds[5] + bounds[4])/2.0;

  double dx = bounds[1] - bounds[0];
  double dy = bounds[3] - bounds[2];
  double dz = bounds[5] - bounds[4];

  double ds = sqrt(dx*dx+dy*dy+dz*dz);

  // Create progress dialog:
  //=========================
  QProgressDialog dialog;
  connect(this, SIGNAL(povrayState(int)), &dialog, SLOT(setValue(int)));
  dialog.show();

  // Headers etc.
  //==============
  text << "#include \"colors.inc\"\n\n"
       << "background {color White}\n\n"
       << "camera {location <1,1,1>*" << ds/1.5 << " look_at 0}\n\n"
       << "light_source {<1,2,3>*" << 5.0*ds << " color White}\n\n";

  // Mesh2:
  //========
  text << "mesh2 {\n";

  // Vertex points:
  //----------------
  dialog.setLabelText("Writing points...");
  dialog.setMaximum(polyData->GetNumberOfPoints());

  text << "  vertex_vectors {\n"
       << "    " << polyData->GetNumberOfPoints() << "\n";
  
  vtkPoints *points = polyData->GetPoints();

  for(int i = 0; i < polyData->GetNumberOfPoints(); i++) {
    emit(povrayState(i));
    double *p = points->GetPoint(i);
    text << "    <" << p[0] << "," << p[1] << "," << p[2] << ">\n";
  }
    
  text << "  }\n"; // vertex_vectors

  // Surface normals:
  //------------------
  dialog.setLabelText("Writing normals...");
  dialog.setMaximum(polyData->GetNumberOfPoints());
 
  text << "  normal_vectors {\n"
       << "    " << polyData->GetNumberOfPoints() << "\n";

  for(int i = 0; i < polyData->GetNumberOfPoints(); i++) {
    emit(povrayState(i));
    double *p = polyData->GetPointData()->GetNormals()->GetTuple(i);
    text << "    <" << p[0] << "," << p[1] << "," << p[2] << ">\n";
  }

  text << "  }\n"; // normal_vectors

  // Face indices:
  //---------------
  dialog.setLabelText("Writing faces...");
  dialog.setMaximum(polyData->GetNumberOfCells());

  text << "  face_indices {\n"
       << "    " << polyData->GetNumberOfCells() << "\n";
  
  for(int i = 0; i < polyData->GetNumberOfCells(); i++) {
    emit(povrayState(i));
    vtkCell *cell = polyData->GetCell(i);
    int n0 = cell->GetPointId(0);
    int n1 = cell->GetPointId(1);
    int n2 = cell->GetPointId(2);
    text << "    <" << n0 << "," << n1 << "," << n2 << ">\n";    
  }

  text << "  }\n"; // face_indices

  // Texture:
  //----------
  text << "  texture {\n"
       << "    pigment {\n"
       << "      rgb <0,1,1>\n"
       << "    }\n"
       << "    finish {\n"
       << "      specular 0.7\n"
       << "      roughness 0.05\n"
       << "    }\n"
       << "  }\n";

  // Translate:
  //------------
  text << "  translate <" << -x0/2.0 << "," << -y0/2.0 << "," << -z0/2.0 << ">\n";

  text << "}\n\n"; // mesh2

  // Floor:
  //--------
  text << "plane {\n"
       << "  y," << -dy/1.5 << "\n"
       << "  pigment {\n"
       << "    checker rgb<0.9,0.5,0.4> rgb<0.9,0.9,1.0> scale " << ds/5.0 << "\n"
       << "  }\n"
       << "  finish {\n"
       << "    reflection 0.0\n"
       << "    ambient 0.4\n"
       << "  }\n"
       << "}\n";
  
  // Finalize:
  //===========
  normals->Delete();
  geometry->Delete();

  file.close();
}

// Read input file (dialog):
//----------------------------------------------------------------------
void VtkPost::readEpFileSlot()
{
  readEpFile->show();
}

// Reload results:
//----------------------------------------------------------------------
void VtkPost::reloadPostSlot()
{
  if(postFileName.isEmpty()) {
    cout << "Unable to open ep-file. File name is empty." << endl;
    return;
  }

  bool surfaceVisible = drawSurfaceAct->isChecked();

  if(!ReadPostFile(postFileName))
    cout << "Reloading results from current ep-file failed." << endl;

  drawSurfaceAct->setChecked(surfaceVisible);
  
  redrawSlot();
}

// Get one line from post text stream:
//----------------------------------------------------------------------
void VtkPost::getPostLineStream(QTextStream* postStream)
{
  // postLine and postLineStream are private for VtkPost
  postLine = postStream->readLine().trimmed();
  while(postLine.isEmpty() || (postLine.at(0) == '#'))
    postLine = postStream->readLine().trimmed();
  postLineStream.setString(&postLine);
}

// Read in data (public slot):
//----------------------------------------------------------------------
bool VtkPost::ReadPostFile(QString postFileName)
{
  // Open the post file:
  //=====================
  this->postFileName = postFileName;
  this->postFileRead = false;

  QFile postFile(postFileName);

  if(!postFile.open(QIODevice::ReadOnly | QIODevice::Text))
    return false;

  cout << "Loading ep-file" << endl;

  readEpFile->ui.applyButton->setEnabled(false);
  readEpFile->ui.cancelButton->setEnabled(false);
  readEpFile->ui.okButton->setEnabled(false);
  readEpFile->setWindowTitle("Reading...");
  readEpFile->repaint();
  
  QTextStream postStream(&postFile);

  // Read in nodes, elements, timesteps, and scalar components:
  //-----------------------------------------------------------
  int nodes, elements, timesteps, components;

  getPostLineStream(&postStream);

  postLineStream >> nodes >> elements >> components >> timesteps;

  cout << "Ep file header says:" << endl;
  cout << "Nodes: " << nodes << endl;
  cout << "Elements: " << elements << endl;
  cout << "Scalar components: " << components << endl;
  cout << "Timesteps: " << timesteps << endl;


  // Read field names & set up menu actions:
  //=========================================
  if(epMesh->epNode) delete [] epMesh->epNode;
  if(epMesh->epElement) delete [] epMesh->epElement;

  for(int i = 0; i < scalarFields; i++ ) {
     ScalarField *sf = &scalarField[i];
#ifdef EG_MATC
     QByteArray nm = sf->name.trimmed().toAscii();
     var_delete( nm.data() );
#else
     if(sf->value) free(sf->value);
#endif
  }

  scalarFields = 0;

  // Add the null field:
  //--------------------
  QString fieldName = "Null";
  ScalarField* nullField = addScalarField(fieldName, nodes*timesteps, NULL);
  nullField->minVal = 0.0;
  nullField->maxVal = 0.0;

  // Add the scalar fields:
  //-----------------------
  for(int i = 0; i < components; i++) {
    QString fieldType, fieldName;
    postLineStream >> fieldType >> fieldName;

    fieldType.replace(":", "");
    fieldType = fieldType.trimmed();
    fieldName = fieldName.trimmed();

    cout << fieldType.toAscii().data() << ": ";
    cout << fieldName.toAscii().data() << endl;

    if(fieldType == "scalar")
      addScalarField(fieldName, nodes*timesteps, NULL);

    if(fieldType == "vector") {
      addVectorField(fieldName, nodes*timesteps);
      i += 2;
    }
  }

  // Nodes:
  //========
  epMesh->epNodes = nodes;
  epMesh->epNode = new EpNode[nodes];
  
  for(int i = 0; i < nodes; i++) {
    EpNode *epn = &epMesh->epNode[i];

    getPostLineStream(&postStream);

    for(int j = 0; j < 3; j++) 
      postLineStream >> epn->x[j];
  }

  // Add nodes to field variables:
  //-------------------------------
  addVectorField("nodes", nodes);
  int index = -1;
  for(int i = 0; i < scalarFields; i++) {
    ScalarField *sf = &scalarField[i];
    if(sf->name == "nodes_x") {
      index = i;
      break;
    }
  }
  ScalarField *sfx = &scalarField[index+0];
  ScalarField *sfy = &scalarField[index+1];
  ScalarField *sfz = &scalarField[index+2];

  for( int i=0; i < nodes; i++ )
  {
    sfx->value[i] = epMesh->epNode[i].x[0];
    sfy->value[i] = epMesh->epNode[i].x[1];
    sfz->value[i] = epMesh->epNode[i].x[2];
  }

  // Elements:
  //==========
  epMesh->epElements = elements;
  epMesh->epElement = new EpElement[elements];

  for(int i = 0; i < elements; i++) {
    EpElement *epe = &epMesh->epElement[i];

    getPostLineStream(&postStream);    

    postLineStream >> epe->groupName >> epe->code;
    
    epe->indexes = epe->code % 100;
    epe->index = new int[epe->indexes];
    
    for(int j = 0; j < epe->indexes; j++) {
      QString tmpString = "";
      postLineStream >> tmpString;
      if(tmpString.isEmpty()) {
	getPostLineStream(&postStream);
        postLineStream >> tmpString;
      }
      epe->index[j] = tmpString.toInt();
    }
  }

  // Data:
  //=======
  int start = readEpFile->ui.start->value() - 1;
  int end = readEpFile->ui.end->value() - 1;

  // skip values before start:
  for(int i = 0; i < nodes * start; i++) {
    if(postStream.atEnd()) break;
    getPostLineStream(&postStream);
  }

  ScalarField *sf;
  int i;
  for(i = 0; i < nodes * (end - start + 1); i++) {
    if(postStream.atEnd()) break;
    getPostLineStream(&postStream);

    for(int j = 0; j < scalarFields-4; j++) { // - 4 = no nodes, no null field
      sf = &scalarField[j+1];                 // + 1 = skip null field
      postLineStream >> sf->value[i];
    }
  }

  int real_timesteps = i/nodes;
  cout << real_timesteps << " timesteps read in." << endl;

  // Initial min & max values:
  //============================
  int ifield=0, size;
  while( ifield<scalarFields )
  {
    ScalarField *sf = &scalarField[ifield];

    int sf_timesteps = sf->values/nodes;
    if ( real_timesteps < sf_timesteps )
    {
      sf->values = real_timesteps*nodes;
#ifdef EG_MATC
      QString name=sf->name;
      int n = sf->name.indexOf("_x");
      if ( n>0 ) 
      {
        name = sf->name.mid(0,n);

        QString cmd = name+"="+name+"(0:2,0:"+QString::number(sf->values-1)+")";
        mtc_domath(cmd.toAscii().data());
        VARIABLE *var = var_check(name.toAscii().data());

        sf = &scalarField[ifield];
        sf->value = &M(var,0,0);
        minMax(sf);

        ifield++;
        sf = &scalarField[ifield];
        sf->value = &M(var,1,0);
        sf->values = real_timesteps*nodes;
        minMax(sf);

        ifield++;
        sf = &scalarField[ifield];
        sf->value = &M(var,2,0);
        sf->values = real_timesteps*nodes;
        minMax(sf);
      } else {
        size=sf->values*sizeof(double);

        VARIABLE *var = var_check(name.toAscii().data());
        sf->value = (double *)ALLOC_PTR(realloc(
              ALLOC_LST(sf->value), ALLOC_SIZE(size)) );
        MATR(var) = sf->value;
        NCOL(var) = sf->values;
        minMax(sf);
      }
#else
      size = sf->values*sizeof(double);
      sf->value = (double *)realloc(sf->value,size);
      minMax(sf);
#endif
    } else {
      minMax(sf);
    }
    ifield++;
  }

  timesteps = real_timesteps;
  timeStep->maxSteps = timesteps;
  timeStep->ui.start->setValue(1);
  timeStep->ui.stop->setValue(timesteps);
  
  postFile.close();

  // Set up the group edit menu:
  //=============================
  groupActionHash.clear();
  editGroupsMenu->clear();

  for(int i = 0; i < elements; i++) {
    EpElement* epe = &epMesh->epElement[i];

    QString groupName = epe->groupName;
    
    if(groupActionHash.contains(groupName))
      continue;

    QAction* groupAction = new QAction(groupName, this);
    groupAction->setCheckable(true);
    groupAction->setChecked(true);
    editGroupsMenu->addAction(groupAction);
    groupActionHash.insert(groupName, groupAction);
  }

  // Populate the widgets in user interface dialogs:
  //-------------------------------------------------
  populateWidgetsSlot();

  this->postFileRead = true;

  groupChangedSlot(NULL);
  connect(editGroupsMenu, SIGNAL(triggered(QAction*)), this, SLOT(groupChangedSlot(QAction*)));
  
  editGroupsMenu->addSeparator();
  editGroupsMenu->addAction(regenerateGridsAct);

  // Set the null field active:
  //---------------------------
  drawSurfaceAct->setChecked(true);

  renderer->ResetCamera();
  
  readEpFile->ui.fileName->setText(postFileName);
  readEpFile->readHeader();
  readEpFile->ui.applyButton->setEnabled(true);
  readEpFile->ui.cancelButton->setEnabled(true);
  readEpFile->ui.okButton->setEnabled(true);
  readEpFile->setWindowTitle("Read input file");
  readEpFile->repaint();

  redrawSlot();

  renderer->GetActiveCamera()->GetPosition(initialCameraPosition);
  initialCameraRoll = renderer->GetActiveCamera()->GetRoll();

  return true;
}

void VtkPost::addVectorField(QString fieldName, int values)
{
   
#ifdef EG_MATC
    QByteArray nm=fieldName.trimmed().toAscii();

    char *name = (char *)malloc( nm.count()+1 );
    strcpy(name,nm.data());

    VARIABLE *var = var_check(name);
    if ( !var || NROW(var) != 3 || NCOL(var) != values )
      var = var_new( name, TYPE_DOUBLE, 3, values );
    free(name);

   addScalarField(fieldName+"_x", values, &M(var,0,0));
   addScalarField(fieldName+"_y", values, &M(var,1,0));
   addScalarField(fieldName+"_z", values, &M(var,2,0));
#else
   addScalarField(fieldName+"_x", values, NULL);
   addScalarField(fieldName+"_y", values, NULL);
   addScalarField(fieldName+"_z", values, NULL);
#endif
}

// Add a scalar field:
//----------------------------------------------------------------------
ScalarField* VtkPost::addScalarField(QString fieldName, int values, double *value)
{
  if(scalarFields >= MAX_SCALARS) {
    cout << "Max. scalar limit exceeded!" << endl;
    return NULL;
  }

  ScalarField *sf = &scalarField[scalarFields++];
  sf->name = fieldName;
  sf->values = values;
  sf->value = value;
 
  if ( !sf->value ) {
#ifdef EG_MATC
    QByteArray nm=fieldName.trimmed().toAscii();

    char *name = (char *)malloc( nm.count()+1 );
    strcpy(name,nm.data());
    VARIABLE *var = var_check(name);
    if ( !var || NROW(var)!=1 || NCOL(var) != values )
      var = var_new( name, TYPE_DOUBLE, 1, values );
    sf->value = MATR(var);
    free(name);
#else
    sf->value = (double *)calloc(values,sizeof(double));
#endif
  }

  sf->minVal = +9.9e99;
  sf->maxVal = -9.9e99;

  return sf;
}

// Close the widget:
//----------------------------------------------------------------------
void VtkPost::exitSlot()
{
  close();
}

// Group selection changed:
//----------------------------------------------------------------------
void VtkPost::regenerateGridsSlot()
{
  groupChangedSlot(NULL);
}

void VtkPost::groupChangedSlot(QAction* groupAction)
{
  // Status of groupAction has changed: regenerate grids
  //-----------------------------------------------------
  volumeGrid->Delete();
  surfaceGrid->Delete();
  lineGrid->Delete();

  volumeGrid = vtkUnstructuredGrid::New();
  surfaceGrid = vtkUnstructuredGrid::New();
  lineGrid = vtkUnstructuredGrid::New();

  // Points:
  //---------
  int index = -1;
  for(int i = 0; i < scalarFields; i++) {
    ScalarField* sf = &scalarField[i];
    if(sf->name == "nodes_x") {
      index = i;
      break;
    }
  }
  
  if((index < 0) || (index + 2 > scalarFields - 1)) return;

  double x[3];
  ScalarField* sfx = &scalarField[index+0];
  ScalarField* sfy = &scalarField[index+1];
  ScalarField* sfz = &scalarField[index+2];
  
  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(epMesh->epNodes);

  for(int i = 0; i < epMesh->epNodes; i++) {
    x[0] = sfx->value[i];
    x[1] = sfy->value[i];
    x[2] = sfz->value[i];
    points->InsertPoint(i, x);
  }
  volumeGrid->SetPoints(points);
  surfaceGrid->SetPoints(points);
  lineGrid->SetPoints(points);
  points->Delete();

  // Volume grid:
  //---------------
  vtkTetra* tetra = vtkTetra::New();
  vtkHexahedron* hexa = vtkHexahedron::New();

  for(int i = 0; i < epMesh->epElements; i++) {
    EpElement* epe = &epMesh->epElement[i];

    if(epe->code == 504) {
      QString groupName = epe->groupName;
      if(groupName.isEmpty()) continue;

      QAction* groupAction = groupActionHash.value(groupName);
      if(groupAction == NULL) continue;
      
      for(int j = 0; j < 4; j++)
	tetra->GetPointIds()->SetId(j, epe->index[j]);
      
      if(groupAction->isChecked())
	volumeGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }
    
    if(epe->code == 808) {
      QString groupName = epe->groupName;
      if(groupName.isEmpty()) continue;
      
      QAction* groupAction = groupActionHash.value(groupName);
      if(groupAction == NULL) continue;
      
      for(int j = 0; j < 8; j++)
	hexa->GetPointIds()->SetId(j, epe->index[j]);
      
      if(groupAction->isChecked())
	volumeGrid->InsertNextCell(hexa->GetCellType(), hexa->GetPointIds());
    }
  }
  tetra->Delete();
  hexa->Delete();

  // Surface grid:
  //---------------
  vtkTriangle* tria = vtkTriangle::New();
  vtkQuad* quad = vtkQuad::New();
  for(int i = 0; i < epMesh->epElements; i++) {
    EpElement* epe = &epMesh->epElement[i];

    if(epe->code == 303) {
      QString groupName = epe->groupName;
      if(groupName.isEmpty()) continue;

      QAction* groupAction = groupActionHash.value(groupName);
      if(groupAction == NULL) continue;
      
      for(int j = 0; j < 3; j++)
	tria->GetPointIds()->SetId(j, epe->index[j]);
      
      if(groupAction->isChecked())
	surfaceGrid->InsertNextCell(tria->GetCellType(), tria->GetPointIds());
    }

    if(epe->code == 404) {
      QString groupName = epe->groupName;
      if(groupName.isEmpty()) continue;

      QAction* groupAction = groupActionHash.value(groupName);
      if(groupAction == NULL) continue;
      
      for(int j = 0; j < 4; j++)
	quad->GetPointIds()->SetId(j, epe->index[j]);
      
      if(groupAction->isChecked())
	surfaceGrid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());
    }

  }
  tria->Delete();
  quad->Delete();

  // Line grid:
  //---------------
  vtkLine* line = vtkLine::New();
  for(int i = 0; i < epMesh->epElements; i++) {
    EpElement* epe = &epMesh->epElement[i];

    if(epe->code == 202) {
      QString groupName = epe->groupName;
      if(groupName.isEmpty()) continue;

      QAction* groupAction = groupActionHash.value(groupName);
      if(groupAction == NULL) continue;
      
      for(int j = 0; j < 3; j++)
	line->GetPointIds()->SetId(j, epe->index[j]);
      
      if(groupAction->isChecked())
	lineGrid->InsertNextCell(line->GetCellType(), line->GetPointIds());
    }
  }
  line->Delete();

  if(timeStep->ui.regenerateBeforeDrawing->isChecked()) return;

  redrawSlot();

  // Place the implicit plane widget:
  //---------------------------------
  double bounds[6];
  GetBounds(bounds);

  double origin[3];
  origin[0] = (bounds[0] + bounds[1]) / 2.0;
  origin[1] = (bounds[2] + bounds[3]) / 2.0;
  origin[2] = (bounds[4] + bounds[5]) / 2.0;

  planeWidget->SetPlaceFactor(1.5);
  planeWidget->PlaceWidget(bounds);
  planeWidget->SetOrigin(origin);
  planeWidget->GetEdgesProperty()->SetColor(0, 0, 0);
  planeWidget->GetPlaneProperty()->SetColor(1, 0, 0);
  planeWidget->GetPlaneProperty()->SetOpacity(0.2);
  planeWidget->GetSelectedPlaneProperty()->SetColor(0, 1, 0);
  planeWidget->GetSelectedPlaneProperty()->SetOpacity(0.1);

  SetClipPlaneOrigin(planeWidget->GetOrigin());
  SetClipPlaneNormal(planeWidget->GetNormal());
}

// Show preferences dialog:
//----------------------------------------------------------------------
void VtkPost::showPreferencesDialogSlot()
{
  if(!postFileRead) return;
  preferences->show();
}

// Maybe redraw:
//----------------------------------------------------------------------
void VtkPost::maybeRedrawSlot(bool value)
{
  if(!value) redrawSlot();
}

// Redraw:
//----------------------------------------------------------------------
void VtkPost::redrawSlot()
{
  if(!postFileRead) return;

#ifdef EG_MATC
   VARIABLE *tvar = var_check((char *)"t");
   if (!tvar) tvar=var_new((char *)"t", TYPE_DOUBLE,1,1 );
   M(tvar,0,0) = (double)timeStep->ui.timeStep->value();

   QString dosome = timeStep->ui.doBefore->text();
   matc->ui.mcEdit->clear();
   matc->ui.mcEdit->insert(dosome);
   matc->domatc(this);
#endif

   if(timeStep->ui.regenerateBeforeDrawing->isChecked())
     regenerateGridsSlot();

  drawMeshPointSlot();
  drawMeshEdgeSlot();
  drawFeatureEdgesSlot();
  drawSurfaceSlot();
  drawVectorSlot();
  drawIsoContourSlot();
  drawIsoSurfaceSlot();
  drawStreamLineSlot();
  drawColorBarSlot();
  drawAxesSlot();
  drawTextSlot();

  vtkRenderWindow *renderWindow = qvtkWidget->GetRenderWindow();
  renderWindow->Render();

  // Check if the "Stop" button of time stepping loop has been pressed:
  QCoreApplication::processEvents();

  emit(canProceedWithNextSignal(renderWindow));

  planeWidget->GetEdgesProperty()->SetColor(0, 0, 0);
  planeWidget->GetOutlineProperty()->SetColor(0, 0, 0);
  planeWidget->GetNormalProperty()->SetColor(1, 0, 0);  
}

// Draw color bar:
//----------------------------------------------------------------------
void VtkPost::showColorBarDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawColorBarAct->isChecked()) {
    colorBar->show();
  } else {
    colorBar->close();
    drawColorBarSlot();
  }
}

void VtkPost::hideColorBarSlot()
{
  drawColorBarAct->setChecked(false);
  drawColorBarSlot();
}

void VtkPost::drawColorBarSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(colorBarActor);
  if(!drawColorBarAct->isChecked()) return;
  colorBar->draw(this);
  renderer->AddActor(colorBarActor);
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw mesh points:
//----------------------------------------------------------------------
void VtkPost::drawMeshPointSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(meshPointActor);
  if(!drawMeshPointAct->isChecked()) return;
  meshPoint->draw(this, preferences);
  renderer->AddActor(meshPointActor);
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw mesh edges:
//----------------------------------------------------------------------
void VtkPost::drawMeshEdgeSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(meshEdgeActor);
  if(!drawMeshEdgeAct->isChecked()) return;
  meshEdge->draw(this, preferences);
  renderer->AddActor(meshEdgeActor);
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw feature edges:
//----------------------------------------------------------------------
void VtkPost::drawFeatureEdgesSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(featureEdgeActor);
  if(!drawFeatureEdgesAct->isChecked()) return;
  featureEdge->draw(this, preferences);
  renderer->AddActor(featureEdgeActor);
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw stream lines:
//----------------------------------------------------------------------
void VtkPost::showStreamLineDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();
  
  if(drawStreamLineAct->isChecked()) {
    streamLine->show();
  } else {
    streamLine->close();
    drawStreamLineSlot();
  }
}

void VtkPost::hideStreamLineSlot()
{
  drawStreamLineAct->setChecked(false);
  drawStreamLineSlot();
}

void VtkPost::drawStreamLineSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(streamLineActor);
  if(!drawStreamLineAct->isChecked()) return;
  streamLine->draw(this, timeStep);
  renderer->AddActor(streamLineActor);
  drawColorBarSlot();
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw vectors:
//----------------------------------------------------------------------
void VtkPost::showVectorDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawVectorAct->isChecked()) {
    vector->show();
  } else {
    vector->close();
    drawVectorSlot();
  }
}

void VtkPost::hideVectorSlot()
{
  drawVectorAct->setChecked(false);
  drawVectorSlot();
}

void VtkPost::drawVectorSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(vectorActor);
  if(!drawVectorAct->isChecked()) return;
  vector->draw(this, timeStep);
  renderer->AddActor(vectorActor);
  drawColorBarSlot();
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw surfaces:
//----------------------------------------------------------------------
void VtkPost::showSurfaceDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawSurfaceAct->isChecked()) {
    surface->show();
  } else {
    surface->close();
    drawSurfaceSlot();
  }
}

void VtkPost::hideSurfaceSlot()
{
  drawSurfaceAct->setChecked(false);
  drawSurfaceSlot();
}

void VtkPost::drawSurfaceSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(surfaceActor);
  if(!drawSurfaceAct->isChecked()) return;
  surface->draw(this, timeStep);
  renderer->AddActor(surfaceActor);
  drawColorBarSlot();
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw iso contours (2D):
//----------------------------------------------------------------------
void VtkPost::showIsoContourDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawIsoContourAct->isChecked()) {
    isoContour->show();
  } else {
    isoContour->close();
    drawIsoContourSlot();
  }
}

void VtkPost::hideIsoContourSlot()
{
  drawIsoContourAct->setChecked(false);
  drawIsoContourSlot();
}

void VtkPost::drawIsoContourSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(isoContourActor);
  if(!drawIsoContourAct->isChecked()) return;
  isoContour->draw(this, timeStep);
  renderer->AddActor(isoContourActor);
  drawColorBarSlot();  
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw isosurfaces (3D):
//----------------------------------------------------------------------
void VtkPost::showIsoSurfaceDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawIsoSurfaceAct->isChecked()) {
    isoSurface->show();
  } else {
    isoSurface->close();
    drawIsoSurfaceSlot();
  }
}

void VtkPost::hideIsoSurfaceSlot()
{
  drawIsoSurfaceAct->setChecked(false);
  drawIsoSurfaceSlot();
}

void VtkPost::drawIsoSurfaceSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(isoSurfaceActor);
  if(!drawIsoSurfaceAct->isChecked()) return;
  isoSurface->draw(this, timeStep);
  renderer->AddActor(isoSurfaceActor);
  drawColorBarSlot();
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw axes:
//----------------------------------------------------------------------
void VtkPost::drawAxesSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor(axesActor);
  renderer->RemoveActor(axesXTextActor);
  renderer->RemoveActor(axesYTextActor);
  renderer->RemoveActor(axesZTextActor);
  if(!drawAxesAct->isChecked()) return;
  axes->draw(this);
  renderer->AddActor(axesActor);
  renderer->AddActor(axesXTextActor);
  renderer->AddActor(axesYTextActor);
  renderer->AddActor(axesZTextActor);
  qvtkWidget->GetRenderWindow()->Render();
}

// Draw text:
//----------------------------------------------------------------------
void VtkPost::showTextDialogSlot()
{
  if(!postFileRead) return;
  qvtkWidget->GetRenderWindow()->Render();

  if(drawTextAct->isChecked()) {
    text->show();
  } else {
    text->close();
    drawTextSlot();
  }
}

void VtkPost::hideTextSlot()
{
  drawTextAct->setChecked(false);
  drawTextSlot();
}

void VtkPost::drawTextSlot()
{
  if(!postFileRead) return;
  renderer->RemoveActor2D(textActor);
  if(!drawTextAct->isChecked()) return;
  text->draw(this);
  renderer->AddActor2D(textActor);
  qvtkWidget->GetRenderWindow()->Render();
}


// Time step control:
//----------------------------------------------------------------------
void VtkPost::showTimeStepDialogSlot()
{
  if(!postFileRead) return;
  timeStep->show();
}

void VtkPost::timeStepChangedSlot()
{
  redrawSlot();
}

// Fit to window:
//----------------------------------------------------------------------
void VtkPost::fitToWindowSlot()
{
  if(!postFileRead) return;
  renderer->ResetCamera();  
}

// Reset model view:
//----------------------------------------------------------------------
void VtkPost::resetModelViewSlot()
{
  if(!postFileRead) return;
  SetInitialCameraPosition();
}

// Clip all -action toggled:
//----------------------------------------------------------------------
void VtkPost::clipAllToggledSlot(bool)
{
  redrawSlot();
}

// Other public methods:
//----------------------------------------------------------------------
QVTKWidget* VtkPost::GetQVTKWidget()
{
  return qvtkWidget;
}

vtkRenderer* VtkPost::GetRenderer()
{
  return renderer;
}

vtkActor* VtkPost::GetSurfaceActor()
{
  return surfaceActor;
}

vtkActor* VtkPost::GetVectorActor()
{
  return vectorActor;
}

vtkActor* VtkPost::GetIsoContourActor()
{
  return isoContourActor;
}

vtkActor* VtkPost::GetIsoSurfaceActor()
{
  return isoSurfaceActor;
}

vtkActor* VtkPost::GetStreamLineActor()
{
  return streamLineActor;
}

vtkScalarBarActor* VtkPost::GetColorBarActor()
{
  return colorBarActor;
}

vtkActor* VtkPost::GetPickedPointActor()
{
  return pickedPointActor;
}

vtkActor* VtkPost::GetAxesActor()
{
  return axesActor;
}

vtkActor* VtkPost::GetFeatureEdgeActor()
{
  return featureEdgeActor;
}

vtkActor* VtkPost::GetMeshPointActor()
{
  return meshPointActor;
}

vtkActor* VtkPost::GetMeshEdgeActor()
{
  return meshEdgeActor;
}

vtkFollower* VtkPost::GetAxesXTextActor()
{
  return axesXTextActor;
}

vtkFollower* VtkPost::GetAxesYTextActor()
{
  return axesYTextActor;
}

vtkFollower* VtkPost::GetAxesZTextActor()
{
  return axesZTextActor;
}

vtkTextActor* VtkPost::GetTextActor()
{
  return textActor;
}

void VtkPost::GetBounds(double* bounds)
{
  volumeGrid->GetBounds(bounds);
}

vtkUnstructuredGrid* VtkPost::GetLineGrid()
{
  return lineGrid;
}

vtkUnstructuredGrid* VtkPost::GetSurfaceGrid()
{
  return surfaceGrid;
}

vtkUnstructuredGrid* VtkPost::GetVolumeGrid()
{
  return volumeGrid;
}

vtkImplicitPlaneWidget* VtkPost::GetPlaneWidget()
{
  return planeWidget;
}

vtkPlane* VtkPost::GetClipPlane()
{
  double px = preferences->ui.clipPointX->text().toDouble();
  double py = preferences->ui.clipPointY->text().toDouble();
  double pz = preferences->ui.clipPointZ->text().toDouble();
  double nx = preferences->ui.clipNormalX->text().toDouble();
  double ny = preferences->ui.clipNormalY->text().toDouble();
  double nz = preferences->ui.clipNormalZ->text().toDouble();

  this->clipPlane->SetOrigin(px, py, pz);
  this->clipPlane->SetNormal(nx, ny, nz);

  return clipPlane;
}

void VtkPost::SetClipPlaneOrigin(double* origin)
{
  clipPlane->SetOrigin(origin);
  preferences->ui.clipPointX->setText(QString::number(origin[0]));
  preferences->ui.clipPointY->setText(QString::number(origin[1]));
  preferences->ui.clipPointZ->setText(QString::number(origin[2]));
}

void VtkPost::SetClipPlaneNormal(double* normal)
{
  clipPlane->SetNormal(normal);
  preferences->ui.clipNormalX->setText(QString::number(normal[0]));
  preferences->ui.clipNormalY->setText(QString::number(normal[1]));
  preferences->ui.clipNormalZ->setText(QString::number(normal[2]));
}

int VtkPost::NofNodes()
{
   return epMesh->epNodes;
}

vtkLookupTable* VtkPost::GetCurrentLut()
{
  return currentLut;
}

QString VtkPost::GetCurrentSurfaceName()
{
  return currentSurfaceName;
}

QString VtkPost::GetCurrentVectorName()
{
  return currentVectorName;
}

QString VtkPost::GetCurrentVectorColorName()
{
  return currentVectorColorName;
}

QString VtkPost::GetCurrentIsoContourName()
{
  return currentIsoContourName;
}

QString VtkPost::GetCurrentIsoContourColorName()
{
  return currentIsoContourColorName;
}

QString VtkPost::GetCurrentIsoSurfaceName()
{
  return currentIsoSurfaceName;
}

QString VtkPost::GetCurrentIsoSurfaceColorName()
{
  return currentIsoContourColorName;
}

QString VtkPost::GetCurrentStreamLineName()
{
  return currentStreamLineName;
}

QString VtkPost::GetCurrentStreamLineColorName()
{
  return currentStreamLineColorName;
}

void VtkPost::SetCurrentSurfaceName(QString name)
{
  currentSurfaceName = name;
}

void VtkPost::SetCurrentVectorName(QString name)
{
  currentVectorName = name;
}

void VtkPost::SetCurrentVectorColorName(QString name)
{
  currentVectorName = name;
}

void VtkPost::SetCurrentIsoContourName(QString name)
{
  currentIsoContourName = name;
}

void VtkPost::SetCurrentIsoContourColorName(QString name)
{
  currentIsoContourColorName = name;
}

void VtkPost::SetCurrentIsoSurfaceName(QString name)
{
  currentIsoSurfaceName = name;
}

void VtkPost::SetCurrentIsoSurfaceColorName(QString name)
{
  currentIsoSurfaceColorName = name;
}

void VtkPost::SetCurrentStreamLineName(QString name)
{
  currentStreamLineName = name;
}

void VtkPost::SetCurrentStreamLineColorName(QString name)
{
  currentStreamLineColorName = name;
}

int VtkPost::GetScalarFields()
{
  return scalarFields;
}

void VtkPost::SetScalarFields(int n)
{
  scalarFields = n;
}

ScalarField* VtkPost::GetScalarField()
{
  return scalarField;
}

EpMesh* VtkPost::GetEpMesh()
{
  return epMesh;
}

double* VtkPost::GetCurrentPickPosition()
{
  return &currentPickPosition[0];
}

void VtkPost::SetCurrentPickPosition(double *p)
{
  currentPickPosition[0] = p[0];
  currentPickPosition[1] = p[1];
  currentPickPosition[2] = p[2];
}

Preferences* VtkPost::GetPreferences()
{
  return preferences;
}

QSize VtkPost::minimumSizeHint() const
{
  return QSize(64, 64);
}

QSize VtkPost::sizeHint() const
{
  return QSize(640, 480);
}

void VtkPost::showHelpSlot()
{
  QMessageBox::about(this, tr("ElmerGUI postprocessor (beta)"),
		     tr("Press 'p' to pick a point\n"
                        "Press 'i' to show/hide the interactive plane widget\n"
			"Press 'w' to show the results in wireframe mode\n"
			"Press 's' to show the results in surface mode"));
}

//------------------------------------------------------------
//                       Public slots:
//------------------------------------------------------------
void VtkPost::SetPostFileStart(int n)
{
  readEpFile->ui.start->setValue(n);
}

void VtkPost::SetPostFileEnd(int n)
{
  readEpFile->ui.end->setValue(n);
}

// ReadPostFile(QString) defined above

//------------------------------------------------------------
void VtkPost::Render()
{
  qvtkWidget->GetRenderWindow()->Render();
}

void VtkPost::Redraw()
{
  this->redrawSlot();
}

void VtkPost::ResetCamera()
{
  renderer->ResetCamera();
}

void VtkPost::ResetAll()
{
  SetInitialCameraPosition();
  SetOrientation(0, 0, 0);
  SetOrigin(0, 0, 0);
  SetScale(1, 1, 1);
}

//------------------------------------------------------------
void VtkPost::SetText(bool b)
{
  drawTextAct->setChecked(b);
  drawTextSlot();
}

void VtkPost::SetSurfaces(bool b)
{
  drawSurfaceAct->setChecked(b);
  drawSurfaceSlot();
}

void VtkPost::SetVectors(bool b)
{
  drawVectorAct->setChecked(b);
  drawVectorSlot();
}

void VtkPost::SetIsoContours(bool b)
{
  drawIsoContourAct->setChecked(b);
  drawIsoContourSlot();
}

void VtkPost::SetIsoSurfaces(bool b)
{
  drawIsoSurfaceAct->setChecked(b);
  drawIsoSurfaceSlot();
}

void VtkPost::SetStreamLines(bool b)
{
  drawStreamLineAct->setChecked(b);
  drawStreamLineSlot();
}

void VtkPost::SetColorBar(bool b)
{
  drawColorBarAct->setChecked(b);
  drawColorBarSlot();
}

void VtkPost::SetMeshPoints(bool b)
{
  drawMeshPointAct->setChecked(b);
  drawMeshPointSlot();
}

void VtkPost::SetMeshEdges(bool b)
{
  drawMeshEdgeAct->setChecked(b);
  drawMeshEdgeSlot();
}

void VtkPost::SetFeatureEdges(bool b)
{
  drawFeatureEdgesAct->setChecked(b);
  drawFeatureEdgesSlot();
}

void VtkPost::SetAxes(bool b)
{
  drawAxesAct->setChecked(b);
  drawAxesSlot();
}

//------------------------------------------------------------
bool VtkPost::GetClipAll()
{
  return clipAllAct->isChecked();
}

void VtkPost::SetClipAll(bool clip)
{
  clipAllAct->setChecked(clip);
}

void VtkPost::SetClipPlaneOx(double x)
{
  preferences->ui.clipPointX->setText(QString::number(x));
  GetClipPlane();
}

void VtkPost::SetClipPlaneOy(double y)
{
  preferences->ui.clipPointY->setText(QString::number(y));
  GetClipPlane();
}

void VtkPost::SetClipPlaneOz(double z)
{
  preferences->ui.clipPointZ->setText(QString::number(z));  
  GetClipPlane();
}

void VtkPost::SetClipPlaneNx(double x)
{
  preferences->ui.clipNormalX->setText(QString::number(x));  
  GetClipPlane();
}

void VtkPost::SetClipPlaneNy(double y)
{
  preferences->ui.clipNormalY->setText(QString::number(y));
  GetClipPlane();
}

void VtkPost::SetClipPlaneNz(double z)
{
  preferences->ui.clipNormalZ->setText(QString::number(z));  
  GetClipPlane();
}

//------------------------------------------------------------
double VtkPost::GetCameraDistance()
{
  return renderer->GetActiveCamera()->GetDistance();
}

void VtkPost::SetCameraDistance(double f)
{
  renderer->GetActiveCamera()->SetDistance(f);
}

double VtkPost::GetCameraPositionX()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  return x;
}

double VtkPost::GetCameraPositionY()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  return y;
}

double VtkPost::GetCameraPositionZ()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  return z;
}

void VtkPost::SetCameraPositionX(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  renderer->GetActiveCamera()->SetPosition(f, y, z);
}

void VtkPost::SetCameraPositionY(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  renderer->GetActiveCamera()->SetPosition(x, f, z);
}

void VtkPost::SetCameraPositionZ(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetPosition(x, y, z);
  renderer->GetActiveCamera()->SetPosition(x, y, f);
}

double VtkPost::GetCameraFocalPointX()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  return x;
}

double VtkPost::GetCameraFocalPointY()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  return y;
}

double VtkPost::GetCameraFocalPointZ()
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  return z;
}

void VtkPost::SetCameraFocalPointX(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  renderer->GetActiveCamera()->SetFocalPoint(f, y, z);
}

void VtkPost::SetCameraFocalPointY(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  renderer->GetActiveCamera()->SetFocalPoint(x, f, z);
}

void VtkPost::SetCameraFocalPointZ(double f)
{
  double x, y, z;
  renderer->GetActiveCamera()->GetFocalPoint(x, y, z);
  renderer->GetActiveCamera()->SetFocalPoint(x, y, f);
}

void VtkPost::CameraDolly(double f)
{
  renderer->GetActiveCamera()->Dolly(f);
}

void VtkPost::CameraRoll(double f)
{
  renderer->GetActiveCamera()->Roll(f);
}

void VtkPost::CameraAzimuth(double f)
{
  renderer->GetActiveCamera()->Azimuth(f);
}

void VtkPost::CameraElevation(double f)
{
  renderer->GetActiveCamera()->Elevation(f);
}

void VtkPost::CameraPitch(double f)
{
  renderer->GetActiveCamera()->Pitch(f);
}

void VtkPost::CameraZoom(double f)
{
  renderer->GetActiveCamera()->Zoom(f);
}

void VtkPost::CameraYaw(double f)
{
  renderer->GetActiveCamera()->Yaw(f);
}

void VtkPost::SetCameraRoll(double f)
{
  renderer->GetActiveCamera()->SetRoll(f);
}

void VtkPost::SetInitialCameraPosition()
{
  renderer->ResetCamera();  
  renderer->GetActiveCamera()->SetPosition(initialCameraPosition[0], 
					   initialCameraPosition[1], 
					   initialCameraPosition[2]);
  renderer->GetActiveCamera()->SetRoll(initialCameraRoll);
}

//------------------------------------------------------------
void VtkPost::RotateX(double f)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->RotateX(f);
  }
}

void VtkPost::RotateY(double f)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->RotateY(f);
  }
}

void VtkPost::RotateZ(double f)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->RotateZ(f);
  }
}

void VtkPost::SetOrientation(double x, double y, double z)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->SetOrientation(x, y, z);
  }
}

void VtkPost::SetPositionX(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetPosition(x);
    actor->SetPosition(f, x[1], x[2]);
  }
}

void VtkPost::SetPositionY(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetPosition(x);
    actor->SetPosition(x[0], f, x[2]);
  }
}

void VtkPost::SetPositionZ(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetPosition(x);
    actor->SetPosition(x[0], x[1], f);
  }
}

void VtkPost::SetPosition(double x, double y, double z)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->SetPosition(x, y, z);
  }
}

void VtkPost::AddPosition(double x, double y, double z)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->AddPosition(x, y, z);
  }
}

void VtkPost::SetOrigin(double x, double y, double z)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->SetOrigin(x, y, z);
  }
}

void VtkPost::SetScaleX(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetScale(x);
    actor->SetScale(f, x[1], x[2]);
  }
}

void VtkPost::SetScaleY(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetScale(x);
    actor->SetScale(x[0], f, x[2]);
  }
}

void VtkPost::SetScaleZ(double f)
{
  double x[3];
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->GetScale(x);
    actor->SetScale(x[0], x[1], f);
  }
}

void VtkPost::SetScale(double x, double y, double z)
{
  vtkActorCollection* actors = renderer->GetActors();
  actors->InitTraversal();
  vtkActor* actor = NULL;
  while((actor = actors->GetNextActor()) != NULL) {
    actor->SetScale(x, y, z);
  }
}

//----------------------------------------------------------------------
double VtkPost::GetLength()
{
  double volumeLength = volumeGrid->GetLength();
  double surfaceLength = surfaceGrid->GetLength();
  double lineLength = lineGrid->GetLength();

  double length = volumeLength;
  if(surfaceLength > length) length = surfaceLength;
  if(lineLength > length) length = lineLength;

  return length;
}

double VtkPost::GetNofNodes()
{
   return epMesh->epNodes;
}

double VtkPost::GetMinX()
{
  return boundingBoxMinX;
}

double VtkPost::GetMaxX()
{
  return boundingBoxMaxX;
}

double VtkPost::GetMinY()
{
  return boundingBoxMinY;
}

double VtkPost::GetMaxY()
{
  return boundingBoxMaxY;
}

double VtkPost::GetMinZ()
{
  return boundingBoxMinZ;
}

double VtkPost::GetMaxZ()
{
  return boundingBoxMaxZ;
}

//----------------------------------------------------------------------
bool VtkPost::SavePngFile(QString fileName)
{ 
  if(fileName.isEmpty()) return false;

  vtkWindowToImageFilter* image = vtkWindowToImageFilter::New();

  image->SetInput(qvtkWidget->GetRenderWindow());
  image->Update();

  vtkPNGWriter* writer = vtkPNGWriter::New();

  writer->SetInputConnection(image->GetOutputPort());
  writer->SetFileName(fileName.toAscii().data());

  qvtkWidget->GetRenderWindow()->Render();
  writer->Write();

  writer->Delete();
  image->Delete();

  return true;
}

//----------------------------------------------------------------------
void VtkPost::SetBgColor(double r, double g, double b)
{
  double R = r;
  double G = g;
  double B = b;

  if(R < 0.0) R = 0.0;
  if(R > 1.0) R = 1.0;
  
  if(G < 0.0) G = 0.0;
  if(G > 1.0) G = 1.0;
  
  if(B < 0.0) B = 0.0;
  if(B > 1.0) B = 1.0;
  
  renderer->SetBackground(R, G, B);
  renderer->GetRenderWindow()->Render();
}

//----------------------------------------------------------------------
bool VtkPost::Execute(QString fileName)
{
  if(fileName.isEmpty()) {
    ecmaConsole->append("Execute: file name must be given");
    return false;
  }

  QFile scriptFile(fileName);

  if(!scriptFile.exists()) {
    ecmaConsole->append("Execute: script file does not exist");
    return false;
  }

  scriptFile.open(QIODevice::ReadOnly);

  if(scriptFile.error()) {
    ecmaConsole->append("Execute: error when opening script file");
    return false;
  }
  
  QByteArray script = scriptFile.readAll();
  scriptFile.close();

  QScriptValue val = engine->evaluate(script);

  if(engine->hasUncaughtException()) {
    int line = engine->uncaughtExceptionLineNumber();
    QString msg = "Uncaught exception at line" + QString::number(line) + ":" + val.toString();
    ecmaConsole->append(msg);
    return false;
  }

  // do not show "undefined"
  if(val.toString().trimmed() != "undefined")
    ecmaConsole->append(val.toString());

  return true;
}


