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
 *  ElmerGUI mainwindow                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter RÂback                   *
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
#include <QFile>
#include <QFont>
#include <QProgressBar>
#include <QAction>
#include <QSystemTrayIcon>
#include <QContextMenuEvent>
#include <QTimeLine>
#include <QFileInfo>
#include <QStringList>
#include <QDir>

#include <iostream>
#include <fstream>

#include <QDebug>

#include "mainwindow.h"

#ifdef EG_VTK
#include "vtkpost/vtkpost.h"
VtkPost *vtkp;
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace std;

#undef MPICH2

// Construct main window...
//-----------------------------------------------------------------------------
MainWindow::MainWindow()
{
#ifdef __APPLE__
  // find "Home directory":
  char executablePath[MAXPATHLENGTH] = {0};
  uint32_t len = MAXPATHLENGTH;
  this->homePath = "";
  if(! _NSGetExecutablePath( (char*) executablePath, &len)){
    // remove executable name from path:
    *(strrchr(executablePath,'/'))='\0';
    // remove last path component name from path:
    *(strrchr(executablePath,'/'))='\0';
    this->homePath = executablePath;
  }
#else
  homePath="";
#endif

  // load ini file:
  egIni = new EgIni(this);

  // splash screen:
  setupSplash();

  // load splash screen:
  updateSplash("Loading images...");
  
  // load tetlib:
  updateSplash("Loading tetlib...");
  tetlibAPI = new TetlibAPI;
  tetlibPresent = tetlibAPI->loadTetlib();
  this->in = tetlibAPI->in;
  this->out = tetlibAPI->out;
  
  // load nglib:
  updateSplash("Loading nglib...");
  nglibAPI = new NglibAPI;
  nglibPresent = true;

  // construct elmergrid:
  updateSplash("Constructing elmergrid...");
  elmergridAPI = new ElmergridAPI;

  // set dynamic limits:
  limit = new Limit;
  setDynamicLimits();
  
  // widgets and utilities:
  updateSplash("ElmerGUI loading...");
  glWidget = new GLWidget(this);
  setCentralWidget(glWidget);
  sifWindow = new SifWindow(this);
  meshControl = new MeshControl(this);
  boundaryDivide = new BoundaryDivide(this);
  meshingThread = new MeshingThread(this);
  meshutils = new Meshutils;
  solverLogWindow = new SifWindow(this);
  solver = new QProcess(this);
  post = new QProcess(this);
  compiler = new QProcess(this);
  meshSplitter = new QProcess(this);
  meshUnifier = new QProcess(this);
  generalSetup = new GeneralSetup(this);
  summaryEditor = new SummaryEditor(this);
  sifGenerator = new SifGenerator;
  sifGenerator->setLimit(this->limit);
  elmerDefs = new QDomDocument;
  edfEditor = new EdfEditor;
  glControl = new GLcontrol(this);
  parallel = new Parallel(this);
  checkMpi = new CheckMpi;
  materialLibrary = new MaterialLibrary(this);
  twodView = new TwodView;
  grabTimeLine = new QTimeLine(1000, this);
  
#ifdef EG_QWT
  convergenceView = new ConvergenceView(limit, this);
#endif

#ifdef EG_VTK
  vtkp = vtkPost = new VtkPost(this);
  vtkPostMeshUnifierRunning = false;
#endif

#ifdef EG_OCC
  cadView = new CadView();
  if(egIni->isPresent("deflection"))
    cadView->setDeflection(egIni->value("deflection").toDouble());
#endif

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();
      
  // Always, when an action from the menu bar has been selected, synchronize menu to state:
  connect(menuBar(), SIGNAL(triggered(QAction*)), this, SLOT(menuBarTriggeredSlot(QAction*)));
  connect(contextMenu, SIGNAL(triggered(QAction*)), this, SLOT(menuBarTriggeredSlot(QAction*)));

  // glWidget emits (list_t*) when a boundary is selected by double clicking:
  connect(glWidget, SIGNAL(signalBoundarySelected(list_t*)), this, SLOT(boundarySelectedSlot(list_t*)));

  // glWidget emits (void) when esc has been pressed:
  connect(glWidget, SIGNAL(escPressed()), this, SLOT(viewNormalModeSlot()));

  // meshingThread emits (void) when the mesh generation has finished or terminated:
  connect(meshingThread, SIGNAL(started()), this, SLOT(meshingStartedSlot()));
  connect(meshingThread, SIGNAL(finished()), this, SLOT(meshingFinishedSlot()));
  connect(meshingThread, SIGNAL(terminated()), this, SLOT(meshingTerminatedSlot()));

  // boundaryDivide emits (double) when "divide button" has been clicked:
  connect(boundaryDivide, SIGNAL(signalDoDivideSurface(double)), this, SLOT(doDivideSurfaceSlot(double)));

  // boundaryDivide emits (double) when "divide button" has been clicked:
  connect(boundaryDivide, SIGNAL(signalDoDivideEdge(double)), this, SLOT(doDivideEdgeSlot(double)));

  // solver emits (int) when finished:
  connect(solver, SIGNAL(finished(int)), this, SLOT(solverFinishedSlot(int))) ;

  // solver emits (void) when there is something to read from stdout:
  connect(solver, SIGNAL(readyReadStandardOutput()), this, SLOT(solverStdoutSlot()));

  // solver emits (void) when there is something to read from stderr:
  connect(solver, SIGNAL(readyReadStandardError()), this, SLOT(solverStderrSlot()));

  // solver emits (QProcess::ProcessError) when error occurs:
  connect(solver, SIGNAL(error(QProcess::ProcessError)), this, SLOT(solverErrorSlot(QProcess::ProcessError)));

  // solver emits (QProcess::ProcessState) when state changed:
  connect(solver, SIGNAL(stateChanged(QProcess::ProcessState)), this, SLOT(solverStateChangedSlot(QProcess::ProcessState)));

  // compiler emits (int) when finished:
  connect(compiler, SIGNAL(finished(int)), this, SLOT(compilerFinishedSlot(int))) ;

  // compiler emits (void) when there is something to read from stdout:
  connect(compiler, SIGNAL(readyReadStandardOutput()), this, SLOT(compilerStdoutSlot()));

  // compiler emits (void) when there is something to read from stderr:
  connect(compiler, SIGNAL(readyReadStandardError()), this, SLOT(compilerStderrSlot()));

  // post emits (int) when finished:
  connect(post, SIGNAL(finished(int)), this, SLOT(postProcessFinishedSlot(int))) ;

  // meshSplitter emits (int) when finished:
  connect(meshSplitter, SIGNAL(finished(int)), this, SLOT(meshSplitterFinishedSlot(int)));

  // meshSplitter emits(void) when there is something to read from stdout:
  connect(meshSplitter, SIGNAL(readyReadStandardOutput()), this, SLOT(meshSplitterStdoutSlot()));
  
  // meshSplitter emits(void) when there is something to read from stderr:
  connect(meshSplitter, SIGNAL(readyReadStandardError()), this, SLOT(meshSplitterStderrSlot()));
  
  // meshUnifier emits (int) when finished:
  connect(meshUnifier, SIGNAL(finished(int)), this, SLOT(meshUnifierFinishedSlot(int)));

  // meshUnifier emits(void) when there is something to read from stdout:
  connect(meshUnifier, SIGNAL(readyReadStandardOutput()), this, SLOT(meshUnifierStdoutSlot()));
  
  // meshUnifier emits(void) when there is something to read from stderr:
  connect(meshUnifier, SIGNAL(readyReadStandardError()), this, SLOT(meshUnifierStderrSlot()));

  // grabTimeLine emits finished() when done:
  connect(grabTimeLine, SIGNAL(finished()), this, SLOT(grabFrameSlot()));
  
  // set initial state:
  operations = 0;
  meshControl->nglibPresent = nglibPresent;
  meshControl->tetlibPresent = tetlibPresent;
  meshControl->defaultControls();
  nglibInputOk = false;
  tetlibInputOk = false;
  activeGenerator = GEN_UNKNOWN;
  bcEditActive = false;
  bodyEditActive = false;
  showConvergence = egIni->isSet("showconvergence");
  geometryInputFileName = "";
  occInputOk = false;

  // background image:
  glWidget->stateUseBgImage = egIni->isSet("bgimage");
  glWidget->stateStretchBgImage = egIni->isSet("bgimagestretch");
  glWidget->stateAlignRightBgImage = egIni->isSet("bgimagealignright");
  glWidget->bgImageFileName = egIni->value("bgimagefile");

  // set font for text editors:
  // QFont sansFont("Courier", 10);
  // sifWindow->getTextEdit()->setCurrentFont(sansFont);
  // solverLogWindow->getTextEdit()->setCurrentFont(sansFont);

  // load definition files:
  updateSplash("Loading definitions...");
  loadDefinitions();

  // initialization ready:
  synchronizeMenuToState();
  setWindowTitle(tr("ElmerGUI"));
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
  finalizeSplash();
  setupSysTrayIcon();

  // default size:
  int defW = egIni->value("width").toInt();
  int defH = egIni->value("height").toInt();
  if(defW <= 200) defW = 200;
  if(defH <= 200) defH = 200;
  this->resize(defW, defH);
}

// dtor...
//-----------------------------------------------------------------------------
MainWindow::~MainWindow()
{
  qApp->closeAllWindows();
}

// Set limits for dynamic editors, materials, bcs, etc...
//-----------------------------------------------------------------------------
void MainWindow::setDynamicLimits()
{
  // Values defined in "edf/egini.xml" that override default limits:

  // Deprecated ** 23/04/09 **
  if(egIni->isPresent("max_boundaries")) {
    limit->setMaxBoundaries(egIni->value("max_boundaries").toInt());
    // cout << "Max boundaries: " << limit->maxBoundaries() << endl;
  }

  // Deprecated ** 23/04/09 **
  if(egIni->isPresent("max_solvers")) {
    limit->setMaxSolvers(egIni->value("max_solvers").toInt());
    // cout << "Max solvers: " << limit->maxSolvers() << endl;
  }

  // Deprecated ** 23/04/09 **
  if(egIni->isPresent("max_bodies")) {
    limit->setMaxBodies(egIni->value("max_bodies").toInt());
    // cout << "Max bodies: " << limit->maxBodies() << endl;
  }

  // Deprecated ** 21/04/09 **
  if(egIni->isPresent("max_equations")) {
    limit->setMaxEquations(egIni->value("max_equations").toInt());
    // cout << "Max equations: " << limit->maxEquations() << endl;
  }

  // Deprecated ** 21/04/09 **
  if(egIni->isPresent("max_materials")) {
    limit->setMaxMaterials(egIni->value("max_materials").toInt());
    // cout << "Max materials: " << limit->maxMaterials() << endl;
  }

  // Deprecated ** 21/04/09 **
  if(egIni->isPresent("max_bodyforces")) {
    limit->setMaxBodyforces(egIni->value("max_bodyforces").toInt());
    // cout << "Max bodyforces: " << limit->maxBodyforces() << endl;
  }

  // Deprecated ** 21/04/09 **
  if(egIni->isPresent("max_initialconditions")) {
    limit->setMaxInitialconditions(egIni->value("max_initialconditions").toInt());
    // cout << "Max initialconditions: " << limit->maxInitialconditions() << endl;
  }

  // Deprecated ** 21/04/09 **  
  if(egIni->isPresent("max_bcs")) {
    limit->setMaxBcs(egIni->value("max_bcs").toInt());
    // cout << "Max bcs: " << limit->maxBcs() << endl;
  }
}


// Always synchronize menu to state when the menubar has been triggered...
//-----------------------------------------------------------------------------
void MainWindow::menuBarTriggeredSlot(QAction *act)
{
  synchronizeMenuToState();
}


// Create actions...
//-----------------------------------------------------------------------------
void MainWindow::createActions()
{
  // File -> Open file
  openAct = new QAction(QIcon(":/icons/document-open.png"), tr("&Open..."), this);
  openAct->setShortcut(tr("Ctrl+O"));
  openAct->setStatusTip(tr("Open geometry input file"));
  connect(openAct, SIGNAL(triggered()), this, SLOT(openSlot()));
  
  // File -> Load mesh...
  loadAct = new QAction(QIcon(":/icons/document-open-folder.png"), tr("&Load mesh..."), this);
  loadAct->setStatusTip(tr("Load Elmer mesh files"));
  connect(loadAct, SIGNAL(triggered()), this, SLOT(loadSlot()));
  
  // File -> Load project...
  loadProjectAct = new QAction(QIcon(":/icons/document-import.png"), tr("Load &project..."), this);
  loadProjectAct->setStatusTip(tr("Load previously saved project"));
  connect(loadProjectAct, SIGNAL(triggered()), this, SLOT(loadProjectSlot()));

  // File -> Definitions...
  editDefinitionsAct = new QAction(QIcon(":/icons/games-config-custom.png"), tr("&Definitions..."), this);
  editDefinitionsAct->setStatusTip(tr("Load and edit Elmer sif definitions file"));
  connect(editDefinitionsAct, SIGNAL(triggered()), this, SLOT(editDefinitionsSlot()));

  // File -> Save...
  saveAct = new QAction(QIcon(":/icons/document-save.png"), tr("&Save..."), this);
  saveAct->setShortcut(tr("Ctrl+S"));
  saveAct->setStatusTip(tr("Save Elmer mesh and sif-files"));
  connect(saveAct, SIGNAL(triggered()), this, SLOT(saveSlot()));

  // File -> Save as...
  saveAsAct = new QAction(QIcon(":/icons/document-save-as.png"), tr("&Save as..."), this);
  saveAsAct->setStatusTip(tr("Save Elmer mesh and sif-files"));
  connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAsSlot()));

  // File -> Save project...
  saveProjectAct = new QAction(QIcon(":/icons/document-export.png"), tr("&Save project..."), this);
  saveProjectAct->setStatusTip(tr("Save current project"));
  connect(saveProjectAct, SIGNAL(triggered()), this, SLOT(saveProjectSlot()));

  // File -> Save picture as...
  savePictureAct = new QAction(QIcon(":/icons/view-preview.png"), tr("&Save picture as..."), this);
  savePictureAct->setStatusTip(tr("Save picture in file"));
  connect(savePictureAct, SIGNAL(triggered()), this, SLOT(savePictureSlot()));

  // File -> Exit
  exitAct = new QAction(QIcon(":/icons/application-exit.png"), tr("E&xit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  exitAct->setStatusTip(tr("Exit"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(closeMainWindowSlot()));

  // Model -> Setup...
  modelSetupAct = new QAction(QIcon(), tr("Setup..."), this);
  modelSetupAct->setStatusTip(tr("Setup simulation environment"));
  connect(modelSetupAct, SIGNAL(triggered()), this, SLOT(modelSetupSlot()));

  // Model -> Equation...
  addEquationAct = new QAction(QIcon(), tr("Add..."), this);
  addEquationAct->setStatusTip(tr("Add a PDE-system to the equation list"));
  connect(addEquationAct, SIGNAL(triggered()), this, SLOT(addEquationSlot()));

  // Model -> Material...
  addMaterialAct = new QAction(QIcon(), tr("Add..."), this);
  addMaterialAct->setStatusTip(tr("Add a material set to the material list"));
  connect(addMaterialAct, SIGNAL(triggered()), this, SLOT(addMaterialSlot()));

  // Model -> Body force...
  addBodyForceAct = new QAction(QIcon(), tr("Add..."), this);
  addBodyForceAct->setStatusTip(tr("Add body forces..."));
  connect(addBodyForceAct, SIGNAL(triggered()), this, SLOT(addBodyForceSlot()));

  // Model -> Initial condition...
  addInitialConditionAct = new QAction(QIcon(), tr("Add..."), this);
  addInitialConditionAct->setStatusTip(tr("Add initial conditions..."));
  connect(addInitialConditionAct, SIGNAL(triggered()), this, SLOT(addInitialConditionSlot()));

  // Model -> Boundary condition...
  addBoundaryConditionAct = new QAction(QIcon(), tr("Add..."), this);
  addBoundaryConditionAct->setStatusTip(tr("Add boundary conditions..."));
  connect(addBoundaryConditionAct, SIGNAL(triggered()), this, SLOT(addBoundaryConditionSlot()));

  // Model -> Set body properties
  bodyEditAct = new QAction(QIcon(), tr("Set body properties"), this);
  bodyEditAct->setStatusTip(tr("Set body properties (equivalent to holding down the SHIFT key)"));
  connect(bodyEditAct, SIGNAL(triggered()), this, SLOT(bodyEditSlot()));
  bodyEditAct->setCheckable(true);

  // Model -> Set boundary conditions
  bcEditAct = new QAction(QIcon(), tr("Set boundary properties"), this);
  bcEditAct->setStatusTip(tr("Set boundary properties (equivalent to holding down the ALT key)"));
  connect(bcEditAct, SIGNAL(triggered()), this, SLOT(bcEditSlot()));
  bcEditAct->setCheckable(true);

  // Model -> Summary...
  modelSummaryAct = new QAction(QIcon(), tr("Summary..."), this);
  modelSummaryAct->setStatusTip(tr("Model summary"));
  connect(modelSummaryAct, SIGNAL(triggered()), this, SLOT(modelSummarySlot()));

  // Model -> Clear
  modelClearAct = new QAction(QIcon(), tr("Clear all"), this);
  modelClearAct->setStatusTip(tr("Clear all model definitions"));
  connect(modelClearAct, SIGNAL(triggered()), this, SLOT(modelClearSlot()));

  // Edit -> Generate sif
  generateSifAct = new QAction(QIcon(""), tr("&Generate"), this);
  generateSifAct->setShortcut(tr("Ctrl+G"));
  generateSifAct->setStatusTip(tr("Genarete solver input file"));
  connect(generateSifAct, SIGNAL(triggered()), this, SLOT(generateSifSlot()));

  // Edit -> Solver input file...
  showsifAct = new QAction(QIcon(":/icons/document-properties.png"), tr("&Edit..."), this);
  showsifAct->setShortcut(tr("Ctrl+S"));
  showsifAct->setStatusTip(tr("Edit solver input file"));
  connect(showsifAct, SIGNAL(triggered()), this, SLOT(showsifSlot()));

  // Mesh -> Control
  meshcontrolAct = new QAction(QIcon(":/icons/configure.png"), tr("&Configure..."), this);
  meshcontrolAct->setShortcut(tr("Ctrl+C"));
  meshcontrolAct->setStatusTip(tr("Configure mesh generators"));
  connect(meshcontrolAct, SIGNAL(triggered()), this, SLOT(meshcontrolSlot()));

  // Mesh -> Remesh
  remeshAct = new QAction(QIcon(":/icons/edit-redo.png"), tr("&Remesh"), this);
  remeshAct->setShortcut(tr("Ctrl+R"));
  remeshAct->setStatusTip(tr("Remesh"));
  connect(remeshAct, SIGNAL(triggered()), this, SLOT(remeshSlot()));

  // Mesh -> Kill generator
  stopMeshingAct = new QAction(QIcon(":/icons/window-close.png"), 
			       tr("&Terminate meshing"), this);
  stopMeshingAct->setStatusTip(tr("Terminate mesh generator"));
  connect(stopMeshingAct, SIGNAL(triggered()), this, SLOT(stopMeshingSlot()));
  stopMeshingAct->setEnabled(false);

  // Mesh -> Divide surface
  surfaceDivideAct = new QAction(QIcon(":/icons/divide.png"), 
				 tr("&Divide surface..."), this);
  surfaceDivideAct->setStatusTip(tr("Divide surface by sharp edges"));
  connect(surfaceDivideAct, SIGNAL(triggered()), this, SLOT(surfaceDivideSlot()));

  // Mesh -> Unify surface
  surfaceUnifyAct = new QAction(QIcon(":/icons/unify.png"), tr("&Unify surface"), this);
  surfaceUnifyAct->setStatusTip(tr("Unify surface (merge selected)"));
  connect(surfaceUnifyAct, SIGNAL(triggered()), this, SLOT(surfaceUnifySlot()));

  // Mesh -> Divide edge
  edgeDivideAct = new QAction(QIcon(":/icons/divide-edge.png"), tr("&Divide edge..."), this);
  edgeDivideAct->setStatusTip(tr("Divide edge by sharp points"));
  connect(edgeDivideAct, SIGNAL(triggered()), this, SLOT(edgeDivideSlot()));

  // Mesh -> Unify edges
  edgeUnifyAct = new QAction(QIcon(":/icons/unify-edge.png"), tr("&Unify edge"), this);
  edgeUnifyAct->setStatusTip(tr("Unify edge (merge selected)"));
  connect(edgeUnifyAct, SIGNAL(triggered()), this, SLOT(edgeUnifySlot()));

  // Mesh -> Clean up
  cleanHangingSharpEdgesAct = new QAction(QIcon(""), tr("Clean up"), this);
  cleanHangingSharpEdgesAct->setStatusTip(tr("Removes hanging/orphan sharp edges (for visualization)"));
  connect(cleanHangingSharpEdgesAct, SIGNAL(triggered()), this, SLOT(cleanHangingSharpEdgesSlot()));

  // View -> Full screen
  viewFullScreenAct = new QAction(QIcon(), tr("Full screen"), this);
  viewFullScreenAct->setShortcut(tr("Ctrl+L"));
  viewFullScreenAct->setStatusTip(tr("Full screen mode"));
  connect(viewFullScreenAct, SIGNAL(triggered()), this, SLOT(viewFullScreenSlot()));
  viewFullScreenAct->setCheckable(true);

  // View -> Show surface mesh
  hidesurfacemeshAct = new QAction(QIcon(), tr("Surface mesh"), this);
  hidesurfacemeshAct->setStatusTip(tr("Show/hide surface mesh "
				      "(do/do not outline surface elements)"));
  connect(hidesurfacemeshAct, SIGNAL(triggered()), this, SLOT(hidesurfacemeshSlot()));
  hidesurfacemeshAct->setCheckable(true);

  // View -> Show volume mesh
  hidevolumemeshAct = new QAction(QIcon(), tr("Volume mesh"), this);
  hidevolumemeshAct->setStatusTip(tr("Show/hide volume mesh "
				      "(do/do not outline volume mesh edges)"));
  connect(hidevolumemeshAct, SIGNAL(triggered()), this, SLOT(hidevolumemeshSlot()));
  hidevolumemeshAct->setCheckable(true);

  // View -> Show sharp edges
  hidesharpedgesAct = new QAction(QIcon(), tr("Sharp edges"), this);
  hidesharpedgesAct->setStatusTip(tr("Show/hide sharp edges"));
  connect(hidesharpedgesAct, SIGNAL(triggered()), this, SLOT(hidesharpedgesSlot()));
  hidesharpedgesAct->setCheckable(true);

  // View -> Compass
  viewCoordinatesAct = new QAction(QIcon(), tr("Compass"), this);
  viewCoordinatesAct->setStatusTip(tr("View coordinates "
				      "(RGB=XYZ modulo translation)"));
  connect(viewCoordinatesAct, SIGNAL(triggered()), this, SLOT(viewCoordinatesSlot()));
  viewCoordinatesAct->setCheckable(true);

  // View -> Select all surfaces
  selectAllSurfacesAct = new QAction(QIcon(), tr("Select all surfaces"), this);
  selectAllSurfacesAct->setStatusTip(tr("Select all surfaces"));
  connect(selectAllSurfacesAct, SIGNAL(triggered()), this, SLOT(selectAllSurfacesSlot()));

  // View -> Select all edges
  selectAllEdgesAct = new QAction(QIcon(), tr("Select all edges"), this);
  selectAllEdgesAct->setStatusTip(tr("Select all edges"));
  connect(selectAllEdgesAct, SIGNAL(triggered()), this, SLOT(selectAllEdgesSlot()));

  // View -> Select defined edges
  selectDefinedEdgesAct = new QAction(QIcon(), tr("Select defined edges"), this);
  selectDefinedEdgesAct->setStatusTip(tr("Select defined edges"));
  connect(selectDefinedEdgesAct, SIGNAL(triggered()), this, SLOT(selectDefinedEdgesSlot()));

  // View -> Select defined surfaces
  selectDefinedSurfacesAct = new QAction(QIcon(), tr("Select defined surfaces"), this);
  selectDefinedSurfacesAct->setStatusTip(tr("Select defined surfaces"));
  connect(selectDefinedSurfacesAct, SIGNAL(triggered()), this, SLOT(selectDefinedSurfacesSlot()));

  // View -> Hide/show selected
  hideselectedAct = new QAction(QIcon(), tr("&Hide/show selected"), this);
  hideselectedAct->setShortcut(tr("Ctrl+H"));
  hideselectedAct->setStatusTip(tr("Show/hide selected objects"));
  connect(hideselectedAct, SIGNAL(triggered()), this, SLOT(hideselectedSlot()));

  // View -> Show surface numbers
  showSurfaceNumbersAct = new QAction(QIcon(), tr("Surface element numbers"), this);
  showSurfaceNumbersAct->setStatusTip(tr("Show surface element numbers "
				      "(Show the surface element numbering)"));
  connect(showSurfaceNumbersAct, SIGNAL(triggered()), this, SLOT(showSurfaceNumbersSlot()));
  showSurfaceNumbersAct->setCheckable(true);

  // View -> Show edge numbers
  showEdgeNumbersAct = new QAction(QIcon(), tr("Edge element numbers"), this);
  showEdgeNumbersAct->setStatusTip(tr("Show edge element numbers "
				      "(Show the node element numbering)"));
  connect(showEdgeNumbersAct, SIGNAL(triggered()), this, SLOT(showEdgeNumbersSlot()));
  showEdgeNumbersAct->setCheckable(true);

  // View -> Show node numbers
  showNodeNumbersAct = new QAction(QIcon(), tr("Node numbers"), this);
  showNodeNumbersAct->setStatusTip(tr("Show node numbers "
				      "(Show the node numbers)"));
  connect(showNodeNumbersAct, SIGNAL(triggered()), this, SLOT(showNodeNumbersSlot()));
  showNodeNumbersAct->setCheckable(true);
  
  // View -> Show boundray index
  showBoundaryIndexAct = new QAction(QIcon(), tr("Boundary index"), this);
  showBoundaryIndexAct->setStatusTip(tr("Show boundary index"));
  connect(showBoundaryIndexAct, SIGNAL(triggered()), this, SLOT(showBoundaryIndexSlot()));
  showBoundaryIndexAct->setCheckable(true);
  
  // View -> Show body index
  showBodyIndexAct = new QAction(QIcon(), tr("Body index"), this);
  showBodyIndexAct->setStatusTip(tr("Show body index"));
  connect(showBodyIndexAct, SIGNAL(triggered()), this, SLOT(showBodyIndexSlot()));
  showBodyIndexAct->setCheckable(true);

  // View -> Colors -> GL controls
  glControlAct = new QAction(QIcon(), tr("GL controls..."), this);
  glControlAct->setStatusTip(tr("Control GL parameters for lights and materials"));
  connect(glControlAct, SIGNAL(triggered()), this, SLOT(glControlSlot()));
  
  // View -> Colors -> Background
  chooseBGColorAct = new QAction(QIcon(), tr("Background..."), this);
  chooseBGColorAct->setStatusTip(tr("Set background color"));
  connect(chooseBGColorAct, SIGNAL(triggered()), this, SLOT(backgroundColorSlot()));
  
  // View -> Colors -> Surface elements
  chooseSurfaceColorAct = new QAction(QIcon(), tr("Surface elements..."), this);
  chooseSurfaceColorAct->setStatusTip(tr("Set surface color"));
  connect(chooseSurfaceColorAct, SIGNAL(triggered()), this, SLOT(surfaceColorSlot()));
  
  // View -> Colors -> Edge elements
  chooseEdgeColorAct = new QAction(QIcon(), tr("Edge elements..."), this);
  chooseEdgeColorAct->setStatusTip(tr("Set edge color"));
  connect(chooseEdgeColorAct, SIGNAL(triggered()), this, SLOT(edgeColorSlot()));
  
  // View -> Colors -> Surface mesh
  chooseSurfaceMeshColorAct = new QAction(QIcon(), tr("Surface mesh..."), this);
  chooseSurfaceMeshColorAct->setStatusTip(tr("Set surface mesh color"));
  connect(chooseSurfaceMeshColorAct, SIGNAL(triggered()), this, SLOT(surfaceMeshColorSlot()));
  
  // View -> Colors -> Sharp edges
  chooseSharpEdgeColorAct = new QAction(QIcon(), tr("Sharp edges..."), this);
  chooseSharpEdgeColorAct->setStatusTip(tr("Set sharp edge color"));
  connect(chooseSharpEdgeColorAct, SIGNAL(triggered()), this, SLOT(sharpEdgeColorSlot()));
  
  // View -> Colors -> Boundaries
  showBoundaryColorAct = new QAction(QIcon(), tr("Boundaries"), this);
  showBoundaryColorAct->setStatusTip(tr("Visualize different boundary parts with color patches"));
  connect(showBoundaryColorAct, SIGNAL(triggered()), this, SLOT(colorizeBoundarySlot()));
  showBoundaryColorAct->setCheckable(true);
  
  // View -> Colors -> Bodies
  showBodyColorAct = new QAction(QIcon(), tr("Bodies"), this);
  showBodyColorAct->setStatusTip(tr("Visualize different body with color patches"));
  connect(showBodyColorAct, SIGNAL(triggered()), this, SLOT(colorizeBodySlot()));
  showBodyColorAct->setCheckable(true);
  
  // View -> Shade model -> Smooth
  smoothShadeAct = new QAction(QIcon(), tr("Smooth"), this);
  smoothShadeAct->setStatusTip(tr("Set shade model to smooth"));
  connect(smoothShadeAct, SIGNAL(triggered()), this, SLOT(smoothShadeSlot()));
  smoothShadeAct->setCheckable(true);

  // View -> Shade model -> Flat
  flatShadeAct = new QAction(QIcon(), tr("Flat"), this);
  flatShadeAct->setStatusTip(tr("Set shade model to flat"));
  connect(flatShadeAct, SIGNAL(triggered()), this, SLOT(flatShadeSlot()));
  flatShadeAct->setCheckable(true);

  // View -> Projection -> Orthogonal
  orthoAct = new QAction(QIcon(), tr("Orthogonal"), this);
  orthoAct->setStatusTip(tr("Set projection to orthogonal"));
  connect(orthoAct, SIGNAL(triggered()), this, SLOT(orthoSlot()));
  orthoAct->setCheckable(true);

  // View -> Projection -> Perspective
  perspectiveAct = new QAction(QIcon(), tr("Perspective"), this);
  perspectiveAct->setStatusTip(tr("Set projection to perspective"));
  connect(perspectiveAct, SIGNAL(triggered()), this, SLOT(perspectiveSlot()));
  perspectiveAct->setCheckable(true);

  // View -> Show all
  showallAct = new QAction(QIcon(), tr("Show all"), this);
  showallAct->setStatusTip(tr("Show all objects"));
  connect(showallAct, SIGNAL(triggered()), this, SLOT(showallSlot()));

  // View -> Reset model view
  resetAct = new QAction(QIcon(), tr("Reset model view"), this);
  resetAct->setStatusTip(tr("Reset model view"));
  connect(resetAct, SIGNAL(triggered()), this, SLOT(resetSlot()));

  // View -> Show cad model
  showCadModelAct = new QAction(QIcon(), tr("Cad model..."), this);
  showCadModelAct->setStatusTip(tr("Displays the cad model in a separate window"));
  connect(showCadModelAct, SIGNAL(triggered()), this, SLOT(showCadModelSlot()));

  // View -> Show 2d view
  showTwodViewAct = new QAction(QIcon(), tr("2D modeler..."), this);
  showTwodViewAct->setStatusTip(tr("Displays the 2d geometry in a separate window"));
  connect(showTwodViewAct, SIGNAL(triggered()), this, SLOT(showTwodViewSlot()));

  // Solver -> Parallel settings
  parallelSettingsAct = new QAction(QIcon(), tr("Parallel settings..."), this);
  parallelSettingsAct->setStatusTip(tr("Choose parameters and methods for parallel solution"));
  connect(parallelSettingsAct, SIGNAL(triggered()), this, SLOT(parallelSettingsSlot()));

  // Solver -> Run solver
  runsolverAct = new QAction(QIcon(":/icons/Solver.png"), tr("Start solver"), this);
  runsolverAct->setStatusTip(tr("Run ElmerSolver"));
  connect(runsolverAct, SIGNAL(triggered()), this, SLOT(runsolverSlot()));

  // Solver -> Kill solver
  killsolverAct = new QAction(QIcon(":/icons/window-close.png"), tr("Kill solver"), this);
  killsolverAct->setStatusTip(tr("Kill ElmerSolver"));
  connect(killsolverAct, SIGNAL(triggered()), this, SLOT(killsolverSlot()));
  killsolverAct->setEnabled(false);

  // Solver -> Show convergence
  showConvergenceAct = new QAction(QIcon(), tr("Show convergence"), this);
  showConvergenceAct->setStatusTip(tr("Show/hide convergence plot"));
  connect(showConvergenceAct, SIGNAL(triggered()), this, SLOT(showConvergenceSlot()));
  showConvergenceAct->setCheckable(true);

  // Solver -> Post process
  resultsAct = new QAction(QIcon(":/icons/Post.png"), tr("Start ElmerPost"), this);
  resultsAct->setStatusTip(tr("Run ElmerPost for visualization"));
  connect(resultsAct, SIGNAL(triggered()), this, SLOT(resultsSlot()));

  // Solver -> Kill post process
  killresultsAct = new QAction(QIcon(":/icons/window-close.png"), tr("Kill ElmerPost"), this);
  killresultsAct->setStatusTip(tr("Kill ElmerPost"));
  connect(killresultsAct, SIGNAL(triggered()), this, SLOT(killresultsSlot()));
  killresultsAct->setEnabled(false);

  // Solver -> Show Vtk postprocessor
  showVtkPostAct = new QAction(QIcon(), tr("Start ElmerVTK"), this);
  showVtkPostAct->setStatusTip(tr("Invokes VTK based ElmerGUI postprocessor"));
  connect(showVtkPostAct, SIGNAL(triggered()), this, SLOT(showVtkPostSlot()));

  // Solver -> Show ParaView postprocessor
  paraviewAct = new QAction(QIcon(), tr("Start ParaView"), this);
  paraviewAct->setStatusTip(tr("Invokes ParaView for visualization"));
  connect(paraviewAct, SIGNAL(triggered()), this, SLOT(showParaViewSlot()));

  // Solver -> Compiler...
  compileSolverAct = new QAction(QIcon(""), tr("Compiler..."), this);
  compileSolverAct->setStatusTip(tr("Compile Elmer specific source code (f90) into a shared library (dll)"));
  connect(compileSolverAct, SIGNAL(triggered()), 
	  this, SLOT(compileSolverSlot()));

  // Help -> About
  aboutAct = new QAction(QIcon(":/icons/help-about.png"), tr("About..."), this);
  aboutAct->setStatusTip(tr("Information about the program"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(showaboutSlot()));

#if WIN32
#else
  compileSolverAct->setEnabled(false);
#endif

  if(egIni->isSet("bgimage"))
    chooseBGColorAct->setEnabled(false);
}


// Create menus...
//-----------------------------------------------------------------------------
void MainWindow::createMenus()
{
  // File menu
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAct);
  fileMenu->addAction(loadAct);
  fileMenu->addAction(loadProjectAct);
  fileMenu->addSeparator();
  fileMenu->addAction(editDefinitionsAct);
  fileMenu->addSeparator();
  fileMenu->addAction(saveAct);
  fileMenu->addAction(saveAsAct);
  fileMenu->addAction(saveProjectAct);
  fileMenu->addSeparator();
  fileMenu->addAction(savePictureAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  // Mesh menu
  meshMenu = menuBar()->addMenu(tr("&Mesh"));
  meshMenu->addAction(meshcontrolAct);
  meshMenu->addAction(remeshAct);
  meshMenu->addAction(stopMeshingAct);
  meshMenu->addSeparator();
  meshMenu->addAction(surfaceDivideAct);
  meshMenu->addAction(surfaceUnifyAct);
  meshMenu->addSeparator();
  meshMenu->addAction(edgeDivideAct);
  meshMenu->addAction(edgeUnifyAct);
  meshMenu->addSeparator();
  meshMenu->addAction(cleanHangingSharpEdgesAct);

  // Model menu
  modelMenu = menuBar()->addMenu(tr("&Model"));

  modelMenu->addAction(modelSetupAct);
  modelMenu->addSeparator();

  equationMenu = modelMenu->addMenu(tr("Equation"));
  equationMenu->addAction(addEquationAct);
  equationMenu->addSeparator();
  connect(equationMenu, SIGNAL(triggered(QAction*)), 
	  this, SLOT(equationSelectedSlot(QAction*)));

  modelMenu->addSeparator();
  materialMenu = modelMenu->addMenu(tr("Material"));
  materialMenu->addAction(addMaterialAct);
  materialMenu->addSeparator();
  connect(materialMenu, SIGNAL(triggered(QAction*)), 
	  this, SLOT(materialSelectedSlot(QAction*)));

  modelMenu->addSeparator();
  bodyForceMenu = modelMenu->addMenu(tr("Body force"));
  bodyForceMenu->addAction(addBodyForceAct);
  bodyForceMenu->addSeparator();
  connect(bodyForceMenu, SIGNAL(triggered(QAction*)), 
	  this, SLOT(bodyForceSelectedSlot(QAction*)));
  
  modelMenu->addSeparator();
  initialConditionMenu = modelMenu->addMenu(tr("Initial condition"));
  initialConditionMenu->addAction(addInitialConditionAct);
  initialConditionMenu->addSeparator();
  connect(initialConditionMenu, SIGNAL(triggered(QAction*)), 
	  this, SLOT(initialConditionSelectedSlot(QAction*)));
  
  modelMenu->addSeparator();
  boundaryConditionMenu = modelMenu->addMenu(tr("Boundary condition"));
  boundaryConditionMenu->addAction(addBoundaryConditionAct);
  boundaryConditionMenu->addSeparator();
  connect(boundaryConditionMenu, SIGNAL(triggered(QAction*)), 
	  this, SLOT(boundaryConditionSelectedSlot(QAction*)));
  
  modelMenu->addSeparator();
  modelMenu->addAction(bodyEditAct);
  modelMenu->addAction(bcEditAct);
  modelMenu->addSeparator();
  modelMenu->addAction(modelSummaryAct);
  modelMenu->addSeparator();
  modelMenu->addAction(modelClearAct);
  modelMenu->addSeparator();

  // View menu
  viewMenu = menuBar()->addMenu(tr("&View"));
  viewMenu->addAction(viewFullScreenAct);
  viewMenu->addSeparator();
  viewMenu->addAction(hidesurfacemeshAct);
  viewMenu->addAction(hidevolumemeshAct);
  viewMenu->addAction(hidesharpedgesAct);
  viewMenu->addAction(viewCoordinatesAct);
  viewMenu->addSeparator();
  viewMenu->addAction(selectAllSurfacesAct);
  viewMenu->addAction(selectAllEdgesAct);
  // Momentarily disabled (see comment *** TODO *** below):
  // viewMenu->addSeparator();
  // viewMenu->addAction(selectDefinedEdgesAct);
  // viewMenu->addAction(selectDefinedSurfacesAct);
  viewMenu->addSeparator();
  viewMenu->addAction(hideselectedAct);
  viewMenu->addSeparator();
  shadeMenu = viewMenu->addMenu(tr("Shade model"));
  shadeMenu->addAction(flatShadeAct);
  shadeMenu->addAction(smoothShadeAct);
  viewMenu->addSeparator();
  projectionMenu = viewMenu->addMenu(tr("Projection"));
  projectionMenu->addAction(orthoAct);
  projectionMenu->addAction(perspectiveAct);
  viewMenu->addSeparator();
  numberingMenu = viewMenu->addMenu(tr("Numbering"));
  numberingMenu->addAction(showSurfaceNumbersAct);
  numberingMenu->addAction(showEdgeNumbersAct);
  numberingMenu->addAction(showNodeNumbersAct);
  numberingMenu->addSeparator();
  numberingMenu->addAction(showBoundaryIndexAct);
  numberingMenu->addAction(showBodyIndexAct);
  viewMenu->addSeparator();
  colorizeMenu = viewMenu->addMenu(tr("Lights and colors"));
  colorizeMenu->addAction(glControlAct);
  colorizeMenu->addSeparator();
  colorizeMenu->addAction(chooseBGColorAct);
  colorizeMenu->addSeparator();
  colorizeMenu->addAction(chooseSurfaceColorAct);
  colorizeMenu->addAction(chooseEdgeColorAct);
  colorizeMenu->addSeparator();
  colorizeMenu->addAction(chooseSurfaceMeshColorAct);
  colorizeMenu->addAction(chooseSharpEdgeColorAct);
  colorizeMenu->addSeparator();
  colorizeMenu->addAction(showBoundaryColorAct);
  colorizeMenu->addAction(showBodyColorAct);
  viewMenu->addSeparator();
  viewMenu->addAction(showallAct);
  viewMenu->addAction(resetAct);
#ifdef EG_OCC
  viewMenu->addSeparator();
  viewMenu->addAction(showCadModelAct);
#endif
  viewMenu->addAction(showTwodViewAct);

  // Edit menu
  editMenu = menuBar()->addMenu(tr("&Sif"));
  editMenu->addAction(generateSifAct);
  editMenu->addSeparator();
  editMenu->addAction(showsifAct);

  //  SolverMenu
  solverMenu = menuBar()->addMenu(tr("&Run"));
  solverMenu->addAction(parallelSettingsAct);
  solverMenu->addSeparator();
  solverMenu->addAction(runsolverAct);
  solverMenu->addAction(killsolverAct);
#ifdef EG_QWT
  solverMenu->addAction(showConvergenceAct);
#endif
  solverMenu->addSeparator();
  solverMenu->addAction(resultsAct);
  solverMenu->addAction(killresultsAct);
#ifdef EG_VTK
  solverMenu->addSeparator();
  solverMenu->addAction(showVtkPostAct);
#endif
#ifdef EG_PARAVIEW
  solverMenu->addSeparator();
  solverMenu->addAction(paraviewAct);
#endif
  solverMenu->addSeparator();
  solverMenu->addAction(compileSolverAct);

  // Help menu
  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutAct);

  // Sys tray menu:
  sysTrayMenu = new QMenu;
  sysTrayMenu->addAction(modelSummaryAct);
  sysTrayMenu->addSeparator();
  sysTrayMenu->addAction(stopMeshingAct);
  sysTrayMenu->addSeparator();
  sysTrayMenu->addAction(killsolverAct);
  sysTrayMenu->addAction(killresultsAct);
  sysTrayMenu->addSeparator();
  sysTrayMenu->addAction(aboutAct);
  sysTrayMenu->addSeparator();
  sysTrayMenu->addAction(exitAct);

  // Context menu:
  contextMenu = new QMenu;
  contextMenu->addMenu(fileMenu);
  contextMenu->addMenu(meshMenu);
  contextMenu->addMenu(modelMenu);
  contextMenu->addMenu(viewMenu);
  contextMenu->addMenu(editMenu);
  contextMenu->addMenu(solverMenu);
  contextMenu->addMenu(helpMenu);

  // Disable unavailable external components:
  //------------------------------------------
  if(!egIni->isSet("checkexternalcomponents"))
     return;

  QProcess testProcess;
  QStringList args;

  cout << "Checking for ElmerSolver... ";
  updateSplash("Checking for ElmerSolver...");
  args << "-v";
  testProcess.start("ElmerSolver", args);
  if(!testProcess.waitForStarted()) {
    logMessage("no - disabling solver features");
    runsolverAct->setEnabled(false);
    showConvergenceAct->setEnabled(false);
    killsolverAct->setEnabled(false);
  } else {
    cout << "yes" << endl;
  }
  testProcess.waitForFinished(2000);

  cout << "Checking for ... ";
  updateSplash("Checking for ElmerPost...");
  args << "-v";
  testProcess.start("ElmerPost", args);
  if(!testProcess.waitForStarted()) {
    logMessage("no - disabling ElmerPost postprocessing features");
    resultsAct->setEnabled(false);
    killresultsAct->setEnabled(false);
  } else {
    cout << "yes" << endl;
  }
  testProcess.waitForFinished(2000);

  cout << "Checking for ElmerGrid... ";
  updateSplash("Checking for ElmerGrid...");
  testProcess.start("ElmerGrid");
  if(!testProcess.waitForStarted()) {
    logMessage("no - disabling parallel features");
    parallelSettingsAct->setEnabled(false);
  } else {
    cout << "yes" << endl;
  }
  testProcess.waitForFinished(2000);

  cout << "Checking for ElmerSolver_mpi... ";
  updateSplash("Checking for ElmerSolver_mpi...");
  args << "-v";
  testProcess.start("ElmerSolver_mpi", args);
  if(!testProcess.waitForStarted()) {
    logMessage("no - disabling parallel features");
    parallelSettingsAct->setEnabled(false);
  } else {
    cout << "yes" << endl;
  }
  testProcess.waitForFinished(2000);

}


// Create tool bars...
//-----------------------------------------------------------------------------
void MainWindow::createToolBars()
{
  // File toolbar
  fileToolBar = addToolBar(tr("&File"));
  fileToolBar->addAction(openAct);
  fileToolBar->addAction(loadAct);
  fileToolBar->addAction(loadProjectAct);
  fileToolBar->addSeparator();
  fileToolBar->addAction(saveAct);
  fileToolBar->addAction(saveAsAct);
  fileToolBar->addAction(saveProjectAct);
  fileToolBar->addSeparator();
  fileToolBar->addAction(savePictureAct);

  // Edit toolbar
  editToolBar = addToolBar(tr("&Edit"));
  editToolBar->addAction(showsifAct);

  // Mesh toolbar
  meshToolBar = addToolBar(tr("&Mesh"));
  meshToolBar->addAction(meshcontrolAct);
  meshToolBar->addAction(remeshAct);
  meshToolBar->addSeparator();
  meshToolBar->addAction(surfaceDivideAct);
  meshToolBar->addAction(surfaceUnifyAct);
  meshToolBar->addSeparator();
  meshToolBar->addAction(edgeDivideAct);
  meshToolBar->addAction(edgeUnifyAct);

  // Solver toolbar
  solverToolBar = addToolBar(tr("&Solver"));
  solverToolBar->addAction(runsolverAct);
  solverToolBar->addAction(resultsAct);

  if(egIni->isSet("hidetoolbars")) {
    fileToolBar->hide();
    editToolBar->hide();
    meshToolBar->hide();
    solverToolBar->hide();
  }
}


// Create status bar...
//-----------------------------------------------------------------------------
void MainWindow::createStatusBar()
{
  progressBar = new QProgressBar;
  progressBar->setMaximumHeight(12);
  progressBar->setMaximumWidth(120);
  progressBar->setTextVisible(false);
  progressBar->hide();

  progressLabel = new QLabel;
  progressLabel->hide();

  statusBar()->addPermanentWidget(progressLabel);
  statusBar()->addPermanentWidget(progressBar);

  statusBar()->showMessage(tr("Ready"));

  connect(grabTimeLine, SIGNAL(frameChanged(int)), progressBar, SLOT(setValue(int)));
}


//*****************************************************************************
//
//                                File MENU
//
//*****************************************************************************

void MainWindow::parseCmdLine()
{
  QStringList args = QCoreApplication::arguments();
  
  if(!args.contains("-nogui"))
    this->show();

  int input = args.indexOf("-i");

  if(input > 0) {
    QString fileName = args.at(input + 1);

    QFileInfo fileInfo(fileName);
    
    if(!fileInfo.exists()) {
#if WITH_QT5
      cout << "Input file \"" << fileName.toLatin1().data() << "\" does not exist" << endl;
#else
      cout << "Input file \"" << fileName.toAscii().data() << "\" does not exist" << endl;
#endif
      QApplication::closeAllWindows();
      exit(0);
    }

    if(fileName.left(1) != "-") {
#if WITH_QT5
      cout << "Reading input file " << fileName.toLatin1().data() << endl;
#else
      cout << "Reading input file " << fileName.toAscii().data() << endl;
#endif
      readInputFile(fileName);
      remeshSlot();
    }
  }
}

// File -> Open...
//-----------------------------------------------------------------------------
void MainWindow::openSlot()
{
  QString defaultDirName = getDefaultDirName();

  QString fileName = QFileDialog::getOpenFileName(this, tr("Open geometry input file"), defaultDirName);

  if (!fileName.isEmpty()) {
    
    QFileInfo fi(fileName);
    QString absolutePath = fi.absolutePath();
    QDir::setCurrent(absolutePath);
    
  } else {
    
    logMessage("Unable to open file: file name is empty");
    return;

  }

  geometryInputFileName = fileName;

  operation_t *p = operation.next;
  operation_t *q = NULL;

  while(p != NULL) {
    if(p->select_set != NULL)
      delete [] p->select_set;

    q = p->next;

    if(p != NULL)
      delete p;

    p = q;
  }

  operations = 0;
  operation.next = NULL;

  saveDirName = "";
  readInputFile(fileName);

  if(egIni->isSet("automesh"))
    remeshSlot();
}


// Read input file and populate mesh generator's input structures:
//-----------------------------------------------------------------------------
void MainWindow::readInputFile(QString fileName)
{
  occInputOk = false;

  char cs[1024];

  QFileInfo fi(fileName);
  QString absolutePath = fi.absolutePath();
  QString baseName = fi.baseName();
  QString fileSuffix = fi.suffix();
  QString baseFileName = absolutePath + "/" + baseName;
#if WITH_QT5
  sprintf(cs, "%s", baseFileName.toLatin1().data());
#else
  sprintf(cs, "%s", baseFileName.toAscii().data());
#endif

  activeGenerator = GEN_UNKNOWN;
  tetlibInputOk = false;
  nglibInputOk = false;
  ngDim = 3;

  // Choose generator according to fileSuffix:
  //------------------------------------------
  if((fileSuffix == "smesh") || 
     (fileSuffix == "poly")) {
    
    if(!tetlibPresent) {
      logMessage("unable to mesh - tetlib unavailable");
      return;
    }

    activeGenerator = GEN_TETLIB;
    cout << "Selected tetlib for smesh/poly-format" << endl;

    in->deinitialize();
    in->initialize();
    in->load_poly(cs);

    tetlibInputOk = true;

  } else if(fileSuffix == "off") {

    if(!tetlibPresent) {
      logMessage("unable to mesh - tetlib unavailable");
      return;
    }

    activeGenerator = GEN_TETLIB;
    cout << "Selected tetlib for off-format" << endl;

    in->deinitialize();
    in->initialize();
    in->load_off(cs);

    tetlibInputOk = true;

  } else if(fileSuffix == "ply") {

    if(!tetlibPresent) {
      logMessage("unable to mesh - tetlib unavailable");
      return;
    }

    activeGenerator = GEN_TETLIB;
    cout << "Selected tetlib for ply-format" << endl;

    in->deinitialize();
    in->initialize();
    in->load_ply(cs);

    tetlibInputOk = true;

  } else if(fileSuffix == "mesh") {

    if(!tetlibPresent) {
      logMessage("unable to mesh - tetlib unavailable");
      return;
    }

    activeGenerator = GEN_TETLIB;
    cout << "Selected tetlib for mesh-format" << endl;

    in->deinitialize();
    in->initialize();
    in->load_medit(cs);
    
    tetlibInputOk = true;
    
  } else if(fileSuffix == "stl") {

    // for stl there are two alternative generators:
    if(meshControl->generatorType == GEN_NGLIB) {
      
      if(!nglibPresent) {
	logMessage("unable to mesh - nglib unavailable");
	return;
      }
      
      activeGenerator = GEN_NGLIB;
      cout << "Selected nglib for stl-format" << endl;

      stlFileName = fileName;
      
      nglibInputOk = true;
      
    } else {

      if(!tetlibPresent) {
	logMessage("unable to mesh - tetlib unavailable");
	return;
      }
      
      activeGenerator = GEN_TETLIB;
      cout << "Selected tetlib for stl-format" << endl;
      
      in->deinitialize();
      in->initialize();
      in->load_stl(cs);
      
      tetlibInputOk = true;
      
    }

  } else if((fileSuffix == "grd") ||
	    (fileSuffix == "FDNEUT") ||
	    (fileSuffix == "msh") ||
	    (fileSuffix == "mphtxt") ||
	    (fileSuffix == "inp") ||    
	    (fileSuffix == "unv") ||
            (fileSuffix == "plt")) {

    activeGenerator = GEN_ELMERGRID;
    cout << "Selected elmergrid" << endl;

#if WITH_QT5
    int errstat = elmergridAPI->loadElmerMeshStructure((const char*)(fileName.toLatin1()));
#else
    int errstat = elmergridAPI->loadElmerMeshStructure((const char*)(fileName.toAscii()));
#endif
    
    if (errstat)
      logMessage("loadElmerMeshStructure failed!");

    return;

#ifdef EG_OCC

  } else if( (fileSuffix.toLower() == "brep") ||
	     (fileSuffix.toLower() == "step") ||
	     (fileSuffix.toLower() == "stp")  || 
	     (fileSuffix.toLower() == "iges")  || 
	     (fileSuffix.toLower() == "igs") ) {

    meshControl->ui.nglibRadioButton->setChecked(true);
    meshControl->generatorType = GEN_NGLIB;
    activeGenerator = meshControl->generatorType;

    if(egIni->isSet("autoview"))
       cadView->show();

    occInputOk = cadView->readFile(fileName);

    ngDim = cadView->getDim();

    if(!occInputOk) {
      logMessage("Cad import: error: Unable to proceed with input file");
      cadView->close();
      return;
    }

    nglibInputOk = true;

#endif

  } else if( (fileSuffix.toLower() == "in2d") ) {
    
    if(!nglibPresent) {
      logMessage("unable to mesh - nglib unavailable");
      return;
    }
    
    activeGenerator = GEN_NGLIB;
    cout << "Selected nglib for in2d-format" << endl;
    
    in2dFileName = fileName;
    
    nglibInputOk = true;

    ngDim = 2;

  } else {

    logMessage("Unable to open file: file type unknown");
    activeGenerator = GEN_UNKNOWN;

    return;
  }
}
  


// Populate elmer's mesh structure and make GL-lists (tetlib):
//-----------------------------------------------------------------------------
void MainWindow::makeElmerMeshFromTetlib()
{
  meshutils->clearMesh(glWidget->getMesh());

  glWidget->setMesh(tetlibAPI->createElmerMeshStructure());

  glWidget->rebuildLists();

  logMessage("Input file processed");
}



// Populate elmer's mesh structure and make GL-lists (nglib):
//-----------------------------------------------------------------------------
void MainWindow::makeElmerMeshFromNglib()
{
  meshutils->clearMesh(glWidget->getMesh());
  nglibAPI->setDim(this->ngDim);
  nglibAPI->setNgmesh(ngmesh);

  glWidget->setMesh(nglibAPI->createElmerMeshStructure());
  glWidget->rebuildLists();

  logMessage("Input file processed");
}


// File -> Load mesh...
//-----------------------------------------------------------------------------
void MainWindow::loadSlot()
{
  QString defaultDirName = getDefaultDirName();

  QString dirName = QFileDialog::getExistingDirectory(this, tr("Open directory"), defaultDirName);

  if (!dirName.isEmpty()) {

    logMessage("Loading from directory " + dirName);

  } else {

    logMessage("Unable to load mesh: directory undefined");
    return;

  }
  
  loadElmerMesh(dirName);
}



// Import mesh files in elmer-format:
//-----------------------------------------------------------------------------
void MainWindow::loadElmerMesh(QString dirName)
{
  logMessage("Loading elmer mesh files");

  if(glWidget->hasMesh()) {
    glWidget->getMesh()->clear();
    glWidget->deleteMesh();
  }

  glWidget->newMesh();

#if WITH_QT5
  bool success = glWidget->getMesh()->load(dirName.toLatin1().data());
#else
  bool success = glWidget->getMesh()->load(dirName.toAscii().data());
#endif

  if(!success) {
    glWidget->getMesh()->clear();
    glWidget->deleteMesh();
    logMessage("Failed loading mesh files");
    return;
  }

  meshutils->findSurfaceElementEdges(glWidget->getMesh());
  meshutils->findSurfaceElementNormals(glWidget->getMesh());
  
  glWidget->rebuildLists();

  QDir::setCurrent(dirName);
  saveDirName = dirName;
  
  logMessage("Ready");
}


// File -> Save...
//-----------------------------------------------------------------------------
void MainWindow::saveSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to save mesh: no data");
    return;
  }

  if(!saveDirName.isEmpty()) {
    logMessage("Output directory " + saveDirName);
  } else {
    saveAsSlot();
    return;
  }

  generateSifSlot();
  saveElmerMesh(saveDirName);
}

// File -> Save as...
//-----------------------------------------------------------------------------
void MainWindow::saveAsSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to save mesh: no data");
    return;
  }

  QString defaultDirName = getDefaultDirName();

  saveDirName = QFileDialog::getExistingDirectory(this, tr("Open directory"), defaultDirName);

  if (!saveDirName.isEmpty()) {
    logMessage("Output directory " + saveDirName);
  } else {
    logMessage("Unable to save: directory undefined");
    return;
  }

  generateSifSlot();
  saveElmerMesh(saveDirName);
}


// File -> Save project...
//-----------------------------------------------------------------------------
void MainWindow::saveProjectSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to save project: no mesh");
    return;
  }

  generateSifSlot();

  QString defaultDirName = getDefaultDirName();

  QString projectDirName = QFileDialog::getExistingDirectory(this, tr("Open directory"), defaultDirName);

  if (!projectDirName.isEmpty()) {
    logMessage("Project directory " + projectDirName);
  } else {
    logMessage("Unable to save project: directory undefined");
    return;
  }

  progressBar->show();
  progressBar->setRange(0, 13);

  progressLabel->setText("Saving");
  progressLabel->show();

  // Create project document:
  //-------------------------
  progressBar->setValue(1);

  QDomDocument projectDoc("egproject");
  QDomElement contents = projectDoc.createElement("contents");
  projectDoc.appendChild(contents);

  //===========================================================================
  //                                  SAVE MESH
  //===========================================================================
  progressBar->setValue(2);
  logMessage("Saving mesh files...");
  saveElmerMesh(projectDirName);

  //===========================================================================
  //                        SAVE GEOMETRY INPUT FILE(S)
  //===========================================================================
  progressBar->setValue(3);

#ifdef Q_OS_LINUX
  QFileInfo fileInfo(geometryInputFileName);
  QString pathName(fileInfo.absolutePath());
  QString baseName(fileInfo.baseName());

  // System copy command:
  QString cmd("cp -f " + pathName + "/" + baseName + ".* "+ projectDirName);

  if(system(cmd.toLatin1().data()))
    logMessage("Geometry input file(s) not copied");

  QDomElement geomInput(projectDoc.createElement("geometryinputfile"));
  QDomText geomInputValue(projectDoc.createTextNode(fileInfo.fileName()));
  geomInput.appendChild(geomInputValue);
  contents.appendChild(geomInput);

#else
  QFileInfo geometryInputFileInfo(geometryInputFileName);
  QString baseName(geometryInputFileInfo.baseName());

  QString srcPathName(geometryInputFileInfo.absolutePath());
  QString dstPathName(QDir(projectDirName).absolutePath());

  // Avoid copying file(s) into it self:

  if( srcPathName  != dstPathName ) {
    QDirIterator srcDirIterator(srcPathName);

    while(srcDirIterator.hasNext()) {
      QString srcFileName(srcDirIterator.next());
      QFileInfo srcFileInfo(srcDirIterator.fileInfo());

      if(srcFileInfo.baseName() == baseName) {
	logMessage("Copying: " + srcFileName);

	QFile src(srcFileName);

	if(!src.open(QFile::ReadOnly)) {
	  logMessage("Unable to read: " + src.fileName());
	  continue;
	}

	QFile dst(dstPathName + "/" + srcFileInfo.fileName());

	if(!dst.open(QFile::WriteOnly)) {
	  logMessage("Unable to write: " + dst.fileName());
	  src.close();
	  continue;
	}

	QTextStream srcStream(&src);
	QTextStream dstStream(&dst);
	dstStream << srcStream.readAll();
	
	dst.close();
	src.close();
      }
    }

  } else {
    logMessage("Geometry input file(s) not copied");
  }

  QDomElement geomInput = projectDoc.createElement("geometryinputfile");
  QDomText geomInputValue = projectDoc.createTextNode(geometryInputFileInfo.fileName());
  geomInput.appendChild(geomInputValue);
  contents.appendChild(geomInput);
#endif

  //===========================================================================
  //                               SAVE OPERATIONS
  //===========================================================================
  progressBar->setValue(4);
  QDomElement ops = projectDoc.createElement("operations");
  contents.appendChild(ops);
  operation.appendToProject(&projectDoc, &ops);
  
  //===========================================================================
  //                              SAVE GENERAL SETUP
  //===========================================================================
  progressBar->setValue(5);
  logMessage("Saving menu contents... ");
  QDomElement gsBlock = projectDoc.createElement("generalsetup");
  projectDoc.documentElement().appendChild(gsBlock);
  generalSetup->appendToProject(&projectDoc, &gsBlock);

  //===========================================================================
  //                            SAVE PARALLEL SETTINGS
  //===========================================================================
  progressBar->setValue(6);
  QDomElement paraBlock = projectDoc.createElement("parallelsettings");
  projectDoc.documentElement().appendChild(paraBlock);
  parallel->appendToProject(&projectDoc, &paraBlock);

  //===========================================================================
  //                            SAVE MESH PARAMETERS
  //===========================================================================
  progressBar->setValue(7);
  QDomElement meshParams = projectDoc.createElement("meshparameters");
  projectDoc.documentElement().appendChild(meshParams);
  meshControl->appendToProject(&projectDoc, &meshParams);

  //===========================================================================
  //                            SAVE SOLVER PARAMETERS
  //===========================================================================
  progressBar->setValue(8);
  QDomElement speBlock = projectDoc.createElement("solverparameters");
  projectDoc.documentElement().appendChild(speBlock);

  for(int index = 0; index < solverParameterEditor.size(); index++) {
    SolverParameterEditor *spe = solverParameterEditor[index];

    if(!spe)
      continue;

    QDomElement item = projectDoc.createElement("item");
    item.setAttribute("index", QString::number(index));
    item.setAttribute("name", spe->solverName);
    speBlock.appendChild(item);
    spe->appendToProject(&projectDoc, &item);
  }

  //===========================================================================
  //                          SAVE DYNAMIC MENU CONTENTS
  //===========================================================================
  progressBar->setValue(9);
  saveProjectContents(projectDoc, "equation", equationEditor);
  saveProjectContents(projectDoc, "material", materialEditor);
  saveProjectContents(projectDoc, "bodyforce", bodyForceEditor);
  saveProjectContents(projectDoc, "initialcondition", initialConditionEditor);
  saveProjectContents(projectDoc, "boundarycondition", boundaryConditionEditor);

  //===========================================================================
  //                          SAVE SOLVER SPECIFIC OPTIONS
  //===========================================================================
  progressBar->setValue(10);
  QDomElement solverOptionsBlock = projectDoc.createElement("solverspecificoptions");
  projectDoc.documentElement().appendChild(solverOptionsBlock);

  for(int index = 0; index < solverParameterEditor.size(); index++) {
    SolverParameterEditor *spe = solverParameterEditor[index];

    if(!spe)
      continue;

    DynamicEditor *dynEdit = spe->generalOptions;
    
    if(!dynEdit)
      continue;

    QDomElement item = projectDoc.createElement("item");
    item.setAttribute("index", QString::number(index));
    item.setAttribute("name", spe->solverName);
    item.setAttribute("id", QString::number(dynEdit->ID));
    solverOptionsBlock.appendChild(item);

    dynEdit->dumpHash(&projectDoc, &item);
  }

  //===========================================================================
  //                            SAVE BODY PROPERTIES
  //===========================================================================
  progressBar->setValue(11);
  QDomElement bodyBlock = projectDoc.createElement("bodyproperties");
  projectDoc.documentElement().appendChild(bodyBlock);

  for(int index = 0; index < bodyPropertyEditor.size(); index++) {
    BodyPropertyEditor *bpe = bodyPropertyEditor[index];

    if(!bpe)
      continue;

    QDomElement item = projectDoc.createElement("item");
    item.setAttribute("index", QString::number(index));
    bodyBlock.appendChild(item);
    bpe->appendToProject(&projectDoc, &item);
  }

  //===========================================================================
  //                          SAVE BOUNDARY PROPERTIES
  //===========================================================================
  progressBar->setValue(12);
  QDomElement boundaryBlock = projectDoc.createElement("boundaryproperties");
  projectDoc.documentElement().appendChild(boundaryBlock);

  for(int index = 0; index < boundaryPropertyEditor.size(); index++) {
    BoundaryPropertyEditor *bpe = boundaryPropertyEditor[index];

    if(!bpe)
      continue;

    QDomElement item = projectDoc.createElement("item");
    item.setAttribute("index", QString::number(index));
    boundaryBlock.appendChild(item);
    bpe->appendToProject(&projectDoc, &item);
  }

  //===========================================================================
  //                             SAVE PROJECT DOCUMENT
  //===========================================================================
  progressBar->setValue(13);
  const int indent = 3;
  QFile projectFile("egproject.xml");
  projectFile.open(QIODevice::WriteOnly);
  QTextStream projectTextStream(&projectFile);
  projectDoc.save(projectTextStream, indent);

  saveDirName = projectDirName;
  logMessage("Ready");

  progressBar->hide();
  progressLabel->hide();
}



// Helper function for saveProject
//-----------------------------------------------------------------------------
void MainWindow::saveProjectContents(QDomDocument projectDoc,
				     QString blockName, 
				     QVector<DynamicEditor*>& editor)
{
  int Nmax = editor.size();

  QDomElement editorBlock = projectDoc.createElement(blockName);
  projectDoc.documentElement().appendChild(editorBlock);

  for(int i = 0; i < Nmax; i++) {
    DynamicEditor *de = editor[i];
    
    if(de->menuAction == NULL)
      continue;

    // Menu item number:
    QDomElement item = projectDoc.createElement("item");
    item.setAttribute("index", QString::number(i));
    editorBlock.appendChild(item);

    // Is active?
    QDomElement itemActive = projectDoc.createElement("active");
    QDomText itemActiveValue = projectDoc.createTextNode(QString::number(de->menuAction != NULL));
    itemActive.appendChild(itemActiveValue);
    item.appendChild(itemActive);
    
    // Name:
    if(de->menuAction != NULL) {
      QDomElement itemName = projectDoc.createElement("name");
      QDomText itemNameValue = projectDoc.createTextNode(de->nameEdit->text().trimmed());
      itemName.appendChild(itemNameValue);
      item.appendChild(itemName);
    }

    de->dumpHash(&projectDoc, &item);
  }
}



// File -> Load project...
//-----------------------------------------------------------------------------
void MainWindow::loadProjectSlot()
{
  QString defaultDirName = getDefaultDirName();

  QString projectDirName = QFileDialog::getExistingDirectory(this, tr("Open directory"), defaultDirName);

  if (!projectDirName.isEmpty()) {
    logMessage("Project directory: " + projectDirName);
  } else {
    logMessage("Unable to load project: directory undefined");
    return;
  }

  QDir::setCurrent(projectDirName);
  saveDirName = projectDirName;

  progressBar->show();
  progressBar->setRange(0, 14);

  progressLabel->setText("Loading");
  progressLabel->show();

  // Clear previous data:
  //----------------------
  progressBar->setValue(1);

  logMessage("Clearing model data");
  modelClearSlot();

  // Load project doc:
  //-------------------
  progressBar->setValue(2);

  logMessage("Loading project document...");
  QDomDocument projectDoc;
  QString errStr;
  int errRow;
  int errCol;
  QFile projectFile("egproject.xml");

  if(!projectFile.exists()) {
    QMessageBox::information(window(), tr("Project loader"),
			     tr("Project file does not exist"));

    progressBar->hide();
    progressLabel->hide();

    return;

  } else {  

    if(!projectDoc.setContent(&projectFile, true, &errStr, &errRow, &errCol)) {
      QMessageBox::information(window(), tr("Project loader"),
			       tr("Parse error at line %1, col %2:\n%3")
			       .arg(errRow).arg(errCol).arg(errStr));
      projectFile.close();

      progressBar->hide();
      progressLabel->hide();

      return;
    }
  }

  projectFile.close();	
  
  if(projectDoc.documentElement().tagName() != "contents") {
    QMessageBox::information(window(), tr("Project loader"),
			     tr("This is not a project file"));

    progressBar->hide();
    progressLabel->hide();

    return;
  }

  QDomElement contents = projectDoc.documentElement();

  //===========================================================================
  //                                 LOAD MESH
  //===========================================================================
  progressBar->setValue(3);
  logMessage("Loading mesh files...");
  loadElmerMesh(projectDirName);
  resetSlot();

  //===========================================================================
  //                          LOAD GEOMETRY INPUT FILE
  //===========================================================================
  progressBar->setValue(4);
  cout << "Loading geometry input file" << endl;
  QDomElement geomInput = contents.firstChildElement("geometryinputfile");
  geometryInputFileName = projectDirName + "/" + geomInput.text().trimmed();
  logMessage("Geometry input file: " + geometryInputFileName);
  readInputFile(geometryInputFileName);

  //===========================================================================
  //                               LOAD OPERATIONS
  //===========================================================================
  progressBar->setValue(5);
  QDomElement ops = contents.firstChildElement("operations");
  operations = operation.readFromProject(&projectDoc, &ops);

  //===========================================================================
  //                            LOAD GENERAL SETUP
  //===========================================================================
  progressBar->setValue(6);
  QDomElement gsBlock = contents.firstChildElement("generalsetup");
  generalSetup->readFromProject(&projectDoc, &gsBlock);

  //===========================================================================
  //                          LOAD PARALLEL SETTINGS
  //===========================================================================
  progressBar->setValue(7);
  QDomElement paraBlock = contents.firstChildElement("parallelsettings");
  parallel->readFromProject(&projectDoc, &paraBlock);

  //===========================================================================
  //                            LOAD MESH PARAMETERS
  //===========================================================================
  progressBar->setValue(8);
  QDomElement meshParams = contents.firstChildElement("meshparameters");
  meshControl->readFromProject(&projectDoc, &meshParams);

  //===========================================================================
  //                          LOAD SOLVER PARAMETERS
  //===========================================================================
  progressBar->setValue(9);
  QDomElement speBlock = contents.firstChildElement("solverparameters");

  QDomElement item = speBlock.firstChildElement("item");
  for( ; !item.isNull(); item = item.nextSiblingElement()) {
    int index = item.attribute("index").toInt();
    QString name = item.attribute("name");

    if(name.trimmed().isEmpty()) continue;

    // Find the real index for the current edf setup:
    int count = 0, realIndex = -1;
    QDomElement root = elmerDefs->documentElement();
    QDomElement elem = root.firstChildElement("PDE");
    while(!elem.isNull()) {
      QDomElement pdeName = elem.firstChildElement("Name");
      if(pdeName.text().trimmed() == name.trimmed()) realIndex = count;
      elem = elem.nextSiblingElement();
      count++;
    }

    if(realIndex < 0) {
      cout << "ERROR: The current edf setup conflicts with the project. Aborting." << endl;

      progressBar->hide();
      progressLabel->hide();

      return;
    }

    index = realIndex - 1;

    if(index < 0) {
      logMessage("Load project: solver parameters: index out of bounds");

      progressBar->hide();
      progressLabel->hide();

      return;
    }

    if(index >= solverParameterEditor.size())
      solverParameterEditor.resize(index + 1);

    if(!solverParameterEditor[index])
      solverParameterEditor[index] = new SolverParameterEditor;

    SolverParameterEditor *spe = solverParameterEditor[index];
    spe->readFromProject(&projectDoc, &item);
  }

#if 0
  // Changed the load order in 19 March 2009 for taking the "use as a body"
  // flags into account. The original boundary property loader is below.
  // 
  // Changed back to original 23 March 2009. Todo...
  //===========================================================================
  //                          LOAD BOUNDARY PROPERTIES
  //===========================================================================
  progressBar->setValue(10);
  QDomElement boundaryBlock = contents.firstChildElement("boundaryproperties");

  item = boundaryBlock.firstChildElement("item");
  for( ; !item.isNull(); item = item.nextSiblingElement()) {
    int index = item.attribute("index").toInt();

    if(index < 0) {
      logMessage("Load project: boundary properties: index out of bounds");
      
      progressBar->hide();
      progressLabel->hide();
      
      return;
    }

    if(index >= boundaryPropertyEditor.size())
      boundaryPropertyEditor.resize(index + 1);

    if(!boundaryPropertyEditor[index])
      boundaryPropertyEditor[index] = new BoundaryPropertyEditor;

    BoundaryPropertyEditor *bpe = boundaryPropertyEditor[index];

    bpe->readFromProject(&projectDoc, &item);

    if(bpe->ui.boundaryAsABody->isChecked()) {
      connect(bpe, SIGNAL(BoundaryAsABodyChanged(BoundaryPropertyEditor*, int)),
	      this, SLOT(boundaryAsABodyChanged(BoundaryPropertyEditor*, int)));

      populateBoundaryComboBoxes(bpe);

      bpe->ui.boundaryAsABody->toggle();
      bpe->ui.boundaryAsABody->toggle();
      bpe->ui.applyButton->click();
    }
  }
#endif

  //===========================================================================
  //                        LOAD DYNAMIC EDITOR CONTENTS
  //===========================================================================
  progressBar->setValue(11);
  QDomElement element = projectDoc.documentElement().firstChildElement("equation");
  loadProjectContents(element, equationEditor, "Equation");
  element = projectDoc.documentElement().firstChildElement("material");
  loadProjectContents(element, materialEditor, "Material");
  element = projectDoc.documentElement().firstChildElement("bodyforce");
  loadProjectContents(element, bodyForceEditor, "BodyForce");
  element = projectDoc.documentElement().firstChildElement("initialcondition");
  loadProjectContents(element, initialConditionEditor, "InitialCondition");
  element = projectDoc.documentElement().firstChildElement("boundarycondition");
  loadProjectContents(element, boundaryConditionEditor, "BoundaryCondition");


  //===========================================================================
  //                          LOAD SOLVER SPECIFIC OPTIONS
  //===========================================================================
  progressBar->setValue(12);
  QDomElement solverOptionsBlock = contents.firstChildElement("solverspecificoptions");

  for(item = solverOptionsBlock.firstChildElement("item"); 
      !item.isNull(); item = item.nextSiblingElement()) {
    
    int index = item.attribute("index").toInt();
    QString name = item.attribute("name");
    int id = item.attribute("id").toInt();

    if(name.trimmed().isEmpty()) continue;

    // Find the real index for the current edf setup:
    int count = 0, realIndex = -1;
    QDomElement root = elmerDefs->documentElement();
    QDomElement elem = root.firstChildElement("PDE");
    while(!elem.isNull()) {
      QDomElement pdeName = elem.firstChildElement("Name");
      if(pdeName.text().trimmed() == name.trimmed()) realIndex = count;
      elem = elem.nextSiblingElement();
      count++;
    }

    if(realIndex < 0) {
      cout << "ERROR: The current edf setup conflicts with the project. Aborting." << endl;

      progressBar->hide();
      progressLabel->hide();
      
      return;
    }

    index = realIndex - 1;

    if(index < 0) {
      logMessage("Load project: solver specific options: index out of bounds");

      progressBar->hide();
      progressLabel->hide();

      return;
    }

    if(index >= solverParameterEditor.size())
      solverParameterEditor.resize(index + 1);

    if(!solverParameterEditor[index])
      solverParameterEditor[index] = new SolverParameterEditor;

    SolverParameterEditor *spe = solverParameterEditor[index];
    spe->solverName = name;

    if(spe->generalOptions == NULL) 
      spe->generalOptions = new DynamicEditor;

    spe->generalOptions->setupTabs(elmerDefs, "Solver", id);
    spe->generalOptions->populateHash(&item);
    spe->ui.solverControlTabs->insertTab(0, spe->generalOptions->tabWidget->widget(id), "Solver specific options");	
  }

  //===========================================================================
  //                           LOAD BODY PROPERTIES
  //===========================================================================
  progressBar->setValue(13);
  QDomElement bodyBlock = contents.firstChildElement("bodyproperties");

  item = bodyBlock.firstChildElement("item");
  for(; !item.isNull(); item = item.nextSiblingElement()) {
    int index = item.attribute("index").toInt();

    if(index < 0) {
      logMessage("Load project: body properties: index out of bounds");

      progressBar->hide();
      progressLabel->hide();

      return;
    }

    if(index >= bodyPropertyEditor.size())
      bodyPropertyEditor.resize(index + 1);

    if(!bodyPropertyEditor[index])
      bodyPropertyEditor[index] = new BodyPropertyEditor;

    BodyPropertyEditor *bpe = bodyPropertyEditor[index];
    bpe->readFromProject(&projectDoc, &item);
  }

  //===========================================================================
  //                          LOAD BOUNDARY PROPERTIES
  //===========================================================================
  progressBar->setValue(13);
  QDomElement boundaryBlock = contents.firstChildElement("boundaryproperties");

  item = boundaryBlock.firstChildElement("item");
  for( ; !item.isNull(); item = item.nextSiblingElement()) {
    int index = item.attribute("index").toInt();

    if(index < 0) {
      logMessage("Load project: boundary properties: index out of bounds");
      
      progressBar->hide();
      progressLabel->hide();
      
      return;
    }

    if(index >= boundaryPropertyEditor.size())
      boundaryPropertyEditor.resize(index + 1);

    if(!boundaryPropertyEditor[index])
      boundaryPropertyEditor[index] = new BoundaryPropertyEditor;

    BoundaryPropertyEditor *bpe = boundaryPropertyEditor[index];
    bpe->readFromProject(&projectDoc, &item);
  }

  //===========================================================================
  //                              REGENERATE SIF
  //===========================================================================
  progressBar->setValue(14);
  if(glWidget->hasMesh()) {
    logMessage("Regenerating and saving the solver input file...");

    generateSifSlot();

    QFile file;
    QString sifName = generalSetup->ui.solverInputFileEdit->text().trimmed();
    file.setFileName(sifName);
    file.open(QIODevice::WriteOnly);
    QTextStream sif(&file);    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    sif << sifWindow->getTextEdit()->toPlainText();
    QApplication::restoreOverrideCursor();
    file.close();
    
    file.setFileName("ELMERSOLVER_STARTINFO");
    file.open(QIODevice::WriteOnly);
    QTextStream startinfo(&file);
#if WITH_QT5
    startinfo << sifName.toLatin1() << "\n1\n";    
#else
    startinfo << sifName.toAscii() << "\n1\n";    
#endif
    file.close();
  }

  logMessage("Ready");

  progressBar->hide();
  progressLabel->hide();
}


// Helper function for load project
//--------------------------------------------------------------------------------------------
void MainWindow::loadProjectContents(QDomElement projectElement, 
				     QVector<DynamicEditor*>& editor,
				     QString Mname)
{
  int Nmax = editor.size();

  QDomElement item = projectElement.firstChildElement("item");

  for(; !item.isNull(); item = item.nextSiblingElement()) {  
    int index = item.attribute("index").toInt();

    if(index < 0) {
      logMessage("Project loader: index out of bounds (dynamic editor)");
      return;
    }

    if(index >= editor.size())
      editor.resize(index+1);

    if(!editor[index])
      editor[index] = new DynamicEditor;
    
    DynamicEditor *de = editor[index];

    bool active = (item.firstChildElement("active").text().toInt() > 0);

    if(!active)
      continue;

    // Set up dynamic editor and connect:
    //------------------------------------
    QString itemName = item.firstChildElement("name").text().trimmed();

    de->setupTabs(elmerDefs, Mname, index);
    de->nameEdit->setText(itemName);
    de->applyButton->setText("Update");
    de->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
    de->discardButton->setText("Remove");
    de->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      
    const QString &tmpName = itemName;
    QAction *act = new QAction(tmpName, this);

    if(Mname == "Equation") {
      connect(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(pdeEditorFinishedSlot(int,int)));
      de->spareButton->setText("Edit Solver Settings");
      de->spareButton->show();
      de->spareButton->setIcon(QIcon(":/icons/tools-wizard.png"));      
      connect(de, SIGNAL(dynamicEditorSpareButtonClicked(int,int)), this, SLOT(editNumericalMethods(int,int)));
      equationMenu->addAction(act);
    }

    if(Mname == "Material") {
      connect(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(matEditorFinishedSlot(int,int)));
      de->spareButton->setText("Material library");
      de->spareButton->show();
      de->spareButton->setIcon(QIcon(":/icons/tools-wizard.png"));      
      connect(de, SIGNAL(dynamicEditorSpareButtonClicked(int,int)), this, SLOT(showMaterialLibrary(int,int)));
      materialMenu->addAction(act);
    }

    if(Mname == "BodyForce") {
      connect(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(bodyForceEditorFinishedSlot(int,int)));
      bodyForceMenu->addAction(act);
    }
    
    if(Mname == "InitialCondition") {
      connect(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(initialConditionEditorFinishedSlot(int,int)));
      initialConditionMenu->addAction(act);
    }
    
    if(Mname == "BoundaryCondition") {
      connect(de, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(boundaryConditionEditorFinishedSlot(int,int)));
      boundaryConditionMenu->addAction(act);
    }
    
    de->menuAction = act;
    
    if(Mname == "Equation") 
      createBodyCheckBoxes(BODY_EQUATION, de);
    
    if(Mname == "Material") 
      createBodyCheckBoxes(BODY_MATERIAL, de);
    
    if(Mname == "BodyForce") 
      createBodyCheckBoxes(BODY_FORCE, de);
    
    if(Mname == "InitialCondition") 
      createBodyCheckBoxes(BODY_INITIAL, de);
    
    if(Mname == "BoundaryCondition") 
      createBoundaryCheckBoxes(de);
    
    de->populateHash(&item);
  }
}


// Export mesh files in elmer-format:
//-----------------------------------------------------------------------------
void MainWindow::saveElmerMesh(QString dirName)
{
  logMessage("Saving elmer mesh files");

  QDir dir(dirName);

  if(!dir.exists())
    dir.mkdir(dirName);

  dir.setCurrent(dirName);

  // Save mesh files:
  //------------------
#if WITH_QT5
  glWidget->getMesh()->save(dirName.toLatin1().data());
#else
  glWidget->getMesh()->save(dirName.toAscii().data());
#endif

  // Save solver input file:
  //-------------------------
  QFile file;
  QString sifName = generalSetup->ui.solverInputFileEdit->text().trimmed();
  file.setFileName(sifName);
  file.open(QIODevice::WriteOnly);
  QTextStream sif(&file);

  QApplication::setOverrideCursor(Qt::WaitCursor);
  sif << sifWindow->getTextEdit()->toPlainText();
  QApplication::restoreOverrideCursor();

  file.close();

  // Save ELMERSOLVER_STARTINFO:
  //-----------------------------
  file.setFileName("ELMERSOLVER_STARTINFO");
  file.open(QIODevice::WriteOnly);
  QTextStream startinfo(&file);

#if WITH_QT5
  startinfo << sifName.toLatin1() << endl << "1" << endl;
#else
  startinfo << sifName.toAscii() << endl << "1" << endl;
#endif

  file.close();

  logMessage("Ready");
}


// File -> Exit
//-----------------------------------------------------------------------------
void MainWindow::closeMainWindowSlot()
{
    saveSlot();
    QApplication::closeAllWindows();
  // close();
}


// File -> Save picture as...
//-----------------------------------------------------------------------------
void MainWindow::savePictureSlot()
{
  QString defaultDirName(getDefaultDirName());

  pictureFileName = QFileDialog::getSaveFileName(this,	tr("Save picture"), defaultDirName, tr("Picture files (*.bmp *.jpg *.png *.pbm *.pgm *.ppm)"));
  
  if(pictureFileName.isEmpty()) {
    logMessage("File name is empty");
    return;
  }

  int delay = egIni->value("screenshotdelay").toInt();

  grabTimeLine->stop();
  grabTimeLine->setDuration(delay);
  grabTimeLine->setCurveShape(QTimeLine::LinearCurve);
  grabTimeLine->setDirection(QTimeLine::Backward);
  grabTimeLine->setFrameRange(0, 10);
  progressLabel->setText("Delay screen shot");
  progressLabel->show();
  progressBar->setRange(0, 10);
  progressBar->show();
  grabTimeLine->start();
}

void MainWindow::grabFrameSlot()
{
  progressLabel->hide();
  progressBar->hide();

  if(pictureFileName.isEmpty()) {
    logMessage("Unable to take screen shot - file name is empty");
    return;
  }

  QFileInfo fi(pictureFileName);
  QString suffix(fi.suffix());
  suffix.toUpper();
  
  int imageQuality(egIni->value("defaultimagequality").toInt());

  bool withAlpha(false);

  glWidget->updateGL();
  glReadBuffer(GL_BACK);

  QImage image(glWidget->grabFrameBuffer(withAlpha));

#if WITH_QT5
  bool success(image.save(pictureFileName, suffix.toLatin1(), imageQuality));
#else
  bool success(image.save(pictureFileName, suffix.toAscii(), imageQuality));
#endif
  
  if(!success)
    logMessage("Failed writing picture file");
}


//*****************************************************************************
//
//                                Model MENU
//
//*****************************************************************************


// Model -> Setup...
//-----------------------------------------------------------------------------
void MainWindow::modelSetupSlot()
{
  generalSetup->show();
}

//-----------------------------------------------------------------------------
void MainWindow::createBodyCheckBoxes(int which, DynamicEditor *pe)
{
  if(!glWidget->hasMesh()) return;

  if ( pe->spareScroll->widget() )
    delete pe->spareScroll->widget();

  QGridLayout *slayout = new QGridLayout;
  QLabel *l = new QLabel(tr("Apply to bodies:"));

  int count = 0, even = 0;

  slayout->addWidget(l,count,0);
  count++;

  for( int i=0; i < glWidget->bodyMap.count(); i++ )
  {
     int n=glWidget->bodyMap.key(i);
     if ( n >= 0 ) {
        int m = glWidget->bodyMap.value(n);

	if(m >= bodyPropertyEditor.size())
	  bodyPropertyEditor.resize(m + 1);

	if(!bodyPropertyEditor[m])
	  bodyPropertyEditor[m] = new BodyPropertyEditor;

	BodyPropertyEditor *body = bodyPropertyEditor[m];
	
        populateBodyComboBoxes(body);

        QString title = body->ui.nameEdit->text().trimmed();
        QCheckBox *a;

        if ( title.isEmpty() )
          a = new QCheckBox("Body " + QString::number(n));
        else
          a = new QCheckBox(title);

        DynamicEditor *p = NULL;

        switch(which) {
          case BODY_MATERIAL:
            p=body->material;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(materialBodyChanged(int)));
          break;
          case BODY_INITIAL:
            p=body->initial;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(initialBodyChanged(int)));
          break;
          case BODY_FORCE:
            p=body->force;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(forceBodyChanged(int)));
          break;
          case BODY_EQUATION:
            p=body->equation;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(equationBodyChanged(int)));
          break;
        }

        a->setProperty( "body", (qulonglong)body );
        a->setProperty( "editor", (qulonglong)pe );

        if ( p==pe )
          a->setChecked(true);
        else if ( p != NULL )
          a->setEnabled(false);
        else
          a->setChecked(false);

        slayout->addWidget(a,count,even);
        even = 1 - even;
        if (!even) count++;
     }
  }

  for( int i = 0; i < boundaryPropertyEditor.size(); i++ )
  {
    BoundaryPropertyEditor *boundary = boundaryPropertyEditor[i];

    if(!boundary)
      continue;

     if ( boundary->bodyProperties ) {
       BodyPropertyEditor *body = boundary->bodyProperties;
        populateBodyComboBoxes(body);

        QString title = body->ui.nameEdit->text().trimmed();
        QCheckBox *a;

        if ( title.isEmpty() )
          a = new QCheckBox("Body{Boundary " + QString::number(i)+ "}");
        else
          a = new QCheckBox(title);

        DynamicEditor *p = NULL;

        switch(which) {
          case BODY_MATERIAL:
            p=body->material;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(materialBodyChanged(int)));
          break;
          case BODY_INITIAL:
            p=body->initial;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(initialBodyChanged(int)));
          break;
          case BODY_FORCE:
            p=body->force;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(forceBodyChanged(int)));
          break;
          case BODY_EQUATION:
            p=body->equation;
            connect(a, SIGNAL(stateChanged(int)), this, SLOT(equationBodyChanged(int)));
          break;
        }

        a->setProperty( "body", (qulonglong)body );
        a->setProperty( "editor", (qulonglong)pe );

        if ( p==pe )
          a->setChecked(true);
        else if ( p != NULL )
          a->setEnabled(false);

        slayout->addWidget(a,count,even);
        even = 1-even;
        if (!even) count++;
     }
  }

  QGroupBox *box = new QGroupBox;
  box->setLayout(slayout);

  pe->spareScroll->setWidget(box);
  pe->spareScroll->setMinimumHeight(80);
  pe->spareScroll->show();
}


//-----------------------------------------------------------------------------

//*****************************************************************************

// Model -> Equation -> Add...
//-----------------------------------------------------------------------------
void MainWindow::addEquationSlot()
{
  DynamicEditor *pe = new DynamicEditor;
  equationEditor.append(pe);
  int current = equationEditor.size() - 1;

  pe->setupTabs(elmerDefs, "Equation", current);

  pe->applyButton->setText("Add");
  pe->applyButton->setIcon(QIcon(":/icons/list-add.png"));
  pe->discardButton->setText("Cancel");
  pe->discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
  pe->show();
  
  connect(pe, SIGNAL(dynamicEditorReady(int,int)),this, SLOT(pdeEditorFinishedSlot(int,int)));

  // Use "spareButton" to invoke solver parameter editor:
  pe->spareButton->setText("Edit Solver Settings");
  pe->spareButton->show();
  pe->spareButton->setIcon(QIcon(":/icons/tools-wizard.png"));
  connect(pe, SIGNAL(dynamicEditorSpareButtonClicked(int, int)), this, SLOT(editNumericalMethods(int, int)));

  // Equation is new - add to menu:
  const QString &equationName = pe->nameEdit->text().trimmed();
  QAction *act = new QAction(equationName, this);
  equationMenu->addAction(act);
  pe->menuAction = act;

  connect( pe->nameEdit, SIGNAL(textChanged(QString)), this,
        SLOT(dynamicEditorNameChange(QString)) );

  createBodyCheckBoxes(BODY_EQUATION,pe);
}


// signal (int, int) emitted by dynamic editor when "spare button" clicked:
//-----------------------------------------------------------------------------
void MainWindow::editNumericalMethods(int current, int id)
{
  QString title="";

  for(int i = 0; i < equationEditor.size(); i++) {
    // ** 23/04/09 **
    if(equationEditor[i]->ID == id) {
      title = equationEditor[i]->tabWidget->tabText(current);
      break;
    }
  }

  if(title == "General") {
    logMessage("No solver controls for 'General' equation options");
    return;
  }

  if(current >= solverParameterEditor.size())
    solverParameterEditor.resize(current + 1);

  if(!solverParameterEditor[current])
    solverParameterEditor[current] = new SolverParameterEditor;

  SolverParameterEditor *spe = solverParameterEditor[current];

  spe->setWindowTitle("Solver control for " + title);

  spe->solverName = title;

  if(spe->generalOptions == NULL) {
    spe->generalOptions = new DynamicEditor(spe);
    spe->generalOptions->setupTabs(elmerDefs, "Solver", current );
    spe->ui.solverControlTabs->insertTab(0, spe->generalOptions->tabWidget->widget(current), "Solver specific options");
    
#if 0
    for( int i=0; i < spe->generalOptions->tabWidget->count(); i++ )
      {
	if ( spe->generalOptions->tabWidget->tabText(i) == title )
	  {
	    spe->ui.solverControlTabs->insertTab(0, spe->generalOptions->tabWidget->widget(i),
						 "Solver specific options");
	    break;
	  }
      }
#endif
  }
  
  spe->show();
  spe->raise();
}


void MainWindow::dynamicEditorNameChange(QString t)
{
  for( int i = 0; i < bodyPropertyEditor.size(); i++ )
    {
      if(!bodyPropertyEditor[i])
	continue;
      
      if ( bodyPropertyEditor[i]->touched )
	populateBodyComboBoxes( bodyPropertyEditor[i] );
    }
  
   for( int i = 0; i < boundaryPropertyEditor.size(); i++ )
     {
       if(!boundaryPropertyEditor[i])
	 continue;
       
       if ( boundaryPropertyEditor[i]->touched )
	 populateBoundaryComboBoxes( boundaryPropertyEditor[i] );
     }
}

// signal (int,int) emitted by equation editor when ready:
//-----------------------------------------------------------------------------
void MainWindow::pdeEditorFinishedSlot(int signal, int id)
{  
  DynamicEditor *pe = equationEditor[id];
  
  const QString &equationName = pe->nameEdit->text().trimmed();

  bool signalOK = signal==MAT_OK || signal==MAT_APPLY;

  if((equationName.isEmpty()) && signalOK ) {
    logMessage("Refusing to add/update equation without name");
    return;
  }
  
  if( signalOK ) {
    if(pe->menuAction != NULL) {
      pe->menuAction->setText(equationName);
      logMessage("Equation updated");
      if ( signal==MAT_OK ) pe->close();
    }
  } else if (signal==MAT_NEW) {
    addEquationSlot();

  } else if(signal == MAT_DELETE) {

    for( int i = 0; i < bodyPropertyEditor.size(); i++ ) {
       BodyPropertyEditor *body = bodyPropertyEditor[i];

       if(!body)
	 continue;
       
       if ( body->equation == pe )
	 body->equation = NULL;
    }

    // Equation is not in menu:
    if(pe->menuAction == NULL) {
      logMessage("Ready");
      pe->close();
      return;
    }

    // Delete from menu:
    delete pe->menuAction;
    pe->menuAction = NULL;
    pe->close();

    int k = equationEditor.indexOf(pe);
    if(k>=0) equationEditor.remove(k);

    logMessage("Equation deleted");
  }
}

// signal (QAction*) emitted by equationMenu when an item has been selected:
//-----------------------------------------------------------------------------
void MainWindow::equationSelectedSlot(QAction* act)
{
  // Edit the selected material:
  for(int i = 0; i < equationEditor.size(); i++) {
    DynamicEditor *pe = equationEditor[i];
    if(pe->menuAction == act) {
      pe->applyButton->setText("Update");
      pe->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      pe->discardButton->setText("Remove");
      pe->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      createBodyCheckBoxes(BODY_EQUATION,pe);
      pe->show(); pe->raise();
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::equationBodyChanged(int state)
{
  QWidget *a = (QWidget *)QObject::sender();
  if(glWidget->getMesh()) {
     BodyPropertyEditor *body = (BodyPropertyEditor *)a->property("body").toULongLong();
     populateBodyComboBoxes(body);
     if ( state ) {
       DynamicEditor *mat  = (DynamicEditor *)a->property("editor").toULongLong();
       QString mat_name = mat->nameEdit->text().trimmed();
       int ind = body->ui.equationCombo->findText(mat_name);
       body->touched = true;
       body->equation = mat;
       body->ui.equationCombo->setCurrentIndex(ind);
     } else {
       body->equation = NULL;
       body->ui.equationCombo->setCurrentIndex(-1);
     }
  }
}


//*****************************************************************************

// Model -> Material -> Add...
//-----------------------------------------------------------------------------
void MainWindow::addMaterialSlot()
{
  DynamicEditor *pe = new DynamicEditor;
  materialEditor.append(pe);
  int current = materialEditor.size() - 1;

  pe->setupTabs(elmerDefs, "Material", current );
  pe->applyButton->setText("Add");
  pe->applyButton->setIcon(QIcon(":/icons/list-add.png"));
  pe->discardButton->setText("Cancel");
  pe->discardButton->setIcon(QIcon(":/icons/dialog-close.png"));

  connect(pe, SIGNAL(dynamicEditorReady(int,int)), this, SLOT(matEditorFinishedSlot(int,int)));

  // Use "spareButton" to invoke material library:
  pe->spareButton->setText("Material library");
  pe->spareButton->show();
  pe->spareButton->setIcon(QIcon(":/icons/tools-wizard.png"));
  connect(pe, SIGNAL(dynamicEditorSpareButtonClicked(int,int)), this, SLOT(showMaterialLibrary(int,int)));

  connect( pe->nameEdit, SIGNAL(textChanged(QString)), this,
        SLOT(dynamicEditorNameChange(QString)) );

  // Material is new - add to menu:
  const QString &materialName = pe->nameEdit->text().trimmed();
  QAction *act = new QAction(materialName, this);
  materialMenu->addAction(act);
  pe->menuAction = act;

  createBodyCheckBoxes(BODY_MATERIAL,pe);
  pe->show(); pe->raise();

}


void MainWindow::showMaterialLibrary(int tab, int ID)
{
  materialLibrary->editor = materialEditor[ID];
  materialLibrary->elmerDefs = this->elmerDefs;
  materialLibrary->show();
}

// signal (int,int) emitted by material editor when ready:
//-----------------------------------------------------------------------------
void MainWindow::matEditorFinishedSlot(int signal, int id)
{
  DynamicEditor *pe = materialEditor[id];
  
  const QString &materialName = pe->nameEdit->text().trimmed();

  bool signalOK = signal==MAT_OK || signal==MAT_APPLY;
  if( materialName.isEmpty() && signalOK ) {
    logMessage("Refusing to add/update material with no name");
    return;
  }
  
  if( signalOK ) {
    if(pe->menuAction != NULL) {
      pe->menuAction->setText(materialName);
      logMessage("Material updated");
      if ( signal == MAT_OK ) pe->close();
      return;
    }
  } else if ( signal==MAT_NEW ) {
    
    addMaterialSlot();

  } else if(signal == MAT_DELETE) {

    for( int i = 0; i < bodyPropertyEditor.size(); i++ ) {
      BodyPropertyEditor *body = bodyPropertyEditor[i];
      
      if(!body)
	continue;

      if ( body->material == pe )
	body->material = NULL;
    }

    // Material is not in menu:
    if(pe->menuAction == NULL) {
      logMessage("Ready");
      pe->close();
      return;
    }

    // Delete from menu:
    delete pe->menuAction;
    pe->menuAction = NULL;
    pe->close();

    int k = materialEditor.indexOf(pe);
    if(k>=0) materialEditor.remove(k);

    logMessage("Material deleted");

  } else {
    cout << "Matedit: unknown signal" << endl;
  }
}

// signal (QAction*) emitted by materialMenu when an item has been selected:
//-----------------------------------------------------------------------------
void MainWindow::materialSelectedSlot(QAction* act)
{
  // Edit the selected material:
  for(int i = 0; i < materialEditor.size(); i++) {
    DynamicEditor *pe = materialEditor[i];

    if(pe->menuAction == act) {
      pe->applyButton->setText("Update");
      pe->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      pe->discardButton->setText("Remove");
      pe->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      createBodyCheckBoxes(BODY_MATERIAL, pe);
      pe->show(); pe->raise();
    }
  }
}


void MainWindow::materialBodyChanged(int state)
{
  QWidget *a = (QWidget *)QObject::sender();
  if(glWidget->hasMesh()) {
     BodyPropertyEditor *body = (BodyPropertyEditor *)a->property("body").toULongLong();
     populateBodyComboBoxes( body);
 
     if ( state > 0 ) {
       DynamicEditor *mat = (DynamicEditor *)a->property("editor").toULongLong();
       QString mat_name = mat->nameEdit->text().trimmed();
       int ind = body->ui.materialCombo->findText(mat_name);

       body->touched = true;
       body->material = mat;
       body->ui.materialCombo->setCurrentIndex(ind);
     } else {
       body->material = NULL;
       body->ui.materialCombo->setCurrentIndex(-1);
     }
  }
}


//*****************************************************************************

// Model -> Body force -> Add...
//-----------------------------------------------------------------------------
void MainWindow::addBodyForceSlot()
{
  DynamicEditor *pe = new DynamicEditor;
  bodyForceEditor.append(pe);
  int current = bodyForceEditor.size() - 1;

  pe->setupTabs(elmerDefs, "BodyForce", current );

  pe->applyButton->setText("Add");
  pe->applyButton->setIcon(QIcon(":/icons/list-add.png"));
  pe->discardButton->setText("Cancel");
  pe->discardButton->setIcon(QIcon(":/icons/dialog-close.png"));

  connect(pe, SIGNAL(dynamicEditorReady(int,int)),
	  this, SLOT(bodyForceEditorFinishedSlot(int,int)));

  // Body force is new - add to menu:
  const QString &bodyForceName = pe->nameEdit->text().trimmed();
  QAction *act = new QAction(bodyForceName, this);
  bodyForceMenu->addAction(act);
  pe->menuAction = act;

  connect( pe->nameEdit, SIGNAL(textChanged(QString)), this,
        SLOT(dynamicEditorNameChange(QString)) );

  createBodyCheckBoxes( BODY_FORCE, pe );
  pe->show(); pe->raise();
}

// signal (int,int) emitted by body force editor when ready:
//-----------------------------------------------------------------------------
void MainWindow::bodyForceEditorFinishedSlot(int signal, int id)
{
  DynamicEditor *pe = bodyForceEditor[id];
  
  const QString &bodyForceName = pe->nameEdit->text().trimmed();

  bool signalOK = signal==MAT_OK || signal==MAT_APPLY;
  
  if((bodyForceName.isEmpty()) && signalOK ) {
    logMessage("Refusing to add/update body force with no name");
    return;
  }
  
  if( signalOK ) {
    if(pe->menuAction != NULL) {
      pe->menuAction->setText(bodyForceName);
      logMessage("Body force updated");
      if ( signal==MAT_OK ) pe->close();
    }

  } else if (signal==MAT_NEW) {
     addBodyForceSlot(); 

  } else if(signal == MAT_DELETE) {
    for( int i = 0; i < bodyPropertyEditor.size(); i++ ) {
      BodyPropertyEditor *body = bodyPropertyEditor[i];

      if(!body)
	continue;

      if ( body->force == pe )
	body->force = NULL;
    }

    if(pe->menuAction == NULL) {
      logMessage("Ready");
      pe->close();
      return;
    }
    
    // Delete from menu:
    delete pe->menuAction;
    pe->menuAction = NULL;
    pe->close();

    int k = bodyForceEditor.indexOf(pe);
    if(k>=0) bodyForceEditor.remove(k);

    logMessage("Body force deleted");
  }
}

// signal (QAction*) emitted by bodyForceMenu when an item has been selected:
//-----------------------------------------------------------------------------
void MainWindow::bodyForceSelectedSlot(QAction* act)
{
  // Edit the selected body force:
  for(int i = 0; i < bodyForceEditor.size(); i++) {
    DynamicEditor *pe = bodyForceEditor[i];
    if(pe->menuAction == act) {
      pe->applyButton->setText("Update");
      pe->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      pe->discardButton->setText("Remove");
      pe->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      createBodyCheckBoxes( BODY_FORCE, pe );
      pe->show(); pe->raise();
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::forceBodyChanged(int state)
{
  QWidget *a = (QWidget *)QObject::sender();
  if(glWidget->hasMesh()) {
     BodyPropertyEditor *body = (BodyPropertyEditor *)a->property("body").toULongLong();
     populateBodyComboBoxes(body);
 
     if ( state ) {
       DynamicEditor *mat  = (DynamicEditor *)a->property("editor").toULongLong();
       QString mat_name = mat->nameEdit->text().trimmed();
       int ind = body->ui.bodyForceCombo->findText(mat_name);

       body->touched = true;
       body->force = mat;
       body->ui.bodyForceCombo->setCurrentIndex(ind);
     } else {
       body->force = NULL;
       body->ui.bodyForceCombo->setCurrentIndex(-1);
     }
  }
}


//*****************************************************************************

// Model -> Initial condition -> Add...
//-----------------------------------------------------------------------------
void MainWindow::addInitialConditionSlot()
{
  DynamicEditor *pe = new DynamicEditor;
  initialConditionEditor.append(pe);
  int current = initialConditionEditor.size() - 1;
  
  pe->setupTabs(elmerDefs, "InitialCondition", current );

  pe->applyButton->setText("Add");
  pe->applyButton->setIcon(QIcon(":/icons/list-add.png"));
  pe->discardButton->setText("Cancel");
  pe->discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
  
  connect(pe, SIGNAL(dynamicEditorReady(int,int)),
	  this, SLOT(initialConditionEditorFinishedSlot(int,int)));

    // Initial condition is new - add to menu:
  const QString &initialConditionName = pe->nameEdit->text().trimmed();
  QAction *act = new QAction(initialConditionName, this);
  initialConditionMenu->addAction(act);
  pe->menuAction = act;

  connect( pe->nameEdit, SIGNAL(textChanged(QString)), this,
        SLOT(dynamicEditorNameChange(QString)) );

  createBodyCheckBoxes( BODY_INITIAL, pe );
  pe->show(); pe->raise();
}

// signal (int,int) emitted by initial condition editor when ready:
//-----------------------------------------------------------------------------
void MainWindow::initialConditionEditorFinishedSlot(int signal, int id)
{
  DynamicEditor *pe = initialConditionEditor[id];
  
  const QString &initialConditionName = pe->nameEdit->text().trimmed();
  
  bool signalOK = signal==MAT_OK || signal==MAT_APPLY;
  if((initialConditionName.isEmpty()) && signalOK ) {
    logMessage("Refusing to add/update initial condition with no name");
    return;
  }
  
  if( signalOK ) {
    if(pe->menuAction != NULL) {
      pe->menuAction->setText(initialConditionName);
      logMessage("Initial condition updated");
      if ( signal==MAT_OK ) pe->close();
    }
  } else if (signal==MAT_NEW ) {
     addInitialConditionSlot();

  } else if(signal == MAT_DELETE) {

    for( int i = 0; i < bodyPropertyEditor.size(); i++ ) {
      BodyPropertyEditor *body = bodyPropertyEditor[i];

      if(!body)
	continue;

      if ( body->initial == pe )
	body->initial = NULL;
    }

    // Initial condition is not in menu:
    if(pe->menuAction == NULL) {
      logMessage("Ready");
      pe->close();
      return;
    }
    
    // Delete from menu:
    delete pe->menuAction;
    pe->menuAction = NULL;
    pe->close();
    
    int k = initialConditionEditor.indexOf(pe);
    if(k>=0) initialConditionEditor.remove(k);

    logMessage("Initial condition deleted");
  }
}

// signal (QAction*) emitted by initialConditionMenu when item selected:
//-----------------------------------------------------------------------------
void MainWindow::initialConditionSelectedSlot(QAction* act)
{
  // Edit the selected initial condition:
  for(int i = 0; i < initialConditionEditor.size(); i++) {
    DynamicEditor *pe = initialConditionEditor[i];
    if(pe->menuAction == act) {
      pe->applyButton->setText("Update");
      pe->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      pe->discardButton->setText("Remove");
      pe->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      createBodyCheckBoxes( BODY_INITIAL, pe );
      pe->show(); pe->raise();
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::initialBodyChanged(int state)
{
  QWidget *a = (QWidget *)QObject::sender();
  if(glWidget->hasMesh()) {
     BodyPropertyEditor *body = (BodyPropertyEditor *)a->property("body").toULongLong();
     populateBodyComboBoxes( body);
 
     if ( state ) {
       DynamicEditor *mat  = (DynamicEditor *)a->property("editor").toULongLong();
       QString mat_name = mat->nameEdit->text().trimmed();
       int ind = body->ui.initialConditionCombo->findText(mat_name);
       body->touched = true;
       body->initial = mat;
       body->ui.initialConditionCombo->setCurrentIndex(ind);
     } else {
       body->initial = NULL;
       body->ui.initialConditionCombo->setCurrentIndex(-1);
     }
  }
}


//*****************************************************************************
//-----------------------------------------------------------------------------
void MainWindow::createBoundaryCheckBoxes(DynamicEditor *pe)
{
  if(!glWidget->hasMesh()) return;

  if ( pe->spareScroll->widget() ) {
    delete pe->spareScroll->widget();
  }

  QGridLayout *slayout = new QGridLayout;
  QLabel *l = new QLabel(tr("Apply to boundaries:"));
  int count=0,even=0;

  slayout->addWidget(l,count,0);
  count++;


  for( int i=0; i<glWidget->boundaryMap.count(); i++ )
  {
     int n=glWidget->boundaryMap.key(i);
     if ( n >= 0 ) {
       int m = glWidget->boundaryMap.value(n);

       if(m >= boundaryPropertyEditor.size())
	 boundaryPropertyEditor.resize(m + 1);

       if(!boundaryPropertyEditor[m])
	 boundaryPropertyEditor[m] = new BoundaryPropertyEditor;

	BoundaryPropertyEditor *boundary = boundaryPropertyEditor[m];
	
        populateBoundaryComboBoxes(boundary);

	// TODO: check this
        QString title =  ""; // boundary->ui.nameEdit->text().trimmed();
        QCheckBox *a;

        if ( title.isEmpty() )
          a = new QCheckBox("Boundary " + QString::number(n));
        else
          a = new QCheckBox(title);

        if (glWidget->stateBcColors) {
          int c[3];
          QPixmap pm(16,16);

          GLWidget::indexColors(c, n);
          pm.fill(qRgb(c[0], c[1], c[2]));
          a->setIcon(QIcon(pm));
        }

        DynamicEditor *p = NULL;

        p=boundary->condition;
        connect(a, SIGNAL(stateChanged(int)), this, SLOT(bcBoundaryChanged(int)));

        a->setProperty( "boundary", (qulonglong)boundary );
        a->setProperty( "condition", (qulonglong)pe );

        if ( p==pe )
          a->setChecked(true);
        else if ( p != NULL )
          a->setEnabled(false);

        slayout->addWidget(a,count,even);
        even = 1-even;
        if (!even) count++;
     }
  }

  QGroupBox *box = new QGroupBox;
  box->setLayout(slayout);

  pe->spareScroll->setWidget(box);
  pe->spareScroll->setMinimumHeight(80);
  pe->spareScroll->show();
}


//-----------------------------------------------------------------------------

// Model -> Boundary condition -> Add...
//-----------------------------------------------------------------------------
void MainWindow::addBoundaryConditionSlot()
{
  DynamicEditor *pe = new DynamicEditor;
  boundaryConditionEditor.append(pe);
  int current = boundaryConditionEditor.size() - 1;

  pe->setupTabs(elmerDefs, "BoundaryCondition", current );
  
  pe->applyButton->setText("Add");
  pe->applyButton->setIcon(QIcon(":/icons/list-add.png"));
  pe->discardButton->setText("Cancel");
  pe->discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
  pe->show();
  
  connect(pe, SIGNAL(dynamicEditorReady(int,int)),
	  this, SLOT(boundaryConditionEditorFinishedSlot(int,int)));

  // Boundary condition is new - add to menu:
  const QString &boundaryConditionName = pe->nameEdit->text().trimmed();
  QAction *act = new QAction(boundaryConditionName, this);
  boundaryConditionMenu->addAction(act);
  pe->menuAction = act;

  connect( pe->nameEdit, SIGNAL(textChanged(QString)), this,
        SLOT(dynamicEditorNameChange(QString)) );

  createBoundaryCheckBoxes(pe);
}

// signal (int,int) emitted by boundary condition editor when ready:
//-----------------------------------------------------------------------------
void MainWindow::boundaryConditionEditorFinishedSlot(int signal, int id)
{
  DynamicEditor *pe = boundaryConditionEditor[id];
  
  const QString &boundaryConditionName = pe->nameEdit->text().trimmed();
  
  bool signalOK = signal==MAT_OK || signal==MAT_APPLY;

  if((boundaryConditionName.isEmpty()) && signalOK ) {
    logMessage("Refusing to add/update boundary condition with no name");
    return;
  }
  
  if( signalOK ) {
    if(pe->menuAction != NULL) {
      pe->menuAction->setText(boundaryConditionName);
      logMessage("Boundary condition updated");
      if ( signal==MAT_OK ) pe->close();
    }
  } else if ( signal==MAT_NEW ) {
    addBoundaryConditionSlot();

  } else if(signal == MAT_DELETE) {

    pe->nameEdit->setText(QString());

    for( int i=0; i < boundaryPropertyEditor.size(); i++ ) {
      BoundaryPropertyEditor *bndry = boundaryPropertyEditor[i];
      
      if(!bndry)
	continue;

       if ( bndry->condition == pe ) {
           bndry->condition=NULL;
       }
    }

    // Boundary condition is not in menu:
    if(pe->menuAction == NULL) {
      logMessage("Ready");
      pe->close();
      return;
    }
    
    // Delete from menu:
    delete pe->menuAction;
    pe->menuAction = NULL;
    pe->close();

    int k = boundaryConditionEditor.indexOf(pe);
    if(k>=0) {
        boundaryConditionEditor.remove(k);
    }

    logMessage("Boundary condition deleted");
  }
}

// signal (QAction*) emitted by boundaryConditionMenu when item selected:
//-----------------------------------------------------------------------------
void MainWindow::boundaryConditionSelectedSlot(QAction* act)
{
  // Edit the selected boundary condition:
  for(int i = 0; i < boundaryConditionEditor.size(); i++) {
    DynamicEditor *pe = boundaryConditionEditor[i];
    if(pe->menuAction == act) {
      pe->applyButton->setText("Update");
      pe->applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      pe->discardButton->setText("Remove");
      pe->discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      createBoundaryCheckBoxes(pe);
      pe->show(); pe->raise();
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::bcBoundaryChanged(int state)
{
  QWidget *a = (QWidget *)QObject::sender();
  if(glWidget->hasMesh()) {
     BoundaryPropertyEditor *boundary = 
           (BoundaryPropertyEditor *)a->property("boundary").toULongLong();
     populateBoundaryComboBoxes(boundary);
 
     if ( state ) {
       DynamicEditor *mat  = (DynamicEditor *)a->property("condition").toULongLong();
       QString mat_name = mat->nameEdit->text().trimmed();
       int ind = boundary->ui.boundaryConditionCombo->findText(mat_name);
       boundary->touched = true;
       boundary->condition = mat;
       boundary->ui.boundaryConditionCombo->setCurrentIndex(ind);
     } else {
       boundary->condition = NULL;
       boundary->ui.boundaryConditionCombo->setCurrentIndex(-1);
     }
  }
}


// Model -> Set body properties
//-----------------------------------------------------------------------------
void MainWindow::bodyEditSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to open body editor - no mesh");
    bodyEditActive = false;
    synchronizeMenuToState();
    return;
  }

  bodyEditActive = !bodyEditActive;
  glWidget->bodyEditActive = bodyEditActive;

  if(bodyEditActive)
    bcEditActive = false;

  synchronizeMenuToState();

  if(bodyEditActive)
    logMessage("Double click a boundary to edit body properties");
}



// Model -> Set boundary conditions
//-----------------------------------------------------------------------------
void MainWindow::bcEditSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to open BC editor - no mesh");
    bcEditActive = false;
    synchronizeMenuToState();
    return;
  }

  bcEditActive = !bcEditActive;

  if(bcEditActive)
    bodyEditActive = false;

  synchronizeMenuToState();

  if(bcEditActive)
    logMessage("Double click a boundary to edit BCs");
}



// Model -> Summary...
//-----------------------------------------------------------------------------
void MainWindow::modelSummarySlot()
{
  mesh_t *mesh = glWidget->getMesh();
  QTextEdit *te = summaryEditor->ui.summaryEdit;
  te->clear();
  summaryEditor->show();

  if(mesh == NULL) {
    te->append("No mesh");
    return;
  }
  
  te->append("FINITE ELEMENT MESH");
  te->append("Mesh dimension: " + QString::number(mesh->getCdim()));
  te->append("Leading element dimension: " + QString::number(mesh->getDim()));
  te->append("Nodes: " + QString::number(mesh->getNodes()));
  te->append("Volume elements: " + QString::number(mesh->getElements()));
  te->append("Surface elements: " + QString::number(mesh->getSurfaces()));
  te->append("Edge elements: " + QString::number(mesh->getEdges()));
  te->append("Point elements: " + QString::number(mesh->getPoints()));
  te->append("");

  // This is almost duplicate info with the above, they might be fused in some way...
  te->append("ELEMENT TYPES");
  int *elementtypes = new int[828];
  for(int i=0;i<=827;i++)
    elementtypes[i] = 0;
  for(int i = 0; i < mesh->getElements(); i++)
    elementtypes[mesh->getElement(i)->getCode()] += 1;
  for(int i = 0; i < mesh->getSurfaces(); i++)
    elementtypes[mesh->getSurface(i)->getCode()] += 1;
  for(int i = 0; i < mesh->getEdges(); i++)
    elementtypes[mesh->getEdge(i)->getCode()] += 1;
  for(int i = 0; i < mesh->getPoints(); i++)
    elementtypes[mesh->getPoint(i)->getCode()] += 1;
  for(int i=827;i>0;i--)
    if(elementtypes[i])  te->append(QString::number(i) + ": " + QString::number(elementtypes[i]));
  te->append("");
  delete [] elementtypes;


  te->append("BOUNDING BOX");
  QString coordnames="XYZ";
  for(int j=0;j<3;j++) {
    double mincoord, maxcoord, coord;
    mincoord = maxcoord = mesh->getNode(0)->getX(j);
    for(int i = 0; i < mesh->getNodes(); i++) {
      coord = mesh->getNode(i)->getX(j);
      if(mincoord > coord) mincoord = coord;
      if(maxcoord < coord) maxcoord = coord;
    }
    te->append(coordnames[j]+"-coordinate: [ " + QString::number(mincoord) + " ,  " + QString::number(maxcoord)+" ]");
  }
  te->append("");

  // Check equations:
  int count = 0;
  for(int i = 0; i < equationEditor.size(); i++) {
    if(equationEditor[i]->menuAction != NULL)
      count++;
  }
  te->append("GENERAL");
  te->append("Equations: " + QString::number(count));

  // Check materials:
  count = 0;
  for(int i = 0; i < materialEditor.size(); i++) {
    if(materialEditor[i]->menuAction != NULL)
      count++;
  }
  te->append("Materials: " + QString::number(count));

  // Check boundary conditions:
  count = 0;
  for(int i = 0; i < boundaryConditionEditor.size(); i++) {
    if( boundaryConditionEditor[i]->touched) count++;
  }
  te->append("Boundary conditions: " + QString::number(count));

  // Check body properties:
  count = 0;
  for(int i = 0; i < bodyPropertyEditor.size(); i++) {

    if(!bodyPropertyEditor[i])
      continue;

    if(bodyPropertyEditor[i]->touched)
      count++;
  }

  te->append("Body properties: " + QString::number(count));
  te->append("");

  // Count volume bodies:
  //---------------------
  int undetermined = 0;
  int *tmp = new int[mesh->getElements()];
  for(int i = 0; i < mesh->getElements(); i++)
    tmp[i] = 0;

  for(int i = 0; i < mesh->getElements(); i++) {
    element_t *e = mesh->getElement(i);
    if(e->getNature() == PDE_BULK) {
      if(e->getIndex() >= 0)
	tmp[e->getIndex()]++;
      else
	undetermined++;
    }
  }

  te->append("VOLUME BODIES");
  count = 0;
  for(int i = 0; i < mesh->getElements(); i++) {
    if( tmp[i]>0 ) {
      count++;
      QString qs = "Body " + QString::number(i) + ": " 
	+ QString::number(tmp[i]) + " volume elements";

      element_t *e = mesh->getElement(i);
      int j = e->getIndex();

      if((j >= 0) && (j < bodyPropertyEditor.size()))
	if(bodyPropertyEditor[j] && bodyPropertyEditor[j]->touched) 
	  qs.append(" (Body property set)");
      
      te->append(qs);
    }
  }
  te->append("Undetermined: " + QString::number(undetermined));
  te->append("Total: " + QString::number(count) + " volume bodies");
  te->append("");

  delete [] tmp;

  // Count surface bodies:
  //---------------------
  undetermined = 0;
  tmp = new int[mesh->getSurfaces()];
  for(int i = 0; i < mesh->getSurfaces(); i++)
    tmp[i] = 0;

  for(int i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *s = mesh->getSurface(i);
    if(s->getNature() == PDE_BULK) {
      if(s->getIndex() >= 0)
	tmp[s->getIndex()]++;
      else
	undetermined++;
    }
  }

  te->append("SURFACE BODIES");
  count = 0;
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    if( tmp[i]>0 ) {
      count++;
      QString qs = "Body " + QString::number(i) + ": " 
	+ QString::number(tmp[i]) + " surface elements";

      surface_t *s = mesh->getSurface(i);
      int j = s->getIndex();

      if((j >= 0) && (j < bodyPropertyEditor.size()))
	if(bodyPropertyEditor[j] && bodyPropertyEditor[j]->touched)
	  qs.append(" (Body property set)");

      te->append(qs);
    }
  }
  te->append("Undetermined: " + QString::number(undetermined));
  te->append("Total: " + QString::number(count) + " surface bodies");
  te->append("");

  delete [] tmp;

  // Count edge bodies:
  //---------------------
  undetermined = 0;
  tmp = new int[mesh->getEdges()];
  for(int i = 0; i < mesh->getEdges(); i++)
    tmp[i] = 0;

  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *e = mesh->getEdge(i);
    if(e->getNature() == PDE_BULK) {
      if(e->getIndex() >= 0)
	tmp[e->getIndex()]++;
      else
	undetermined++;
    }
  }

  te->append("EDGE BODIES");
  count = 0;
  for(int i = 0; i < mesh->getEdges(); i++) {
    if( tmp[i]>0 ) {
      count++;
      QString qs = "Body " + QString::number(i) + ": " 
	+ QString::number(tmp[i]) + " edge elements";

      edge_t *e = mesh->getEdge(i);
      int j = e->getIndex();

      if((j >= 0) && (j < bodyPropertyEditor.size()))
	if(bodyPropertyEditor[j] && bodyPropertyEditor[j]->touched) 
	  qs.append(" (Body property set)");

      te->append(qs);
    }
  }
  te->append("Undetermined: " + QString::number(undetermined));
  te->append("Total: " + QString::number(count) + " edge bodies");
  te->append("");

  delete [] tmp;

  // Count surface boundaries:
  //--------------------------
  undetermined = 0;
  tmp = new int[mesh->getSurfaces()];
  for(int i = 0; i < mesh->getSurfaces(); i++)
    tmp[i] = 0;
  
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    surface_t *s = mesh->getSurface(i);
    if(s->getNature() == PDE_BOUNDARY) {
      if(s->getIndex() >= 0)
	tmp[s->getIndex()]++;
      else
	undetermined++;
    }
  }

  te->append("SURFACE BOUNDARIES");
  count = 0;
  for(int i = 0; i < mesh->getSurfaces(); i++) {
    if( tmp[i]>0 ) {
      count++;
      QString qs = "Boundary " + QString::number(i) + ": " 
	+ QString::number(tmp[i]) + " surface elements";
      
      surface_t *s = mesh->getSurface(i);     
      int j = s->getIndex();
      if((j >= 0) &&(j < boundaryConditionEditor.size()))
	if(boundaryConditionEditor[j]->touched)
	  qs.append(" (BC set)");

      te->append(qs);
    }
  }
  te->append("Undetermined: " + QString::number(undetermined));
  te->append("Total: " + QString::number(count) + " surface boundaries");
  te->append("");

  delete [] tmp;

  // Count edge boundaries:
  //--------------------------
  undetermined = 0;
  tmp = new int[mesh->getEdges()];
  for(int i = 0; i < mesh->getEdges(); i++)
    tmp[i] = 0;
  
  for(int i = 0; i < mesh->getEdges(); i++) {
    edge_t *e = mesh->getEdge(i);
    if(e->getNature() == PDE_BOUNDARY) {
      if(e->getIndex() >= 0)
	tmp[e->getIndex()]++;
      else
	undetermined++;
    }
  }

  te->append("EDGE BOUNDARIES");
  count = 0;
  for(int i = 0; i < mesh->getEdges(); i++) {
    if( tmp[i]>0 ) {
      count++;
      QString qs = "Boundary " + QString::number(i) + ": " 
	+ QString::number(tmp[i]) + " edge elements";

      edge_t *e = mesh->getEdge(i);
      int j = e->getIndex();
      if((j >= 0) && (j < boundaryConditionEditor.size()))
	if( boundaryConditionEditor[j]->touched)
	  qs.append(" (BC set)");

      te->append(qs);
    }
  }
  te->append("Undetermined: " + QString::number(undetermined));
  te->append("Total: " + QString::number(count) + " edge boundaries");
  te->append("");

  delete [] tmp;
}


// Model -> Clear
//-----------------------------------------------------------------------------
void MainWindow::modelClearSlot()
{
  // clear equations:
  for(int i = 0; i < equationEditor.size(); i++) {
    DynamicEditor *pe = equationEditor[i];
    if(pe->menuAction != NULL)
      delete pe->menuAction;
  }

  for(int i = 0; i < equationEditor.size(); i++)
    delete equationEditor[i];

  equationEditor.clear();

  // clear materials:
  for(int i = 0; i < materialEditor.size(); i++) {
    DynamicEditor *de = materialEditor[i];
    if(de->menuAction != NULL)
      delete de->menuAction;
  }

  for(int i = 0; i < materialEditor.size(); i++)
    delete materialEditor[i];

  materialEditor.clear();

  // clear body forces:
  for(int i = 0; i < bodyForceEditor.size(); i++) {
    DynamicEditor *de = bodyForceEditor[i];
    if(de->menuAction != NULL)
      delete de->menuAction;
  }
  
  for(int i = 0; i < bodyForceEditor.size(); i++)
    delete bodyForceEditor[i];

  bodyForceEditor.clear();

  // clear initial conditions:
  for(int i = 0; i < initialConditionEditor.size(); i++) {
    DynamicEditor *de = initialConditionEditor[i];
    if(de->menuAction != NULL)
      delete de->menuAction;
  }
  
  for(int i = 0; i < initialConditionEditor.size(); i++)
    delete initialConditionEditor[i];

  initialConditionEditor.clear();

  // clear boundary conditions:
  for(int i = 0; i < boundaryConditionEditor.size(); i++) {
    DynamicEditor *de = boundaryConditionEditor[i];
    if(de->menuAction != NULL)
      delete de->menuAction;
  }

  for(int i = 0; i < boundaryConditionEditor.size(); i++)
    if(boundaryConditionEditor[i])
      delete boundaryConditionEditor[i];

  boundaryConditionEditor.clear();

  // clear boundary setting:
  for(int i = 0; i < boundaryPropertyEditor.size(); i++)
    if(boundaryPropertyEditor[i])
      delete boundaryPropertyEditor[i];
  
  boundaryPropertyEditor.clear();

  // clear body settings:
  for(int i = 0; i < bodyPropertyEditor.size(); i++)
    if(bodyPropertyEditor[i])
      delete bodyPropertyEditor[i];
  
  bodyPropertyEditor.clear();
}


//*****************************************************************************
//
//                                View MENU
//
//*****************************************************************************


// View -> Full screen
//-----------------------------------------------------------------------------
void MainWindow::viewFullScreenSlot()
{
  if(!isFullScreen()) {
    cout << "Switching to full screen mode" << endl;
    menuBar()->hide();
    statusBar()->hide();
    fileToolBar->hide();
    editToolBar->hide();
    meshToolBar->hide();
    solverToolBar->hide();
    this->showFullScreen();
  } else {
    viewNormalModeSlot();
  }
  synchronizeMenuToState();
}

// Return to normal mode (GLWidget emits (void) when esc is pressed)...
//-----------------------------------------------------------------------------
void MainWindow::viewNormalModeSlot()
{
  if(isFullScreen()) {
    cout << "Switching to normal window mode" << endl;
    this->showNormal();
    menuBar()->show();
    statusBar()->show();
    if(!egIni->isSet("hidetoolbars")) {
      fileToolBar->show();
      editToolBar->show();
      meshToolBar->show();
      solverToolBar->show();
    }
  }
  synchronizeMenuToState();
  statusBar()->showMessage(tr("Ready"));
}


// Context menu event (usually mouse has been right clicked)...
//-----------------------------------------------------------------------------
void MainWindow::contextMenuEvent(QContextMenuEvent *event)
{
  // if(isFullScreen())
  contextMenu->popup(event->globalPos());
}


// View -> Surface mesh
//-----------------------------------------------------------------------------
void MainWindow::hidesurfacemeshSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There is no surface mesh to hide/show");
    return;
  }
  
  glWidget->stateDrawSurfaceMesh = !glWidget->stateDrawSurfaceMesh;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == SURFACEMESHLIST) 
    {
      l->setVisible(glWidget->stateDrawSurfaceMesh);

      // do not set visible if the parent surface list is hidden
      int p = l->getParent();
      if(p >= 0) {
	list_t *lp = glWidget->getList(p);
	if(!lp->isVisible())
	  l->setVisible(false);
      }
    }
  }

  synchronizeMenuToState();  

  if(!glWidget->stateDrawSurfaceMesh) 
    logMessage("Surface mesh hidden");
  else
    logMessage("Surface mesh shown");
}


// View -> Volume mesh
//-----------------------------------------------------------------------------
void MainWindow::hidevolumemeshSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There is no volume mesh to hide/show");
    return;
  }

  glWidget->stateDrawVolumeMesh = !glWidget->stateDrawVolumeMesh;


  for(int i = 0; i < lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == VOLUMEMESHLIST)
      l->setVisible(glWidget->stateDrawVolumeMesh);
  }


  synchronizeMenuToState();  

  if(!glWidget->stateDrawVolumeMesh) 
    logMessage("Volume mesh hidden");
  else
    logMessage("Volume mesh shown");

}


// View -> Sharp edges
//-----------------------------------------------------------------------------
void MainWindow::hidesharpedgesSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There are no sharp edges to hide/show");
    return;
  }

  glWidget->stateDrawSharpEdges = !glWidget->stateDrawSharpEdges;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == SHARPEDGELIST)  
      l->setVisible(glWidget->stateDrawSharpEdges);
  }

  
  synchronizeMenuToState();

  if ( !glWidget->stateDrawSharpEdges ) 
    logMessage("Sharp edges hidden");
  else 
    logMessage("Sharp edges shown");
}


// View -> Coordinates
//-----------------------------------------------------------------------------
void MainWindow::viewCoordinatesSlot()
{
  if( glWidget->toggleCoordinates() )
    logMessage("Coordinates shown");
  else 
    logMessage("Cordinates hidden");

  synchronizeMenuToState();
}



// View -> Select defined edges
//-----------------------------------------------------------------------------
void MainWindow::selectDefinedEdgesSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There are no entities from which to select");
    return;
  }

  // At the moment only edges are included in search:
  int nmax = 0;
  for( int i=0; i<glWidget->boundaryMap.count(); i++ ) {
    int n = glWidget->boundaryMap.key(i);
    if(n > nmax) nmax = n;
  }

  bool *activeboundary = new bool[nmax+1];
  for (int i=0;i<=nmax;i++)
    activeboundary[i] = false;

  for( int i=0; i<glWidget->boundaryMap.count(); i++ ) {
    int n=glWidget->boundaryMap.key(i);
    if ( n >= 0 ) {
      int m = glWidget->boundaryMap.value(n);
      
      if(m >= boundaryPropertyEditor.size())
	boundaryPropertyEditor.resize(m + 1);

      if(!boundaryPropertyEditor[m])
	boundaryPropertyEditor[m] = new BoundaryPropertyEditor;

      BoundaryPropertyEditor *boundary = boundaryPropertyEditor[m];
      activeboundary[n] = boundary->condition;
    }
  }

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == EDGELIST) {
      int j = l->getIndex();
      if( j < 0 ) continue;
      
      // *** TODO ***
      //
      // This is wrong: Comparing body indices with boundary indices
      if( activeboundary[j] ) l->setSelected(true);
    }
  }

  for( int i = 0; i < mesh->getEdges(); i++ ) {
    edge_t *edge = mesh->getEdge(i);
    if( edge->getNature() == PDE_BOUNDARY ) { 
      int j = edge->getIndex();
      if( j < 0) continue;
      if( activeboundary[j] ) edge->setSelected(true);
    }
  }
  delete [] activeboundary;

  glWidget->rebuildEdgeLists();
  glWidget->updateGL();
  
  logMessage("Defined edges selected");
}


// View -> Select defined surfaces
//-----------------------------------------------------------------------------
void MainWindow::selectDefinedSurfacesSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There are no entities from which to select");
    return;
  }

  // At the moment only surfaces are included in search:
  int nmax = 0;
  for( int i=0; i<glWidget->bodyMap.count(); i++ ) {
    int n = glWidget->bodyMap.key(i);
    if(n > nmax) nmax = n;
  }

  bool *activebody = new bool[nmax+1];
  for (int i=0;i<=nmax;i++)
    activebody[i] = false;

  for( int i=0; i<glWidget->bodyMap.count(); i++ ) {
    int n=glWidget->bodyMap.key(i);
    if ( n >= 0 ) {
      int m = glWidget->bodyMap.value(n);

      BodyPropertyEditor *body = bodyPropertyEditor[m];

      if(!body) {
	cout << "MainWindow: Body index out of bounds" << endl;
	continue;
      }

      activebody[n] = body->material && body->equation;
    }
  }

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == SURFACELIST) {
      int j = l->getIndex();
      if( j < 0 ) continue;

      // *** TODO ***
      //
      // This is wrong: Comparing body indices with boundary indexes
      if( activebody[j] ) l->setSelected(true);
    }
  }

  for( int i = 0; i < mesh->getSurfaces(); i++ ) {
    surface_t *surface = mesh->getSurface(i);
    if( surface->getNature() == PDE_BULK ) { 
      int j = surface->getIndex();
      if( j < 0) continue;
      if( activebody[j] ) surface->setSelected(true);
    }
  }
  delete [] activebody;

  glWidget->rebuildSurfaceLists();
  glWidget->updateGL();
  
  logMessage("Defined surfaces selected");
}



// View -> Select all surfaces
//-----------------------------------------------------------------------------
void MainWindow::selectAllSurfacesSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There are no surfaces to select");
    return;
  }

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->getType() == SURFACELIST)
    {
      l->setSelected(true);
      for( int j=0; j < mesh->getSurfaces(); j++ ) {
        surface_t *surf = mesh->getSurface(j);
        if( l->getIndex() == surf->getIndex() )
          surf->setSelected(l->isSelected());
      }
    }
  }

  glWidget->rebuildSurfaceLists();
  glWidget->updateGL();
  
  logMessage("All surfaces selected");
}



// View -> Select all edges
//-----------------------------------------------------------------------------
void MainWindow::selectAllEdgesSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There are no edges to select");
    return;
  }

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);

    if(l->getType() == EDGELIST)
      l->setSelected(true);

    for(int j = 0; j < mesh->getEdges(); j++ ) {
      edge_t *edge = mesh->getEdge(j);
      if( l->getIndex() == edge->getIndex() )
	edge->setSelected(l->isSelected());
    }
  }

  glWidget->rebuildEdgeLists();
  glWidget->updateGL();
  
  logMessage("All edges selected");
}



// View -> Hide/Show selected
//-----------------------------------------------------------------------------
void MainWindow::hideselectedSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("There is nothing to hide/show");
    return;
  }

  bool something_selected = false;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    something_selected |= l->isSelected();
  }

  if(!something_selected) {
    logMessage("Nothing selected");
    return;
  }

  bool vis = false;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected()) {
      l->setVisible(!l->isVisible());
      if(l->isVisible())
	vis = true;

      // hide the child surface edge list if parent is hidden
      int c = l->getChild();
      if(c >= 0) {
	list_t *lc = glWidget->getList(c);
	lc->setVisible(l->isVisible());
	if(!glWidget->stateDrawSurfaceMesh)
	  lc->setVisible(false);
      }
    }
  }
  glWidget->updateGL();
  
  if( !vis )
    logMessage("Selected objects hidden");
  else 
    logMessage("Selected objects shown");
}



// View -> Show all
//-----------------------------------------------------------------------------
void MainWindow::showallSlot()
{
  int lists = glWidget->getLists();
  
  glWidget->stateDrawSurfaceMesh = true;
  glWidget->stateDrawSharpEdges = true;
  glWidget->stateDrawSurfaceElements = true;
  glWidget->stateDrawEdgeElements = true;

  synchronizeMenuToState();

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    l->setVisible(true);
  }

  logMessage("All objects visible");
}



// View -> Reset model view
//-----------------------------------------------------------------------------
void MainWindow::resetSlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();
  
  if(mesh == NULL) {
    logMessage("There is nothing to reset");
    return;
  }

  glWidget->stateFlatShade = true;
  glWidget->stateDrawSurfaceMesh = true;
  glWidget->stateDrawSharpEdges = true;
  glWidget->stateDrawSurfaceElements = true;
  glWidget->stateDrawEdgeElements = true;
  glWidget->stateDrawSurfaceNumbers = false;
  glWidget->stateDrawEdgeNumbers = false;
  glWidget->stateDrawNodeNumbers = false;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    l->setVisible(true);
    l->setSelected(false);

    for( int j=0; j < mesh->getSurfaces(); j++ ) {
      surface_t *surf = mesh->getSurface(j);
      if( l->getIndex() == surf->getIndex() )
        surf->setSelected(l->isSelected());
    }
    for( int j=0; j<mesh->getEdges(); j++ ) {
      edge_t *edge = mesh->getEdge(j);
      if( l->getIndex() == edge->getIndex() )
        edge->setSelected(l->isSelected());
    }
  }

  glWidget->stateBcColors = false;
  glWidget->stateBodyColors = false;

  glLoadIdentity();
  glWidget->rebuildLists();
  glWidget->updateGL();

  synchronizeMenuToState();
  logMessage("Reset model view");
}



// View -> Shade model -> Flat
//-----------------------------------------------------------------------------
void MainWindow::flatShadeSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to change shade model when mesh is empty");
    return;
  }

  glWidget->stateFlatShade = true;
  glWidget->rebuildSurfaceLists();
  glWidget->updateGL();

  synchronizeMenuToState();
  logMessage("Shade model: flat");
}


// View -> Shade model -> Smooth
//-----------------------------------------------------------------------------
void MainWindow::smoothShadeSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to change shade model when mesh is empty");
    return;
  }

  glWidget->stateFlatShade = false;
  glWidget->rebuildSurfaceLists();
  glWidget->updateGL();

  synchronizeMenuToState();
  logMessage("Shade model: smooth");
}

// View -> Projection -> Orthogonal
//-----------------------------------------------------------------------------
void MainWindow::orthoSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to change projection when mesh is empty");
    return;
  }

  glWidget->stateOrtho = true;
  glWidget->changeProjection();
  glWidget->updateGL();

  synchronizeMenuToState();
  logMessage("Projection: orthogonal");
}

// View -> Projection -> Perspective
//-----------------------------------------------------------------------------
void MainWindow::perspectiveSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to change projection when mesh is empty");
    return;
  }

  glWidget->stateOrtho = false;
  glWidget->changeProjection();
  glWidget->updateGL();

  synchronizeMenuToState();
  logMessage("Projection: perspective");
}


// View -> Show numbering -> Surface numbering
//-----------------------------------------------------------------------------
void MainWindow::showSurfaceNumbersSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to show surface element numbering when mesh is empty");
    return;
  }
  glWidget->stateDrawSurfaceNumbers = !glWidget->stateDrawSurfaceNumbers;
  glWidget->updateGL();
  synchronizeMenuToState();

  if(glWidget->stateDrawSurfaceNumbers) 
    logMessage("Surface element numbering turned on");
  else
    logMessage("Surface element numbering turned off");   
}

// View -> Show numbering -> Edge numbering
//-----------------------------------------------------------------------------
void MainWindow::showEdgeNumbersSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to show edge element numbering when mesh is empty");
    return;
  }
  glWidget->stateDrawEdgeNumbers = !glWidget->stateDrawEdgeNumbers;
  glWidget->updateGL();
  synchronizeMenuToState();

  if(glWidget->stateDrawEdgeNumbers) 
    logMessage("Edge element numbering turned on");
  else
    logMessage("Edge element numbering turned off");   
}

// View -> Numbering -> Node numbers
//-----------------------------------------------------------------------------
void MainWindow::showNodeNumbersSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to show node numbering when mesh is empty");
    return;
  }
  glWidget->stateDrawNodeNumbers = !glWidget->stateDrawNodeNumbers;
  glWidget->updateGL();
  synchronizeMenuToState();
  
  if(glWidget->stateDrawNodeNumbers) 
    logMessage("Node numbering turned on");
  else 
    logMessage("Node numbering turned off");    
}

// View -> Numbering -> Boundary index
//-----------------------------------------------------------------------------
void MainWindow::showBoundaryIndexSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to show boundary indices when mesh is empty");
    return;
  }
  glWidget->stateDrawBoundaryIndex = !glWidget->stateDrawBoundaryIndex;
  glWidget->updateGL();
  synchronizeMenuToState();
  
  if(glWidget->stateDrawBoundaryIndex) 
    logMessage("Boundary indices visible");
  else 
    logMessage("Boundary indices hidden");    
}


// View -> Numbering -> Body index
//-----------------------------------------------------------------------------
void MainWindow::showBodyIndexSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Refusing to show body indices when mesh is empty");
    return;
  }

  glWidget->stateDrawBodyIndex = !glWidget->stateDrawBodyIndex;
  glWidget->updateGL();
  synchronizeMenuToState();
  
  if(glWidget->stateDrawBodyIndex) 
    logMessage("Body indices visible");
  else 
    logMessage("Body indices hidden");    
}


// View -> Colors -> GL controls
//-----------------------------------------------------------------------------
void MainWindow::glControlSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("No mesh - unable to set GL parameters when the mesh is empty");
    return;
  }

  glControl->glWidget = this->glWidget;
  glControl->show();
}


// View -> Colors -> Boundaries
//-----------------------------------------------------------------------------
void MainWindow::colorizeBoundarySlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("No mesh - unable to colorize boundaries");
    return;
  }

  glWidget->stateBcColors = !glWidget->stateBcColors;

  if(glWidget->stateBcColors)
    glWidget->stateBodyColors = false;

  glWidget->rebuildLists();
  synchronizeMenuToState();
}

// View -> Colors -> Bodies
//-----------------------------------------------------------------------------
void MainWindow::colorizeBodySlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("No mesh - unable to colorize bodies");
    return;
  }

  glWidget->stateBodyColors = !glWidget->stateBodyColors;

  if(glWidget->stateBodyColors)
    glWidget->stateBcColors = false;

  glWidget->rebuildLists();
  synchronizeMenuToState();
}


// View -> Colors -> Background
//-----------------------------------------------------------------------------
void MainWindow::backgroundColorSlot()
{
  QColor newColor = QColorDialog::getColor(glWidget->backgroundColor, this);
  glWidget->qglClearColor(newColor);
  glWidget->backgroundColor = newColor;
}



// View -> Colors -> Surface
//-----------------------------------------------------------------------------
void MainWindow::surfaceColorSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to change surface color when the mesh is empty");
    return;
  }

  QColor newColor = QColorDialog::getColor(glWidget->surfaceColor, this);
  glWidget->surfaceColor = newColor;
  glWidget->rebuildLists();
}


// View -> Colors -> Edge
//-----------------------------------------------------------------------------
void MainWindow::edgeColorSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to change edge color when the mesh is empty");
    return;
  }

  QColor newColor = QColorDialog::getColor(glWidget->edgeColor, this);
  glWidget->edgeColor = newColor;
  glWidget->rebuildLists();
}


// View -> Colors -> Surface mesh
//-----------------------------------------------------------------------------
void MainWindow::surfaceMeshColorSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to change surface mesh color when the mesh is empty");
    return;
  }

  QColor newColor = QColorDialog::getColor(glWidget->surfaceMeshColor, this);
  glWidget->surfaceMeshColor = newColor;
  glWidget->rebuildLists();
}


// View -> Colors -> Sharp edges
//-----------------------------------------------------------------------------
void MainWindow::sharpEdgeColorSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("Unable to change sharp edge colors when the mesh is empty");
    return;
  }

  QColor newColor = QColorDialog::getColor(glWidget->sharpEdgeColor, this);
  glWidget->sharpEdgeColor = newColor;
  glWidget->rebuildLists();
}


// View -> Cad model...
//-----------------------------------------------------------------------------
void MainWindow::showCadModelSlot()
{
#ifdef EG_OCC
  cadView->show();
#endif
}


// View -> Twod model...
//-----------------------------------------------------------------------------
void MainWindow::showTwodViewSlot()
{
  twodView->show();
}


// View -> VTK post...
//-----------------------------------------------------------------------------
void MainWindow::showVtkPostSlot()
{
#ifdef EG_VTK
  QString postFileName = saveDirName + "/"  + generalSetup->ui.postFileEdit->text().trimmed();

  // Parallel solution:
  //====================
  Ui::parallelDialog ui = parallel->ui;
  bool parallelActive = ui.parallelActiveCheckBox->isChecked();

  if(parallelActive) {
    
    // unify mesh:
    if(meshUnifier->state() == QProcess::Running) {
      logMessage("Mesh unifier is already running - aborted");
      return;
    }
    
    if(saveDirName.isEmpty()) {
      logMessage("saveDirName is empty - unable to locate result files");
      return;
    }

    // Set up log window:
    solverLogWindow->setWindowTitle(tr("ElmerGrid log"));
    solverLogWindow->getTextEdit()->clear();
    solverLogWindow->setFound(false);
    solverLogWindow->show();

    QString postName = generalSetup->ui.postFileEdit->text().trimmed();
    QStringList postNameSplitted = postName.split(".");
    int nofProcessors = ui.nofProcessorsSpinBox->value();

    QString unifyingCommand = ui.mergeLineEdit->text().trimmed();
    unifyingCommand.replace(QString("%ep"), postNameSplitted.at(0).trimmed());
    unifyingCommand.replace(QString("%n"), QString::number(nofProcessors));
    
    logMessage("Executing: " + unifyingCommand);
    
    meshUnifier->start(unifyingCommand);
    
    if(!meshUnifier->waitForStarted()) {
      solverLogWindow->getTextEdit()->append("Unable to start ElmerGrid for mesh unification - aborted");
      logMessage("Unable to start ElmerGrid for mesh unification - aborted");
      vtkPostMeshUnifierRunning = false;
      return;
    }
    
    // The rest is done in meshUnifierFinishedSlot:
    vtkPostMeshUnifierRunning = true;
    return;
  }

  // Scalar solution:
  //-----------------
  vtkPost->show();

  vtkPost->ReadPostFile(postFileName);
  // if(!vtkPost->ReadPostFile(postFileName)) vtkPost->readEpFileSlot();
#endif
}


// View -> Paraview
//-----------------------------------------------------------------------------
void MainWindow::showParaViewSlot()
{
#ifdef EG_PARAVIEW
  QString postFileName = generalSetup->ui.postFileEdit->text().trimmed();
  QString pvFileName;
  QFileInfo pvFile(postFileName);

  Ui::parallelDialog ui = parallel->ui;
  bool parallelActive = ui.parallelActiveCheckBox->isChecked();

  // Serial solution
  //================
  if(!parallelActive) {
    pvFileName = pvFile.baseName() + "????.vtu";
  }

  // Parallel solution
  //==================
  if(parallelActive) {
     pvFileName = pvFile.baseName() + "????.pvtu";
  }

  // Launch ParaView
  //================
  QDir currentDir;

  currentDir = QDir(saveDirName);
  QStringList pvFiles;

  pvFiles = currentDir.entryList(QStringList(pvFileName), QDir::Files | QDir::NoSymLinks);

  post->start("paraview", pvFiles);
#endif
}


//*****************************************************************************
//
//                                Mesh MENU
//
//*****************************************************************************


// Mesh -> Control...
//-----------------------------------------------------------------------------
void MainWindow::meshcontrolSlot()
{
  meshControl->tetlibPresent = this->tetlibPresent;
  meshControl->nglibPresent = this->nglibPresent;

  if(!tetlibPresent) {
    meshControl->tetlibPresent = false;
    meshControl->ui.nglibRadioButton->setChecked(true);
    meshControl->ui.tetlibRadioButton->setEnabled(false);
    meshControl->ui.tetlibStringEdit->setEnabled(false);
  }

  if(!nglibPresent) {
    meshControl->nglibPresent = false;
    meshControl->ui.tetlibRadioButton->setChecked(true);
    meshControl->ui.nglibRadioButton->setEnabled(false);
    meshControl->ui.nglibMaxHEdit->setEnabled(false);
    meshControl->ui.nglibFinenessEdit->setEnabled(false);
    meshControl->ui.nglibBgmeshEdit->setEnabled(false);
  }

  if(!tetlibPresent && !nglibPresent) 
    meshControl->ui.elmerGridRadioButton->setChecked(true);  

  meshControl->show();
}



// Mesh -> Remesh
//-----------------------------------------------------------------------------
void MainWindow::remeshSlot()
{
  if(activeGenerator == GEN_UNKNOWN) {
    logMessage("Unable to (re)mesh: no input data or mesh generator (please make sure that your input file suffix is in lower case)");
    return;
  }  
    
  // ***** ELMERGRID *****

  if(activeGenerator == GEN_ELMERGRID) {

    meshutils->clearMesh(glWidget->getMesh());
    glWidget->newMesh();
    mesh_t *mesh = glWidget->getMesh();
    
#if WITH_QT5
    elmergridAPI->createElmerMeshStructure(mesh, meshControl->elmerGridControlString.toLatin1());
#else
    elmergridAPI->createElmerMeshStructure(mesh, meshControl->elmerGridControlString.toAscii());
#endif
    
    if(mesh->getSurfaces() == 0) meshutils->findSurfaceElements(mesh);
    
    for(int i = 0; i < mesh->getSurfaces(); i++ ) {
      surface_t *surface = mesh->getSurface(i);
      
      surface->setEdges((int)(surface->getCode() / 100));
      surface->newEdgeIndexes(surface->getEdges());
      for(int j=0; j<surface->getEdges(); j++)
	surface->setEdgeIndex(j, -1);
    }

    meshutils->findSurfaceElementEdges(mesh);
    meshutils->findSurfaceElementNormals(mesh);

    glWidget->rebuildLists();
    applyOperations();
    
    return;
  }
  
  // ***** Threaded generators *****

  if(!remeshAct->isEnabled()) {
    logMessage("Meshing thread is already running - aborting");
    return;
  }

  if(activeGenerator == GEN_TETLIB) {

    if(!tetlibPresent) {
      logMessage("tetlib functionality unavailable");
      return;
    }
    
    if(!tetlibInputOk) {
      logMessage("Remesh: error: no input data for tetlib");
      return;
    }

    // Usually "J" should be included in the control string: 
    tetlibControlString = meshControl->tetlibControlString;

  } else if(activeGenerator == GEN_NGLIB) {

    if(!nglibPresent) {
      logMessage("nglib functionality unavailable");
      return;
    }

    if(!nglibInputOk) {
      logMessage("Remesh: error: no input data for nglib");
      return;
    }

    // Init & set mesh params.:
    //--------------------------    
    cout << "Initializing nglib" << endl;
    nglib::Ng_Init();

    char backgroundmesh[1024];
#if WITH_QT5
    sprintf(backgroundmesh, "%s", meshControl->nglibBackgroundmesh.toLatin1().data());
#else
    sprintf(backgroundmesh, "%s", meshControl->nglibBackgroundmesh.toAscii().data());
#endif

    mp.maxh = meshControl->nglibMaxH.toDouble();
    mp.fineness = meshControl->nglibFineness.toDouble();
    mp.secondorder = 0;
    mp.meshsize_filename = backgroundmesh;

    if(ngDim == 3) {

      // STL (3D):
      //-----------
      cout << "Start meshing..." << endl;

      nggeom = nglib::Ng_STL_NewGeometry();

      ngmesh = nglib::Ng_NewMesh();
      
      if(!occInputOk) {
	
	// STL: regenerate structures for nglib:
	//--------------------------------------
#if WITH_QT5
	nggeom = nglib::Ng_STL_LoadGeometry(stlFileName.toLatin1().data(), 0);
#else
	nggeom = nglib::Ng_STL_LoadGeometry(stlFileName.toAscii().data(), 0);
#endif
	
	if(!nggeom) {
	  logMessage("Ng_STL_LoadGeometry failed");
	  return;
	}
	
	nglib::Ng_STL_InitSTLGeometry(nggeom);

	nglib::Ng_STL_MakeEdges(nggeom, ngmesh, &mp);
	
	double maxMeshSize = mp.maxh;

	if(maxMeshSize <= 0) maxMeshSize = 10000000;

	nglib::Ng_RestrictMeshSizeGlobal(ngmesh, maxMeshSize);      
	
#ifdef EG_OCC
      } else {
	
	// OCC: (re)generate STL for nglib:
	//----------------------------------
	cadView->setMesh(ngmesh);
	cadView->setGeom(nggeom);
	cadView->setMp(&mp);
	cadView->generateSTL();
#endif
      }
      
    } else if(ngDim == 2) {

      // IN2D (2D):
      //------------
      cout << "Start 2D meshing..." << endl;

      if(!occInputOk) {
	
	// Native 2D geometry input for Ng:
	//----------------------------------
	if(in2dFileName.isEmpty()) {
	  logMessage("File name is empty - aborting");
	  return;
	}
	
	ngmesh = nglib::Ng_NewMesh();
	
#if WITH_QT5
	nggeom2d = nglib::Ng_LoadGeometry_2D(in2dFileName.toLatin1().data());
#else
	nggeom2d = nglib::Ng_LoadGeometry_2D(in2dFileName.toAscii().data());
#endif
	
	if(!nggeom2d) {
	  logMessage("Ng_LoadGeometry_2D failed");
	  return;
	}
	
	nglibAPI->setNggeom2D(nggeom2d);
	
	double maxMeshSize = mp.maxh;
	
	if(maxMeshSize <= 0) maxMeshSize = 10000000;
	
	nglib::Ng_RestrictMeshSizeGlobal(ngmesh, maxMeshSize);

#ifdef EG_OCC
      } else {

	// Model originates from a 2D cad file:
	//--------------------------------------
	cadView->generateIn2dFile();

	ngmesh = nglib::Ng_NewMesh();
	
	nggeom2d = nglib::Ng_LoadGeometry_2D("iges2ng.in2d");
	
	if(!nggeom2d) {
	  logMessage("Ng_LoadGeometry_2D failed");
	  return;
	}
	
	nglibAPI->setNggeom2D(nggeom2d);
	
	double maxMeshSize = mp.maxh;
	
	if(maxMeshSize <= 0) maxMeshSize = 10000000;
	
	nglib::Ng_RestrictMeshSizeGlobal(ngmesh, maxMeshSize);
#endif
      }

    } else {
      
      // Unknown spatial dimension:
      //----------------------------
      cout << "Unknown spatial dimension" << endl;
      return;

    }
    
  } else {

    logMessage("Remesh: uknown generator type");
    return;

  }

  // ***** Start meshing thread *****

  logMessage("Sending start request to mesh generator...");

  meshutils->clearMesh(glWidget->getMesh());
  glWidget->newMesh();
  mesh_t *mesh = glWidget->getMesh();

  // Re-enable when finished() or terminated() signal is received:
  remeshAct->setEnabled(false);
  stopMeshingAct->setEnabled(true);

  if(activeGenerator == GEN_NGLIB) 
    stopMeshingAct->setEnabled(false);

  meshingThread->generate(activeGenerator, tetlibControlString,
			  tetlibAPI, ngmesh, nggeom, nggeom2d,
			  ngDim, &mp);
}



// Mesh -> Kill generator
//-----------------------------------------------------------------------------
void MainWindow::stopMeshingSlot()
{
  if(remeshAct->isEnabled()) {
    logMessage("Mesh generator is not running");
    return;
  }

  logMessage("Sending termination request to mesh generator...");
  meshingThread->stopMeshing();
}



// Meshing has started (signaled by meshingThread):
//-----------------------------------------------------------------------------
void MainWindow::meshingStartedSlot()
{
  logMessage("Mesh generator started");

  updateSysTrayIcon("Mesh generator started",
		    "Use Mesh->Terminate to stop processing");

  statusBar()->showMessage(tr("Mesh generator started"));
  
  progressBar->show();
  progressBar->setRange(0, 0);

  progressLabel->show();
  progressLabel->setText("Meshing");
}


// Meshing has been terminated (signaled by meshingThread):
//-----------------------------------------------------------------------------
void MainWindow::meshingTerminatedSlot()
{
  logMessage("Mesh generator terminated");

  progressBar->hide();
  progressBar->setRange(0, 100);

  progressLabel->hide();

  stopMeshingAct->setEnabled(true);

  updateSysTrayIcon("Mesh generator terminated",
		    "Use Mesh->Remesh to restart");

  statusBar()->showMessage(tr("Ready"));
  
  // clean up:
  if(activeGenerator == GEN_TETLIB) {
    cout << "Cleaning up...";
    out->deinitialize();
    cout << "done" << endl;
    cout.flush();
  }

  if(activeGenerator == GEN_NGLIB) {
    nglib::Ng_DeleteMesh(ngmesh);
    nglib::Ng_Exit();
  }

  remeshAct->setEnabled(true);
  stopMeshingAct->setEnabled(false);
}

// Mesh is ready (signaled by meshingThread):
//-----------------------------------------------------------------------------
void MainWindow::meshingFinishedSlot()
{
  logMessage("Mesh generation ready");

  progressBar->hide();
  progressBar->setRange(0, 100);

  progressLabel->hide();

  if(activeGenerator == GEN_TETLIB) {

    makeElmerMeshFromTetlib();

  } else if(activeGenerator == GEN_NGLIB) {

    this->ngmesh = meshingThread->getNgMesh();

    makeElmerMeshFromNglib();

    nglib::Ng_DeleteMesh(ngmesh);
    nglib::Ng_Exit();

  } else {
    
    logMessage("MeshOk: error: unknown mesh generator");

  }

  applyOperations();

  statusBar()->showMessage(tr("Ready"));
  
  updateSysTrayIcon("Mesh generator has finished",
		    "Select Model->Summary for statistics");

  remeshAct->setEnabled(true);
  stopMeshingAct->setEnabled(false);

  // Check cmd line arguments:
  //---------------------------
  QStringList args = QCoreApplication::arguments();

  int output = args.indexOf("-o");

  if(output > 0) {
    QString dirName = args.at(output + 1);
    
    if(dirName.left(1) != "-") {
      cout << "Saving mesh files" << endl;
      saveElmerMesh(dirName);
    }
  }

  if(args.contains("-e") || args.contains("-nogui")) {
    cout << "Exiting" << endl;
    QApplication::closeAllWindows();
    exit(0);
  }

  resetSlot();
}


// Mesh -> Divide surface...
//-----------------------------------------------------------------------------
void MainWindow::surfaceDivideSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("There is nothing to divide - mesh is empty");
    return;
  }

  boundaryDivide->target = TARGET_SURFACES;
  boundaryDivide->show();
}



// Make surface division by sharp edges (signalled by boundaryDivide)...
//-----------------------------------------------------------------------------
void MainWindow::doDivideSurfaceSlot(double angle)
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("No mesh to divide");
    return;
  }
  
  operations++;
  operation_t *p = new operation_t;
  operation_t *q = NULL;

  for( q=&operation; q->next; q=q->next );
  q->next = p;
  p->next = NULL;
  
  p->type = OP_DIVIDE_SURFACE;
  p->angle = angle;

  int selected=0;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);

    if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY))
      selected++;
  }
  p->selected = selected;
  p->select_set = new int[selected];
  selected = 0;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY))
      p->select_set[selected++] = i;
  }
  

  meshutils->findSharpEdges(mesh, angle);
  int parts = meshutils->divideSurfaceBySharpEdges(mesh);

  QString qs = "Surface divided into " + QString::number(parts) + " parts";
  statusBar()->showMessage(qs);
  
  synchronizeMenuToState();
  glWidget->rebuildLists();
  glWidget->updateGL();

  // Added 05 September 2009
  boundaryPropertyEditor.clear();
  boundaryPropertyEditor.resize(parts);
  for(int i = 0; i < parts; i++)
    boundaryPropertyEditor[i] = new BoundaryPropertyEditor;
  
}



// Mesh -> Unify surface
//-----------------------------------------------------------------------------
void MainWindow::surfaceUnifySlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("No surfaces to unify");
    return;
  }
  
  int targetindex = -1, selected=0;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY)) {
      selected++;
      if(targetindex < 0) targetindex = l->getIndex();
    }
  }
  
  if(targetindex < 0) {
    logMessage("No surfaces selected");
    return;
  }


  operations++;
  operation_t *p = new operation_t, *q;
  for( q=&operation; q->next; q=q->next );
  q->next = p;
  p->next = NULL;
  p->type = OP_UNIFY_SURFACE;
  p->selected = selected;
  p->select_set = new int[selected]; 
  
  selected = 0;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY)) {
      p->select_set[selected++] = i;
      for(int j=0; j < mesh->getSurfaces(); j++) {
	surface_t *s = mesh->getSurface(j);
	if((s->getIndex() == l->getIndex()) && (s->getNature() == PDE_BOUNDARY)) 
	  s->setIndex(targetindex);
      }
    }
  }
  
  cout << "Selected surfaces marked with index " << targetindex << endl;
  cout.flush();

  glWidget->rebuildLists();

  logMessage("Selected surfaces unified");
}


void MainWindow::applyOperations()
{
  mesh_t *mesh = glWidget->getMesh();

  cout << "Apply " << operations << " operations" << endl;
  cout.flush();

  operation_t *p = operation.next;
  for( ; p; p=p->next )
  {
    int lists = glWidget->getLists();

    for( int i=0; i<lists; i++ )
      glWidget->getList(i)->setSelected(false);

    for( int j=0; j<mesh->getSurfaces(); j++ )
      mesh->getSurface(j)->setSelected(false);

    for( int j=0; j<mesh->getEdges(); j++ )
      mesh->getEdge(j)->setSelected(false);

    for(int i=0; i < p->selected; i++) {
      list_t *l = glWidget->getList(p->select_set[i]);

      l->setSelected(true);
      if ( p->type < OP_UNIFY_EDGE ) {
        for( int j=0; j<mesh->getSurfaces(); j++ ) {
          surface_t *surf = mesh->getSurface(j);
          if( l->getIndex() == surf->getIndex() )
            surf->setSelected(l->isSelected());
        }
      } else {
        for( int j=0; j<mesh->getEdges(); j++ ) {
          edge_t *edge = mesh->getEdge(j);
          if( l->getIndex() == edge->getIndex() )
            edge->setSelected(l->isSelected());
        }
      }
    }

    if ( p->type == OP_DIVIDE_SURFACE ) {
      meshutils->findSharpEdges(mesh, p->angle);
      int parts = meshutils->divideSurfaceBySharpEdges(mesh);
      QString qs = "Surface divided into " + QString::number(parts) + " parts";
      statusBar()->showMessage(qs);

    } else if ( p->type == OP_DIVIDE_EDGE ) {
      meshutils->findEdgeElementPoints(mesh);
      meshutils->findSharpPoints(mesh, p->angle);
      int parts = meshutils->divideEdgeBySharpPoints(mesh);
      QString qs = "Edges divided into " + QString::number(parts) + " parts";
      statusBar()->showMessage(qs);

    } else if (p->type == OP_UNIFY_SURFACE ) {
      int targetindex = -1;

      for(int i=0; i<lists; i++) {
        list_t *l = glWidget->getList(i);
        if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY)) {
          if(targetindex < 0) {
            targetindex = l->getIndex();
            break;
          }
        }
      }
      for(int i=0; i<lists; i++) {
        list_t *l = glWidget->getList(i);
        if(l->isSelected() && (l->getType() == SURFACELIST) && (l->getNature() == PDE_BOUNDARY)) {
          for(int j=0; j < mesh->getSurfaces(); j++) {
            surface_t *s = mesh->getSurface(j);
            if((s->getIndex() == l->getIndex()) && (s->getNature() == PDE_BOUNDARY)) 
              s->setIndex(targetindex);
          }
        }
      }
      cout << "Selected surfaces marked with index " << targetindex << endl;
      cout.flush();

    } else if (p->type == OP_UNIFY_EDGE ) {
      int targetindex = -1;
      for(int i=0; i<lists; i++) {
        list_t *l = glWidget->getList(i);
        if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY) {
          if(targetindex < 0) {
            targetindex = l->getIndex();
            break;
          }
        }
      }
      for(int i=0; i<lists; i++) {
        list_t *l = glWidget->getList(i);
        if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY) {
          for(int j=0; j < mesh->getEdges(); j++) {
            edge_t *e = mesh->getEdge(j);
            if(e->getIndex() == l->getIndex() && e->getNature() == PDE_BOUNDARY)
              e->setIndex(targetindex);
          }
        }
      }
      cout << "Selected edges marked with index " << targetindex << endl;
      cout.flush();
    }
    glWidget->rebuildLists();
  }
  
  synchronizeMenuToState();
  glWidget->updateGL();
  
  // Added 05 September 2009
  boundaryPropertyEditor.clear();
  int parts = glWidget->getLists();
  boundaryPropertyEditor.resize(parts);
  for(int i = 0; i < parts; i++)
    boundaryPropertyEditor[i] = new BoundaryPropertyEditor;
}



// Mesh -> Divide edge...
//-----------------------------------------------------------------------------
void MainWindow::edgeDivideSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("There is nothing to divide - mesh is empty");
    return;
  }

  boundaryDivide->target = TARGET_EDGES;
  boundaryDivide->show();
}



// Make edge division by sharp points (signalled by boundaryDivide)...
//-----------------------------------------------------------------------------
void MainWindow::doDivideEdgeSlot(double angle)
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("No mesh to divide");
    return;
  }

  operations++;
  operation_t *p = new operation_t, *q;
  for( q=&operation; q->next; q=q->next );
  q->next = p;
  p->next = NULL;

  p->type = OP_DIVIDE_EDGE;
  p->angle = angle;

  int selected=0;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY)
      selected++;
  }
  p->selected = selected;
  p->select_set = new int[selected];
  selected = 0;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY)
      p->select_set[selected++] = i;
  }
  

  meshutils->findEdgeElementPoints(mesh);
  meshutils->findSharpPoints(mesh, angle);
  int parts = meshutils->divideEdgeBySharpPoints(mesh);
  
  QString qs = "Edge divided into " + QString::number(parts) + " parts";
  statusBar()->showMessage(qs);

  synchronizeMenuToState();
  glWidget->rebuildLists();
  glWidget->updateGL();

  // Added 05 September 2009
  boundaryPropertyEditor.clear();
  boundaryPropertyEditor.resize(parts);
  for(int i = 0; i < parts; i++)
    boundaryPropertyEditor[i] = new BoundaryPropertyEditor;

}



// Mesh -> Unify edge
//-----------------------------------------------------------------------------
void MainWindow::edgeUnifySlot()
{
  mesh_t *mesh = glWidget->getMesh();
  int lists = glWidget->getLists();

  if(mesh == NULL) {
    logMessage("No edges to unify");
    return;
  }
  
  int targetindex = -1, selected=0;
  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY) {
      selected++;
      if(targetindex < 0) targetindex = l->getIndex();
    }
  }
  

  if(targetindex < 0) {
    logMessage("No edges selected");
    return;
  }

  operations++;
  operation_t *p = new operation_t, *q;
  for( q=&operation; q->next; q=q->next );
  q->next = p;
  p->next = NULL;
  p->type = OP_UNIFY_EDGE;
  p->selected = selected;
  p->select_set = new int[selected]; 
  
  selected = 0;

  for(int i=0; i<lists; i++) {
    list_t *l = glWidget->getList(i);
    if(l->isSelected() && l->getType() == EDGELIST && l->getNature() == PDE_BOUNDARY) {
      p->select_set[selected++] = i;
      for(int j=0; j < mesh->getEdges(); j++) {
	edge_t *e = mesh->getEdge(j);
	if(e->getIndex() == l->getIndex() && e->getNature() == PDE_BOUNDARY) 
	  e->setIndex(targetindex);
      }
    }
  }
  
  cout << "Selected edges marked with index " << targetindex << endl;
  cout.flush();

  glWidget->rebuildLists();

  logMessage("Selected edges unified");
}


// Mesh -> Clean up
//-----------------------------------------------------------------------------
void MainWindow::cleanHangingSharpEdgesSlot()
{
  mesh_t *mesh = glWidget->getMesh();

  if(mesh == NULL)
    return;

  int count = meshutils->cleanHangingSharpEdges(mesh);

  cout << "Removed " << count << " hanging sharp edges" << endl;
  cout.flush();

  glWidget->rebuildLists();
}


//*****************************************************************************
//
//                                Edit MENU
//
//*****************************************************************************


// Edit -> Sif...
//-----------------------------------------------------------------------------
void MainWindow::showsifSlot()
{
  // QFont sansFont("Courier", 10);
  // sifWindow->getTextEdit()->setCurrentFont(sansFont);
  sifWindow->show();
}


// Edit -> Generate sif
//-----------------------------------------------------------------------------
void MainWindow::generateSifSlot()
{
  mesh_t *mesh = glWidget->getMesh();

  if(mesh == NULL) {
    logMessage("Unable to create SIF: no mesh");
    return;
  }
  
  if((mesh->getDim() < 1) || (mesh->getCdim() < 1)) {
    logMessage("Model dimension inconsistent with SIF syntax");
    return;
  }

  // Clear SIF text editor:
  //------------------------
  sifWindow->getTextEdit()->clear();
  sifWindow->setFirstTime(true);
  sifWindow->setFound(false);
  // QFont sansFont("Courier", 10);
  // sifWindow->getTextEdit()->setCurrentFont(sansFont);

  // Set up SIF generator:
  //-----------------------
  sifGenerator->setMesh(mesh);
  sifGenerator->setTextEdit(sifWindow->getTextEdit());
  sifGenerator->setDim(mesh->getDim());
  sifGenerator->setCdim(mesh->getCdim());
  sifGenerator->setGeneralSetup(generalSetup);

  sifGenerator->setEquationEditor(equationEditor);
  sifGenerator->setMaterialEditor(materialEditor);
  sifGenerator->setBodyForceEditor(bodyForceEditor);
  sifGenerator->setInitialConditionEditor(initialConditionEditor);
  sifGenerator->setBoundaryConditionEditor(boundaryConditionEditor);
  sifGenerator->setSolverParameterEditor(solverParameterEditor);
  sifGenerator->setBoundaryPropertyEditor(boundaryPropertyEditor);
  sifGenerator->setBodyPropertyEditor(bodyPropertyEditor);
  sifGenerator->setMeshControl(meshControl);
  sifGenerator->setElmerDefs(elmerDefs);
  sifGenerator->bodyMap = glWidget->bodyMap;
  sifGenerator->boundaryMap = glWidget->boundaryMap;

  // Make SIF:
  //----------
  sifGenerator->makeHeaderBlock();
  sifGenerator->makeSimulationBlock();
  sifGenerator->makeConstantsBlock();
  sifGenerator->makeBodyBlocks();
  sifGenerator->makeEquationBlocks();
  sifGenerator->makeMaterialBlocks();
  sifGenerator->makeBodyForceBlocks();
  sifGenerator->makeInitialConditionBlocks();
  sifGenerator->makeBoundaryBlocks();
}


// Boundary selected by double clicking (signaled by glWidget::select):
//-----------------------------------------------------------------------------
void MainWindow::boundarySelectedSlot(list_t *l)
{
  QString qs;

  if(l->getIndex() < 0) {
    statusBar()->showMessage("Ready");    
    return;
  }

  if(l->isSelected()) {
    if(l->getType() == SURFACELIST) {
      qs = "Selected surface " + QString::number(l->getIndex());
    } else if(l->getType() == EDGELIST) {
      qs = "Selected edge " + QString::number(l->getIndex());
    } else {
      qs = "Selected object " + QString::number(l->getIndex()) + " (type unknown)";
    }
  } else {
    if(l->getType() == SURFACELIST) {
      qs = "Unselected surface " + QString::number(l->getIndex());
    } else if(l->getType() == EDGELIST) {
      qs = "Unselected edge " + QString::number(l->getIndex());
    } else {
      qs = "Unselected object " + QString::number(l->getIndex()) + " (type unknown)";
    }
  }

  logMessage(qs);
  

  // Open bc property sheet for selected boundary:
  //-----------------------------------------------
  if(l->isSelected() && (glWidget->altPressed || bcEditActive)) {
    glWidget->ctrlPressed = false;
    glWidget->shiftPressed = false;
    glWidget->altPressed = false;

    // renumbering:
    int n = glWidget->boundaryMap.value(l->getIndex());

    if(n >= boundaryPropertyEditor.size()) {
      logMessage("Error: Boundary index mismatch");
      return;
    }
    
    BoundaryPropertyEditor *boundaryEdit = boundaryPropertyEditor[n];
    populateBoundaryComboBoxes(boundaryEdit);
    
    connect( boundaryEdit, SIGNAL(BoundaryAsABodyChanged(BoundaryPropertyEditor *,int)),
	     this, SLOT(boundaryAsABodyChanged(BoundaryPropertyEditor *,int)) );
    
    if(boundaryEdit->touched) {
      boundaryEdit->ui.applyButton->setText("Update");
      // boundaryEdit->ui.discardButton->setText("Remove");
      boundaryEdit->ui.discardButton->setText("Cancel");
      boundaryEdit->ui.applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      // boundaryEdit->ui.discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      boundaryEdit->ui.discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
    } else {
      boundaryEdit->ui.applyButton->setText("Add");
      boundaryEdit->ui.discardButton->setText("Cancel");
      boundaryEdit->ui.applyButton->setIcon(QIcon(":/icons/list-add.png"));
      boundaryEdit->ui.discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
    }

    boundaryEdit->setWindowTitle("Properties for boundary " + QString::number(l->getIndex()));
    boundaryEdit->show();
  }

  BodyPropertyEditor *bodyEdit = NULL;
  int current = -1, n =-1;

  // boundary as a body treatment
  // ----------------------------
  if(l->isSelected() && glWidget->ctrlPressed ) {

    // renumbering:
    int n = glWidget->boundaryMap.value(l->getIndex());

    if(n >= boundaryPropertyEditor.size()) {
      logMessage("Error: Boundary index mismatch");
      return;
    }

    BoundaryPropertyEditor *boundaryEdit = boundaryPropertyEditor[n];

    if(!boundaryEdit) {
      cout << "MainWindow: Boundary index out of bounds" << endl;
      return;
    }
      

    bodyEdit = boundaryEdit->bodyProperties;

    if ( bodyEdit ) {
      glWidget->ctrlPressed = false;
      glWidget->shiftPressed = false;
      glWidget->altPressed = false;

      bodyEdit->setWindowTitle("Properties for body " + QString::number(current));

      // if(bodyEdit->ui.nameEdit->text().trimmed().isEmpty())
      bodyEdit->ui.nameEdit->setText("Body Property{Boundary " + QString::number(n+1) +  "}");
    }
  }

  // Open body property sheet for selected body:
  //---------------------------------------------
  if( (glWidget->currentlySelectedBody >= 0) &&
      (glWidget->shiftPressed || bodyEditActive) ) {
    
    glWidget->ctrlPressed = false;
    glWidget->shiftPressed = false;
    glWidget->altPressed = false;
    
    current = glWidget->currentlySelectedBody;
    
    cout << "Current selection uniquely determines body: " << current << endl;
    cout.flush();
 
    // renumbering:
    n = glWidget->bodyMap.value(current);

    if(n >= bodyPropertyEditor.size()) {
      logMessage("MainWindow: Body index out of bounds)");
      return;
    }
     
    bodyEdit = bodyPropertyEditor[n];

    if(!bodyEdit)
      cout << "MainWindow: Undetermined body index" << endl;

    bodyEdit->setWindowTitle("Properties for body " + QString::number(current));

    if(bodyEdit->ui.nameEdit->text().trimmed().isEmpty())
      bodyEdit->ui.nameEdit->setText("Body Property " + QString::number(n+1));

  }

  if ( bodyEdit ) {
    populateBodyComboBoxes(bodyEdit);
    
    if(bodyEdit->touched) {
      bodyEdit->ui.applyButton->setText("Update");
      // bodyEdit->ui.discardButton->setText("Remove");
      bodyEdit->ui.discardButton->setText("Cancel");
      bodyEdit->ui.applyButton->setIcon(QIcon(":/icons/dialog-ok-apply.png"));
      // bodyEdit->ui.discardButton->setIcon(QIcon(":/icons/list-remove.png"));
      bodyEdit->ui.discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
    } else {
      bodyEdit->ui.applyButton->setText("Add");
      bodyEdit->ui.discardButton->setText("Cancel");
      bodyEdit->ui.applyButton->setIcon(QIcon(":/icons/list-add.png"));
      bodyEdit->ui.discardButton->setIcon(QIcon(":/icons/dialog-close.png"));
    }

    bodyEdit->show();
  }
}


// Populate boundary editor's comboboxes:
//---------------------------------------
void MainWindow::populateBoundaryComboBoxes(BoundaryPropertyEditor *boundary)
{
  boundary->disconnect(SIGNAL(BoundaryComboChanged(BoundaryPropertyEditor *,QString)));
  while(boundary && boundary->ui.boundaryConditionCombo && boundary->ui.boundaryConditionCombo->count() > 0) 
    boundary->ui.boundaryConditionCombo->removeItem(0);
    
  int takethis = 1; //-1;
  int count = 1;
  boundary->ui.boundaryConditionCombo->insertItem(count++, "");

  for(int i = 0; i < boundaryConditionEditor.size(); i++) {
    DynamicEditor *bcEdit = boundaryConditionEditor[i];
    if(bcEdit->menuAction != NULL) {
      const QString &name = bcEdit->nameEdit->text().trimmed();
      boundary->ui.boundaryConditionCombo->insertItem(count, name);
      if ( boundary->condition == bcEdit ) takethis = count;
      count++;
    }
  }
  connect( boundary,SIGNAL(BoundaryComboChanged(BoundaryPropertyEditor *,QString)), 
        this, SLOT(boundaryComboChanged(BoundaryPropertyEditor *,QString)) );

  boundary->ui.boundaryConditionCombo->setCurrentIndex(takethis-1);
}

//-----------------------------------------------------------------------------
void MainWindow::boundaryComboChanged(BoundaryPropertyEditor *b, QString text)
{
  b->condition = 0;
  b->touched = false;

  for( int i=0; i < boundaryConditionEditor.size(); i++ )
  {
    DynamicEditor *bc = boundaryConditionEditor[i];
    if ( bc->ID >= 0 ) {
       if ( bc->nameEdit->text().trimmed() == text ) {
         b->condition = bc; 
         b->touched = true;
         break;
       }
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::boundaryAsABodyChanged(BoundaryPropertyEditor *b, int status)
{
  int indx=glWidget->bodyMap.count();

  if ( status ) {
    for( int i = 0; i < boundaryPropertyEditor.size(); i++ )
      if ( boundaryPropertyEditor[i]
	   && boundaryPropertyEditor[i]->bodyProperties ) indx++;

    if(indx >= bodyPropertyEditor.size()) {
      cout << "MainWindow: Body index out of bounds" << endl;
      return;
    }

    if(bodyPropertyEditor[indx])
      b->bodyProperties = bodyPropertyEditor[indx];

  } else {

    b->bodyProperties = NULL;
  }
}

    
// Populate body editor's comboboxes:
//-----------------------------------
void MainWindow::populateBodyComboBoxes(BodyPropertyEditor *bodyEdit)
{
  // Equation:
  // =========
  bodyEdit->disconnect(SIGNAL(BodyEquationComboChanged(BodyPropertyEditor *,QString)));
  bodyEdit->ui.equationCombo->clear();

  int count = 1;
  int takethis = 1; // -1
  bodyEdit->ui.equationCombo->insertItem(count++, "");

  for(int i = 0; i < equationEditor.size(); i++) {
    DynamicEditor *eqEdit = equationEditor[i];
    if(eqEdit->menuAction != NULL) {
      const QString &name = eqEdit->nameEdit->text().trimmed();
      bodyEdit->ui.equationCombo->insertItem(count, name);
      if ( bodyEdit->equation == eqEdit ) takethis = count;
      count++;
    }
  }
  connect( bodyEdit,SIGNAL(BodyEquationComboChanged(BodyPropertyEditor *,QString)), 
        this, SLOT(equationComboChanged(BodyPropertyEditor *,QString)) );
  bodyEdit->ui.equationCombo->setCurrentIndex(takethis-1);
    
  // Material
  // =========
  bodyEdit->disconnect(SIGNAL(BodyMaterialComboChanged(BodyPropertyEditor *,QString)));
  bodyEdit->ui.materialCombo->clear();

  count = 1;
  takethis = 1;
  bodyEdit->ui.materialCombo->insertItem(count, "");
  count++;

  for(int i = 0; i < materialEditor.size(); i++) {
    DynamicEditor *matEdit = materialEditor[i];

    if(matEdit->menuAction != NULL) {
      const QString &name = matEdit->nameEdit->text().trimmed();
      bodyEdit->ui.materialCombo->insertItem(count, name);
      if ( bodyEdit->material==matEdit ) takethis = count;
      count++;
    }
  }

  connect( bodyEdit,SIGNAL(BodyMaterialComboChanged(BodyPropertyEditor *,QString)), 
        this, SLOT(materialComboChanged(BodyPropertyEditor *,QString)) );


  bodyEdit->ui.materialCombo->setCurrentIndex(takethis-1);


  // Bodyforce:
  //===========
  bodyEdit->disconnect(SIGNAL(BodyForceComboChanged(BodyPropertyEditor *,QString)));
  bodyEdit->ui.bodyForceCombo->clear();

  count = 1;
  takethis = 1; // -1
  bodyEdit->ui.bodyForceCombo->insertItem(count++, "");

  for(int i = 0; i < bodyForceEditor.size(); i++) {
    DynamicEditor *bodyForceEdit = bodyForceEditor[i];
    if(bodyForceEdit->menuAction != NULL) {
      const QString &name = bodyForceEdit->nameEdit->text().trimmed();
      bodyEdit->ui.bodyForceCombo->insertItem(count, name);
      if ( bodyEdit->force == bodyForceEdit ) takethis = count;
      count++;
    }
  }
  connect( bodyEdit,SIGNAL(BodyForceComboChanged(BodyPropertyEditor *,QString)), 
        this, SLOT(forceComboChanged(BodyPropertyEditor *,QString)) );
  bodyEdit->ui.bodyForceCombo->setCurrentIndex(takethis-1);
    
  // Initial Condition:
  //====================
  bodyEdit->disconnect(SIGNAL(BodyInitialComboChanged(BodyPropertyEditor *,QString)));
  bodyEdit->ui.initialConditionCombo->clear();

  count = 1;
  takethis = 1; // -1
  bodyEdit->ui.initialConditionCombo->insertItem(count++, "");

  for(int i = 0; i < initialConditionEditor.size(); i++) {
    DynamicEditor *initialConditionEdit = initialConditionEditor[i];
    if(initialConditionEdit->menuAction != NULL) {
      const QString &name = initialConditionEdit->nameEdit->text().trimmed();
      bodyEdit->ui.initialConditionCombo->insertItem(count, name);
      if ( bodyEdit->initial == initialConditionEdit ) takethis = count;
      count++;
    }
  }
  connect( bodyEdit,SIGNAL(BodyInitialComboChanged(BodyPropertyEditor *,QString)), 
        this, SLOT(initialComboChanged(BodyPropertyEditor *,QString)) );
  bodyEdit->ui.initialConditionCombo->setCurrentIndex(takethis-1);
}

//-----------------------------------------------------------------------------
void MainWindow::materialComboChanged(BodyPropertyEditor *b, QString text)
{
  b->material = 0;
  b->touched = false;
  for(int i=0; i < materialEditor.size(); i++)
  {
    DynamicEditor *mat = materialEditor[i];

    if ( mat->ID >= 0 ) {
       if ( mat->nameEdit->text().trimmed()==text ) {
         b->material = mat; 
         b->touched = true;
         break;
       }
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::initialComboChanged(BodyPropertyEditor *b, QString text)
{
  b->initial = 0;
  b->touched = false;

  for( int i=0; i < initialConditionEditor.size(); i++ )
  {
    DynamicEditor *ic = initialConditionEditor[i];
    if ( ic->ID >= 0 ) {
       if ( ic->nameEdit->text().trimmed()==text ) {
         b->initial = ic; 
         b->touched = true;
         break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::forceComboChanged(BodyPropertyEditor *b, QString text)
{
  b->force = 0;
  b->touched = false;

  for( int i=0; i < bodyForceEditor.size(); i++ ) {
    DynamicEditor *bf = bodyForceEditor[i];
    if ( bf->ID >= 0 ) {
       if ( bf->nameEdit->text().trimmed()==text ) {
         b->force = bf; 
         b->touched = true;
         break;
       }
    }
  }
}

//-----------------------------------------------------------------------------
void MainWindow::equationComboChanged(BodyPropertyEditor *b, QString text)
{
  b->equation = 0;
  b->touched = false;

  for( int i=0; i < equationEditor.size(); i++ )
  {
    DynamicEditor *equ = equationEditor[i];
    if ( equ->ID >= 0 ) {
       if ( equ->nameEdit->text().trimmed() == text ) {
         b->equation = equ; 
         b->touched = true;
         break;
       }
    }
  }
}


// Edit -> Definitions...
//-----------------------------------------------------------------------------
void MainWindow::editDefinitionsSlot()
{
  if(elmerDefs == NULL)
    return;

  edfEditor->show();
}


//*****************************************************************************
//
//                                Solver MENU
//
//*****************************************************************************


// Solver -> Parallel settings
//-----------------------------------------------------------------------------
void MainWindow::parallelSettingsSlot()
{
  parallel->show();
}


// Solver -> Run solver
//-----------------------------------------------------------------------------
void MainWindow::runsolverSlot()
{
  if(!glWidget->hasMesh()) {
    logMessage("No mesh - unable to start solver");
    return;
  }
  
  if(solver->state() == QProcess::Running) {
    logMessage("Solver is already running - returning");
    return;
  }

  // Parallel solution:
  //====================
  Ui::parallelDialog ui = parallel->ui;
  bool parallelActive = ui.parallelActiveCheckBox->isChecked();
  bool partitioningActive = !ui.skipPartitioningCheckBox->isChecked();
  int nofProcessors = ui.nofProcessorsSpinBox->value();

  if(parallelActive) {

    // Set up log window:
    solverLogWindow->setWindowTitle(tr("Solver log"));
    solverLogWindow->getTextEdit()->clear();
    solverLogWindow->setFound(false);
    solverLogWindow->show();
    
    if(!partitioningActive) {
      
      // skip splitting:
      meshSplitterFinishedSlot(0);

    } else {

      // split mesh:
      if(meshSplitter->state() == QProcess::Running) {
	logMessage("Mesh partitioner is already running - aborted");
	return;
      }
      
      if (saveDirName.isEmpty()) {
	logMessage("Please save the mesh before running the parallel solver - aborted");
	return;
      }
      
      QString partitioningCommand = ui.divideLineEdit->text().trimmed();
      partitioningCommand.replace(QString("%msh"), saveDirName);
      partitioningCommand.replace(QString("%n"), QString::number(nofProcessors));
 
      logMessage("Executing: " + partitioningCommand);
      
      meshSplitter->start(partitioningCommand);
      
      if(!meshSplitter->waitForStarted()) {
	logMessage("Unable to start ElmerGrid for mesh paritioning - aborted");
	return;
      }
    }
    
    // the rest is done in meshSplitterFinishedSlot:
    return;
  }

  // Scalar solution:
  //==================
  solver->start("ElmerSolver");
  killsolverAct->setEnabled(true);

  if(!solver->waitForStarted()) {
    logMessage("Unable to start solver");
    return;
  }
  
  solverLogWindow->setWindowTitle(tr("Solver log"));
  solverLogWindow->getTextEdit()->clear();
  solverLogWindow->setFound(false);
  solverLogWindow->show();

  // convergence plot:
#ifdef EG_QWT
  convergenceView->removeData();
  convergenceView->title = "Convergence history";
#endif

  logMessage("Solver started");

  runsolverAct->setIcon(QIcon(":/icons/Solver-red.png"));

  updateSysTrayIcon("ElmerSolver started",
		    "Use Run->Kill solver to stop processing");
}


// meshSplitter emits (int) when ready...
//-----------------------------------------------------------------------------
void MainWindow::meshSplitterFinishedSlot(int exitCode)
{
  if(exitCode != 0) {
    solverLogWindow->getTextEdit()->append("MeshSplitter failed - aborting");
    logMessage("MeshSplitter failed - aborting");
    return;
  }

  logMessage("MeshSplitter ready");

  // Prepare for parallel solution:
  Ui::parallelDialog ui = parallel->ui;

  int nofProcessors = ui.nofProcessorsSpinBox->value();

  QString parallelExec = ui.parallelExecLineEdit->text().trimmed();
#ifdef WIN32
  parallelExec = "\"" + parallelExec + "\"";
#endif

  QString parallelArgs = ui.parallelArgsLineEdit->text().trimmed();
  parallelArgs.replace(QString("%n"), QString::number(nofProcessors));

  QString parallelCmd = parallelExec + " " + parallelArgs;

  logMessage("Executing: " + parallelCmd);

  solver->start(parallelCmd);
  killsolverAct->setEnabled(true);

  if(!solver->waitForStarted()) {
    solverLogWindow->getTextEdit()->append("Unable to start parallel solver");
    logMessage("Unable to start parallel solver");
    return;
  }
  
  // Set up convergence plot:
#ifdef EG_QWT
  convergenceView->removeData();
  convergenceView->title = "Convergence history";
#endif
  logMessage("Parallel solver started");
  runsolverAct->setIcon(QIcon(":/icons/Solver-red.png"));

  updateSysTrayIcon("ElmerSolver started",
		    "Use Run->Kill solver to stop processing");
}


// meshSplitter emits (void) when there is something to read from stdout:
//-----------------------------------------------------------------------------
void MainWindow::meshSplitterStdoutSlot()
{
  QString qs = meshSplitter->readAllStandardOutput();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}

// meshSplitter emits (void) when there is something to read from stderr:
//-----------------------------------------------------------------------------
void MainWindow::meshSplitterStderrSlot()
{
  QString qs = meshSplitter->readAllStandardError();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}

// meshUnifier emits (int) when ready...
//-----------------------------------------------------------------------------
void MainWindow::meshUnifierFinishedSlot(int exitCode)
{
  QStringList args;

  if(exitCode != 0) {
    solverLogWindow->getTextEdit()->append("MeshUnifier failed - aborting");
    logMessage("MeshUnifier failed - aborting");
    vtkPostMeshUnifierRunning = false;
    return;
  }

  logMessage("MeshUnifier ready");

  // Prepare for post processing parallel reults:
  //----------------------------------------------
  QString postName = generalSetup->ui.postFileEdit->text().trimmed();

  // VtkPost:
  //---------
#ifdef EG_VTK
  if(vtkPostMeshUnifierRunning) {
    vtkPost->show();
    
    vtkPost->ReadPostFile(postName);
    // if(!vtkPost->ReadPostFile(postName)) vtkPost->readEpFileSlot();
    
    vtkPostMeshUnifierRunning = false;
    return;
  }
#endif

  QFile file(postName);
  if(!file.exists()) {
    solverLogWindow->getTextEdit()->append("Elmerpost input file does not exist.");
    logMessage("Elmerpost input file does not exist.");
    vtkPostMeshUnifierRunning = false;
    return;
  }

  file.open(QIODevice::ReadOnly);
  QTextStream header(&file);

  int nn, ne, nt, nf;
  QString type, name, tstep;

  header >> nn >> ne >> nf >> nt >> type >> name;
  if ( type == "vector:" )
    name = name + "_abs";

  file.close();

  QString  simtype=generalSetup->ui.simulationTypeCombo->currentText().trimmed();
  if ( simtype.toLower() == "transient" ) {
    tstep = QString::number(1);
  } else {
    tstep = QString::number(nt);
  }

  args << "readfile " + postName + " " + tstep + 
      " " + tstep + "1; "
    "set ColorScaleY -0.85; "
    "set ColorScaleEntries  4;"
    "set ColorScaleDecimals 2;"
    "set ColorScaleColor " + name + ";"
    "set DisplayStyle(ColorScale) 1; "
    "set MeshStyle 1; "
    "set MeshColor " + name + ";"
    "set DisplayStyle(ColorMesh) 1; "
    "translate -y 0.2; "
    "UpdateObject; ";

  post->start("ElmerPost", args);
  killresultsAct->setEnabled(true);
  
  if(!post->waitForStarted()) {
    logMessage("Unable to start post processor");
    return;
  }
  
  resultsAct->setIcon(QIcon(":/icons/Post-red.png"));
  
  logMessage("Post processor started");

  updateSysTrayIcon("ElmerPost started",
		    "Use Run->Kill ElmerPost to stop processing");
}


// meshUnifier emits (void) when there is something to read from stdout:
//-----------------------------------------------------------------------------
void MainWindow::meshUnifierStdoutSlot()
{
  QString qs = meshUnifier->readAllStandardOutput();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}

// meshUnifier emits (void) when there is something to read from stderr:
//-----------------------------------------------------------------------------
void MainWindow::meshUnifierStderrSlot()
{
  QString qs = meshUnifier->readAllStandardError();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}


// solver process emits (void) when there is something to read from stdout:
//-----------------------------------------------------------------------------
void MainWindow::solverStdoutSlot()
{
  static QString qs_save = "";

  QString qs = qs_save + solver->readAllStandardOutput();

  int n = qs.lastIndexOf('\n');

  if((n > 0) && (n < qs.size()-1)) {
    qs_save = qs.mid(n+1);
    qs = qs.mid(0, n);

  } else if(n == 0) {
    if(qs.size() == 1) {
      qs_save = "";
      return;
    }
    qs_save = qs.mid(1);
    return;

  } else if(n < 0) {
      qs_save = qs;
      return;

  } else qs_save = "";

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  if(qs.isEmpty())
    return;

  solverLogWindow->getTextEdit()->append(qs);

#ifdef EG_QWT
  if(!showConvergence) {

    // hide convergence plot
    //----------------------
    convergenceView->hide();

  } else {

    // show convergence plot
    //----------------------
    if(!convergenceView->isVisible())
      convergenceView->show();
  }
#endif
    
  QStringList qsl = qs.split("\n");
  for(int i = 0; i < qsl.count(); i++) {
    QString tmp = qsl.at(i).trimmed();
    
    if(tmp.contains("Time:")) {
      QStringList tmpSplitted = tmp.split(" ");
      int last = tmpSplitted.count() - 1;
      QString timeString = tmpSplitted.at(last);
      double timeDouble = timeString.toDouble();
#ifdef EG_QWT
      convergenceView->title = "Convergence history (time="
	+ QString::number(timeDouble) + ")";
#endif
    }   
    
    if(tmp.contains("ComputeChange")) { // && tmp.contains("NS")) {
      QString copyOfTmp = tmp;
      
      // check solver name:
      QStringList tmpSplitted = tmp.split(":");
      int last = tmpSplitted.count() - 1;
      QString name = tmpSplitted.at(last).trimmed();
      
      // parse rest of the line:
      double res1 = 0.0;
      double res2 = 0.0;
      int n = tmp.indexOf("NRM,RELC");
      tmp = tmp.mid(n);
      tmpSplitted = tmp.split("(");
      
      if(tmpSplitted.count() >= 2) {
	QString tmp2 = tmpSplitted.at(1).trimmed();
	QStringList tmp2Splitted = tmp2.split(" ");
	QString qs1 = tmp2Splitted.at(0).trimmed();
	res1 = qs1.toDouble();
	int pos = 1;
	// while(tmp2Splitted.at(pos).trimmed() == "") {
	while(tmp2Splitted.at(pos).trimmed().isEmpty()) {
	  pos++;
	  if(pos > tmp2Splitted.count())
	    break;
	}
	QString qs2 = tmp2Splitted.at(pos).trimmed();
	res2 = max( qs2.toDouble(), 1.0e-16 );
      }
      
      // res1 = norm, res2 = relative change
#ifdef EG_QWT
      if(copyOfTmp.contains("NS"))	
	convergenceView->appendData(res2, "NS/" + name);
      
      if(copyOfTmp.contains("SS"))
	convergenceView->appendData(res2, "SS/" + name);
#endif
    }
  }
}


// solver process emits (void) when there is something to read from stderr:
//-----------------------------------------------------------------------------
void MainWindow::solverStderrSlot()
{
  QString qs = solver->readAllStandardError();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}



// solver process emits (int) when ready...
//-----------------------------------------------------------------------------
void MainWindow::solverFinishedSlot(int)
{
  logMessage("Solver ready");
  runsolverAct->setIcon(QIcon(":/icons/Solver.png"));
  updateSysTrayIcon("ElmerSolver has finished",
            "Use Run->Start ElmerPost, ElmerVTK or Paraview to view results");
  killsolverAct->setEnabled(false);
}


// solver process emits (QProcess::ProcessError) when error occurs...
//-----------------------------------------------------------------------------
void MainWindow::solverErrorSlot(QProcess::ProcessError error)
{
  logMessage("Solver emitted error signal: " + QString::number(error));
  solver->kill();
  runsolverAct->setIcon(QIcon(":/icons/Solver.png"));
}

// solver process emits (QProcess::ProcessState) when state changed...
//-----------------------------------------------------------------------------
void MainWindow::solverStateChangedSlot(QProcess::ProcessState state)
{
  logMessage("Solver emitted signal: QProcess::ProcessState: " + QString::number(state));
  // solver->kill();
  // runsolverAct->setIcon(QIcon(":/icons/Solver.png"));
}


// Solver -> Kill solver
//-----------------------------------------------------------------------------
void MainWindow::killsolverSlot()
{
  solver->kill();

  logMessage("Solver killed");
  runsolverAct->setIcon(QIcon(":/icons/Solver.png"));
}

// Solver -> Show convergence
//-----------------------------------------------------------------------------
void MainWindow::showConvergenceSlot()
{
  showConvergence = !showConvergence;
  synchronizeMenuToState();
}



// Solver -> Run post process
//-----------------------------------------------------------------------------
void MainWindow::resultsSlot()
{
  QStringList args;
  
  if(post->state() == QProcess::Running) {
    logMessage("Post processor is already running");
    return;
  }

  // Parallel solution:
  //====================
  Ui::parallelDialog ui = parallel->ui;
  bool parallelActive = ui.parallelActiveCheckBox->isChecked();

  if(parallelActive) {
    
    // unify mesh:
    if(meshUnifier->state() == QProcess::Running) {
      logMessage("Mesh unifier is already running - aborted");
      return;
    }
    
    if(saveDirName.isEmpty()) {
      logMessage("saveDirName is empty - unable to locate result files");
      return;
    }

    // Set up log window:
    solverLogWindow->setWindowTitle(tr("ElmerGrid log"));
    solverLogWindow->getTextEdit()->clear();
    solverLogWindow->setFound(false);
    solverLogWindow->show();

    QString postName = generalSetup->ui.postFileEdit->text().trimmed();
    QStringList postNameSplitted = postName.split(".");
    int nofProcessors = ui.nofProcessorsSpinBox->value();

    QString unifyingCommand = ui.mergeLineEdit->text().trimmed();
    unifyingCommand.replace(QString("%ep"), postNameSplitted.at(0).trimmed());
    unifyingCommand.replace(QString("%n"), QString::number(nofProcessors));
    
    logMessage("Executing: " + unifyingCommand);
    
    meshUnifier->start(unifyingCommand);
    
    if(!meshUnifier->waitForStarted()) {
      solverLogWindow->getTextEdit()->append("Unable to start ElmerGrid for mesh unification - aborted");
      logMessage("Unable to start ElmerGrid for mesh unification - aborted");
      return;
    }
    
    // The rest is done in meshUnifierFinishedSlot:
    return;
  }
   
  // Scalar solution:
  //==================
  QString postName = generalSetup->ui.postFileEdit->text().trimmed();
  QFile file(postName);
  if(!file.exists()) {
    logMessage("Elmerpost input file does not exist.");
    return;
  }

  file.open(QIODevice::ReadOnly);
  QTextStream header(&file);

  int nn, ne, nt, nf;
  QString type, name, tstep;

  header >> nn >> ne >> nf >> nt >> type >> name;
  if ( type == "vector:" )
    name = name + "_abs";

  file.close();

  QString  simtype=generalSetup->ui.simulationTypeCombo->currentText().trimmed();
  if ( simtype.toLower() == "transient" ) {
    tstep = QString::number(1);
  } else {
    tstep = QString::number(nt);
  }

  args << "readfile " + postName + " " + tstep + 
    " " + tstep + " 1; "
    "set ColorScaleY -0.85; "
    "set ColorScaleEntries  4;"
    "set ColorScaleDecimals 2;"
    "set ColorScaleColor " + name + ";"
    "set DisplayStyle(ColorScale) 1; "
    "set MeshStyle 1; "
    "set MeshColor " + name + ";"
    "set DisplayStyle(ColorMesh) 1; "
    "translate -y 0.2; "
    "UpdateObject; ";

  post->start("ElmerPost", args);
  killresultsAct->setEnabled(true);
  
  if(!post->waitForStarted()) {
    logMessage("Unable to start ElmerPost");
    return;
  }
  
  resultsAct->setIcon(QIcon(":/icons/Post-red.png"));
  
  logMessage("ElmerPost started");

  updateSysTrayIcon("ElmerPost started",
		    "Use Run->Kill ElmerPost to stop processing");
}


// Signal (int) emitted by postProcess when finished:
//-----------------------------------------------------------------------------
void MainWindow::postProcessFinishedSlot(int)
{
  logMessage("ElmerPost finished");
  resultsAct->setIcon(QIcon(":/icons/Post.png"));
  updateSysTrayIcon("ElmerPost has finished",
		    "Use Run->Start ElmerPost to restart");
  killresultsAct->setEnabled(false);
}


// Solver -> Kill post process
//-----------------------------------------------------------------------------
void MainWindow::killresultsSlot()
{
  post->kill();

  logMessage("Post process killed");
  resultsAct->setIcon(QIcon(":/icons/Post.png"));
}

// Solver -> Compile...
//-----------------------------------------------------------------------------
void MainWindow::compileSolverSlot()
{
  QString defaultDirName = getDefaultDirName();

  QString fileName = QFileDialog::getOpenFileName(this,
       tr("Open source file"), defaultDirName, tr("F90 files (*.f90)"));

  if (!fileName.isEmpty()) {
    QFileInfo fi(fileName);
    QString absolutePath = fi.absolutePath();
    QDir::setCurrent(absolutePath);
  } else {
    logMessage("Unable to open file: file name is empty");
    return;
  }

  if(compiler->state() == QProcess::Running) {
    logMessage("Compiler is currently running");
    return;
  }

  QStringList args;
#ifdef WIN32
  args << "/C";
  args << "" + QString(qgetenv("ELMER_HOME")) + "\\bin\\elmerf90.bat";
  QString dllFileName;
  dllFileName = fileName.left(fileName.lastIndexOf(".")) + ".dll";
  args << "-o";
  args << dllFileName;
#endif
  args << fileName;

#ifdef WIN32
  compiler->start("cmd.exe", args);
#else
  logMessage("Run->compiler is currently not implemented on this platform");
  return;
#endif
  
  if(!compiler->waitForStarted()) {
    logMessage("Unable to start compiler");
    return;
  }

  solverLogWindow->setWindowTitle(tr("Compiler log"));
  solverLogWindow->getTextEdit()->clear();
  solverLogWindow->setFound(false);
  solverLogWindow->show();
  solverLogWindow->statusBar()->showMessage("Compiling...");

  logMessage("Compiling...");
}


// compiler process emits (void) when there is something to read from stdout:
//-----------------------------------------------------------------------------
void MainWindow::compilerStdoutSlot()
{
  QString qs = compiler->readAllStandardOutput();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}


// compiler process emits (void) when there is something to read from stderr:
//-----------------------------------------------------------------------------
void MainWindow::compilerStderrSlot()
{
  QString qs = compiler->readAllStandardError();

  while(qs.at(qs.size()-1).unicode() == '\n')
    qs.chop(1);

  solverLogWindow->getTextEdit()->append(qs);
}

// Signal (int) emitted by compiler when finished:
//-----------------------------------------------------------------------------
void MainWindow::compilerFinishedSlot(int)
{
  logMessage("Ready");
  solverLogWindow->statusBar()->showMessage("Ready");
  solverLogWindow->getTextEdit()->append("Ready");
}



//*****************************************************************************
//
//                                Help MENU
//
//*****************************************************************************


// About dialog...
//-----------------------------------------------------------------------------
void MainWindow::showaboutSlot()
{
  QMessageBox::about(this, tr("Information about ElmerGUI"),
		     tr("ElmerGUI is a preprocessor for two and "
			"three dimensional modeling with Elmer "
			"finite element software. The program "
			"uses elmergrid, nglib, and optionally tetlib, "
			"as finite element mesh generators:\n\n"
			"http://www.csc.fi/elmer/\n"
			"http://www.hpfem.jku.at/netgen/\n"
			"http://tetgen.berlios.de/\n\n"
			"ElmerGUI uses the Qt4 Cross-Platform "
			"Application Framework by Qtsoftware:\n\n"
			"http://www.qtsoftware.com/\n\n"
#ifdef EG_VTK
			"This version of ElmerGUI contains a built-in "
			"postprocessor based on the Visualization Toolkit "
			"(VTK):\n\n"
			"http://www.vtk.org/\n\n"
#endif

#ifdef EG_PARAVIEW
            "This version of ElmerGUI has been linked "
            "against ParaView visualization software."
            "\n\n"
            "http://www.paraview.org\n\n"
#endif

#ifdef EG_OCC
			"This version of ElmerGUI has been compiled with "
			"the OpenCascade solids modeling library:\n\n"
			"http://www.opencascade.org/\n\n"
#endif

#ifdef MPICH2
			"The parallel solver of this package has been linked "
			"against the MPICH2 library v. 1.0.7 from Argonne "
			"national laboratory. In order to use the parallel "
			"solver, the MPICH2 runtime environment should be "
			"installed and configured on your system. For more "
			"details, see:\n\n"
			"http://www.mcs.anl.gov/research/projects/mpich2/\n\n"
#endif
			"The GPL-licensed source code of ElmerGUI is available "
			"from the SVN repository at Sourceforge.net\n\n"
			"http://sourceforge.net/projects/elmerfem/\n\n"
            "Written by Mikko Lyly, Juha Ruokolainen, "
            "Peter RÂback and Sampo Sillanp‰‰ 2008-2014"));
}



//*****************************************************************************
//
//                           Auxiliary non-menu items
//
//*****************************************************************************


// Log message...
//-----------------------------------------------------------------------------
void MainWindow::logMessage(QString message)
{
#if WITH_QT5
  cout << string(message.toLatin1()) << endl;
#else
  cout << string(message.toAscii()) << endl;
#endif
  statusBar()->showMessage(message, 0);
  cout.flush();
}



// Synchronize menu to GL glwidget state variables:
//-----------------------------------------------------------------------------
void MainWindow::synchronizeMenuToState()
{
  // glwidget state variables:
  if(glWidget->stateDrawSurfaceMesh)
    hidesurfacemeshAct->setChecked(true);
  else
    hidesurfacemeshAct->setChecked(false);
  
  if(glWidget->stateDrawVolumeMesh)
    hidevolumemeshAct->setChecked(true);
  else
    hidevolumemeshAct->setChecked(false);
  
  if(glWidget->stateDrawSharpEdges)
    hidesharpedgesAct->setChecked(true);
  else
    hidesharpedgesAct->setChecked(false);
  
  if(glWidget->stateFlatShade) {
    flatShadeAct->setChecked(true);
    smoothShadeAct->setChecked(false);
  } else {
    flatShadeAct->setChecked(false);
    smoothShadeAct->setChecked(true);
  }

  if(glWidget->stateOrtho) {
    orthoAct->setChecked(true);
    perspectiveAct->setChecked(false);
  } else {
    orthoAct->setChecked(false);
    perspectiveAct->setChecked(true);
  }

  if(glWidget->stateDrawSurfaceNumbers) 
    showSurfaceNumbersAct->setChecked(true);
  else 
    showSurfaceNumbersAct->setChecked(false);

  if(glWidget->stateDrawEdgeNumbers) 
    showEdgeNumbersAct->setChecked(true);
  else 
    showEdgeNumbersAct->setChecked(false);

  if(glWidget->stateDrawNodeNumbers) 
    showNodeNumbersAct->setChecked(true);
  else 
    showNodeNumbersAct->setChecked(false);

  if(glWidget->stateDrawBoundaryIndex)
    showBoundaryIndexAct->setChecked(true);
  else 
    showBoundaryIndexAct->setChecked(false);

  if(glWidget->stateDrawBodyIndex)
    showBodyIndexAct->setChecked(true);
  else 
    showBodyIndexAct->setChecked(false);

  if(glWidget->stateDrawCoordinates) 
    viewCoordinatesAct->setChecked(true);
  else 
    viewCoordinatesAct->setChecked(false);

  if(bodyEditActive)
    bodyEditAct->setChecked(true);
  else
    bodyEditAct->setChecked(false);

  if(bcEditActive)
    bcEditAct->setChecked(true);
  else
    bcEditAct->setChecked(false);

  if(showConvergence)
    showConvergenceAct->setChecked(true);
  else
    showConvergenceAct->setChecked(false);

  if(glWidget->stateBcColors)
    showBoundaryColorAct->setChecked(true);
  else
    showBoundaryColorAct->setChecked(false);

  if(glWidget->stateBodyColors)
    showBodyColorAct->setChecked(true);
  else
    showBodyColorAct->setChecked(false);

  if(isFullScreen())
    viewFullScreenAct->setChecked(true);
  else
    viewFullScreenAct->setChecked(false);
}


// Load definitions...
//-----------------------------------------------------------------------------
void MainWindow::loadDefinitions()
{
  // Determine edf-file location and name:
  //--------------------------------------
  QString elmerGuiHome;

#ifdef __APPLE__DONTGO_HERE_TODO
  QString generalDefs = this->homePath +  "/edf/edf.xml";            
#else
  QString generalDefs = QCoreApplication::applicationDirPath() + "/../share/ElmerGUI/edf/edf.xml";  // @TODO: fix path to share/ElmerGUI/edf

  elmerGuiHome = QString(getenv("ELMERGUI_HOME"));
  
  if(!elmerGuiHome.isEmpty())
    generalDefs = elmerGuiHome + "/edf/edf.xml";  

  // ML 5. August 2010
  generalDefs.replace('\\', '/');
#endif

  // Load general definitions file:
  //--------------------------------
#if WITH_QT5
  cout << "Load " << string(generalDefs.toLatin1()) << "... ";
#else
  cout << "Load " << string(generalDefs.toAscii()) << "... ";
#endif
  cout.flush();
  updateSplash("Loading general definitions...");

  QFile file(generalDefs);
  
  QString errStr;
  int errRow;
  int errCol;
  
  if(!file.exists()) {

    elmerDefs = NULL;
    QMessageBox::information(window(), tr("Edf loader: ") + generalDefs,
			     tr("Definitions file does not exist"));
    return;

  } else {  

    if(!elmerDefs->setContent(&file, true, &errStr, &errRow, &errCol)) {
      QMessageBox::information(window(), tr("Edf loader: ") + generalDefs,
			       tr("Parse error at line %1, col %2:\n%3")
			       .arg(errRow).arg(errCol).arg(errStr));
      file.close();
      return;

    } else {

      if(elmerDefs->documentElement().tagName() != "edf") {
	QMessageBox::information(window(), tr("Edf loader: ") + generalDefs,
				 tr("This is not an edf file"));
	delete elmerDefs;
	file.close();	
	return;

      }
    }
  }

  edfEditor->setupEditor(elmerDefs);
  file.close();

  cout << "done" << endl;
  cout.flush();

  // load additional definitions:
  //-----------------------------
#ifdef __APPLE__DONTGO_HERE_TODO
  QDirIterator iterator( homePath+"/edf", QDirIterator::Subdirectories);  
#else
  QString additionalEdfs = QCoreApplication::applicationDirPath() + "/../share/ElmerGUI/edf";  // @TODO: fix path to share/ElmerGUI/edf

  if(!elmerGuiHome.isEmpty()) 
    additionalEdfs = elmerGuiHome + "/edf";  

  QDirIterator iterator(additionalEdfs, QDirIterator::Subdirectories);
#endif

  while (iterator.hasNext()) {
    QString fileName = iterator.next();

    // ML 5. August 2010
    fileName.replace('\\', '/');

    QFileInfo fileInfo(fileName);
    QString fileSuffix = fileInfo.suffix();

    // The names "egini" and "egmaterials" are reserved, skip them:
    if(fileInfo.completeBaseName() == "egini")
      continue;

    if(fileInfo.completeBaseName() == "egmaterials")
      continue;

    if((fileSuffix == "xml") && (fileName != generalDefs)) {

#if WITH_QT5
      cout << "Load " << string(fileName.toLatin1()) << "... ";
#else
      cout << "Load " << string(fileName.toAscii()) << "... ";
#endif
      cout.flush();

      updateSplash("Loading " + fileName + "...");

      file.setFileName(fileName);

      QDomDocument tmpDoc;
      tmpDoc.clear();

      if(!tmpDoc.setContent(&file, true, &errStr, &errRow, &errCol)) {
	QMessageBox::information(window(), tr("Edf loader: ") + fileName,
				 tr("Parse error at line %1, col %2:\n%3")
				 .arg(errRow).arg(errCol).arg(errStr));
	file.close();
	return;

      } else {

	if(tmpDoc.documentElement().tagName() != "edf") {
	  QMessageBox::information(window(), tr("Edf loader: ") + fileName,
				   tr("This is not an edf file"));
	  file.close();
	  return;      
	}
      }
      
      // add new elements to the document
      QDomElement root = elmerDefs->documentElement();
      QDomElement tmpRoot = tmpDoc.documentElement();
      QDomElement element = tmpRoot.firstChildElement();

      while(!element.isNull()) {
	root.appendChild(element);
	element = tmpRoot.firstChildElement();
      }
      
      edfEditor->setupEditor(elmerDefs);

      file.close();

      cout << "done" << endl;
      cout.flush();
    }
  }

  // Load qss:
  //-----------
  QString qssFileName = QCoreApplication::applicationDirPath() + "/elmergui.qss"; // @TODO: fix path to share/ElmerGUI

#ifdef __APPLE__
  qssFileName = homePath + "/elmergui.qss";
#else
  if(!elmerGuiHome.isEmpty()) 
    qssFileName = elmerGuiHome + "/elmergui.qss";
#endif

  QFile qssFile(qssFileName);
  
  if(qssFile.exists()) {
    cout << "Loading QSS style sheet... ";
    qssFile.open(QFile::ReadOnly);
    QString styleSheet = QLatin1String(qssFile.readAll());
    qssFile.close();
    qApp->setStyleSheet(styleSheet);    
    cout << "done" << endl;
  }
}

// Setup splash...
//-----------------------------------------------------------------------------
void MainWindow::setupSplash()
{
  QStringList args = QCoreApplication::arguments();
  
  if(args.contains("-nogui"))
    return;

  if(egIni->isSet("splashscreen")) {
    pixmap.load(":/images/splash.png");
    splash.setPixmap(pixmap);
    splash.show();
    qApp->processEvents();
  }
}


// Update splash...
//-----------------------------------------------------------------------------
void MainWindow::updateSplash(QString text)
{
  QStringList args = QCoreApplication::arguments();
  
  if(args.contains("-nogui"))
    return;

  if(!egIni->isSet("splashscreen"))
    return;

  if(splash.isVisible()) {
    splash.showMessage(text, Qt::AlignBottom);
    qApp->processEvents();
  }
}

// Finalize splash...
//-----------------------------------------------------------------------------
void MainWindow::finalizeSplash()
{
  if(!egIni->isSet("splashscreen"))
    return;

  if(splash.isVisible())
    splash.finish(this);
}

// Setup system tray icon...
//-----------------------------------------------------------------------------
void MainWindow::setupSysTrayIcon()
{
  sysTrayIcon = NULL;

  QStringList args = QCoreApplication::arguments();
  
  if(args.contains("-nogui"))
    return;

  if(!egIni->isSet("systrayicon"))
    return;

  if(QSystemTrayIcon::isSystemTrayAvailable()) {
    sysTrayIcon = new QSystemTrayIcon(this);
    sysTrayIcon->setIcon(QIcon(":/icons/Mesh3D.png"));
    sysTrayIcon->setVisible(true);
    sysTrayIcon->setContextMenu(sysTrayMenu);
  }
}

// Update system tray icon...
//-----------------------------------------------------------------------------
void MainWindow::updateSysTrayIcon(QString label, QString msg)
{
  int duration = 3000;

  QStringList args = QCoreApplication::arguments();
  
  if(args.contains("-nogui"))
    return;

  if(!sysTrayIcon)
    return;

  if(!egIni->isSet("systraymessages"))
    return;

  if(isFullScreen())
    return;
  
  if(egIni->isPresent("systraymsgduration"))
    duration = egIni->value("systraymsgduration").toInt();

  if(sysTrayIcon->supportsMessages())
    sysTrayIcon->showMessage(label, msg, QSystemTrayIcon::Information, duration);

}

// Finalize system tray icon...
//-----------------------------------------------------------------------------
void MainWindow::finalizeSysTrayIcon()
{
}

// Get default open/save directory
//-----------------------------------------------------------------------------
QString MainWindow::getDefaultDirName()
{
  QString defaultDirName = "";

#ifdef WIN32
  defaultDirName = egIni->value("win32defaultdir");
#else
#ifdef __APPLE__
  defaultDirName = egIni->value("macxdefaultdir");
#else
  defaultDirName = egIni->value("unixdefaultdir");
#endif
#endif

  if(!saveDirName.isEmpty())
    defaultDirName = saveDirName;

  return defaultDirName;
}
