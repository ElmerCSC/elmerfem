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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDomDocument>
#include <QMainWindow>
#include <QProcess>

#include "bodypropertyeditor.h"
#include "boundarydivision.h"
#include "boundarypropertyeditor.h"
#include "checkmpi.h"
#include "dynamiceditor.h"
#include "edfeditor.h"
#include "egini.h"
#include "generalsetup.h"
#include "glcontrol.h"
#include "glwidget.h"
#include "materiallibrary.h"
#include "maxlimits.h"
#include "meshcontrol.h"
#include "meshingthread.h"
#include "meshutils.h"
#include "objectbrowser.h"
#include "operation.h"
#include "parallel.h"
#include "plugins/elmergrid_api.h"
#include "plugins/nglib_api.h"
#include "plugins/tetlib_api.h"
#include "sifgenerator.h"
#include "sifwindow.h"
#include "solverlogwindow.h"
#include "summaryeditor.h"
#include "twod/twodview.h"

#ifdef EG_QWT
#include "convergenceview.h"
#endif

#if __APPLE__
// #ifndef EG_OCC
// #define EG_OCC
// #endif
#endif

#ifdef EG_OCC
#include "cad/cadview.h"
#endif

class QAction;
class QMenu;
class GLWidget;
class QProgressBar;
class QSystemTrayIcon;
class QContextMenuEvent;

#ifdef EG_VTK
class VtkPost;
#endif

#define MAXPATHLENGTH 600

class MainWindow : public QMainWindow {
  Q_OBJECT

  friend class ObjectBrowser;

public:
  MainWindow();
  ~MainWindow();

  void parseCmdLine();

  static QVariant settings_value(const QString &key,
                          const QVariant &defaultValue = QVariant());
  static void settings_setValue(const QString &key, const QVariant &value);
  void saveAndRun(bool generateSif);
  void showContextMenu(QPoint);

protected:
  void contextMenuEvent(QContextMenuEvent *event);
  void closeEvent(QCloseEvent *event);
  
private slots:
  // menu slots:
  void openSlot();                      // File -> Open...
  void loadSlot();                      // File -> Load...
  void loadProjectSlot();               // File -> Load project...
  void saveSlot();                      // File -> Save...
  void saveAsSlot();                    // File -> Save As...
  void saveProjectSlot();               // File -> Save project
  void saveProjectAsSlot();             // File -> Save project as...
  void newProjectSlot();                // File -> New project...
  void savePictureSlot();               // File -> Save picture as...
  void grabFrameSlot();                 // utility slot
  void closeMainWindowSlot();           // File -> exit
  void modelSetupSlot();                // Model -> Setup...
  void addEquationSlot();               // Model -> Equation...
  void addMaterialSlot();               // Model -> Material...
  void addBodyForceSlot();              // Model -> Body force...
  void addInitialConditionSlot();       // Model -> Initial condition...
  void addBoundaryConditionSlot();      // Model -> Boundary condition...
  void bodyEditSlot();                  // Model -> Set body properties
  void bcEditSlot();                    // Model -> Set boundary conditions
  void modelSummarySlot();              // Model -> Summary...
  void modelClearSlot();                // Model -> Clear
  void generateSifSlot();               // Edit -> Generate sif
  void showsifSlot();                   // Edit -> Solver input file...
  void editDefinitionsSlot();           // Edit -> Definitions...
  void meshcontrolSlot();               // Mesh -> Control...
  void remeshSlot();                    // Mesh -> Remesh
  void stopMeshingSlot();               // Mesh -> Kill generator
  void surfaceDivideSlot();             // Mesh -> Divide surface...
  void surfaceUnifySlot();              // Mesh -> Unify surface
  void edgeUnifySlot();                 // Mesh -> Unify edge
  void edgeDivideSlot();                // Mesh -> Divide edge...
  void cleanHangingSharpEdgesSlot();    // Mesh -> Clean up
  void viewFullScreenSlot();            // View -> Full screen
  void hidesurfacemeshSlot();           // View -> Surface mesh
  void hidevolumemeshSlot();            // View -> Volume mesh
  void hidesharpedgesSlot();            // View -> Sharp edges
  void viewCoordinatesSlot();           // View -> Coordinates
  void selectAllSurfacesSlot();         // View -> Select all surfaces
  void selectAllEdgesSlot();            // View -> Select all edges
  void selectDefinedEdgesSlot();        // View -> Select defined edges
  void showSurfaceNumbersSlot();    // View -> Show numbering -> surface numbering
  void showEdgeNumbersSlot();       // View -> Show numbering -> edge numbering
  void showNodeNumbersSlot();       // View -> Show numbering -> node numbering
  void showBoundaryIndexSlot();     // View -> Show numbering -> boundary index
  void showBodyIndexSlot();         // View -> Show numbering -> body index
  void glControlSlot();             // View -> Colors -> GL controls
  void backgroundColorSlot();       // View -> Colors -> Background
  void surfaceColorSlot();          // View -> Colors -> Surfaces
  void edgeColorSlot();             // View -> Colors -> Edges
  void surfaceMeshColorSlot();      // View -> Colors -> Surface mesh
  void sharpEdgeColorSlot();        // View -> Colors -> Sharp edges
  void selectionColorSlot();        // View -> Colors -> Selection
  void colorizeBoundarySlot();      // View -> Colors -> Boundaries
  void colorizeBodySlot();          // View -> Colors -> Bodies
  void selectDefinedSurfacesSlot(); // View -> Select defined surfaces
  void hideselectedSlot();          // View -> Hide/show selected
  void showallSlot();               // View -> Show all
  void resetSlot();                 // View -> Reset model view
  void flatShadeSlot();             // View -> Shade model -> flat
  void smoothShadeSlot();           // View -> Shade model -> smooth
  void orthoSlot();                 // View -> Projection -> Ortho
  void perspectiveSlot();           // View -> Projection -> Perspective
  void showCadModelSlot();          // View -> Show cad model...
  void showTwodViewSlot();          // View -> Show 2D view...
  void showVtkPostSlot();           // View -> Show VTK post processor...
  void showParaViewSlot();          // View -> Use ParaView for postprocessing
  void showObjectBrowserSlot();     // view -> Show Object Browser
  void parallelSettingsSlot();      // Solver -> Parallel settings
  void runsolverSlot();             // Solver -> Run solver
  void killsolverSlot();            // Solver -> Kill solver
  void showConvergenceSlot();       // Solver -> Show convergence...
  void resultsSlot();               // Solver -> Post process
  void killresultsSlot();           // Solver -> Kill post process
  void compileSolverSlot();         // Solver -> Compile...
  void showaboutSlot();             // Help -> About...
  void getStartedSlot();            // Help -> Get Started...
  void generateAndSaveAndRunSlot();

  // other private slots:
  void meshingStartedSlot();    // signal emitted by meshingThread
  void meshingTerminatedSlot(); // signal emitted by meshingThread
  void meshingFinishedSlot();   // signal emitted by meshingThread

  void boundarySelectedSlot(list_t *, Qt::KeyboardModifiers); // signal emitted by glWidget
  void doDivideSurfaceSlot(double);    // signal emitted by boundaryDivide
  void doDivideEdgeSlot(double);       // signal emitted by boundaryDivide

  void postProcessFinishedSlot(int); // signal emitted by postProcess
  void paraviewProcessFinishedSlot(int); // signal emitted by paraview
  
  void solverStdoutSlot();      // solver's stdout redirection
  void solverStderrSlot();      // solver's stderr redirection
  void solverFinishedSlot(int); // signal emitted by solver process
  void solverErrorSlot(QProcess::ProcessError);        // solver error signal
  void solverStateChangedSlot(QProcess::ProcessState); // state changed

  void compilerStdoutSlot();      // compiler's stdout redirection
  void compilerStderrSlot();      // compiler's stderr redirection
  void compilerFinishedSlot(int); // signal emitted by compiler

  void meshSplitterStdoutSlot();      // meshSplitter's stdout redirection
  void meshSplitterStderrSlot();      // meshSplitter's stderr redirection
  void meshSplitterFinishedSlot(int); // signal emitted by meshSplitter

  void meshUnifierStdoutSlot();      // meshUnifier's stdout redirection
  void meshUnifierStderrSlot();      // meshUnifier's stderr redirection
  void meshUnifierFinishedSlot(int); // signal emitted by meshUnifier

  void pdeEditorFinishedSlot(int, int);       // signal emitted by pde editor
  void matEditorFinishedSlot(int, int);       // signal emitted by mat editor
  void bodyForceEditorFinishedSlot(int, int); // signal emitted by bf editor
  void initialConditionEditorFinishedSlot(int, int);  // emitted by ic editor
  void boundaryConditionEditorFinishedSlot(int, int); // emitted by bc editor

  void equationSelectedSlot(QAction *);  // signal emitted by Equation menu
  void materialSelectedSlot(QAction *);  // signal emitted by Material menu
  void bodyForceSelectedSlot(QAction *); // signal emitted by BodyForce menu
  void initialConditionSelectedSlot(QAction *);  // emitted by ic menu
  void boundaryConditionSelectedSlot(QAction *); // emitted by bc menu

  void materialComboChanged(BodyPropertyEditor *, QString);
  void initialComboChanged(BodyPropertyEditor *, QString);
  void forceComboChanged(BodyPropertyEditor *, QString);
  void equationComboChanged(BodyPropertyEditor *, QString);
  void boundaryAsABodyChanged(BoundaryPropertyEditor *, int);
  void boundaryComboChanged(BoundaryPropertyEditor *, QString);

  void dynamicEditorNameChange(QString);

  void editNumericalMethods(int, int); // signal emitted by dynamic editor
  void showMaterialLibrary(int, int);  // signal emitted by dynamic editor
  void materialBodyChanged(int);
  void initialBodyChanged(int);
  void forceBodyChanged(int);
  void bcBoundaryChanged(int);
  void equationBodyChanged(int);

  void viewNormalModeSlot();

  void menuBarTriggeredSlot(QAction *);

  void loadRecentProject0Slot();
  void loadRecentProject1Slot();
  void loadRecentProject2Slot();
  void loadRecentProject3Slot();
  void loadRecentProject4Slot();
  void loadRecentProject5Slot();
  void loadRecentProject6Slot();
  void loadRecentProject7Slot();
  void loadRecentProject8Slot();
  void loadRecentProject9Slot();

  void selectElmerPostSlot();
  void selectVtkPostSlot();  
  void selectParaViewSlot();
  
private:
  // widgets and helpers:
  GLWidget *glWidget;               // central gl widget
  SifWindow *sifWindow;             // sif text editor
  MeshControl *meshControl;         // mesh generator control
  BoundaryDivide *boundaryDivide;   // boundary division control
  Meshutils *meshutils;             // mesh manipulation utilities
  MeshingThread *meshingThread;     // meshing thread
  SolverLogWindow *solverLogWindow; // Solver log
  SifGenerator *sifGenerator;       // SIF generator
  EdfEditor *edfEditor;             // Edf editor
#ifdef EG_QWT
  ConvergenceView *convergenceView; // Convergence plotter
#endif

  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
  void applyOperations();
  void populateBodyComboBoxes(BodyPropertyEditor *);
  void populateBoundaryComboBoxes(BoundaryPropertyEditor *);
  void saveProjectContents(QDomDocument, QString, QVector<DynamicEditor *> &);
  void loadProjectContents(QDomElement, QVector<DynamicEditor *> &, QString);
  QString getDefaultDirName();
  void loadProject(QString);
  bool saveProject(QString);
  void loadSettings();
  void saveSettings();
  bool loadExtraSolver(
      QString); // load the solver with specified solver name, Nov 2019 by TS
  void checkAndLoadExtraSolvers(QFile *);

  QMenu *fileMenu;              // File menu
  QMenu *recentProjectsMenu;    // File -> Recent projects menu
  QMenu *modelMenu;             // Model menu
  QMenu *equationMenu;          // Model -> Equation menu
  QMenu *materialMenu;          // Model -> Material menu
  QMenu *bodyForceMenu;         // Model -> Body force...
  QMenu *initialConditionMenu;  // Model -> Initial condition...
  QMenu *boundaryConditionMenu; // Model -> Boundary condition...
  QMenu *editMenu;              // Edit menu
  QMenu *viewMenu;              // View menu
  QMenu *shadeMenu;             // View -> Shade model menu
  QMenu *projectionMenu;        // View -> Projection menu
  QMenu *numberingMenu;         // View -> Show numbering menu
  QMenu *colorizeMenu;          // View -> Colors menu
  QMenu *meshMenu;              // Mesh menu
  QMenu *solverMenu;            // Solver menu
  QMenu *helpMenu;              // Help menu
  QMenu *sysTrayMenu;           // System tray menu
  QMenu *contextMenu;           // Context menu
  QMenu *selectPostMenu;

  QToolBar *fileToolBar;   // File toolbar
  QToolBar *editToolBar;   // Edit toolbar
  QToolBar *meshToolBar;   // Mesh toolbar
  QToolBar *solverToolBar; // Solver toolbar

  QAction *openAct;        // File -> Open...
  QAction *loadAct;        // File -> Load...
  QAction *loadProjectAct; // File -> Load project....
  QAction *newProjectAct;  // File -> New Project...
  QAction *recentProject0Act;
  QAction *recentProject1Act;
  QAction *recentProject2Act;
  QAction *recentProject3Act;
  QAction *recentProject4Act;
  QAction *recentProject5Act;
  QAction *recentProject6Act;
  QAction *recentProject7Act;
  QAction *recentProject8Act;
  QAction *recentProject9Act;
  QAction *saveAct;                      // File -> Save...
  QAction *saveAsAct;                    // File -> Save As...
  QAction *saveProjectAct;               // File -> Save project
  QAction *saveProjectAsAct;             // File -> Save project as...
  QAction *savePictureAct;               // File -> Save picture as...
  QAction *exitAct;                      // File -> Exit
  QAction *modelSetupAct;                // Model -> Setup...
  QAction *addEquationAct;               // Model -> Equation...
  QAction *addMaterialAct;               // Model -> Material...
  QAction *addBodyForceAct;              // Model -> Body force...
  QAction *addInitialConditionAct;       // Model -> Initial condition...
  QAction *addBoundaryConditionAct;      // Model -> Boundary condition...
  QAction *bodyEditAct;                  // Model -> Set body properties
  QAction *bcEditAct;                    // Model -> Set boundary conditions
  QAction *modelSummaryAct;              // Model -> Summary...
  QAction *modelClearAct;                // Model -> Clear
  QAction *generateSifAct;               // Edit -> Generate sif
  QAction *showsifAct;                   // Edit -> Edit SIF...
  QAction *editDefinitionsAct;           // Edit -> Edit SIF...
  QAction *viewFullScreenAct;            // View -> Full screen
  QAction *hidesurfacemeshAct;           // View -> Show surface mesh
  QAction *hidevolumemeshAct;            // View -> Show volume mesh
  QAction *hidesharpedgesAct;            // View -> Show sharp edges
  QAction *viewCoordinatesAct;           // View -> Show sharp edges
  QAction *selectAllSurfacesAct;         // View -> Select all surfaces
  QAction *selectAllEdgesAct;            // View -> Select all edges
  QAction *selectDefinedEdgesAct;        // View -> Select defined edges
  QAction *selectDefinedSurfacesAct;     // View -> Select defined surfaces
  QAction *showSurfaceNumbersAct; // View -> Show numbering -> element numbers
  QAction *showEdgeNumbersAct;    // View -> Show numbering -> edge numbers
  QAction *showNodeNumbersAct;    // View -> Show numbering -> node numbers
  QAction *showBoundaryIndexAct;  // View -> Show numbering -> boundary index
  QAction *showBodyIndexAct;      // View -> Show numbering -> body index
  QAction *glControlAct;          // View -> Colors -> GL controls
  QAction *chooseBGColorAct;      // View -> Colors -> Background color
  QAction *chooseSurfaceColorAct; // View -> Colors -> Surface color
  QAction *chooseSurfaceMeshColorAct; // View -> Colors -> Surface mesh
  QAction *chooseSharpEdgeColorAct;   // View -> Colors -> Sharp edges
  QAction *chooseEdgeColorAct;        // View -> Colors -> Edge color
  QAction *chooseSelectionColorAct;   // View -> Colors -> Selection
  QAction *showBoundaryColorAct;      // View -> Colors -> Boundaries
  QAction *showBodyColorAct;          // View -> Colors -> Body
  QAction *hideselectedAct;           // View -> Show selected
  QAction *flatShadeAct;              // View -> Shade model -> Flat
  QAction *smoothShadeAct;            // View -> Shade model -> Smooth
  QAction *orthoAct;                  // View -> Projection -> Ortho
  QAction *perspectiveAct;            // View -> Projection -> Perspective
  QAction *showallAct;                // View -> Show all
  QAction *resetAct;                  // View -> Reset model view
  QAction *showCadModelAct;           // View -> Show cad model...
  QAction *showTwodViewAct;           // View -> Show 2d view...
  QAction *showVtkPostAct;            // View -> Show VTK post processor...
  QAction *showObjectBrowserAct;      // View -> Show Object Browser
  QAction *meshcontrolAct;            // Mesh -> Control...
  QAction *remeshAct;                 // Mesh -> Remesh
  QAction *stopMeshingAct;            // Mesh -> Kill generator
  QAction *surfaceDivideAct;          // Mesh -> Divide surface...
  QAction *surfaceUnifyAct;           // Mesh -> Unify surface
  QAction *edgeDivideAct;             // Mesh -> Divide edges...
  QAction *edgeUnifyAct;              // Mesh -> Unify edge
  QAction *cleanHangingSharpEdgesAct; // Mesh -> Clean up
  QAction *parallelSettingsAct;       // Solver -> Parallel settings
  QAction *runsolverAct;              // Solver -> Run solver
  QAction *killsolverAct;             // Solver -> Kill solver
  QAction *showConvergenceAct;        // Solver -> Show convergence...
  QAction *resultsAct;                // Solver -> Post process
  QAction *killresultsAct;            // Solver -> Kill post process
  QAction *paraviewAct;               // Solver -> Launch Paraview
  QAction *compileSolverAct;          // Solver -> Compile...
  QAction *aboutAct;                  // Help -> About...
  QAction *getStartedAct;             // Help -> Get Started...
  QAction *generateAndSaveAndRunAct;
  QAction *runPostProcessorAct;  
  QAction *selectElmerPostAct;
  QAction *selectVtkPostAct;
  QAction *selectParaViewAct;  

  // property editors etc:
  GeneralSetup *generalSetup;

  QVector<DynamicEditor *> equationEditor;
  QVector<DynamicEditor *> materialEditor;
  QVector<DynamicEditor *> bodyForceEditor;
  QVector<DynamicEditor *> initialConditionEditor;
  QVector<DynamicEditor *> boundaryConditionEditor;
  QVector<BoundaryPropertyEditor *> boundaryPropertyEditor;
  QVector<BodyPropertyEditor *> bodyPropertyEditor;
  QVector<SolverParameterEditor *> solverParameterEditor;

  SummaryEditor *summaryEditor;
  GLcontrol *glControl;
  Parallel *parallel;
  CheckMpi *checkMpi;
  MaterialLibrary *materialLibrary;

#ifdef EG_OCC
  CadView *cadView;
#endif

  TwodView *twodView;

#ifdef EG_VTK
  VtkPost *vtkPost;
#endif

  // elmer definitions:
  QDomDocument *elmerDefs;

  // tetlib:
  bool tetlibPresent;
  TetlibAPI *tetlibAPI;
  tetgenio *in;
  tetgenio *out;
  QString tetlibControlString;
  bool tetlibInputOk;

  // nglib:
  bool nglibPresent;
  NglibAPI *nglibAPI;
  nglib::Ng_Mesh *ngmesh;
  nglib::Ng_STL_Geometry *nggeom;
  nglib::Ng_Geometry_2D *nggeom2d;
  nglib::Ng_Meshing_Parameters mp;
  int ngDim;
  QString stlFileName;
  QString in2dFileName;
  bool nglibInputOk;

  // occ:
  bool occInputOk;

  // vtkPost:
  bool vtkPostMeshUnifierRunning;

  // elmergrid:
  ElmergridAPI *elmergridAPI;

  // solver, post processor, and other processes:
  QProcess *solver;
  QProcess *post;
  QProcess *paraview;
  QProcess *compiler;
  QProcess *meshSplitter;
  QProcess *meshUnifier;

  // utility functions:
  void readInputFile(QString);
  void loadElmerMesh(QString);
  void saveElmerMesh(QString);
  void makeElmerMeshFromTetlib();
  void makeElmerMeshFromNglib();
  void logMessage(QString);
  void loadDefinitions();
  void createBoundaryCheckBoxes(DynamicEditor *);
  void createBodyCheckBoxes(int, DynamicEditor *);
  void synchronizeMenuToState();

  // state variables:
  int activeGenerator;
  bool bcEditActive;
  bool bodyEditActive;
  bool showConvergence;
  QString saveDirName;
  QString geometryInputFileName;

  // splash screen:
  QPixmap pixmap;
  QSplashScreen splash;
  void setupSplash();
  void updateSplash(QString);
  void finalizeSplash();

  // sys tray icon:
  QSystemTrayIcon *sysTrayIcon;
  void setupSysTrayIcon();
  void updateSysTrayIcon(QString, QString);
  void finalizeSysTrayIcon();

  // initialization:
  EgIni *egIni;

  // limits:
  void setDynamicLimits();
  Limit *limit;

  // operations:
  int operations;
  operation_t operation;

  // progress bar:
  QProgressBar *progressBar;
  QLabel *progressLabel;

  // screen shot:
  QTimeLine *grabTimeLine;
  QString pictureFileName;

  //  #ifdef __APPLE__
  //  This is only needed for Mac OS X, but it's easier to include in all
  //  architectures and it's small so there's no marked adverse effects
  QString homePath;
  //  #endif

  // variables and functions for "Recent projects..." menu
  QStringList recentProject;
  void addRecentProject(QString, bool);

  // String to store current project dir for "generate, save and run" button
  QString currentProjectDirName;


  ObjectBrowser *objectBrowser;

public:  
  /*
  rebuildGLLists() is assumed to be called from ObjectBrowser to avoid a problem of 3D surface 
  mesh not shown correctly when project loading (typically, TemperatureGeneric sample)
  */
  void rebuildGLLists();
};

#endif // MAINWINDOW_H
