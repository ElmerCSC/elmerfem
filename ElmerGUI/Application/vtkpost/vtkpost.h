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

#ifndef VTKPOST_H
#define VTKPOST_H

#define MAX_SCALARS 100

#include <QMainWindow>
#include <QHash>
#include <QTextStream>

#ifdef EG_PYTHONQT
#include <PythonQt.h>
#include <gui/PythonQtScriptingConsole.h>
#endif

class EpMesh;
class ScalarField;
class QVTKWidget;
class vtkRenderer;
class vtkRenderWindow;
class vtkActor;
class vtkTextActor;
class vtkFollower;
class vtkScalarBarActor;
class vtkDataSetMapper;
class vtkUnstructuredGrid;
class vtkLookupTable;
class vtkPlane;
class vtkAxes;
class vtkImplicitPlaneWidget;
class vtkCamera;
class IsoSurface;
class IsoContour;
class ColorBar;
class Surface;
class Preferences;
class Vector;
class ReadEpFile;
class StreamLine;
class TimeStep;
class Axes;
class Text;
class FeatureEdge;
class MeshPoint;
class MeshEdge;
class Matc;
class EcmaConsole;
class QScriptEngine;

class VtkPost : public QMainWindow
{
  Q_OBJECT

public:
  VtkPost(QWidget *parent = 0);
  ~VtkPost();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  QVTKWidget* GetQVTKWidget();
  vtkRenderer* GetRenderer();
  vtkActor* GetSurfaceActor();
  vtkActor* GetVectorActor();
  vtkActor* GetIsoContourActor();
  vtkActor* GetIsoSurfaceActor();
  vtkActor* GetStreamLineActor();
  vtkActor* GetAxesActor();
  vtkFollower* GetAxesXTextActor();
  vtkFollower* GetAxesYTextActor();
  vtkFollower* GetAxesZTextActor();
  vtkScalarBarActor* GetColorBarActor();
  vtkActor* GetFeatureEdgeActor();
  vtkActor* GetPickedPointActor();
  vtkActor* GetMeshPointActor();
  vtkActor* GetMeshEdgeActor();
  vtkUnstructuredGrid* GetLineGrid();
  vtkUnstructuredGrid* GetSurfaceGrid();
  vtkUnstructuredGrid* GetVolumeGrid();
  vtkPlane* GetClipPlane();
  vtkImplicitPlaneWidget* GetPlaneWidget();
  vtkLookupTable* GetCurrentLut();
  ScalarField* GetScalarField();
  EpMesh* GetEpMesh();
  Preferences* GetPreferences();
  vtkTextActor* GetTextActor();

  void minMax(ScalarField *);
  ScalarField* addScalarField(QString, int, double *);
  void SetClipPlaneOrigin(double*);
  void SetClipPlaneNormal(double*);
  void GetBounds(double*);
  double* GetCurrentPickPosition();
  void SetCurrentPickPosition(double*);
  int GetScalarFields();
  void SetScalarFields(int);
  QString GetCurrentSurfaceName();
  void SetCurrentSurfaceName(QString);
  QString GetCurrentVectorName();
  QString GetCurrentVectorColorName();
  void SetCurrentVectorName(QString);
  void SetCurrentVectorColorName(QString);
  QString GetCurrentIsoContourName();
  QString GetCurrentIsoContourColorName();
  void SetCurrentIsoContourName(QString);
  void SetCurrentIsoContourColorName(QString);
  QString GetCurrentIsoSurfaceName();
  QString GetCurrentIsoSurfaceColorName();
  void SetCurrentIsoSurfaceName(QString);
  void SetCurrentIsoSurfaceColorName(QString);
  QString GetCurrentStreamLineName();
  QString GetCurrentStreamLineColorName();
  void SetCurrentStreamLineName(QString);
  void SetCurrentStreamLineColorName(QString);
  int NofNodes();

signals:
  void canProceedWithNextSignal(vtkRenderWindow*);
  void povrayState(int value);

public slots:
  void redrawSlot();                                // redraw all actors

#ifdef EG_MATC
  QString MatcCmd(QString);                         // evaluate matc cmd
  QString domatcSlot();                             // flush matc console
  void matcOpenSlot();                              // open matc console
  void matcCutPasteSlot();                          // handle cut-paste
#endif

  void SetPostFileStart(int);                       // first time step
  void SetPostFileEnd(int);                         // last time step
  bool ReadPostFile(QString);                       // read result file

  void Render();                                    // render
  void ResetCamera();                               // reset camera
  void Redraw();                                    // redraw actors
  void ResetAll();                                  // reset view

  void SetSurfaces(bool);                           // show/hide surfaces
  void SetVectors(bool);                            // show/hide vectors
  void SetIsoContours(bool);                        // show/hide isocontours
  void SetIsoSurfaces(bool);                        // show/hide isosurfaces
  void SetStreamLines(bool);                        // show/hide streamlines
  void SetColorBar(bool);                           // show/hide colorbar
  void SetMeshPoints(bool);                         // show/hide/nodes
  void SetMeshEdges(bool);                          // show/hide edges
  void SetFeatureEdges(bool);                       // show/hide f-edges
  void SetAxes(bool);                               // show/hide axes
  void SetText(bool);                               // show/hide text

  bool GetClipAll();                                // is clipping on?
  void SetClipAll(bool);                            // clipping on/off
  void SetClipPlaneOx(double);                      // clip plane origin
  void SetClipPlaneOy(double);                      // clip plane origin
  void SetClipPlaneOz(double);                      // clip plane origin
  void SetClipPlaneNx(double);                      // clip plane normal
  void SetClipPlaneNy(double);                      // clip plane normal
  void SetClipPlaneNz(double);                      // clip plane normal

  double GetCameraDistance();                       // get camera distance
  void SetCameraDistance(double);                   // set camera distance
  double GetCameraPositionX();                      // get camera position
  double GetCameraPositionY();                      // get camera position
  double GetCameraPositionZ();                      // get camera position
  void SetCameraPositionX(double);                  // set camera position
  void SetCameraPositionY(double);                  // set camera position
  void SetCameraPositionZ(double);                  // set camera position
  double GetCameraFocalPointX();                    // get focal point
  double GetCameraFocalPointY();                    // get focal point
  double GetCameraFocalPointZ();                    // get focal point
  void SetCameraFocalPointX(double);                // set focal point
  void SetCameraFocalPointY(double);                // set focal point
  void SetCameraFocalPointZ(double);                // set focal point
  void CameraDolly(double);                         // dolly
  void CameraRoll(double);                          // roll
  void CameraAzimuth(double);                       // azimuth
  void CameraYaw(double);                           // yaw
  void CameraElevation(double);                     // elevation
  void CameraPitch(double);                         // pitch
  void CameraZoom(double);                          // zoom
  void SetCameraRoll(double);                       // set roll
  void SetInitialCameraPosition();                  // set initial position

  void RotateX(double);                             // rotate visible actors
  void RotateY(double);                             // rotate visible actors
  void RotateZ(double);                             // rotate visible actors
  void SetOrientation(double, double, double);      // set orientation
  void SetPositionX(double);                        // set position
  void SetPositionY(double);                        // set position
  void SetPositionZ(double);                        // set position
  void SetPosition(double, double, double);         // set position
  void AddPosition(double, double, double);         // add position
  void SetOrigin(double, double, double);           // set origin
  void SetScaleX(double);                           // set scale
  void SetScaleY(double);                           // set scale
  void SetScaleZ(double);                           // set scale
  void SetScale(double, double, double);            // set scale

  double GetLength();                               // get model size
  double GetNofNodes();                             // get nof nodes
  double GetMinX();                                 // bounding box
  double GetMaxX();                                 // bounding box
  double GetMinY();                                 // bounding box
  double GetMaxY();                                 // bounding box
  double GetMinZ();                                 // bounding box
  double GetMaxZ();                                 // bounding box

  bool SavePngFile(QString);                        // save image file

  bool Execute(QString);                            // Execute ECMA script

  void SetBgColor(double, double, double);          // bg color (rgb)

private slots:
  void exitSlot();
  void readEpFileSlot();
  void showSurfaceDialogSlot();
  void showVectorDialogSlot();
  void showIsoContourDialogSlot();
  void showIsoSurfaceDialogSlot();
  void showColorBarDialogSlot();
  void showStreamLineDialogSlot();
  void showTimeStepDialogSlot();
  void showPreferencesDialogSlot();
  void showTextDialogSlot();

  void drawMeshPointSlot();
  void drawMeshEdgeSlot();
  void drawFeatureEdgesSlot();
  void drawSurfaceSlot();
  void drawVectorSlot();
  void drawIsoContourSlot();
  void drawIsoSurfaceSlot();
  void drawColorBarSlot();
  void drawStreamLineSlot();
  void drawAxesSlot();
  void drawTextSlot();

  void hideSurfaceSlot();
  void hideVectorSlot();
  void hideIsoContourSlot();
  void hideIsoSurfaceSlot();
  void hideColorBarSlot();
  void hideStreamLineSlot();
  void hideTextSlot();

  void groupChangedSlot(QAction*);
  void regenerateGridsSlot();
  void maybeRedrawSlot(bool);
  void populateWidgetsSlot();
  void fitToWindowSlot();
  void resetModelViewSlot();
  void clipAllToggledSlot(bool);

  void savePictureSlot();
  void savePovraySlot();
  void timeStepChangedSlot();
  void reloadPostSlot();

  void showHelpSlot();

#ifdef EG_PYTHONQT
  void showPythonQtConsoleSlot();
#endif

  void showECMAScriptConsoleSlot();
  void evaluateECMAScriptSlot(QString);

private:
  QMenu *fileMenu;
  QMenu *editMenu;
  QMenu *editGroupsMenu;
  QMenu *viewMenu;
  QMenu *helpMenu;

  QToolBar *viewToolBar;

  QAction *regenerateGridsAct;
  QAction *matcAct;
  QAction *exitAct;
  QAction *redrawAct;
  QAction *savePictureAct;
  QAction *savePovrayAct;
  QAction *preferencesAct;
  QAction *drawMeshPointAct;
  QAction *drawMeshEdgeAct;
  QAction *drawFeatureEdgesAct;
  QAction *drawSurfaceAct;
  QAction *drawVectorAct;
  QAction *drawIsoContourAct;
  QAction *drawIsoSurfaceAct;
  QAction *drawColorBarAct;
  QAction *drawStreamLineAct;
  QAction *timeStepAct;
  QAction *fitToWindowAct;
  QAction *resetModelViewAct;
  QAction *drawAxesAct;
  QAction *drawTextAct;
  QAction *reloadPostAct;
  QAction *readEpFileAct;
  QAction *clipAllAct;
  QAction *showHelpAct;

  void createActions();
  void createMenus();
  void createToolbars();
  void createStatusBar();

  EpMesh* epMesh;
  QString postFileName;
  bool postFileRead;
  int scalarFields;
  ScalarField* scalarField;

  void addVectorField(QString, int);
  void getPostLineStream(QTextStream*);

  QHash<QString, QAction*> groupActionHash;

  QVTKWidget* qvtkWidget;

  vtkRenderer* renderer;
  vtkUnstructuredGrid* volumeGrid;
  vtkUnstructuredGrid* surfaceGrid;
  vtkUnstructuredGrid* lineGrid;
  vtkLookupTable *currentLut;
  vtkPlane* clipPlane;
  vtkActor* meshPointActor;
  vtkActor* meshEdgeActor;
  vtkActor* featureEdgeActor;
  vtkActor* surfaceActor;
  vtkActor* vectorActor;
  vtkActor* isoContourActor;
  vtkActor* isoSurfaceActor;
  vtkActor* streamLineActor;
  vtkActor* axesActor;
  vtkActor* pickedPointActor;
  vtkFollower* axesXTextActor;
  vtkFollower* axesYTextActor;
  vtkFollower* axesZTextActor;
  vtkScalarBarActor* colorBarActor;
  vtkImplicitPlaneWidget* planeWidget;
  vtkTextActor* textActor;
  double initialCameraPosition[3];
  double initialCameraRoll;

  Surface* surface;         // ui
  Vector* vector;           // ui
  IsoContour* isoContour;   // ui
  IsoSurface* isoSurface;   // ui
  StreamLine* streamLine;   // ui
  ColorBar* colorBar;       // ui
  Preferences* preferences; // ui
  Matc* matc;               // ui
  TimeStep* timeStep;       // ui
  Axes* axes;               // ui
  Text* text;               // ui
  FeatureEdge* featureEdge; // ui
  MeshPoint* meshPoint;     // ui
  MeshEdge* meshEdge;       // ui
  ReadEpFile* readEpFile;   // ui

  QString currentSurfaceName;
  QString currentVectorName;
  QString currentVectorColorName;
  QString currentIsoContourName;
  QString currentIsoContourColorName;
  QString currentIsoSurfaceName;
  QString currentIsoSurfaceColorName;
  QString currentStreamLineName;
  QString currentStreamLineColorName;
  double currentPickPosition[3];

  // Post file input:
  QString postLine;
  QTextStream postLineStream;

  // bounding box:
  double boundingBoxMinX;
  double boundingBoxMaxX;
  double boundingBoxMinY;
  double boundingBoxMaxY;
  double boundingBoxMinZ;
  double boundingBoxMaxZ;

#ifdef EG_PYTHONQT
  QAction *showPythonQtConsoleAct;
  PythonQtObjectPtr mainContext;
  PythonQtScriptingConsole *console;
#endif

  // ECMAScript
  QAction* showECMAScriptConsoleAct;
  EcmaConsole* ecmaConsole;
  QScriptEngine* engine;
};

#endif // VTKPOST_H
