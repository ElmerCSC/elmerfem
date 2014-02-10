/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland   *
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
 *  ElmerGUI sifgenerator                                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                  *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 15 Mar 2008                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef SIFGENERATOR_H
#define SIFGENERATOR_H

#include <QTextEdit>
#include <QHash>
#include <QScriptEngine>

#include "meshtype.h"
#include "maxlimits.h"
#include "generalsetup.h"
#include "boundarypropertyeditor.h"
#include "bodypropertyeditor.h"
#include "solverparameters.h"
#include "meshcontrol.h"
#include "dynamiceditor.h"

enum EquationTypes {
  HEAT_EQUATION,
  LINEAR_ELASTICITY,
  NAVIER_STOKES,
  ADVECTION_DIFFUSION,
  HELMHOLTZ_EQUATION 
};

class SifGenerator  {
 public:
  SifGenerator();
  ~SifGenerator();

  void setMesh(mesh_t *mesh);
  void setTextEdit(QTextEdit *textEdit);
  void setDim(int dim);
  void setCdim(int cdim);
  void setElmerDefs(QDomDocument *doc);
  void setGeneralSetup(GeneralSetup *setup);
  void setEquationEditor(const QVector<DynamicEditor*> &d);
  void setMaterialEditor(const QVector<DynamicEditor*> &d);
  void setBodyForceEditor(const QVector<DynamicEditor*> &d);
  void setInitialConditionEditor(const QVector<DynamicEditor*> &d);
  void setBoundaryConditionEditor(const QVector<DynamicEditor*> &d);
  void setSolverParameterEditor(const QVector<SolverParameterEditor*> &d);
  void setBoundaryPropertyEditor(const QVector<BoundaryPropertyEditor*> &d);
  void setBodyPropertyEditor(const QVector<BodyPropertyEditor*> &d);
  void setMeshControl(MeshControl *meshControl);
  void setLimit(Limit *limit);

  void makeHeaderBlock();
  void makeSimulationBlock();
  void makeConstantsBlock();
  void makeBodyBlocks();
  void makeEquationBlocks();
  void makeSolverBlocks(const QString &name);
  void makeMaterialBlocks();
  void makeBodyForceBlocks();
  void makeInitialConditionBlocks();
  void makeBoundaryBlocks();
  
  QHash<int, int> bodyMap;
  QHash<int, int> boundaryMap;

 private:
  mesh_t* mesh;
  QTextEdit* te;
  int dim, cdim;
  QDomDocument* elmerDefs;
  GeneralSetup* generalSetup;
  QVector<DynamicEditor*> equationEditor;
  QVector<DynamicEditor*> materialEditor;
  QVector<DynamicEditor*> bodyForceEditor;
  QVector<DynamicEditor*> initialConditionEditor;
  QVector<DynamicEditor*> boundaryConditionEditor;
  QVector<SolverParameterEditor*> solverParameterEditor;
  QVector<BoundaryPropertyEditor*> boundaryPropertyEditor;
  QVector<BodyPropertyEditor*> bodyPropertyEditor;
  MeshControl *meshControl;
  Limit *limit;

  int  findHashValue(DynamicEditor*, const QString&, const QString&);
  bool parseSolverSpecificTab(DynamicEditor *, const QString&);
  void parseExecSolverTab(Ui::solverParameterEditor);
  void parseNumericalTechniquesTab(Ui::solverParameterEditor);
  void parseSteadyStateTab(Ui::solverParameterEditor);
  void parseNonlinearSystemTab(Ui::solverParameterEditor);
  void parseLinearSystemTab(Ui::solverParameterEditor);
  void parseParallelTab(Ui::solverParameterEditor);
  void addSifLine(const QString&, const QString&);
  void addSifLineBool(const QString&, bool);
  void handleBCLineEdit(const QDomElement&, QWidget*, const QMap<int, int>&);
  void handleLineEdit(const QDomElement&, QWidget*);
  void handleComboBox(const QDomElement&, QWidget*);
  void handleCheckBox(const QDomElement&, QWidget*);
  void handleTextEdit(const QDomElement&, QWidget*);
};

#endif // SIFGENERATOR_H
