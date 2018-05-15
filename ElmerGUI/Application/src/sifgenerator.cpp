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
 *  ElmerGUI sifgenerator                                                    *
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

#include "sifgenerator.h"
#include <iostream>
#include <QDebug>

using namespace std;

SifGenerator::SifGenerator()
{
}

SifGenerator::~SifGenerator()
{
}

void SifGenerator::setMesh(mesh_t* m)
{
  this->mesh = m;
}

void SifGenerator::setTextEdit(QTextEdit* textEdit)
{
  this->te = textEdit;
}

void SifGenerator::setDim(int n)
{
  this->dim = n;
}

void SifGenerator::setCdim(int n)
{
  this->cdim = n;
}

void SifGenerator::setElmerDefs(QDomDocument* d)
{
  this->elmerDefs = d;
}

void SifGenerator::setGeneralSetup(GeneralSetup* g)
{
  this->generalSetup = g;
}

void SifGenerator::setEquationEditor(const QVector<DynamicEditor*>& d)
{
  this->equationEditor = d;
}

void SifGenerator::setMaterialEditor(const QVector<DynamicEditor*>& d)
{
  this->materialEditor = d;
}

void SifGenerator::setBodyForceEditor(const QVector<DynamicEditor*>& d)
{
  this->bodyForceEditor = d;
}

void SifGenerator::setInitialConditionEditor(const QVector<DynamicEditor*>& d)
{
  this->initialConditionEditor = d;
}

void SifGenerator::setBoundaryConditionEditor(const QVector<DynamicEditor*>& d)
{
  this->boundaryConditionEditor = d;
}

void SifGenerator::setSolverParameterEditor(const QVector<SolverParameterEditor*>& s)
{
  this->solverParameterEditor = s;
}

void SifGenerator::setBoundaryPropertyEditor(const QVector<BoundaryPropertyEditor*>& b)
{
  this->boundaryPropertyEditor = b;
}

void SifGenerator::setBodyPropertyEditor(const QVector<BodyPropertyEditor*>& b)
{
  this->bodyPropertyEditor = b;
}

void SifGenerator::setMeshControl(MeshControl* m)
{
  this->meshControl = m;
}

void SifGenerator::setLimit(Limit* l)
{
  this->limit = l;
}

// Make Header-block:
//-----------------------------------------------------------------------------
void SifGenerator::makeHeaderBlock()
{
  Ui::setupDialog ui = generalSetup->ui;

  te->append("Header");
  
  if(ui.checkKeywordsWarn->isChecked())
    te->append("  CHECK KEYWORDS Warn");

  QString qs1 = ui.meshDBEdit1->text().trimmed();
  QString qs2 = ui.meshDBEdit2->text().trimmed(); 
  QString qs3 = ui.includePathEdit->text().trimmed(); 
  QString qs4 = ui.resultsDirectoryEdit->text().trimmed(); 
  QString qs5 = ui.headerFreeTextEdit->toPlainText();

  te->append("  Mesh DB \"" +  qs1 + "\" \"" + qs2 + "\"");
  te->append("  Include Path \"" + qs3 + "\"");
  te->append("  Results Directory \"" + qs4 + "\"");

  if(!qs5.isEmpty())
    te->append(qs5);

  te->append("End\n");
}

// Make Simulation-block:
//-----------------------------------------------------------------------------
void SifGenerator::makeSimulationBlock()
{
  Ui::setupDialog ui = generalSetup->ui;

  te->append("Simulation");

  addSifLine("  Max Output Level = ", 
	     ui.maxOutputLevelCombo->currentText().trimmed());
  addSifLine("  Coordinate System = ",
	     ui.coordinateSystemCombo->currentText().trimmed());
  addSifLine("  Coordinate Mapping(3) = ",
	     ui.coordinateMappingEdit->text().trimmed());
  addSifLine("  Simulation Type = ", 
	     ui.simulationTypeCombo->currentText().trimmed());
  addSifLine("  Steady State Max Iterations = ",
	     ui.steadyStateMaxIterEdit->text().trimmed());
  addSifLine("  Output Intervals = ",
	     ui.outputIntervalsEdit->text().trimmed());
  addSifLine("  Timestepping Method = ",
	     ui.timesteppingMethodCombo->currentText().trimmed());
  addSifLine("  BDF Order = ",
	     ui.bdfOrderCombo->currentText().trimmed());
  addSifLine("  Timestep intervals = ",
	     ui.timeStepIntervalsEdit->text().trimmed());
  addSifLine("  Timestep Sizes = ",
	     ui.timestepSizesEdit->text().trimmed());

  addSifLine("  Coordinate Scaling = ",
	     ui.coordinateScalingEdit->text().trimmed());
  addSifLine("  Angular Frequency = ",
	     ui.angularFrequencyEdit->text().trimmed());
    
  addSifLine("  Solver Input File = ", 
	     ui.solverInputFileEdit->text().trimmed());
  addSifLine("  Post File = ", 
	     ui.postFileEdit->text().trimmed());

  QString qs = ui.simulationFreeTextEdit->toPlainText();

  if(!qs.isEmpty())
    te->append(qs);

  te->append("End\n");
}

// Make Constants-block:
//-----------------------------------------------------------------------------
void SifGenerator::makeConstantsBlock()
{
  Ui::setupDialog ui = generalSetup->ui;

  te->append("Constants");
  
  addSifLine("  Gravity(4) = ",
	     ui.gravityEdit->text().trimmed());
  addSifLine("  Stefan Boltzmann = ",
	     ui.stefanBoltzmannEdit->text().trimmed());
  addSifLine("  Permittivity of Vacuum = ",
	     ui.vacuumPermittivityEdit->text().trimmed());
  addSifLine("  Boltzmann Constant = ",
	     ui.boltzmannEdit->text().trimmed());
  addSifLine("  Unit Charge = ",
	     ui.unitChargeEdit->text().trimmed());

  QString qs = ui.constantsFreeTextEdit->toPlainText();

  if(!qs.isEmpty())
    te->append(qs);

  te->append("End\n");
}


// Make Body-blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeBodyBlocks()
{
  int i;

  int sifIndex = 0, maxOriginalIndex=-1;

  for(int index = 0; index < bodyMap.count(); index++) {

    if(index >= bodyPropertyEditor.size()) {
      cout << "SifGenerator: Body index out of bounds" << endl;
      continue;
    }

    BodyPropertyEditor *bodyEdit = bodyPropertyEditor[index];

    if(!bodyEdit)
      continue;
    
    int originalIndex = bodyMap.key(index);
    maxOriginalIndex = max(maxOriginalIndex, originalIndex );

    if(bodyEdit->touched) {
      te->append("Body " + QString::number(++sifIndex));

      te->append("  Target Bodies(1) = " + QString::number(originalIndex));

      if ( bodyEdit->ui.nameEdit->text().trimmed() == "" )
        te->append("  Name = \"Body " + QString::number(sifIndex) + "\"");
      else
        te->append("  Name = \"" + bodyEdit->ui.nameEdit->text().trimmed() + "\"");

      i = bodyEdit->ui.equationCombo->currentIndex();
      if(i > 0)
	te->append("  Equation = " + QString::number(i));
      
      i = bodyEdit->ui.materialCombo->currentIndex();
      if(i > 0)
	te->append("  Material = " + QString::number(i));
      
      i = bodyEdit->ui.bodyForceCombo->currentIndex();
      if(i > 0)
	te->append("  Body Force = " + QString::number(i));
      
      i = bodyEdit->ui.initialConditionCombo->currentIndex();
      if(i > 0)
	te->append("  Initial condition = " + QString::number(i));
      
      te->append("End\n");      
    }
  }

  for( int index = 0; index < boundaryPropertyEditor.size(); index++ )
  {
    if(!boundaryPropertyEditor[index])
      continue;

    BodyPropertyEditor *bodyEdit = boundaryPropertyEditor[index]->bodyProperties;

    if(!bodyEdit)
      continue;

    if(bodyEdit && bodyEdit->touched ) {
      te->append("Body " + QString::number(++sifIndex));

      boundaryPropertyEditor[index]->bodyID = ++maxOriginalIndex;

      te->append("  Target Bodies(1) = " + QString::number(maxOriginalIndex));

      if ( bodyEdit->ui.nameEdit->text().trimmed() == "" )
        te->append("  Name = \"Body " + QString::number(sifIndex) + "\"");
      else
        te->append("  Name = \"" + bodyEdit->ui.nameEdit->text().trimmed() + "\"");

      i = bodyEdit->ui.equationCombo->currentIndex();
      if(i > 0)
	te->append("  Equation = " + QString::number(i));
      
      i = bodyEdit->ui.materialCombo->currentIndex();
      if(i > 0)
	te->append("  Material = " + QString::number(i));
      
      i = bodyEdit->ui.bodyForceCombo->currentIndex();
      if(i > 0)
	te->append("  Body Force = " + QString::number(i));
      
      i = bodyEdit->ui.initialConditionCombo->currentIndex();
      if(i > 0)
	te->append("  Initial condition = " + QString::number(i));
      
      te->append("End\n");      
    }
  }
}


int SifGenerator::findHashValue(DynamicEditor *de, const QString &sname, const QString &name)
{
    for(int i = 0; i < de->hash.count(); i++) {
      hash_entry_t entry = de->hash.values().at(i); 

      QWidget *widget = entry.widget;
      if (widget->isEnabled()) {
        QString key = de->hash.keys().at(i);
        QStringList keySplitted = key.split("/");
        QString solverName = keySplitted.at(1).trimmed();
        QString labelName = keySplitted.at(3).trimmed();

        QDomElement elem = entry.elem;
        if ( solverName==sname && labelName==name &&
	     elem.attribute("Widget","")=="Edit" ) {
	  QLineEdit *line = static_cast<QLineEdit*>(widget);
	  QString str = line->text().trimmed();
	  if ( str.isEmpty() ) return 0;
	  return str.toInt();
	}
      }
    }
    return 0;
}

// Make Equation/Solver -blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeEquationBlocks()
{
  // Enumerate solvers && write solver blocks:
  //-------------------------------------------
  QMap<QString, int> numberForSolver;
  numberForSolver.clear();

  int solverNumber = 0;

  for(int index = 0; index < equationEditor.size(); index++) {
    DynamicEditor *eqEditor = equationEditor[index];

    if(eqEditor->menuAction != NULL) {
      for(int i = 0; i < eqEditor->hash.count(); i++) {
	hash_entry_t entry = eqEditor->hash.values().at(i); 

	QWidget *widget = entry.widget;
        if (widget->isEnabled()) {
          QDomElement elem = entry.elem;

	  QString key = eqEditor->hash.keys().at(i);
	  QStringList keySplitted = key.split("/");
	  QString labelName  = keySplitted.at(3).trimmed();
	  QString solverName = keySplitted.at(1).trimmed();

	  if(labelName=="Active" && elem.attribute("Widget", "")=="CheckBox") {
	    QCheckBox *checkBox = static_cast<QCheckBox*>(widget);
	    if(checkBox->isChecked()) {
	      if(!numberForSolver.contains(solverName)) {
                int pri = findHashValue( eqEditor, solverName, "Priority");
                numberForSolver.insert(solverName, pri);
	      }
	    }
	  }
	}
      }
    }
  }

  // Sort and enumerate solvers according to their priority:
  //---------------------------------------------------------
  QList<QPair<int, QString> > tmpList;

  foreach(const QString &key, numberForSolver.keys()) {
    int value = numberForSolver.value(key);
    tmpList << qMakePair(value, key);
  }

  qSort(tmpList);

  numberForSolver.clear();

  int n = tmpList.count();

  for(int i = 0; i < n; i++) {
    const QPair<int, QString> &pair = tmpList[i];
    const QString &key = pair.second;
    numberForSolver.insert(key, n-i);
  }

  // Generate solver blocks:
  //-------------------------
  QMap<int, int> handled;

  for(int index = 0; index < equationEditor.size(); index++) {
    DynamicEditor *eqEditor = equationEditor[index];
    if(eqEditor->menuAction != NULL) {
      for(int i = 0; i < eqEditor->hash.count(); i++) {
	hash_entry_t entry = eqEditor->hash.values().at(i); 
	
	QWidget *widget = entry.widget;
        if (widget->isEnabled()) {
	  QString key = eqEditor->hash.keys().at(i);
	  QStringList keySplitted = key.split("/");
	  QString solverName = keySplitted.at(1).trimmed();
	  QString labelName = keySplitted.at(3).trimmed();
          QDomElement elem = entry.elem;

	  if(labelName=="Active" && elem.attribute("Widget", "")=="CheckBox") {
	    QCheckBox *checkBox = static_cast<QCheckBox*>(widget);
	    if(checkBox->isChecked()) {
	      solverNumber = numberForSolver.value(solverName);
	      if((solverNumber>0) && (handled[solverNumber]==0)) {
		handled[solverNumber] = 1; 
		te->append("Solver " + QString::number(solverNumber));
		te->append("  Equation = " + solverName);
		makeSolverBlocks(solverName);
		te->append("End");
		te->append("");
	      }
	    }
	  }
	}
      }
    }
  }
  
  // Generate equation blocks:
  //---------------------------
  QMap<int, bool> solverActive;

  int sifIndex = 0;
  for(int index = 0; index < equationEditor.size(); index++) {
    DynamicEditor *eqEditor = equationEditor[index];

    if(eqEditor->menuAction != NULL) {
      te->append("Equation " + QString::number(++sifIndex));
      
      QString name = eqEditor->nameEdit->text().trimmed();
      te->append("  Name = \"" + name + "\"");

      QString solverString = "";
      int nofSolvers = 0;

      for( int i=0; i < solverParameterEditor.size(); i++ )
        solverActive[i] = false;

      for(int i = 0; i < eqEditor->hash.count(); i++) {
	hash_entry_t entry = eqEditor->hash.values().at(i); 
	QWidget *widget = entry.widget;

        if(widget->isEnabled()) {
          QDomElement elem = entry.elem;
	  QString key = eqEditor->hash.keys().at(i);
	  QStringList keySplitted = key.split("/");	  
	  QString solverName = keySplitted.at(1).trimmed();
	  QString labelName = keySplitted.at(3).trimmed();

	  // solver active?
	  if((labelName == "Active") && (elem.attribute("Widget", "") == "CheckBox")) {
	    QCheckBox *checkBox = static_cast<QCheckBox*>(widget);
	    if(checkBox->isChecked()) {
	      nofSolvers++;
	      solverNumber = numberForSolver.value(solverName);
              solverActive[solverNumber] = true;
              solverString += " " + QString::number(solverNumber);
	    }
	  }
        }
      }

      for(int i = 0; i < eqEditor->hash.count(); i++) {
	hash_entry_t entry = eqEditor->hash.values().at(i); 
	QWidget *widget = entry.widget;

        if(widget->isEnabled()) {
          QDomElement elem = entry.elem;
	  QString key = eqEditor->hash.keys().at(i);
	  QStringList keySplitted = key.split("/");	  
	  QString solverName = keySplitted.at(1).trimmed();
	  QString labelName = keySplitted.at(3).trimmed();

          solverNumber = numberForSolver.value(solverName);
          if ( !solverActive[solverNumber] ) continue;

          if((elem.attribute("Widget", "") == "CheckBox") &&
	     (labelName != "Active")) 
	    handleCheckBox(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Edit" &&
             labelName != "Priority")
	    handleLineEdit(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Combo")
	    handleComboBox(elem, widget);

	  if(elem.attribute("Widget", "") == "TextEdit")
	    handleTextEdit(elem, widget);
        }
      }

      if(nofSolvers > 0)
	te->append("  Active Solvers(" 
		   + QString::number(nofSolvers) 
		   + ") =" + solverString);
      
      te->append("End\n");
    }
  }
}

//-------------------------------------------------------------------------
void SifGenerator::makeSolverBlocks(const QString &solverName)
{
  SolverParameterEditor *spe, *tmp;
  Ui::solverParameterEditor ui;
  
  bool found = false;
  int current=-1;

  for(int i = 0; i < solverParameterEditor.size(); i++) {
    spe = solverParameterEditor[i];

    if(!spe)
      continue;

    QString currentName = spe->solverName.trimmed();
    if(currentName == solverName) {
      found = true;
      current = i;
      break;
    }
  }

  if(!found) {
    tmp = new SolverParameterEditor;
  } else {
    tmp = spe;
  }

  if ( !tmp->generalOptions ) {
    tmp->generalOptions = new DynamicEditor;
    tmp->generalOptions->setupTabs(elmerDefs, "Solver", current );
  }

  bool hasMatrix = parseSolverSpecificTab(tmp->generalOptions, solverName);

  ui = tmp->ui; 

  // Parse the exec solver also for non-PDE solvers
  parseExecSolverTab(ui);

  if(hasMatrix) {
    parseNumericalTechniquesTab(ui);
    parseSteadyStateTab(ui);
    parseNonlinearSystemTab(ui);
    parseLinearSystemTab(ui);
    parseParallelTab(ui);
    // todo: add adaptivity & multigrid
  }

  if(!found)
  {
    delete tmp->generalOptions;
    delete tmp;
  }
}

// Make Material-blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeMaterialBlocks()
{
  int sifIndex = 0;

  for(int index = 0; index < materialEditor.size(); index++) {
    DynamicEditor *matEditor = materialEditor[index];
    
    if(matEditor->menuAction != NULL) {      
      te->append("Material " + QString::number(++sifIndex));
      
      QString name = matEditor->nameEdit->text().trimmed();
      te->append("  Name = \"" + name + "\"");
      
      for(int i = 0; i < matEditor->hash.count(); i++) {
	hash_entry_t entry = matEditor->hash.values().at(i); 
	
	QWidget *widget = entry.widget;

	QDomElement elem;
        if ( widget->isEnabled() ) {
          elem = entry.elem;
	  
          if(elem.attribute("Widget", "") == "CheckBox") 
	   handleCheckBox(elem, widget);
	
	 if(elem.attribute("Widget", "") == "Edit")
	   handleLineEdit(elem, widget);
	
	 if(elem.attribute("Widget", "") == "Combo")
	   handleComboBox(elem, widget);

	 if(elem.attribute("Widget", "") == "TextEdit")
	   handleTextEdit(elem, widget);
        }
      }
      te->append("End\n");
    }
  }
}


// Make body force blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeBodyForceBlocks()
{
  int sifIndex = 0;

  for(int index = 0; index < bodyForceEditor.size(); index++) {
    DynamicEditor *bfEdit = bodyForceEditor[index];
    
    if(bfEdit->menuAction != NULL) { 
      te->append("Body Force " + QString::number(++sifIndex));
      
      QString name = bfEdit->nameEdit->text().trimmed();
      te->append("  Name = \"" + name + "\"");
      
      for(int i = 0; i < bfEdit->hash.count(); i++) {
	hash_entry_t entry = bfEdit->hash.values().at(i); 
	
	QWidget *widget = entry.widget;

        if ( widget->isEnabled() ) {
          QDomElement elem = entry.elem;
	  
          if(elem.attribute("Widget", "") == "CheckBox") 
	    handleCheckBox(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Edit")
	    handleLineEdit(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Combo")
	    handleComboBox(elem, widget);

	  if(elem.attribute("Widget", "") == "TextEdit")
	    handleTextEdit(elem, widget);
        }
      }
      te->append("End\n");
    }
  }
}


// Make initial condition blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeInitialConditionBlocks()
{
  int sifIndex = 0;
  
  for(int index = 0; index < initialConditionEditor.size(); index++) {
    DynamicEditor *icEdit = initialConditionEditor[index];
    
    if(icEdit->menuAction != NULL) { 
      te->append("Initial Condition " + QString::number(++sifIndex));
      
      QString name = icEdit->nameEdit->text().trimmed();
      te->append("  Name = \"" + name + "\"");
      
      for(int i = 0; i < icEdit->hash.count(); i++) {
	hash_entry_t entry = icEdit->hash.values().at(i); 
	
	QWidget *widget = entry.widget;
	
        if ( widget->isEnabled() ) {
          QDomElement elem = entry.elem;
	  
          if(elem.attribute("Widget", "") == "CheckBox") 
	    handleCheckBox(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Edit")
	    handleLineEdit(elem, widget);
	  
	  if(elem.attribute("Widget", "") == "Combo")
	    handleComboBox(elem, widget);

	  if(elem.attribute("Widget", "") == "TextEdit")
	    handleTextEdit(elem, widget);
        }
      }
      te->append("End\n");
    }
  }
}



// Make boundary blocks:
//-----------------------------------------------------------------------------
void SifGenerator::makeBoundaryBlocks()
{
    int sifIndex = 0;
    int bcnum = 0;
    int diff = 0;
    QMap<int, int> boundaryBC;
    QMap<QString, int> boundaryList; /*(Boundary condition, edge) value pairs */
    QList<QString> boundaryConditions; /* List of different boundary conditions */
    QList<int> boundaryEdges; /* List of edges relating to some specific boundary condition */
    QString tmp;

    boundaryBC.clear();
    boundaryList.clear();
    boundaryConditions.clear();
    boundaryEdges.clear();
    tmp.clear();

    //Find the available boundary conditions
    for (int index = 0; index < boundaryConditionEditor.count(); index++) {
            DynamicEditor *bc = boundaryConditionEditor[index];
            boundaryConditions.append(bc->nameEdit->text().trimmed());
    }

    //Find the boundary conditions and edges related to them.
    for (int index = 0; index < boundaryConditions.count(); index++) {
        for(int k = 0; k < boundaryMap.count(); k++) {
            BoundaryPropertyEditor *bEdit = boundaryPropertyEditor[k];
            if(bEdit->touched) {
                boundaryBC[index] = ++sifIndex;
                int originalIndex = boundaryMap.key(k);
                const QString bcname = bEdit->ui.boundaryConditionCombo->currentText().trimmed();
                if (boundaryConditions.value(index) == bcname)
                    boundaryList.insertMulti(bcname, originalIndex);
            }
        }
    }
    //qDebug() << "boundaryMap: " << boundaryMap;
    //qDebug() << "boundaryList: " << boundaryList;
    qDebug() << "boundaryConditions: " << boundaryConditions;

    //Arrange and sort boundary conditions
    for(int index = 0; index < boundaryConditions.count(); index++) {
        tmp.clear();
        boundaryEdges.clear();
        const QString name = boundaryConditions[index];
        BoundaryPropertyEditor *bEdit = boundaryPropertyEditor[index];
        DynamicEditor *bc = boundaryConditionEditor[index];

        if(boundaryList.contains(name)) {
            bcnum++;
            te->append("Boundary Condition " + QString::number(bcnum));
            if(boundaryConditions.count() > 1) {
                QMap <QString,int>::ConstIterator l = boundaryList.find(boundaryConditions[index]);
                while (l != boundaryList.end() && l.key()==boundaryConditions[index]){
                    boundaryEdges.append(l.value());
                    l++;
                    }
                while ((l--) != boundaryList.begin() && l.key()==boundaryConditions[index]){
                    tmp.append(QString::number(l.value()));
                    tmp.append(" ");
                    }
            }
            if(boundaryConditions.count() <= 1) {
                QMap <QString,int>::ConstIterator l = boundaryList.begin();
                while (l != boundaryList.end()) {
                    boundaryEdges.append(l.value());
                    l++;
                    }
                while ((l--) != boundaryList.begin()){
                    tmp.append(QString::number(l.value()));
                    tmp.append(" ");
                    }
                }

            te->append("  Target Boundaries("
                + QString::number(boundaryEdges.count())
                + ") = " + tmp);

            if ( bEdit->bodyProperties ) {
                te->append("  Body id = " + QString::number(bEdit->bodyID) );
                }

            te->append("  Name = \"" + name + "\"");

            // check which one of the dynamic editors has "name" typed in nameEdit:
            for(int j = 0; j < boundaryConditionEditor.size(); j++) {
                DynamicEditor *bc = boundaryConditionEditor[j];
                if(bc->menuAction != NULL) {
                    if(bc->nameEdit->text().trimmed() == name && name != NULL) {

                    // go through the hash of this dynamic editor:
                    //--------------------------------------------
                    for(int i = 0; i < bc->hash.count(); i++) {
                        hash_entry_t entry = bc->hash.values().at(i);
	      
                        QWidget *widget = entry.widget;
	      
                        QDomElement elem;
                        if ( widget->isEnabled() ) {
                            elem = entry.elem;

                        if(elem.attribute("Widget", "") == "CheckBox")
                            handleCheckBox(elem, widget);
		
                        if(elem.attribute("Widget", "") == "Edit")
                            handleBCLineEdit(elem, widget, boundaryBC);

                        if(elem.attribute("Widget", "") == "Combo")
                            handleComboBox(elem, widget);

                        if(elem.attribute("Widget", "") == "TextEdit")
                            handleTextEdit(elem, widget);
                            }
                        }
                    }
                }
            }
            te->append("End\n");
        }
    }
}

// Parse "Solver specific tab"
//-----------------------------------------------------------------------------
bool SifGenerator::parseSolverSpecificTab(DynamicEditor *solEditor, const QString &solverName)
{
  // Returns true if there is a matrix involved. otherwise returns false.
  if ( !solEditor ) return false;

  bool hasMatrix = true;

  QScriptEngine engine; 

  QScriptValue dim_QSV  = QScriptValue(&engine,dim);
  engine.globalObject().setProperty( "dim", dim_QSV );

  QScriptValue cdim_QSV = QScriptValue(&engine,cdim);
  engine.globalObject().setProperty( "cdim", cdim_QSV );

  for(int i = 0; i < solEditor->hash.count(); i++) {
    hash_entry_t entry = solEditor->hash.values().at(i);

    QString key = solEditor->hash.keys().at(i);
    QStringList keySplitted = key.split("/");	  
    QString tabName   = keySplitted.at(1).trimmed();
    QString labelName = keySplitted.at(3).trimmed();

    if ( tabName != solverName ) continue;

    // Has matrix?
    if(labelName == "No Matrix Equation") {
      if(entry.elem.attribute("Widget", "") == "CheckBox") {
	QCheckBox *cb = static_cast<QCheckBox*>(entry.widget);
	hasMatrix = !cb->isChecked();
      }
    }

    // variable names handled separately...
    // ------------------------------------
    if ( labelName=="Variable" || labelName.mid(0,17)=="Exported Variable" ) {
      if( entry.elem.attribute("Widget", "") != "Edit") continue;

      QLineEdit *l = static_cast<QLineEdit*>(entry.widget);
      QString varName = l->text().trimmed();

      if ( varName == "" ) continue;

      int dofs=1;
      QStringList dofsplit = varName.split("[");
      if ( dofsplit.count()>1 ) {
        varName = dofsplit.at(0).trimmed() + "[";
        QString dof = dofsplit.at(1).trimmed();
        dof = dof.split("]").at(0).trimmed();

        dofsplit = dof.split(":");
        QString subVarName = dofsplit.at(0).trimmed();
        for( int i=1; i<dofsplit.count(); i++)
        {
          dof = dofsplit.at(i).trimmed();

          QStringList subDofSplit = dof.split(" ");
          QString subDof = subDofSplit.at(0).trimmed();

          dofs = engine.evaluate(subDof).toInt32();
          if (i>1) varName = varName + " ";
          varName = varName + subVarName + ":" + QString::number(dofs);

          if ( subDofSplit.count() > 1 )
            subVarName = subDofSplit.at(1).trimmed();
        }
        varName = varName + "]";
        addSifLine( "  " + labelName + " = ", varName );
      } else {
        dofsplit = varName.split("(");
        if ( dofsplit.count()>1 ) {
          varName = dofsplit.at(0).trimmed();
          QString dof = dofsplit.at(1).trimmed();
          dofsplit = dof.split(")");
          dof = dofsplit.at(0).trimmed();
          dofs = engine.evaluate(dof).toInt32();
        }
	// Don't write the the trivial dof==1 case as this leaves possibility to define the number of
	// dofs internally within the solver. 
        if ( dofs <= 1 ) {
	  addSifLine( "  "+labelName+" = ", varName );
	  dofs = 1;
	}
	else {
	  addSifLine( "  "+labelName+" = -dofs ",  QString::number(dofs) + " " + varName );
	}
      }
      continue;
    }

    QWidget *widget = entry.widget;

    QDomElement elem;
    if ( widget->isEnabled() ) {
      elem = entry.elem;

      if(elem.attribute("Widget", "") == "CheckBox")
       handleCheckBox(elem, widget);

     if(elem.attribute("Widget", "") == "Edit")
       handleLineEdit(elem, widget);

     if(elem.attribute("Widget", "") == "Combo")
       handleComboBox(elem, widget);

     if(elem.attribute("Widget", "") == "TextEdit")
       handleTextEdit(elem, widget);
    }
  }

  return hasMatrix;
}

// Parse "Exec Solver" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseExecSolverTab(Ui::solverParameterEditor ui)
{
  if(ui.execAlways->isChecked())
    te->append("  Exec Solver = Always");
  
  if(ui.execBeforeSimulation->isChecked())
    te->append("  Exec Solver = Before Simulation");
  
  if(ui.execAfterSimulation->isChecked())
    te->append("  Exec Solver = After Simulation");
  
  if(ui.execBeforeTimestep->isChecked())
    te->append("  Exec Solver = Before Timestep");
  
  if(ui.execAfterTimestep->isChecked())
    te->append("  Exec Solver = After Timestep");

  if(ui.execBeforeSaving->isChecked())
    te->append("  Exec Solver = Before Saving");

  if(ui.execAfterSaving->isChecked())
    te->append("  Exec Solver = After Saving");
  
  if(ui.execNever->isChecked())
    te->append("  Exec Solver = Never");
}

// Parse "Numerical Techniques" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseNumericalTechniquesTab(Ui::solverParameterEditor ui)
{
  addSifLineBool("  Stabilize = ", ui.stabilizeCheck->isChecked());
  addSifLineBool("  Bubbles = ", ui.bubblesCheck->isChecked());
  addSifLineBool("  Lumped Mass Matrix = ", ui.lumpedMassCheck->isChecked());
  addSifLineBool("  Optimize Bandwidth = ", ui.optimizeBandwidthCheck->isChecked());
}


// Parse "Steady state" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseSteadyStateTab(Ui::solverParameterEditor ui)
{
  if(ui.steadyStateConvergenceToleranceEdit->text() == "") {
    cout << "Steady state convergence tolerance is undefined - aborting" << endl;
    return;
  }
  
  addSifLine("  Steady State Convergence Tolerance = ",
	      ui.steadyStateConvergenceToleranceEdit->text());

  if( ui.steadyStateConvergenceMeasureCombo->currentText().trimmed() != "Norm") 
    addSifLine("  Steady State Convergence Measure = ",
	       ui.steadyStateConvergenceMeasureCombo->currentText().trimmed());
}


// Parse "Nonlinear system" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseNonlinearSystemTab(Ui::solverParameterEditor ui)
{
  addSifLine("  Nonlinear System Convergence Tolerance = ",
	      ui.nonlinSystemConvergenceToleranceEdit->text());
  
  addSifLine("  Nonlinear System Max Iterations = ", 
	      ui.nonlinSystemMaxIterationEdit->text());
  
  addSifLine("  Nonlinear System Newton After Iterations = ",
	      ui.nonlinSystemNewtonAfterIterEdit->text());
  
  addSifLine("  Nonlinear System Newton After Tolerance = ", 
	      ui.nonlinSystemNewtonAfterTolEdit->text());
  
  addSifLine("  Nonlinear System Relaxation Factor = ", 
	      ui.nonlinSystemRelaxationFactorEdit->text());

  if( ui.nonlinSystemConvergenceMeasureCombo->currentText().trimmed() != "Norm") 
    addSifLine("  Nonlinear System Convergence Measure = ",
	       ui.nonlinSystemConvergenceMeasureCombo->currentText().trimmed());

}


// Parse "Parallel" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseParallelTab(Ui::solverParameterEditor ui)
{
  if(ui.useHypre->isChecked()) {
    addSifLine("  Linear System Use HYPRE = ", "True");
    
    if(ui.useParasails->isChecked()) {
      addSifLine("  Linear System Preconditioning = ", "ParaSails");
      
      addSifLine("  ParaSails Threshold = ",
		 ui.thresholdEdit->text().trimmed());
      
      addSifLine("  ParaSails Filter = ",
		 ui.filterEdit->text().trimmed());
      
      addSifLine("  ParaSails MaxLevel = ",
		 ui.maxLevelEdit->text().trimmed());
      
      addSifLine("  ParaSails Symmetry = ",
		 ui.symmetryEdit->text().trimmed());
    }
    
    if(ui.useBoomerAMG->isChecked()) {
      addSifLine("  Linear System Preconditioning = ", "BoomerAMG");

      addSifLine("  BoomerAMG Relax Type = ",
		 QString::number(ui.boomerRelaxation->currentIndex()));

      addSifLine("  BoomerAMG Coarsen Type = ",
		 QString::number(ui.boomerCoarsening->currentIndex()));

      addSifLine("  BoomerAMG Num Sweeps = ",
		 ui.boomerSweeps->text().trimmed());

      addSifLine("  BoomerAMG Max Levels = ",
		 ui.boomerMaxLevels->text().trimmed());

      addSifLine("  BoomerAMG Interpolation = ",
		 QString::number(ui.boomerInterpolation->currentIndex()));

      addSifLine("  BoomerAMG Smooth Type = ",
		 QString::number(ui.boomerSmoother->currentIndex()));

      addSifLine("  BoomerAMG Cycle Type = ",
		 QString::number(ui.boomerCycle->currentIndex()));

    }
  }
}


// Parse "Linear system" tab from ui to sif:
//-----------------------------------------------------------------------------
void SifGenerator::parseLinearSystemTab(Ui::solverParameterEditor ui)
{
  bool hyprePreconditioning 
    = ui.useParasails->isChecked() | ui.useBoomerAMG->isChecked();

  if(ui.linearSystemSolverDirect->isChecked()) {
    
    addSifLine("  Linear System Solver = ", "Direct");
    
    addSifLine("  Linear System Direct Method = ",
		ui.linearSystemDirectMethod->currentText());
    
  } else if(ui.linearSystemSolverIterative->isChecked()) {
    
    addSifLine("  Linear System Solver = ", "Iterative");
    
    addSifLine("  Linear System Iterative Method = ",
		ui.linearSystemIterativeMethod->currentText());
    
    addSifLine("  Linear System Max Iterations = ", 
		ui.linearSystemMaxIterationsEdit->text());
    
    addSifLine("  Linear System Convergence Tolerance = ",
		ui.linearSystemConvergenceToleranceEdit->text());

    addSifLine("  BiCGstabl polynomial degree = ",
               ui.linearSystemBiCGstablPolDeg->text());

    if(!hyprePreconditioning)
      addSifLine("  Linear System Preconditioning = ",
		 ui.linearSystemPreconditioning->currentText());
    
    addSifLine("  Linear System ILUT Tolerance = ",
		ui.linearSystemILUTToleranceEdit->text());
    
    addSifLineBool("  Linear System Abort Not Converged = ",
		ui.linearSystemAbortWhenNotConvergedCheck->isChecked());
    
    addSifLine("  Linear System Residual Output = ",
		ui.linearSystemResiduaOutputEdit->text());
    
    addSifLine("  Linear System Precondition Recompute = ",
		ui.linearSystemPreconditionRecomputeEdit->text());
    
  } else if(ui.linearSystemSolverMultigrid->isChecked()) {
    
    addSifLine("  Linear System Solver = ", "Multigrid");
    
    // TODO: rest of the less common params etc.
  }
}

//------------------------------------------------------------------------
//
//                         COMMON UTILITY FUNCTIONS
//
//------------------------------------------------------------------------

void SifGenerator::addSifLine(const QString &var, const QString &val)
{
  if(val != "")
    te->append(var + val);
}

void SifGenerator::addSifLineBool(const QString &var, bool val)
{
  if(val == true)
    te->append(var + "True");
  else
    te->append(var + "False");
}


void SifGenerator::handleBCLineEdit(const QDomElement &elem, QWidget *widget, const QMap<int, int> &boundaryBC)
{
  QString name = elem.firstChildElement("SifName").text().trimmed();
  if( name == "" )
    name= elem.firstChildElement("Name").text().trimmed();

  QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
  QString value = lineEdit->text().trimmed();

  if ( name=="Periodic BC" && value != "" ) 
  {
     int val = value.toInt(); 
     val = boundaryMap.value(val);
     value=QString::number(boundaryBC[val]);
  }

  addSifLine("  " + name + " = ", value);
}

void SifGenerator::handleLineEdit(const QDomElement &elem, QWidget *widget)
{
  QString name = elem.firstChildElement("SifName").text().trimmed();
  if( name == "" )
    name= elem.firstChildElement("Name").text().trimmed();

  QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
  QString value = lineEdit->text().trimmed();
  addSifLine("  " + name + " = ", value);
}

void SifGenerator::handleTextEdit(const QDomElement &elem, QWidget *widget)
{
  QTextEdit *textEdit = static_cast<QTextEdit*>(widget);
  QString value = textEdit->toPlainText();
  if(!value.isEmpty()) te->append(value);
}

void SifGenerator::handleComboBox(const QDomElement &elem, QWidget *widget)
{  
  QString name = elem.firstChildElement("SifName").text().trimmed();
  if( name == "" )
    name= elem.firstChildElement("Name").text().trimmed();

  QComboBox *comboBox = static_cast<QComboBox*>(widget);
  QString value = comboBox->currentText().trimmed();

  if(value != "None")
    addSifLine("  " + name + " = ", value);
}

void SifGenerator::handleCheckBox(const QDomElement &elem, QWidget *widget)
{
  QString name = elem.firstChildElement("SifName").text().trimmed();
  if( name == "" )
    name = elem.firstChildElement("Name").text().trimmed();
  
  QString def_val = elem.firstChildElement("DefaultValue").text().trimmed();
  if ( def_val == "" )
    def_val = "False";

  QCheckBox *checkBox = static_cast<QCheckBox*>(widget);
  
  if(checkBox->isChecked()) {
    if ( def_val != "True" )
      te->append("  " + name + " = True");
  } else {
    if ( def_val != "False" )
      te->append("  " + name + " = False");
  }
}
