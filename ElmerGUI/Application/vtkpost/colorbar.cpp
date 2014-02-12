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
 *  ElmerGUI colorbar                                                        *
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
#include <iostream>
#include "vtkpost.h"
#include "colorbar.h"

#include <vtkActor.h>
#include <vtkMapper.h>
#include <vtkScalarBarActor.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkActor.h>

using namespace std;

ColorBar::ColorBar(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));
  connect(ui.colorCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(colorSelectionChanged(int)));

  setWindowTitle("Colorbar");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

ColorBar::~ColorBar()
{
}

void ColorBar::cancelButtonClicked()
{
  emit(hideColorBarSignal());
  close();
}

void ColorBar::okButtonClicked()
{
  applyButtonClicked();
  close();
}

void ColorBar::applyButtonClicked()
{
  emit(drawColorBarSignal());
}

void ColorBar::colorSelectionChanged(int newIndex)
{
}

void ColorBar::populateWidgets(VtkPost* vtkPost)
{
  ui.colorCombo->clear();
  ui.colorCombo->addItem("Surface");
  ui.colorCombo->addItem("Vector");
  ui.colorCombo->addItem("Isocontour");
  ui.colorCombo->addItem("Isosurface");
  ui.colorCombo->addItem("Streamline");
}

void ColorBar::draw(VtkPost* vtkPost)
{
  vtkScalarBarActor* colorBarActor = vtkPost->GetColorBarActor();

  vtkTextMapper* tMapper = vtkTextMapper::New();
  colorBarActor->SetMapper(tMapper);
  tMapper->Delete();

  QString actorName = ui.colorCombo->currentText().trimmed();

  if(actorName.isEmpty()) return;

  QString fieldName = "";

  vtkScalarsToColors* lut = NULL;

  if(actorName == "Surface") {
    fieldName = vtkPost->GetCurrentSurfaceName();
    if(fieldName.isEmpty()) return;
    lut = vtkPost->GetSurfaceActor()->GetMapper()->GetLookupTable();
  }

  if(actorName == "Vector") {
    fieldName = vtkPost->GetCurrentVectorName();
    if(fieldName.isEmpty()) return;
    lut = vtkPost->GetVectorActor()->GetMapper()->GetLookupTable();
  }

  if(actorName == "Isocontour") {
    fieldName = vtkPost->GetCurrentIsoContourName();
    if(fieldName.isEmpty()) return;
    lut = vtkPost->GetIsoContourActor()->GetMapper()->GetLookupTable();
  }

  if(actorName == "Isosurface") {
    fieldName = vtkPost->GetCurrentIsoSurfaceName();
    if(fieldName.isEmpty()) return;
    lut = vtkPost->GetIsoSurfaceActor()->GetMapper()->GetLookupTable();
  }
  
  if(actorName == "Streamline") {
    fieldName = vtkPost->GetCurrentStreamLineName();
    if(fieldName.isEmpty()) return;
    lut = vtkPost->GetStreamLineActor()->GetMapper()->GetLookupTable();
  }

  if(!lut) return;

  colorBarActor->SetLookupTable(lut);

  bool horizontal = ui.horizontalRButton->isChecked();
  bool annotate = ui.annotateBox->isChecked();
  int labels = ui.labelsSpin->value();
  double width = ui.widthEdit->text().toDouble();
  double height = ui.heightEdit->text().toDouble();

  if(width < 0.01) width = 0.01;
  if(width > 1.00) width = 1.00;
  if(height < 0.01) height = 0.01;
  if(height > 1.00) height = 1.00;

  colorBarActor->SetPosition(0.05, 0.05);

  if(horizontal) {
    colorBarActor->SetOrientationToHorizontal();
    colorBarActor->SetWidth(height);
    colorBarActor->SetHeight(width);
  } else {
    colorBarActor->SetOrientationToVertical();
    colorBarActor->SetWidth(width);
    colorBarActor->SetHeight(height);
  }
  
  colorBarActor->SetNumberOfLabels(labels);

  colorBarActor->GetLabelTextProperty()->SetFontSize(16);
  colorBarActor->GetLabelTextProperty()->SetFontFamilyToArial();
  colorBarActor->GetLabelTextProperty()->BoldOn();
  colorBarActor->GetLabelTextProperty()->ItalicOn();
  colorBarActor->GetLabelTextProperty()->SetColor(0, 0, 1);
  
  colorBarActor->GetTitleTextProperty()->SetFontSize(16);
  colorBarActor->GetTitleTextProperty()->SetFontFamilyToArial();
  colorBarActor->GetTitleTextProperty()->BoldOn();
  colorBarActor->GetTitleTextProperty()->ItalicOn();
  colorBarActor->GetTitleTextProperty()->SetColor(0, 0, 1);
  
  if(annotate) {
    colorBarActor->SetTitle(fieldName.toAscii().data());
  } else {
    colorBarActor->SetTitle("");
  }
}

// Public slots:
//---------------
bool ColorBar::SetFieldName(QString name)
{
  for(int i = 0; i < ui.colorCombo->count(); i++) {
    if(ui.colorCombo->itemText(i) == name) {
      ui.colorCombo->setCurrentIndex(i);
      return true;
    }
  }
  return false;
}

void ColorBar::UseHorizontalLayout(bool b)
{
  ui.horizontalRButton->setChecked(b);
}

void ColorBar::UseVerticalLayout(bool b)
{
  ui.verticalRButton->setChecked(b);
}

void ColorBar::AnnotateFieldName(bool b)
{
  ui.annotateBox->setChecked(b);
}

void ColorBar::SetLabels(int n)
{
  ui.labelsSpin->setValue(n);
}

void ColorBar::SetLineWidth(double f)
{
  ui.widthEdit->setText(QString::number(f));
}

void ColorBar::SetLength(double f)
{
  ui.heightEdit->setText(QString::number(f));
}
