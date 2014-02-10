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
 *  ElmerGUI text                                                            *
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
#include "text.h"

#include <vtkTextActor.h>
#include <vtkTextProperty.h>

using namespace std;

Text::Text(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);
  setWindowTitle("Text");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));

  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClicked()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));

  red = 0.0;
  green = 0.0;
  blue = 0.0;
}

Text::~Text()
{
}

void Text::applyButtonClicked()
{
  emit(drawTextSignal());
}

void Text::cancelButtonClicked()
{
  emit(hideTextSignal());
  close();
}

void Text::okButtonClicked()
{
  applyButtonClicked();
  close();
}

void Text::draw(VtkPost* vtkPost)
{
  QString message = ui.textEdit->text();
  int posX = ui.posxEdit->text().toInt();
  int posY = ui.posyEdit->text().toInt();
  int size = ui.size->value();
  bool left = ui.left->isChecked();
  bool centered = ui.centered->isChecked();
  bool right = ui.right->isChecked();
  bool bold = ui.bold->isChecked();
  bool italic = ui.italic->isChecked();
  bool shadow = ui.shadow->isChecked();

  vtkTextActor* textActor = vtkPost->GetTextActor();
  if(textActor == NULL) return;

  textActor->SetDisplayPosition(posX, posY);
  textActor->SetInput(message.toAscii().data());

  vtkTextProperty* tprop = textActor->GetTextProperty();
  if(tprop == NULL) return;

  tprop->SetFontFamilyToArial();
  tprop->SetFontSize(size);
  
  if(left) tprop->SetJustificationToLeft();
  if(centered) tprop->SetJustificationToCentered();
  if(right) tprop->SetJustificationToRight();

  if(bold) {
    tprop->BoldOn();
  } else {
    tprop->BoldOff();
  }

  if(italic) {
    tprop->ItalicOn();
  } else {
    tprop->ItalicOff();
  }

  if(shadow) {
    tprop->ShadowOn();
  } else {
    tprop->ShadowOff();
  }

  tprop->SetColor(red, green, blue);
}

void Text::SetMessage(QString message)
{
  ui.textEdit->setText(message);
}

void Text::SetPosX(int x)
{
  ui.posxEdit->setText(QString::number(x));
}

void Text::SetPosY(int y)
{
  ui.posyEdit->setText(QString::number(y));
}

void Text::SetLeft()
{
  ui.left->setChecked(true);
}

void Text::SetCentered()
{
  ui.centered->setChecked(true);
}

void Text::SetRight()
{
  ui.right->setChecked(true);
}

void Text::SetSize(int n)
{
  ui.size->setValue(n);
}

void Text::SetBold(bool b)
{
  ui.bold->setChecked(b);
}

void Text::SetItalic(bool b)
{
  ui.italic->setChecked(b);
}

void Text::SetShadow(bool b)
{
  ui.shadow->setChecked(b);
}

void Text::SetRed(double d)
{
  red = d;
  if(red < 0.0) red = 0.0;
  if(red > 1.0) red = 1.0;
  
}

void Text::SetGreen(double d)
{
  green = d;
  if(green < 0.0) green = 0.0;
  if(green > 1.0) green = 1.0;
  
}

void Text::SetBlue(double d)
{
  blue = d;
  if(blue < 0.0) blue = 0.0;
  if(blue > 1.0) blue = 1.0; 
}

void Text::SetRGB(double r, double g, double b)
{
  this->SetRed(r);
  this->SetGreen(g);
  this->SetBlue(b);
}
