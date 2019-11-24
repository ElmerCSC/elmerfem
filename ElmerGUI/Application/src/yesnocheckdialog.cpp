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
 *  ElmerGUI yesnotcheckdialog.h                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Takayuki Saeki                                                  *
 *  Original Date: 10 Nov 2019                                               *
 *                                                                           *
 *****************************************************************************/

#include <QtGui>
#include <iostream>
#include "yesnocheckdialog.h"

using namespace std;

YesNoCheckDialog::YesNoCheckDialog(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);
  
  connect(ui.yesButton, SIGNAL(clicked()), this, SLOT(yesSlot()));
  connect(ui.noButton, SIGNAL(clicked()), this, SLOT(noSlot()));
		   
  ui.noButton->setDefault(true);
}

YesNoCheckDialog::~YesNoCheckDialog()
{
}

Qt::CheckState YesNoCheckDialog::checkState()
{
  return ui.checkBox->checkState();
}

void YesNoCheckDialog::yesSlot()
{
  accept();
}

void YesNoCheckDialog::noSlot()
{
  reject();
}