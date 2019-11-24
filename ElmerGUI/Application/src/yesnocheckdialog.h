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

#ifndef YESNOCHECKDIALOG_H
#define YESNOCHECKDIALOG_H

#include <QWidget>
#include <QDomDocument>
#include "projectio.h"
#include "dynamiceditor.h"
#include "yesnocheckdialog.h"
#include "ui_yesnocheckdialog.h"

class YesNoCheckDialog : public QDialog
{
  Q_OBJECT

public:
  YesNoCheckDialog(QWidget *parent = 0);
  ~YesNoCheckDialog();
  Qt::CheckState checkState();

  Ui::yesNoCheckDialog ui;

private:
  ProjectIO projectIO;
  
private slots:
  void yesSlot();
  void noSlot();

};

#endif // YESNOCHECKDIALOG_H