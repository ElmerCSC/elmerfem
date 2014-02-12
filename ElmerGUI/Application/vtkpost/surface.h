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
 *  ElmerGUI surface                                                         *
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

#ifndef SURFACE_H
#define SURFACE_H

#include <QWidget>
#include "ui_surface.h"

class ScalarField;
class VtkPost;
class TimeStep;

class Surface : public QDialog
{
  Q_OBJECT

public:
  Surface(QWidget *parent = 0);
  ~Surface();

  Ui::surfaceDialog ui;

  void populateWidgets(VtkPost*);
  void draw(VtkPost*, TimeStep*);

signals:
  void drawSurfaceSignal();
  void hideSurfaceSignal();

public slots:
  QString GetFieldName();                           // get field name
  bool SetFieldName(QString);                       // set field name
  void SetMinVal(double);                           // set minimum
  void SetMaxVal(double);                           // set maximum
  void KeepLimits(bool);                            // keep limits
  void SetComputeNormals(bool);                     // shade model
  void SetFeatureAngle(int);                        // feature angle
  void SetOpacity(int);                             // set opacity
  void SetClipPlane(bool);                          // set clipping

private slots:
  void cancelButtonClicked();
  void okButtonClicked();
  void applyButtonClicked();
  void surfaceSelectionChanged(int);
  void keepLimitsSlot(int);

private:
  ScalarField *scalarField;
  int scalarFields;

};

#endif // SURFACE_H
