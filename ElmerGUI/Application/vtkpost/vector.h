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
 *  ElmerGUI vector                                                          *
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

#ifndef VECTOR_H
#define VECTOR_H

#include <QWidget>
#include "ui_vector.h"

class VtkPost;
class TimeStep;
class ScalarField;

class Vector : public QDialog
{
  Q_OBJECT

public:
  Vector(QWidget *parent = 0);
  ~Vector();

  Ui::vectorDialog ui;

  void populateWidgets(VtkPost*);
  void draw(VtkPost*, TimeStep*);

public slots:
  QString GetFieldName();                           // get field name
  QString GetColorName();                           // get color name
  QString GetThresholdName();                       // get threshold name
  bool SetFieldName(QString);                       // set field name
  bool SetColorName(QString);                       // set color name
  bool SetThresholdName(QString);                   // set threshold name
  void SetLength(int);                              // set arrow length
  void SetQuality(int);                             // set arrow quality
  void SetEveryNth(int);                            // reduce arrows
  void SetClipPlane(bool);                          // set clipping on/off
  void ComputeNormals(bool);                        // draw better arrows
  void SetRandomMode(bool);                         // randomly reduce data
  void ScaleByMagnitude(bool);                      // scale arrows
  void SetMinColorVal(double);                      // set color min
  void SetMaxColorVal(double);                      // set color max
  void KeepColorLimits(bool);                       // keep color limits
  void SetMinThresholdVal(double);                  // set threshold min
  void SetMaxThresholdVal(double);                  // set threshold max
  void UseThreshold(bool);                          // use threshold
  void KeepThresholdLimits(bool);                   // keep thrshld limits

signals:
  void drawVectorSignal();
  void hideVectorSignal();

private slots:
  void cancelButtonClicked();
  void okButtonClicked();
  void applyButtonClicked();
  void colorSelectionChanged(int);
  void thresholdSelectionChanged(int);
  void keepLimitsSlot(int);
  void keepThresholdLimitsSlot(int);

private:
  ScalarField *scalarField;
  int scalarFields;

};

#endif // VECTOR_H
