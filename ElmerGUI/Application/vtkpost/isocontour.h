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
 *  ElmerGUI isocontour                                                      *
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

#ifndef ISOCONTOUR_H
#define ISOCONTOUR_H

#include <QWidget>
#include "ui_isocontour.h"

class ScalarField;
class VtkPost;
class TimeStep;

class IsoContour : public QDialog
{
  Q_OBJECT

public:
  IsoContour(QWidget *parent = 0);
  ~IsoContour();

  Ui::isoContourDialog ui;

  void populateWidgets(VtkPost*);
  void draw(VtkPost*, TimeStep*);

signals:
  void drawIsoContourSignal();
  void hideIsoContourSignal();

public slots:
  QString GetFieldName();                           // get field name    
  QString GetColorName();                           // get color name
  bool SetFieldName(QString);                       // set field name
  bool SetColorName(QString);                       // set color name
  void SetMinFieldVal(double);                      // set min value
  void SetMaxFieldVal(double);                      // set max value
  void SetContours(int);                            // set nof contours
  void SetContourValues(QString);                   // set contour values
  void KeepFieldLimits(bool);                       // keep limits
  void SetMinColorVal(double);                      // set color min
  void SetMaxColorVal(double);                      // set color max
  void KeepColorLimits(bool);                       // keep color limits
  void UseTubeFilter(bool);                         // use tube filter
  void UseClipPlane(bool);                          // set clipping on/off
  void SetLineWidth(int);                           // set line width
  void SetTubeQuality(int);                         // set tube quality
  void SetTubeRadius(int);                          // set tube radius

private slots:
  void cancelButtonClicked();
  void okButtonClicked();
  void applyButtonClicked();
  void contoursSelectionChanged(int);
  void colorSelectionChanged(int);
  void keepContourLimitsSlot(int);
  void keepColorLimitsSlot(int);

private:
  ScalarField *scalarField;
  int scalarFields;

};

#endif // ISOCONTOUR_H
