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
 *  ElmerGUI streamline                                                      *
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

#ifndef STREAMLINE_H
#define STREAMLINE_H

#include <QWidget>
#include "ui_streamline.h"

class ScalarField;
class VtkPost;
class TimeStep;

class StreamLine : public QDialog
{
  Q_OBJECT

public:
  StreamLine(QWidget *parent = 0);
  ~StreamLine();

  Ui::streamLineDialog ui;

  void populateWidgets(VtkPost*);
  void draw(VtkPost*, TimeStep*);

signals:
  void drawStreamLineSignal();
  void hideStreamLineSignal();

public slots:
  QString GetFieldName();                           // get field name
  QString GetColorName();                           // get color name
  bool SetFieldName(QString);                       // set field name
  bool SetColorName(QString);                       // set color name
  void SetMaxTime(double);                          // max time
  void SetStepLength(double);                       // step length
  void SetThreads(int);                             // nof threads
  void SetIntegStepLength(double);                  // integ. step length
  void UseSurfaceMesh(bool);                        // use forface mesh
  void UseVolumeMesh(bool);                         // use volume mesh
  void IntegrateForwards(bool);                     // integrate forwards
  void IntegrateBackwards(bool);                    // integrate backwards
  void SetMinColorVal(double);                      // color min val
  void SetMaxColorVal(double);                      // color max value
  void KeepColorLimits(bool);                       // keep color limits
  void DrawLines(bool);                             // draw using lines
  void DrawRibbons(bool);                           // draw using ribbons
  void SetLineWidth(int);                           // line width
  void SetRibbonWidth(int);                         // ribbon width
  void UseSphereSource(bool);                       // use sphere source
  void UseLineSource(bool);                         // use line source
  void UsePointSource(bool);                        // use point source
  void SetSphereSourceX(double);                    // sphere origin
  void SetSphereSourceY(double);                    // sphere origin
  void SetSphereSourceZ(double);                    // sphere origin
  void SetSphereSourceRadius(double);               // sphere radius
  void SetSphereSourcePoints(int);                  // nof pts in sphere
  void SetLineSourceStartX(double);                 // line start point
  void SetLineSourceStartY(double);                 // line start point
  void SetLineSourceStartZ(double);                 // line start point
  void SetLineSourceEndX(double);                   // line end point
  void SetLineSourceEndY(double);                   // line end point
  void SetLineSourceEndZ(double);                   // line end point
  void SetLineSourcePoints(int);                    // nof pts on line

private slots:
  void cancelButtonClicked();
  void okButtonClicked();
  void applyButtonClicked();
  void colorSelectionChanged(int);
  void keepLimitsSlot(int);

private:
  ScalarField *scalarField;
  int scalarFields;

};

#endif // STREAMLINE_H
