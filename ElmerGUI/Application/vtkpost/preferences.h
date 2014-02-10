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
 *  ElmerGUI preferences                                                     *
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

#ifndef PREFERENCES_H
#define PREFERENCES_H

#include <QWidget>
#include "ui_preferences.h"

class Preferences : public QDialog
{
  Q_OBJECT

public:
  Preferences(QWidget *parent = 0);
  ~Preferences();

  Ui::preferencesDialog ui;

signals:
  void redrawSignal();

public slots:
  void UseSurfaceMeshForPoints(bool);              // nodes of surface mesh
  void UseVolumeMeshForPoints(bool);               // nodes of volume mesh
  void SetPointSize(int);                          // set node point size
  void SetPointQuality(int);                       // set node point quality
  void UseClipPlaneForPoints(bool);                // clip nodes
  void UseSurfaceMeshForEdges(bool);               // edges of surface mesh
  void UseVolumeMeshForEdges(bool);                // edges of volume mesh
  void UseTubeFilterForEdges(bool);                // use tube filter
  void UseClipPlaneForEdges(bool);                 // clip edges
  void SetLineWidthForEdges(int);                  // edge line width
  void SetTubeQualityForEdges(int);                // edge tube quality
  void SetTubeRadiusForEdges(int);                 // edge tube radius
  void UseSurfaceMeshForFeatureEdges(bool);        // surface mesh: f-edges
  void UseVolumeMeshForFeatureEdges(bool);         // volume mesh: f-edges
  void UseTubeFilterForFeatureEdges(bool);         // use tube filter
  void UseClipPlaneForFeatureEdges(bool);          // clip f-edges
  void DrawBoundaryEdges(bool);                    // draw boundary edges
  int GetFeatureAngle();                           // get feature angle
  void SetFeatureAngle(int);                       // set feature angle
  void SetLineWidthForFeatureEdges(int);           // f-edge line width
  void SetTubeQualityForFeatureEdges(int);         // f-edge tube quality
  void SetTubeRadiusForFeatureEdges(int);          // f-edge tube radius
  void SetClipPlaneOx(double);                     // clip plane origin
  void SetClipPlaneOy(double);                     // clip plane origin
  void SetClipPlaneOz(double);                     // clip plane origin
  void SetClipPlaneNx(double);                     // clip plane normal
  void SetClipPlaneNy(double);                     // clip plane normal
  void SetClipPlaneNz(double);                     // clip plane normal

private slots:
  void okButtonClicked();
  void cancelButtonClicked();
  void applyButtonClicked();

private:

};

#endif // PREFERENCES_H
