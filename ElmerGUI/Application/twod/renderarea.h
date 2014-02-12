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
 *  ElmerGUI RenderArea                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                   *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 15 Mar 2008                                               *
 *                                                                           *
 *****************************************************************************/
#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QWidget>
#include <QHash>

class CurveEditor;
class QTableWidget;

class Spline
{
 public:
  int out;
  int in;
  int np;
  int p[3];
};

class RenderArea : public QWidget
{
 Q_OBJECT

 public:
  RenderArea(QWidget *parent = 0);
  ~RenderArea();

  void setCurveEditor(CurveEditor *curveEditor);
  void modifyPoint(int idx, double x, double y);
  void modifyCurve(int idx, int in, int out, int np, int p1, int p2, int p3);
  void updatePoints(QTableWidget *table);
  void updateCurves(QTableWidget *table);

 public slots:
  void fitSlot();
  void readSlot(QString fileName);
  void saveSlot(QString fileName);
  void drawPointsSlot(bool state);
  void drawSplinesSlot(bool state);
  void drawTangentsSlot(bool state);
  void drawPointNumbersSlot(bool state);
  void drawSplineNumbersSlot(bool state);
  void drawMaterialNumbersSlot(bool state);

 signals:
  void statusMessage(QString message);

 protected:
  void paintEvent(QPaintEvent *event);  
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);

 private:
  QHash<int, QPointF> points;
  QHash<int, Spline> splines;
  QVector<int> bodies;
  QRectF viewport;
  QRectF renderport;
  int selectedPoint;  
  int pointRadius;
  QPointF mapToViewport(QPointF point) const;
  QPointF mapToRenderport(QPointF point) const;
  QPointF quadNurbs(double u, QPointF P0, QPointF P1, QPointF P2) const;
  QPoint lastPos;
  bool drawPoints;
  bool drawSplines;
  bool drawTangents;
  bool drawPointNumbers;
  bool drawSplineNumbers;
  bool drawMaterialNumbers;
  CurveEditor *curveEditor;
  bool reading;
};

#endif //  RENDERAREA_H
