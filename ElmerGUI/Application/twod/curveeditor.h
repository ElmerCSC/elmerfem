/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland   *
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
 *  ElmerGUI CurveEditor                                                     *
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
#ifndef CURVEEDITOR_H
#define CURVEEDITOR_H

#include <QTabWidget>

class QTableWidget;
class RenderArea;

class CurveEditor : public QTabWidget
{
Q_OBJECT

public:
  CurveEditor(QWidget *parent = 0);
  ~CurveEditor();

  void setRenderArea(RenderArea *renderArea);
  void addPoint(int idx, double x, double y);
  void addCurve(int in, int out, int pts, int *p);
  void modifyPoint(int idx, double x, double y);
  void clearAll();
  void addPoint();
  void addCurve();
  void deletePoint();
  void deleteCurve();

signals:
  void statusMessage(QString message);

private slots:
  void cCellChanged(int row, int col);
  void pCellChanged(int row, int col);

private:
  QTableWidget *pTable;
  QTableWidget *cTable;
  RenderArea *renderArea;

};

#endif // CURVEEDITOR_H
