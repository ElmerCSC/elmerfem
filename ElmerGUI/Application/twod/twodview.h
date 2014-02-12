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
 *  ElmerGUI TwodView                                                        *
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
#ifndef TWODVIEW_H
#define TWODVIEW_H

#include <QMainWindow>

class RenderArea;
class CurveEditor;
class QAction;
class QMenu;

class TwodView : public QMainWindow
{
Q_OBJECT

public:
  TwodView(QWidget *parent = 0);
  ~TwodView();

public slots:
  void statusMessage(QString message);
  void openSlot();
  void saveSlot();
  void helpSlot();
  void addPointSlot();
  void addCurveSlot();
  void deletePointSlot();
  void deleteCurveSlot();

private:
  void createActions();
  void createMenus();
  void createStatusBar();

  RenderArea *renderArea;
  CurveEditor *curveEditor;
  QAction *openAction;
  QAction *saveAction;
  QAction *quitAction;
  QAction *addPointAction;
  QAction *addCurveAction;
  QAction *deletePointAction;
  QAction *deleteCurveAction;
  QAction *fitAction;
  QAction *drawPointsAction;
  QAction *drawSplinesAction;
  QAction *drawTangentsAction;
  QAction *drawPointNumbersAction;
  QAction *drawSplineNumbersAction;
  QAction *drawMaterialNumbersAction;
  QAction *helpAction;
  QMenu *fileMenu;
  QMenu *editMenu;
  QMenu *viewMenu;
  QMenu *helpMenu;
};

#endif // TWODVIEW_H
