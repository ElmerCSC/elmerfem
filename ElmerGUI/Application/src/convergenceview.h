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
 *  ElmerGUI convergenceview                                                 *
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

#ifndef CONVERGENCEVIEW_H
#define CONVERGENCEVIEW_H

#include <QMainWindow>
#include <QHash>
#include <QIcon>
#include "maxlimits.h"

#include <qwt_plot.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
/*#include <qwt_data.h> <-- deprecated in Qwt6, using qwt_compat.h instead*/
#include <qwt_compat.h>
#include <qwt_text.h>
#include <qwt_scale_engine.h>

#define MAX_NOF_PENS 25

class CurveData
{
public:
  CurveData();
  
  void append(double*, double*, int);
  
  int count() const;
  int size() const;
  const double *x() const;
  const double *y() const;
  
private:
  int d_count;
  QwtArray<double> d_x;
  QwtArray<double> d_y;
};

class Curve
{
public:
  CurveData *d_data;
  QwtPlotCurve *d_curve;
};

class ConvergenceView : public QMainWindow
{
  Q_OBJECT

public:
  ConvergenceView(Limit *limit, QWidget *parent = 0);
  ~ConvergenceView();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  void appendData(double, QString);
  void appendData(double*, int, QString);
  void removeData();

  QString title;

private slots:
  void savePictureSlot();
  void showGridSlot();
  void showLegendSlot();
  void showNSHistorySlot();
  void showSSHistorySlot();
  void clearHistorySlot();

private:
  QwtPlot *plot;
  QwtPlotGrid *grid;
  QwtLegend *legend;
  QwtLog10ScaleEngine *scaleEngine;

  QHash<QString, Curve*> curveList;
  QPen *pen;

  QAction *savePictureAct;
  QAction *exitAct;
  QAction *showGridAct;
  QAction *showLegendAct;
  QAction *showNSHistoryAct;
  QAction *showSSHistoryAct;
  QAction *clearHistoryAct;

  QMenu *fileMenu;
  QMenu *viewMenu;

  QToolBar *fileToolBar;
  QToolBar *viewToolBar;

  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();  

  bool showGrid;
  bool showLegend;
  bool showNSHistory;
  bool showSSHistory;

  QIcon iconChecked;
};

#endif // CONVERGENCEVIEW_H
