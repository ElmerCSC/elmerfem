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

#include <QtGui>
#include <iostream>
#include "convergenceview.h"

using namespace std;

//-----------------------------------------------------------------------------
CurveData::CurveData()
  : d_count(0)
{
}

void CurveData::append(double *x, double *y, int count)
{
  int newSize = ((d_count + count) / 1000 + 1) * 1000;
  if(newSize > size()) {
    d_x.resize(newSize);
    d_y.resize(newSize);
  }
  
  for(register int i = 0; i < count; i++) {
    d_x[d_count + i] = x[i];
    d_y[d_count + i] = y[i];
  }
  
  d_count += count;
}

int CurveData::count() const
{
  return d_count;
}

int CurveData::size() const
{
  return d_x.size();
}

const double *CurveData::x() const
{
  return d_x.data();
}

const double *CurveData::y() const
{
  return d_y.data();
}

//-----------------------------------------------------------------------------
ConvergenceView::ConvergenceView(Limit *limit, QWidget *parent)
  : QMainWindow(parent)
{
  iconChecked = QIcon(":/icons/dialog-ok.png");

  plot = new QwtPlot(this);
  plot->resize(200, 200);

  plot->setAutoReplot(true);
  plot->setCanvasBackground(QColor(Qt::white));
  title = "Nonlinear system convergence";
  plot->setTitle(title);

  // legend
  legend = new QwtLegend;
  legend->setFrameStyle(QFrame::Box|QFrame::Sunken);
  plot->insertLegend(legend, QwtPlot::RightLegend);
  
  // grid
  grid = new QwtPlotGrid;
  grid->enableXMin(true);
#if QWT_VERSION >= 0x060100
  grid->setMajorPen(QPen(Qt::black, 0, Qt::DotLine));
  grid->setMinorPen(QPen(Qt::gray, 0 , Qt::DotLine));
#else
  grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
  grid->setMinPen(QPen(Qt::gray, 0 , Qt::DotLine));
#endif 
  grid->attach(plot);

  // axis
  plot->setAxisTitle(QwtPlot::xBottom, "Iteration step");
  plot->setAxisTitle(QwtPlot::yLeft, "Relative change");
  plot->setAxisMaxMajor(QwtPlot::xBottom, 20);
  plot->setAxisMaxMinor(QwtPlot::xBottom, 1);
  plot->setAxisMaxMajor(QwtPlot::yLeft, 10);
  plot->setAxisMaxMinor(QwtPlot::yLeft, 10);

  // scale engine
  scaleEngine = new QwtLog10ScaleEngine;
  plot->setAxisScaleEngine(QwtPlot::yLeft, scaleEngine);

  // default pen
  pen = new QPen[MAX_NOF_PENS];
  for(int i = 0; i < MAX_NOF_PENS; i++)
    pen[i] = QPen(Qt::black);

  // available colors
  pen[0] = QPen(Qt::red);
  pen[1] = QPen(Qt::green);
  pen[2] = QPen(Qt::blue);
  pen[3] = QPen(Qt::black);
  pen[4] = QPen(Qt::cyan);
  pen[5] = QPen(Qt::yellow);

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  setCentralWidget(plot);
  setWindowTitle("Convergence monitor");

  showLegend = true;
  showGrid = true;
  showNSHistory = true;
  showSSHistory = false;

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

ConvergenceView::~ConvergenceView()
{
#if 0
  curveList.clear();
  delete plot;
#endif
}


void ConvergenceView::appendData(double y, QString name)
{
  appendData(&y, 1, name);
}

void ConvergenceView::appendData(double *y, int size, QString name)
{
  QStringList nameSplitted = name.trimmed().split("/");

  bool NS = false;
  if(nameSplitted.at(0) == "NS")
    NS = true;

  bool SS = false;
  if(nameSplitted.at(0) == "SS")
    SS = true;

  Curve *curve = curveList.value(name, NULL);
  
  if(curve == NULL) {
    curve = new Curve;
    curve->d_data = new CurveData;    
    curve->d_curve = new QwtPlotCurve(name);
//    curve->d_curve->setRenderHint(QwtPlotItem::RenderAntialiased);
    QPen currentPen = pen[curveList.count()];
    currentPen.setWidth(2);
    curve->d_curve->setPen(currentPen);

    if(NS && showNSHistory)
      curve->d_curve->attach(plot);

    if(SS && showSSHistory)
      curve->d_curve->attach(plot);

    curveList.insert(name, curve);
  }

  double x = (double)(curve->d_data->count());
  curve->d_data->append(&x, y, size);
  /*curve->d_curve->setRawData(curve->d_data->x(), 
			     curve->d_data->y(), 
			     curve->d_data->count());*/
  curve->d_curve->setRawSamples(curve->d_data->x(), 
			     curve->d_data->y(), 
			     curve->d_data->count());
  plot->setTitle(title);
}

void ConvergenceView::removeData()
{
  for( int i = 0; i < curveList.count(); i++) {
    Curve *curve = curveList.values().at(i);
    delete curve->d_data;
    curve->d_data = NULL;
    curve->d_curve->detach();
    delete curve->d_curve;
    curve->d_curve = NULL;
  }
  curveList.clear();  
  title = "Convergence history";
  /*  plot->clear(); */
  plot->detachItems();
  plot->replot();
}

QSize ConvergenceView::minimumSizeHint() const
{
  return QSize(64, 64);
}


QSize ConvergenceView::sizeHint() const
{
  return QSize(480, 640);
}

void ConvergenceView::createActions()
{
  savePictureAct = new QAction(QIcon(""), tr("&Save picture as..."), this);
  savePictureAct->setStatusTip("Save picture in jpg-format");
  connect(savePictureAct, SIGNAL(triggered()), this, SLOT(savePictureSlot())); 

  exitAct = new QAction(QIcon(":/icons/application-exit.png"), tr("&Quit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  exitAct->setStatusTip("Quit convergence monitor");
  connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

  showGridAct = new QAction(QIcon(iconChecked), tr("&Grid"), this);
  showGridAct->setStatusTip("Show grid");
  connect(showGridAct, SIGNAL(triggered()), this, SLOT(showGridSlot()));  

  showLegendAct = new QAction(QIcon(iconChecked), tr("&Legend"), this);
  showLegendAct->setStatusTip("Show legend");
  connect(showLegendAct, SIGNAL(triggered()), this, SLOT(showLegendSlot())); 

  showNSHistoryAct = new QAction(QIcon(iconChecked), 
				 tr("&Nonlinear system"), this);
  showNSHistoryAct->setStatusTip("Show nonlinear system convergence history");
  connect(showNSHistoryAct, SIGNAL(triggered()), this, SLOT(showNSHistorySlot())); 

  showSSHistoryAct = new QAction(QIcon(""), tr("&Steady state"), this);
  showSSHistoryAct->setStatusTip("Show steady state convergence history");
  connect(showSSHistoryAct, SIGNAL(triggered()), this, SLOT(showSSHistorySlot())); 

  clearHistoryAct = new QAction(QIcon(""), tr("&Clear history"), this);
  clearHistoryAct->setStatusTip("Clear current convergence history");
  connect(clearHistoryAct, SIGNAL(triggered()), this, SLOT(clearHistorySlot())); 
}

void ConvergenceView::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(savePictureAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);
  
  viewMenu = menuBar()->addMenu(tr("&View"));
  viewMenu->addAction(showGridAct);
  viewMenu->addAction(showLegendAct);
  viewMenu->addSeparator();
  viewMenu->addAction(showNSHistoryAct);
  viewMenu->addAction(showSSHistoryAct);
  viewMenu->addSeparator();
  viewMenu->addAction(clearHistoryAct);
}

void ConvergenceView::createToolBars()
{
#if 0
  viewToolBar = addToolBar(tr("&View"));
  viewToolBar->addAction(showGridAct);
  viewToolBar->addAction(showLegendAct);
#endif
}

void ConvergenceView::createStatusBar()
{
  statusBar()->showMessage("Ready");
}

void ConvergenceView::savePictureSlot()
{
  QPixmap pixmap = QPixmap::grabWidget(plot);
  
  QString fileName = QFileDialog::getSaveFileName(this,
	  tr("Save picture"), "", tr("Picture files (*.bmp *.jpg *.png *.pbm *.pgm *.ppm)"));

  if(fileName == "")
    return;

  QFileInfo fi(fileName);
  QString suffix = fi.suffix();
  suffix.toUpper();

#if WITH_QT5
  bool ok = pixmap.save(fileName, suffix.toLatin1(), 95); // fixed quality
#else
  bool ok = pixmap.save(fileName, suffix.toAscii(), 95); // fixed quality
#endif

  if(!ok) {
    cout << "Failed writing picture" << endl;
    cout.flush();
  }
}

void ConvergenceView::showGridSlot()
{
  showGrid = !showGrid;
  
  if(showGrid) {
    grid->attach(plot);
    showGridAct->setIcon(QIcon(iconChecked));
  } else {
    grid->detach();
    showGridAct->setIcon(QIcon(""));
  }

  plot->replot();
}

void ConvergenceView::showLegendSlot()
{
  showLegend = !showLegend;
  
  if(showLegend) {
    legend = new QwtLegend;
    legend->setFrameStyle(QFrame::Box|QFrame::Sunken);
    plot->insertLegend(legend, QwtPlot::RightLegend);
    showLegendAct->setIcon(QIcon(iconChecked));
  } else {
    delete legend;
    showLegendAct->setIcon(QIcon(""));
  }

  plot->replot();
}

void ConvergenceView::showNSHistorySlot()
{
  showNSHistory = !showNSHistory;
  
  if(showNSHistory) {
    showNSHistoryAct->setIcon(QIcon(iconChecked));
  } else {
    showNSHistoryAct->setIcon(QIcon(""));
  }

  // attach/detach curves:
  for(int i = 0; i < curveList.count(); i++) {
    QString key = curveList.keys().at(i);
    Curve *curve = curveList.values().at(i);
    QStringList keySplitted = key.trimmed().split("/");
    if(keySplitted.at(0) == "NS") {
      if(showNSHistory)
	curve->d_curve->attach(plot);
      else
	curve->d_curve->detach();
    }
  }

  plot->replot();
}

void ConvergenceView::showSSHistorySlot()
{
  showSSHistory = !showSSHistory;
  
  if(showSSHistory) {
    showSSHistoryAct->setIcon(QIcon(iconChecked));
  } else {
    showSSHistoryAct->setIcon(QIcon(""));
  }

  // attach/detach curves:
  for(int i = 0; i < curveList.count(); i++) {
    QString key = curveList.keys().at(i);
    Curve *curve = curveList.values().at(i);
    QStringList keySplitted = key.trimmed().split("/");
    if(keySplitted.at(0) == "SS") {
      if(showSSHistory)
	curve->d_curve->attach(plot);
      else
	curve->d_curve->detach();
    }
  }

  plot->replot();
}

void ConvergenceView::clearHistorySlot()
{
  removeData();
}
