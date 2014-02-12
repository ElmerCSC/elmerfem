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
 *  ElmerClips                                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly                                                      *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 14 Nov 2010                                               *
 *                                                                           *
 *****************************************************************************/
#include "preview.h"

Preview::Preview(QWidget *parent) : QLabel(parent), currentProgress(0)
{
  setWindowTitle("ElmerClips");

  setWindowIcon(QIcon(":/img/ElmerClips.ico"));

  setAlignment(Qt::AlignCenter);

  setMinimumSize(400, 400);

  connect(&encoder, SIGNAL(information(const QString &)),
	  this, SLOT(information(const QString &)),
	  Qt::BlockingQueuedConnection);

  connect(&encoder, SIGNAL(progress(int)),
	  this, SLOT(progress(int)),
	  Qt::BlockingQueuedConnection);

  setupContextMenu();
}

void Preview::setupContextMenu()
{
  smallAction = new QAction(QIcon(""), "Small (width 640 pixels)", this);
  smallAction->setCheckable(true);
  smallAction->setChecked(true);

  mediumAction = new QAction(QIcon(""), "Medium (width 720 pixels)", this);
  mediumAction->setCheckable(true);
  mediumAction->setChecked(true);

  bigAction = new QAction(QIcon(""), "Big (width 1280 pixels)", this);
  bigAction->setCheckable(true);
  bigAction->setChecked(true);

  hugeAction = new QAction(QIcon(""), "Huge (width 1920 pixels)", this);
  hugeAction->setCheckable(true);
  hugeAction->setChecked(true);

  resolutionActionGroup = new QActionGroup(this);
  resolutionActionGroup->addAction(smallAction);
  resolutionActionGroup->addAction(mediumAction);
  resolutionActionGroup->addAction(bigAction);
  resolutionActionGroup->addAction(hugeAction);
  resolutionActionGroup->setExclusive(false);;  

  lowQualityAction = new QAction(QIcon(""), "Low (qmax = 16)", this);
  lowQualityAction->setCheckable(true);
  lowQualityAction->setChecked(false);

  mediumQualityAction = new QAction(QIcon(""), "Medium (qmax = 8)", this);
  mediumQualityAction->setCheckable(true);
  mediumQualityAction->setChecked(false);

  highQualityAction = new QAction(QIcon(""), "High (qmax = 4)", this);
  highQualityAction->setCheckable(true);
  highQualityAction->setChecked(true);

  bestQualityAction = new QAction(QIcon(""), "Best (qmax = 2)", this);
  bestQualityAction->setCheckable(true);
  bestQualityAction->setChecked(false);

  qualityActionGroup = new QActionGroup(this);
  qualityActionGroup->addAction(lowQualityAction);
  qualityActionGroup->addAction(mediumQualityAction);
  qualityActionGroup->addAction(highQualityAction);
  qualityActionGroup->addAction(bestQualityAction);
  qualityActionGroup->setExclusive(true);

  quitAction = new QAction(QIcon(""), "Quit", this);
  connect(quitAction, SIGNAL(triggered()), this, SLOT(quitSlot()));

  resolutionMenu = new QMenu("Resolution", this);
  resolutionMenu->addAction(smallAction);
  resolutionMenu->addAction(mediumAction);
  resolutionMenu->addAction(bigAction);
  resolutionMenu->addAction(hugeAction);

  qualityMenu = new QMenu("Quality", this);
  qualityMenu->addAction(lowQualityAction);
  qualityMenu->addAction(mediumQualityAction);
  qualityMenu->addAction(highQualityAction);
  qualityMenu->addAction(bestQualityAction);

  aboutAction = new QAction(QIcon(""), "About...", this);
  connect(aboutAction, SIGNAL(triggered()), this, SLOT(aboutSlot()));

  contextMenu = new QMenu(this);
  contextMenu->addMenu(resolutionMenu);
  contextMenu->addMenu(qualityMenu);
  contextMenu->addAction(aboutAction);
  contextMenu->addAction(quitAction);
}

void Preview::checkCommandLine()
{
  QList<QUrl> urls;

  if(qApp->arguments().count() < 2) {
    showInfo();

  } else {
    resolutionMenu->setEnabled(false);
    qualityMenu->setEnabled(false);
    setAcceptDrops(false);

    foreach(const QString &arg, qApp->arguments())
      urls << QUrl(arg);

    encoder.setUrls(urls);
    encoder.setResolutions(getResolutions());
    encoder.setQuality(getQuality());
    setWindowTitle("Starting...");
    encoder.start();
  }
}

void Preview::dragEnterEvent(QDragEnterEvent *event)
{
  if(event->mimeData()->hasUrls())
    event->acceptProposedAction();
}

void Preview::dropEvent(QDropEvent *event)
{
  if(!encoder.isRunning()) {
    resolutionMenu->setEnabled(false);
    qualityMenu->setEnabled(false);
    setAcceptDrops(false);

    encoder.setUrls(event->mimeData()->urls());
    encoder.setResolutions(getResolutions());
    encoder.setQuality(getQuality());
    setWindowTitle("Starting...");
    encoder.start();
  }  

  event->acceptProposedAction();
}

void Preview::closeEvent(QCloseEvent *event)
{
  Q_UNUSED(event);

  quitSlot();
}

void Preview::contextMenuEvent(QContextMenuEvent *event)
{
  contextMenu->popup(event->globalPos());
}

void Preview::information(const QString &fileName)
{
  if(fileName.startsWith("ERROR")) {
    QString text = fileName;
    text.replace("ERROR: ", "");
    setWindowTitle(text);
    return;
  }

  if(fileName.startsWith("FILE")) {
    QString text = fileName;
    text.replace("FILE: ", "");
    setWindowTitle(text);
    return;
  }

  if(fileName.startsWith("DONE")) {
    showInfo();

    if(qApp->arguments().count() > 1)
      quitSlot();

    return;
  }

  QPixmap background(size());
  background.fill(Qt::transparent);

  QPainter painter(&background);

  QPixmap pixmap(fileName);
  QPixmap scaled = pixmap.scaled(size(), Qt::KeepAspectRatio,
				 Qt::SmoothTransformation);
  painter.drawPixmap((background.width()-scaled.width())/2,
		     (background.height()-scaled.height())/2,
		     scaled);

  QPixmap overlay = sub(QString::number(currentProgress) + "%");
  painter.drawPixmap((background.width()-overlay.width())/2,
		     (background.height()-overlay.height())/2,
		     overlay);

  setPixmap(background);
}

void Preview::showInfo()
{
  clear();

  QPixmap background(400, 400);
  background.fill(Qt::transparent);

  QPainter painter(&background);

  QPixmap overlay(":/img/500px-Crystal_Clear_mimetype_video.svg.png");
  int overlayHeight = 250;
  overlay = overlay.scaledToHeight(overlayHeight, Qt::SmoothTransformation);
  int overlayWidth = overlay.width();
  QRect targetRect((400-overlayWidth)/2, 25, overlayWidth, overlayHeight);
  painter.drawPixmap(targetRect, overlay, overlay.rect());

  QFont defaultFont = font();
  QFont boldFont = defaultFont;
  boldFont.setBold(true);
  painter.setFont(boldFont);

  painter.drawText(QRect(0, overlayHeight+40, 400, 20), Qt::AlignCenter,
		   "Drag and drop image files/folders here");

  painter.setFont(defaultFont);

  painter.drawText(QRect(0, overlayHeight+70, 400, 20), Qt::AlignCenter,
		   "Supported formats: png, jpg (jpeg), tiff, gif");

  painter.drawText(QRect(0, overlayHeight+90, 400, 20), Qt::AlignCenter,
		   "Sorting: first integer in file name");


  painter.drawText(QRect(0, overlayHeight+110, 400, 20), Qt::AlignCenter,
		   "Right-click for preferences");

  setPixmap(background);

  resolutionMenu->setEnabled(true);
  qualityMenu->setEnabled(true);
  setAcceptDrops(true);
}

void Preview::quitSlot()
{
  exit(0);
}

void Preview::aboutSlot()
{
  QMessageBox msg;
  msg.setText("ElmerClips is a small utility program for making "
	      "video clips\nfrom image files produced e.g. by ElmerPost.\n\n"
	      "ElmerClips is open source software (GPL v2):\n"
	      "http://sourceforge.net/projects/elmerfem/\n\n"
	      "It uses the mpeg1video encoder from the ffmpeg-project:\n"
	      "http://www.ffmpeg.org/");
  msg.exec();
}

QList<int> Preview::getResolutions() const
{
  QList<int> resolutions;

  if(smallAction->isChecked())
    resolutions << 640;

  if(mediumAction->isChecked())
    resolutions << 720;

  if(bigAction->isChecked())
    resolutions << 1280;

  if(hugeAction->isChecked())
    resolutions << 1920;

  return resolutions;
}

int Preview::getQuality() const
{
  if(lowQualityAction->isChecked())
    return 16;

  if(mediumQualityAction->isChecked())
    return 8;

  if(highQualityAction->isChecked())
    return 4;

  return 2; // best
}

void Preview::progress(int value)
{
  currentProgress = value;
}

QPixmap Preview::sub(const QString &text) const
{
  QFont defaultFont = font();
  defaultFont.setPixelSize(60);
  defaultFont.setBold(true);

  QFontMetrics fontMetrics(defaultFont);
  int w = fontMetrics.width(text);
  int h = fontMetrics.height();

  QPixmap background(w+4, h+4);
  background.fill(Qt::transparent);

  QPainter painter(&background);
  painter.setRenderHint(QPainter::Antialiasing);
  painter.setFont(defaultFont);
  painter.setPen(QPen(Qt::black));

  for(int x = -2; x <= 2; x += 2) {
    for(int y = -2; y <= 2; y += 2)
      painter.drawText(QRect(x+2, y+2, w, h), text);
  }
  
  painter.setPen(QPen(Qt::white));
  painter.drawText(QRect(2, 2, w, h), text);
  
  return background;
}
