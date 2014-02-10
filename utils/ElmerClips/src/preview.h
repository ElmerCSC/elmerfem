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
#ifndef PREVIEW_H
#define PREVIEW_H

#include <QtGui>
#include "encoder.h"

class Preview : public QLabel
{
 Q_OBJECT

 public:
  Preview(QWidget *parent = 0);
  void checkCommandLine();

 public slots:
  void information(const QString &fileName);
  void progress(int value);

 protected:
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);
  void contextMenuEvent(QContextMenuEvent *event);

 private slots:
  void quitSlot();
  void aboutSlot();

 private:
  void setupContextMenu();
  void showInfo();
  QList<int> getResolutions() const;
  int getQuality() const;
  QPixmap sub(const QString &text) const;

  Encoder encoder;
  QMenu *resolutionMenu;
  QMenu *qualityMenu;
  QMenu *contextMenu;
  QAction *smallAction;
  QAction *mediumAction;
  QAction *bigAction;
  QAction *hugeAction;
  QAction *quitAction;
  QAction *lowQualityAction;
  QAction *mediumQualityAction;
  QAction *highQualityAction;
  QAction *bestQualityAction;
  QAction *aboutAction;
  QActionGroup *resolutionActionGroup;
  QActionGroup *qualityActionGroup;
  int currentProgress;
};

#endif // PREVIEW_H
