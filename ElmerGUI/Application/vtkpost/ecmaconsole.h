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
 *  ElmerGUI ecmaconsole                                                     *
 *                                                                           *
 *  Modified from the PythonQt console by Florian Link / MeVis Research      *
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

#ifndef ECMACONSOLE_H
#define ECMACONSOLE_H

#include <QTextEdit>
#include <QString>
#include <QStringList>
#include <QHash>

class QWidget;
class QKeyEvent;
class QMouseEvent;
class QMetaObject;
class QCompleter;

class EcmaConsole : public QTextEdit
{
  Q_OBJECT

public:
  EcmaConsole(QWidget* parent = 0);
  ~EcmaConsole();
  void clearHistory();
  void addNames(QString, const QMetaObject*);
  void initCompleter();

public slots:
  void keyPressEvent(QKeyEvent*);
  void mousePressEvent(QMouseEvent*);
  void mouseReleaseEvent(QMouseEvent*);
  void mouseDoubleClickEvent(QMouseEvent*);
  void insertCompletion(const QString&);

signals:
  void cmd(QString);

private:
  QString prompt;
  int getPromptPos();
  void execLine();
  void scanHistory();
  void handleTabCompletion();
  int historyPtr;
  QStringList history;
  QHash<QString, QStringList> names;
  QCompleter* completer;
};

#endif
