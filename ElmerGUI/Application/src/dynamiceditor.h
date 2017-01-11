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
 *  ElmerGUI dynamiceditor                                                   *
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

#ifndef DYNAMICEDITOR_H
#define DYNAMICEDITOR_H

enum BodyTypes {
  BODY_MATERIAL,
  BODY_INITIAL,
  BODY_FORCE,
  BODY_EQUATION
};

enum MatTypes {
  MAT_APPLY,
  MAT_OK,
  MAT_NEW,
  MAT_DELETE
};

// #ifndef WITH_QT5
#include <QWidget>
#include <QtGui>
#include <QIcon>
#include <QDomDocument>
#include <QLayout>
// #else
// #include <QtGui>
// #include <QIcon>
// #include <QDomDocument>
// #include <QLayout>
// #endif

#ifdef WITH_QT5
#include <QtWidgets>
#endif

class QTabWidget;
class QPushButton;

class hash_entry_t
{
 public:
  QWidget *widget,*label;
  QDomElement elem;
};

class DynLineEdit : public QWidget
{
  Q_OBJECT

public:
  DynLineEdit(QWidget *parent=0);
 ~DynLineEdit();
  QString name;
  QLineEdit *lineEdit;

private:
  QTextEdit *textEdit;
  QLabel *label;
  QFrame *frame;
  QLayout *layout;
  QPushButton *closeButton;

private slots:
  void editSlot();
  void lineEditClose();
};

class DynamicEditor : public QWidget
{
  Q_OBJECT

public:
  DynamicEditor(QWidget *parent = 0);
  ~DynamicEditor();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  void setupTabs(QDomDocument *doc, const QString &Section, int ID);

  QTabWidget *tabWidget;
  QLineEdit *nameEdit;
  int tabs;

  QGroupBox *spareBox;
  QScrollArea *spareScroll;

  QPushButton *okButton;
  QPushButton *newButton;
  QPushButton *applyButton;
  QPushButton *spareButton;
  QPushButton *discardButton;

  QAction *menuAction;  // action for menu item
  int ID;               // id in propertyarray

  bool touched;

  QHash<QString, hash_entry_t> hash;

  void dumpHash(QDomDocument*, QDomElement*);
  void populateHash(QDomElement*);

signals:
  void dynamicEditorReady(int, int);
  void dynamicEditorSpareButtonClicked(int, int);

private slots:
  void okButtonClicked();
  void newButtonClicked();
  void applyButtonClicked();
  void discardButtonClicked();
  void spareButtonClicked();
  void lSlot(int);
  void comboSlot(QString);
  void textChangedSlot(QString);

private:
  hash_entry_t h;

  QIcon newIcon;
  QIcon addIcon;
  QIcon okIcon;
  QIcon removeIcon;

  QDomElement root;
  QDomElement all_stuff;
  QDomElement element;
  QDomElement name;
  QDomElement section;
  QDomElement param;
};

#endif // DYNAMICEDITOR_H
