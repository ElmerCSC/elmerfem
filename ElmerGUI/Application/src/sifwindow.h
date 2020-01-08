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
 *  ElmerGUI sifwindow                                                       *
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

#ifndef SIFWINDOW_H
#define SIFWINDOW_H

#define SIF_HIGHLIGHTING_NONE	0
#define SIF_HIGHLIGHTING_LIGHT	1
#define SIF_HIGHLIGHTING_DARK	2

#include <QMainWindow>
#include <QSyntaxHighlighter>

class QTextEdit;
class QLineEdit;

class SifHighlighter : public QSyntaxHighlighter
{
    Q_OBJECT

public:
    SifHighlighter(int type, QTextDocument *parent = 0);

protected:
    void highlightBlock(const QString &text);

private:
    struct HighlightingRule
    {
        QRegExp pattern;
        QTextCharFormat format;
    };
    QVector<HighlightingRule> highlightingRules;

    QRegExp commentStartExpression;
    QRegExp commentEndExpression;

    QTextCharFormat blockFormat;
    QTextCharFormat classFormat;
    QTextCharFormat commentFormat;
    QTextCharFormat quotationFormat;
    QTextCharFormat keywordFormat;
    QTextCharFormat simulationTypeFormat;
    QTextCharFormat valueFormat;  
		QTextCharFormat suffixFormat;
};

class SifWindow : public QMainWindow
{
  Q_OBJECT

public:
  SifWindow(QWidget *parent = 0);
  ~SifWindow();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  QTextEdit* getTextEdit(void);
  SifHighlighter *highlighter; 
  void setFirstTime(bool);
  void setFound(bool);

private slots:
  void newSlot();
  void openSlot();
  void saveSlot();
  void printSlot();
  void findSlot();
	void fontSlot();
	void highlightingNoneSlot();
	void highlightingLightSlot();
	void highlightingDarkSlot();
	void saveAndRunSlot();

private:
  QTextEdit *textEdit;
  bool firstTime;
  bool found;

  QLineEdit *lineEdit;

  QAction *newAct;
  QAction *openAct;
  QAction *saveAct;
  QAction *printAct;
  QAction *exitAct;
  QAction *cutAct;
  QAction *copyAct;
  QAction *pasteAct;
  QAction *findAct;
	QAction *fontAct;
	QAction *highlightingNoneAct;
	QAction *highlightingLightAct;
	QAction *highlightingDarkAct;
	QAction *saveAndRunAct;

  QMenu *fileMenu;
  QMenu *editMenu;
	QMenu *preferenceMenu;
	QMenu *highlightingMenu;

  QToolBar *fileToolBar;
  QToolBar *editToolBar;

  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
};

#endif
