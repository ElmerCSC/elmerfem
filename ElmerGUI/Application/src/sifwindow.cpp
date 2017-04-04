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

#include <QtGui>
#include <iostream>
#include "sifwindow.h"

#if WITH_QT5
#include <QtWidgets>
#include <QPrinter>
#include <QPrintDialog>
#endif

using namespace std;

SifWindow::SifWindow(QWidget *parent)
  : QMainWindow(parent)
{
  setWindowFlags(Qt::Window);

  textEdit = new QTextEdit;
  textEdit->setLineWrapMode(QTextEdit::NoWrap);

  setCentralWidget(textEdit);

  lineEdit = new QLineEdit;
  connect(lineEdit, SIGNAL(returnPressed()), this, SLOT(findSlot()));

  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  firstTime = true;
  found = false;

  setWindowTitle(tr("Solver Input File"));
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

SifWindow::~SifWindow()
{
}

QTextEdit* SifWindow::getTextEdit(void)
{
  return this->textEdit;
}

void SifWindow::setFirstTime(bool b)
{
  this->firstTime = b;
}

void SifWindow::setFound(bool b)
{
  this->found = b;
}

QSize SifWindow::minimumSizeHint() const
{
  return QSize(64, 64);
}


QSize SifWindow::sizeHint() const
{
  return QSize(640, 640);
}

void SifWindow::createActions()
{
  newAct = new QAction(QIcon(":/icons/document-new.png"), tr("&New"), this);
  newAct->setShortcut(tr("Ctrl+N"));
  newAct->setStatusTip(tr("New text document"));
  connect(newAct, SIGNAL(triggered()), this, SLOT(newSlot()));

  openAct = new QAction(QIcon(":/icons/document-open.png"), tr("&Open..."), this);
  openAct->setShortcut(tr("Ctrl+O"));
  openAct->setStatusTip(tr("Open text file"));
  connect(openAct, SIGNAL(triggered()), this, SLOT(openSlot()));

  saveAct = new QAction(QIcon(":/icons/document-save.png"), tr("&Save as..."), this);
  saveAct->setShortcut(tr("Ctrl+S"));
  saveAct->setStatusTip(tr("Save text file"));
  connect(saveAct, SIGNAL(triggered()), this, SLOT(saveSlot()));

  printAct = new QAction(QIcon(":/icons/document-print.png"), tr("&Print..."), this);
  printAct->setShortcut(tr("Ctrl+P"));
  printAct->setStatusTip(tr("Print document"));
  connect(printAct, SIGNAL(triggered()), this, SLOT(printSlot()));

  exitAct = new QAction(QIcon(":/icons/application-exit.png"), tr("&Quit"), this);
  exitAct->setShortcut(tr("Ctrl+Q"));
  exitAct->setStatusTip(tr("Quit editor"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

  cutAct = new QAction(QIcon(":/icons/edit-cut.png"), tr("Cu&t"), this);
  cutAct->setShortcut(tr("Ctrl+X"));
  cutAct->setStatusTip(tr("Cut the current selection to clipboard"));
  connect(cutAct, SIGNAL(triggered()), this->textEdit, SLOT(cut()));

  copyAct = new QAction(QIcon(":/icons/edit-copy.png"), tr("&Copy"), this);
  copyAct->setShortcut(tr("Ctrl+C"));
  copyAct->setStatusTip(tr("Copy the current selection to clipboard"));
  connect(copyAct, SIGNAL(triggered()), this->textEdit, SLOT(copy()));

  pasteAct = new QAction(QIcon(":/icons/edit-paste.png"), tr("&Paste"), this);
  pasteAct->setShortcut(tr("Ctrl+V"));
  pasteAct->setStatusTip(tr("Paste clipboard into the current selection"));
  connect(pasteAct, SIGNAL(triggered()), this->textEdit, SLOT(paste()));

  findAct = new QAction(QIcon(":/icons/edit-find.png"), tr("&Find"), this);
  findAct->setShortcut(tr("Ctrl+F"));
  findAct->setStatusTip(tr("Find text in document"));
  connect(findAct, SIGNAL(triggered()), this, SLOT(findSlot()));
}

void SifWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newAct);
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);
  fileMenu->addSeparator();
  fileMenu->addAction(printAct);
  fileMenu->addSeparator();
  fileMenu->addAction(exitAct);

  editMenu = menuBar()->addMenu(tr("&Edit"));
  editMenu->addAction(cutAct);
  editMenu->addAction(copyAct);
  editMenu->addAction(pasteAct);
  editMenu->addSeparator();
  editMenu->addAction(findAct);
}

void SifWindow::createToolBars()
{
  fileToolBar = addToolBar(tr("&File"));
  fileToolBar->addAction(newAct);
  fileToolBar->addAction(openAct);
  fileToolBar->addAction(saveAct);
  fileToolBar->addAction(printAct);

  editToolBar = addToolBar(tr("&Edit"));
  editToolBar->addAction(cutAct);
  editToolBar->addAction(copyAct);
  editToolBar->addAction(pasteAct);
  editToolBar->addSeparator();
  editToolBar->addWidget(lineEdit);
  editToolBar->addAction(findAct);
}

void SifWindow::createStatusBar()
{
  statusBar()->showMessage(tr("Ready"));
}

void SifWindow::newSlot()
{
  textEdit->clear();

  firstTime = true;
  found = false;

  statusBar()->showMessage(tr("Ready"));
}

void SifWindow::openSlot()
{
  QString fileName;
  
  fileName = QFileDialog::getOpenFileName(this, tr("Open text file"));

  if(fileName.isEmpty())
    return;

  QFile file;
  file.setFileName(fileName);
  if(!file.open(QIODevice::ReadOnly))
    return;
  
  QTextStream inputStream(&file);

  statusBar()->showMessage(tr("Opening file..."));

  textEdit->clear();

  QString line = inputStream.readAll();

  file.close();

  textEdit->append(line);

  firstTime = true;
  found = false;
  
  statusBar()->showMessage(tr("Ready"));
}

void SifWindow::saveSlot()
{
  QString fileName;
  
  fileName = QFileDialog::getSaveFileName(this, tr("Save text file"));
  
  if(fileName.isEmpty())
    return;
  
  QFile file;
  file.setFileName(fileName);
  if(!file.open(QIODevice::WriteOnly))
    return;
  
  QTextStream outputStream(&file);

  statusBar()->showMessage(tr("Saving file..."));

  outputStream << textEdit->toPlainText();

  file.close();

  statusBar()->showMessage(tr("Ready"));
}

void SifWindow::printSlot()
{
  QTextDocument *document = textEdit->document();
  QPrinter printer;

  QPrintDialog *printDialog = new QPrintDialog(&printer, this);
  if (printDialog->exec() != QDialog::Accepted)
    return;
  
  statusBar()->showMessage(tr("Printing..."));

  document->print(&printer);
  
  statusBar()->showMessage(tr("Ready"));
}

void SifWindow::findSlot()
{
  QString searchString = lineEdit->text().trimmed();
  QTextDocument *document = textEdit->document();
  
  if(!firstTime && found)
    document->undo();

  found = false;
  
  if(searchString == "") {
    QMessageBox::information(this,
			     tr("Empty string"),
			     "Please enter a string in the "
			     "line edit box in the tool bar");
  } else {
    
    QTextCursor highlightCursor(document);  
    QTextCursor cursor(document);
    
    cursor.beginEditBlock();
    
    QTextCharFormat plainFormat(highlightCursor.charFormat());
    QTextCharFormat colorFormat = plainFormat;
    colorFormat.setForeground(Qt::red);
    colorFormat.setFontWeight(QFont::Bold);
    
    while(!highlightCursor.isNull() && !highlightCursor.atEnd()) {
      highlightCursor = document->find(searchString, highlightCursor);
      
      if(!highlightCursor.isNull()) {
	found = true;
	highlightCursor.mergeCharFormat(colorFormat);
      }
    }
    
    cursor.endEditBlock();
    firstTime = false;
    
    if(!found)
      QMessageBox::information(this, tr("String not found"),
			"The string was not found in the document");

  }

  statusBar()->showMessage(tr("Ready"));
}
