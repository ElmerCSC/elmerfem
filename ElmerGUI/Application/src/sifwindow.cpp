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
 *  ElmerGUI sifwindow                                                       *
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

#include <QtGui>
#include <iostream>
#include "sifwindow.h"
#include "mainwindow.h"

#if WITH_QT5
#include <QtWidgets>
#include <QPrinter>
#include <QPrintDialog>
#endif

using namespace std;

SifHighlighter::SifHighlighter(int type, QTextDocument *parent)
    : QSyntaxHighlighter(parent)
{
  if(type != SIF_HIGHLIGHTING_LIGHT && type != SIF_HIGHLIGHTING_DARK) return;
  
	QColor yellow  = QColor( 161,  97,   0);
	QColor orange  = QColor( 183,  55,  22);
	QColor red     = QColor( 100,  30,  27);
	QColor magenta = QColor( 191,  34, 110);
	QColor violet  = QColor(  88,  93, 176);
	QColor blue    = QColor(  18,  99, 190);
	QColor cyan    = QColor(  22, 141, 132);
	QColor green   = QColor( 113, 113,   0);

  if(type == SIF_HIGHLIGHTING_DARK)
  {
     yellow  = QColor( 201, 157,  20);
     orange  = QColor( 223,  95,  42);
     red     = QColor( 240,  70,  67);
     magenta = QColor( 231,  74, 150);
     violet  = QColor( 128, 133, 216);
     blue    = QColor(  58, 159, 250);
     cyan    = QColor(  62, 181, 172);
     green   = QColor( 153, 173,   0);  
    
  }
	
	QColor cBlock 	= blue;
	QColor cKeyword = violet;
	QColor cValue 	= cyan;
	QColor cQuotation = orange;
	QColor cComment	= green;
	QColor cSuffix = magenta;
	
    HighlightingRule rule;
    
    blockFormat.setForeground(cBlock);
    QStringList blockPatterns;
    blockPatterns << "^\\s*Equation\\s+[1-9][0-9]*\\s*$" << "^\\s*Material\\s+[1-9][0-9]*\\s*$" 
                  << "^\\s*Body\\s+Force\\s+[1-9][0-9]*\\s*$" << "^\\s*Initial\\s+Condition\\s+[1-9][0-9]*\\s*$"
                  << "^\\s*Boundary\\s+Condition\\s+[1-9][0-9]*\\s*$" << "^\\s*Component\\s+[1-9][0-9]*\\s*$"
                  << "^\\s*Solver\\s+[1-9][0-9]*\\s*$" << "^\\s*Body\\s+[1-9][0-9]*\\s*$"
                  << "^\\s*Header\\s*$" << "^\\s*Simulation\\s*$" << "^\\s*Constants\\s*$" << "^\\s*End\\s*$"
                  
                  << "^\\s*Equation\\s+[1-9][0-9]*\\s*!" << "^\\s*Material\\s+[1-9][0-9]*\\s*!" 
                  << "^\\s*Body\\s+Force\\s+[1-9][0-9]*\\s*!" << "^\\s*Initial\\s+Condition\\s+[1-9][0-9]*\\s*!"
                  << "^\\s*Boundary\\s+Condition\\s+[1-9][0-9]*\\s*!" << "^\\s*Component\\s+[1-9][0-9]*\\s*!"
                  << "^\\s*Solver\\s+[1-9][0-9]*\\s*!" << "^\\s*Body\\s+[1-9][0-9]*\\s*!"
                  << "^\\s*Header\\s*!" << "^\\s*Simulation\\s*!" << "^\\s*Constants\\s*!" << "^\\s*End\\s*!"
                  ;
    foreach (const QString &pattern, blockPatterns) {
        rule.pattern = QRegExp(pattern, Qt::CaseInsensitive);
        rule.format = blockFormat;
        highlightingRules.append(rule);
    }
    
    keywordFormat.setForeground(cKeyword);
    rule.pattern = QRegExp("\\bInteger\\b", Qt::CaseInsensitive);
    rule.format = keywordFormat;
    highlightingRules.append(rule);  
    rule.pattern = QRegExp("\\bReal\\b", Qt::CaseInsensitive);
    rule.format = keywordFormat;
    highlightingRules.append(rule);
    rule.pattern = QRegExp("\\bLogical\\b", Qt::CaseInsensitive);
    rule.format = keywordFormat;
    highlightingRules.append(rule);
    rule.pattern = QRegExp("\\bMATC\\b", Qt::CaseInsensitive);
    rule.format = keywordFormat;
    highlightingRules.append(rule); 
    rule.pattern = QRegExp("\\bEquals\\b", Qt::CaseInsensitive);
    rule.format = keywordFormat;
    highlightingRules.append(rule);
    
    valueFormat.setForeground(cValue);
    rule.pattern = QRegExp("\\bTrue\\b", Qt::CaseInsensitive);
    rule.format = valueFormat;
    highlightingRules.append(rule);  
    rule.pattern = QRegExp("\\bFalse\\b", Qt::CaseInsensitive);
    rule.format = valueFormat;
    highlightingRules.append(rule);
         
    quotationFormat.setForeground(cQuotation);
    rule.pattern = QRegExp("\".*\"");
    rule.format = quotationFormat;
    highlightingRules.append(rule); 
		
	  suffixFormat.setForeground(cSuffix);
    rule.pattern = QRegExp("\\([1-9][0-9]*\\)", Qt::CaseInsensitive);
    rule.format = suffixFormat;
    highlightingRules.append(rule); 
    
    commentFormat.setForeground(cComment);
		commentFormat.setFontWeight(QFont::Normal); 	
    rule.pattern = QRegExp("![^\n]*");
    rule.format = commentFormat;
    highlightingRules.append(rule); 
}

void SifHighlighter::highlightBlock(const QString &text)
{
    foreach (const HighlightingRule &rule, highlightingRules) {
        QRegExp expression(rule.pattern);
        int index = expression.indexIn(text);
        while (index >= 0) {
            int length = expression.matchedLength();
            setFormat(index, length, rule.format);
            index = expression.indexIn(text, index + length);
        }
    }
}

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
  
  highlighter = NULL;

	QString strFont = ((MainWindow*) parent)->settings_value(QString("sifWindow/font"), QString("")).toString();
	QFont font;
	if(!strFont.isEmpty() && font.fromString(strFont)){
		font.setFixedPitch(true);
		textEdit->setFont(font);
		lineEdit->setFont(font);
	}
	
	int syntaxHighlighting = ((MainWindow*) parent)->settings_value(QString("sifWindow/highlighting"), SIF_HIGHLIGHTING_NONE).toInt();
	if(syntaxHighlighting == SIF_HIGHLIGHTING_NONE){
		highlightingNoneSlot();
	}else	if(syntaxHighlighting == SIF_HIGHLIGHTING_LIGHT){
		highlightingLightSlot();
	}else if(syntaxHighlighting == SIF_HIGHLIGHTING_DARK){
		highlightingDarkSlot();
	}

  int x,y,w,h;  
  x = ((MainWindow*) parent)->settings_value("sifWindow/x", -10000).toInt();
  y = ((MainWindow*) parent)->settings_value("sifWindow/y", -10000).toInt();
  w = ((MainWindow*) parent)->settings_value("sifWindow/width", -10000).toInt(); 
  h = ((MainWindow*) parent)->settings_value("sifWindow/height", -10000).toInt();
  if(x != -10000 && y != -10000 && w != -10000 && h != -10000 && x < QApplication::desktop()->width()-100 && y < QApplication::desktop()->height()-100){
    move(x,y);
    if(w > QApplication::desktop()->width()) w = QApplication::desktop()->width();
    if(h > QApplication::desktop()->height()) h = QApplication::desktop()->height();    
    resize(w,h);
  }
  if(((MainWindow*) parent)->settings_value("sifWindow/maximized", false).toBool()){
    setWindowState(windowState() ^ Qt::WindowMaximized);
  }
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
	
	fontAct = new QAction(QIcon(""), tr("&Font"), this);
  findAct->setStatusTip(tr("Select font"));
  connect(fontAct, SIGNAL(triggered()), this, SLOT(fontSlot()));
	
	highlightingNoneAct = new QAction(QIcon(""), tr("&None"), this);
  highlightingNoneAct->setStatusTip(tr("No highlighting"));
  highlightingNoneAct->setCheckable(true);
  connect(highlightingNoneAct, SIGNAL(triggered()), this, SLOT(highlightingNoneSlot()));

	highlightingLightAct = new QAction(QIcon(""), tr(" &Light"), this);
  highlightingLightAct->setStatusTip(tr("Highlight in light theme"));
  highlightingLightAct->setCheckable(true);
  connect(highlightingLightAct, SIGNAL(triggered()), this, SLOT(highlightingLightSlot()));

	highlightingDarkAct = new QAction(QIcon(""), tr(" &Dark"), this);
  highlightingDarkAct->setStatusTip(tr("Highlight in dark theme"));
  highlightingDarkAct->setCheckable(true);
  connect(highlightingDarkAct, SIGNAL(triggered()), this, SLOT(highlightingDarkSlot()));
	
	saveAndRunAct = new QAction(QIcon(":/icons/arrow-right.png"), tr("&Save and run"), this);
	saveAndRunAct->setStatusTip(tr("Save the sif and project (without auto-generating sif),  then run solver "));
  connect(saveAndRunAct, SIGNAL(triggered()), this, SLOT(saveAndRunSlot()));	
}

void SifWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newAct);
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAndRunAct);
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
	
	preferenceMenu = menuBar()->addMenu(tr("&Preference"));
	preferenceMenu->addAction(fontAct);
	highlightingMenu = preferenceMenu->addMenu(tr("&Syntax highlighting"));
	highlightingMenu->addAction(highlightingNoneAct);
  highlightingMenu->addAction(highlightingLightAct);
  highlightingMenu->addAction(highlightingDarkAct);	
}

void SifWindow::createToolBars()
{
  fileToolBar = addToolBar(tr("&File"));
  fileToolBar->addAction(newAct);
  fileToolBar->addAction(openAct);
  fileToolBar->addAction(saveAct);
	fileToolBar->addAction(saveAndRunAct);
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

void SifWindow::fontSlot()
{
	bool ok;
	QFont font = QFontDialog::getFont(&ok, textEdit->font());
	if (ok) {
	  font.setFixedPitch(true);
		textEdit->setFont(font);
		lineEdit->setFont(font);
	}
	((MainWindow*) parent())->settings_setValue(QString("sifWindow/font"), font.toString());
}
void SifWindow::highlightingNoneSlot()
{
	QString style = "";
	setStyleSheet(style);
	delete highlighter;
	highlighter = NULL;
	((MainWindow*) parent())->settings_setValue(QString("sifWindow/highlighting"), SIF_HIGHLIGHTING_NONE);
	highlightingNoneAct->setChecked(true);	
	highlightingLightAct->setChecked(false);
	highlightingDarkAct->setChecked(false);
}

void SifWindow::highlightingLightSlot()
{
	QString style = "QTextEdit { color: #384e55; background: #fffdf6}";	
	setStyleSheet(style);
	delete highlighter;
  highlighter = new SifHighlighter(SIF_HIGHLIGHTING_LIGHT, textEdit->document());
	((MainWindow*) parent())->settings_setValue(QString("sifWindow/highlighting"), SIF_HIGHLIGHTING_LIGHT);
	highlightingNoneAct->setChecked(false);	
	highlightingLightAct->setChecked(true);
	highlightingDarkAct->setChecked(false);
}

void SifWindow::highlightingDarkSlot()
{
	QString style = "QTextEdit { color: #c3b1b1; background: #000814}";
	setStyleSheet(style);
	delete highlighter;
  highlighter = new SifHighlighter(SIF_HIGHLIGHTING_DARK,textEdit->document());
	((MainWindow*) parent())->settings_setValue(QString("sifWindow/highlighting"), SIF_HIGHLIGHTING_DARK);	
	highlightingNoneAct->setChecked(false);	
	highlightingLightAct->setChecked(false);
	highlightingDarkAct->setChecked(true);
}

void SifWindow::saveAndRunSlot()
{
  ((MainWindow*) parent())->saveAndRun(false);
}