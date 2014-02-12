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

#include "ecmaconsole.h"

#include <QWidget>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QTextCursor>
#include <QTextBlock>
#include <QMetaObject>
#include <QMetaMethod>
#include <QCompleter>
#include <QStringListModel>
#include <QScrollBar>

#include <iostream>
using namespace std;

EcmaConsole::EcmaConsole(QWidget* parent)
  : QTextEdit(parent)
{
  prompt = "qs> ";
  this->clearHistory();
}

EcmaConsole::~EcmaConsole()
{
}

void EcmaConsole::mouseDoubleClickEvent(QMouseEvent* event)
{
  event->ignore();
}

void EcmaConsole::mousePressEvent(QMouseEvent* event)
{
  event->ignore();
}

void EcmaConsole::mouseReleaseEvent(QMouseEvent* event)
{
  event->ignore();
}

void EcmaConsole::keyPressEvent(QKeyEvent* event)
{
  if(completer && completer->popup()->isVisible()) {
    switch(event->key()) {
    case Qt::Key_Return:
      if(!completer->popup()->currentIndex().isValid()) {
        insertCompletion(completer->currentCompletion());
        completer->popup()->hide();
        event->accept();
	return;
      }
      event->ignore();
      return;
    case Qt::Key_Enter:
    case Qt::Key_Escape:
    case Qt::Key_Tab:
    case Qt::Key_Backtab:
      event->ignore();
      return;
    default:
      break;
    }
  }
  
  bool eventHandled = false;
  
  switch(event->key()) {
  case Qt::Key_Return:
    execLine();
    eventHandled = true;
    break;
    
  case Qt::Key_Up:
    if(historyPtr > 0) {
      historyPtr--;
      scanHistory();
    }
    eventHandled = true;
    break;
    
  case Qt::Key_Down:
    if(historyPtr < history.count()-1) {
      historyPtr++;
      scanHistory();
    }
    eventHandled = true;
    break;
    
  case Qt::Key_Left:
  case Qt::Key_Backspace:
  case Qt::Key_Backtab:
    if(this->textCursor().position() <= getPromptPos())
      eventHandled = true;
    break;
    
  default:
    break;
  }
  
  if(eventHandled) {
    completer->popup()->hide();
    event->ignore();
  } else {
    QTextEdit::keyPressEvent(event);
    QString text = event->text();
    if(!text.isEmpty()) {
      handleTabCompletion();
    } else {
      completer->popup()->hide();
    }
  }
}

int EcmaConsole::getPromptPos()
{
  QTextCursor textCursor(this->textCursor());
  textCursor.movePosition(QTextCursor::End);
  int position = textCursor.block().position() + prompt.length();
  return position;
}

void EcmaConsole::execLine()
{
  QTextCursor textCursor = this->textCursor();
  textCursor.movePosition(QTextCursor::End);
  textCursor.setPosition(getPromptPos());
  textCursor.movePosition(QTextCursor::End, QTextCursor::KeepAnchor);
  QString line = textCursor.selectedText().trimmed();

  if(!line.isEmpty()) {
    emit(cmd(line));
    history << line;
    historyPtr = history.count();
  }

  this->append(prompt);
  textCursor = this->textCursor();
  textCursor.movePosition(QTextCursor::End);
  setTextCursor(textCursor);
}

void EcmaConsole::scanHistory()
{
  QTextCursor textCursor = this->textCursor();
  textCursor.movePosition(QTextCursor::End);
  textCursor.setPosition(getPromptPos(), QTextCursor::KeepAnchor);
  textCursor.insertText(history.value(historyPtr));
  textCursor.movePosition(QTextCursor::End);
  setTextCursor(textCursor);
}

void EcmaConsole::clearHistory()
{
  this->clear();
  history.clear();
  historyPtr = 0;
  this->append(prompt);
}

void EcmaConsole::addNames(QString className, const QMetaObject* metaObject)
{
  QStringList publicSlots;
  int methodCount = metaObject->methodCount();
  for(int i = 0; i < methodCount; i++) {
    QMetaMethod method = metaObject->method(i);
    QMetaMethod::Access access = method.access();
    QMetaMethod::MethodType methodType = method.methodType();
    if((access == QMetaMethod::Public) && (methodType == QMetaMethod::Slot)) {
      QString signature = method.signature();
      int j = signature.indexOf("(");
      QString slotName = signature.left(j);
      publicSlots << slotName;
    }
  }
  publicSlots.sort();
  names.insert(className, publicSlots);
}

void EcmaConsole::initCompleter()
{
  completer = new QCompleter(this);
  completer->setWidget(this);
  connect(completer, SIGNAL(activated(const QString&)), this, SLOT(insertCompletion(const QString&)));
}


void EcmaConsole::handleTabCompletion()
{
  QTextCursor textCursor = this->textCursor();
  int pos = textCursor.position();
  textCursor.setPosition(getPromptPos());
  textCursor.movePosition(QTextCursor::End, QTextCursor::KeepAnchor);
  int startPos = textCursor.selectionStart();

  int offset = pos-startPos;
  QString text = textCursor.selectedText();

  QString textToComplete;
  int cur = offset;

  while(cur--) {
    QChar c = text.at(cur);
    if(c.isLetterOrNumber() || (c == '.') || (c == '_')) {
      textToComplete.prepend(c);
    } else {
      break;
    }
  }

  QString lookup;
  QString compareText = textToComplete;
  int dot = compareText.lastIndexOf('.');

  if(dot != -1) {
    lookup = compareText.mid(0, dot);
    compareText = compareText.mid(dot+1, offset);
  }

  if(!lookup.isEmpty() || !compareText.isEmpty()) {
    compareText = compareText.toLower();
    QStringList found;

    QStringList list;
    if(lookup.isEmpty()) {
      // all class names
      list = names.keys();
    } else {
      // methods for a class
      list = names.value(lookup);
    }

    foreach(QString name, list) {
      if(name.toLower().startsWith(compareText))
	found << name;
    }
    
    if(!found.isEmpty()) {
      completer->setCompletionPrefix(compareText);
      completer->setCompletionMode(QCompleter::PopupCompletion);
      completer->setModel(new QStringListModel(found, completer));
      completer->setCaseSensitivity(Qt::CaseInsensitive);
      QTextCursor c = this->textCursor();
      c.movePosition(QTextCursor::StartOfWord);
      QRect cr = cursorRect(c);
      cr.setWidth(completer->popup()->sizeHintForColumn(0)
        + completer->popup()->verticalScrollBar()->sizeHint().width());
      cr.translate(0,8);
      completer->complete(cr);
    } else {
      completer->popup()->hide();
    }
  } else {
    completer->popup()->hide();
  }
}

void EcmaConsole::insertCompletion(const QString& completion)
{
  QTextCursor tc = textCursor();
  tc.movePosition(QTextCursor::Left, QTextCursor::KeepAnchor);
  if (tc.selectedText() == ".") {
    tc.insertText(QString(".") + completion);
  } else {
    tc = textCursor();
    tc.movePosition(QTextCursor::StartOfWord, QTextCursor::MoveAnchor);
    tc.movePosition(QTextCursor::EndOfWord, QTextCursor::KeepAnchor);
    tc.insertText(completion);
    setTextCursor(tc);
  }
}
