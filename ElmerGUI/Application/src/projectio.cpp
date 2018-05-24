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
 *  ElmerGUI projectio                                                       *
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
#include "projectio.h"

#if WITH_QT5
#include <QtWidgets>
#endif

using namespace std;

ProjectIO::ProjectIO(QWidget *parent)
  : QDialog(parent)
{
  parentWidget = parent;
}

ProjectIO::~ProjectIO()
{
}

void ProjectIO::appendToProject(QDomDocument *projectDoc, QDomElement *item)
{
  // Radio buttons:
  //----------------
  QList<QRadioButton *> allRadioButtons = parentWidget->findChildren<QRadioButton *>(); 
  
  for(int i = 0; i < allRadioButtons.size(); i++) {
    QRadioButton *rb = allRadioButtons.at(i);

    if(!rb)
      continue;

    QString rbObjectName = rb->objectName();
    QString rbValue = QString::number(rb->isChecked());

    if(rbObjectName.isEmpty())
      continue;

    QDomElement widget = projectDoc->createElement("widget");
    widget.setAttribute("type", "RadioButton");
    item->appendChild(widget);
    
    QDomElement objectName = projectDoc->createElement("objectName");
    QDomText objectNameValue = projectDoc->createTextNode(rbObjectName);
    objectName.appendChild(objectNameValue);
    widget.appendChild(objectName);

    QDomElement isChecked = projectDoc->createElement("isChecked");
    QDomText isCheckedValue = projectDoc->createTextNode(rbValue);
    isChecked.appendChild(isCheckedValue);
    widget.appendChild(isChecked);
  }

  // Check boxes:
  //--------------
  QList<QCheckBox *> allCheckBoxes = parentWidget->findChildren<QCheckBox *>(); 
  
  for(int i = 0; i < allCheckBoxes.size(); i++) {
    QCheckBox *cb = allCheckBoxes.at(i);
    
    if(!cb)
      continue;

    QString cbObjectName = cb->objectName();
    QString cbValue = QString::number(cb->isChecked());

    if(cbObjectName.isEmpty())
      continue;

    QDomElement widget = projectDoc->createElement("widget");
    widget.setAttribute("type", "CheckBox");
    item->appendChild(widget);
    
    QDomElement objectName = projectDoc->createElement("objectName");
    QDomText objectNameValue = projectDoc->createTextNode(cbObjectName);
    objectName.appendChild(objectNameValue);
    widget.appendChild(objectName);

    QDomElement isChecked = projectDoc->createElement("isChecked");
    QDomText isCheckedValue = projectDoc->createTextNode(cbValue);
    isChecked.appendChild(isCheckedValue);
    widget.appendChild(isChecked);
  }

  // Line edits:
  //-------------
  QList<QLineEdit *> allLineEdits = parentWidget->findChildren<QLineEdit *>();
  
  for(int i = 0; i < allLineEdits.size(); i++) {
    QLineEdit *le = allLineEdits.at(i);

    if(!le)
      continue;

    QString leObjectName = le->objectName();
    QString leValue = le->text().trimmed();

    if(leObjectName.isEmpty())
      continue;

    QDomElement widget = projectDoc->createElement("widget");
    widget.setAttribute("type", "LineEdit");
    item->appendChild(widget);
    
    QDomElement objectName = projectDoc->createElement("objectName");
    QDomText objectNameValue = projectDoc->createTextNode(leObjectName);
    objectName.appendChild(objectNameValue);
    widget.appendChild(objectName);

    QDomElement text = projectDoc->createElement("text");
    QDomText textValue = projectDoc->createTextNode(leValue);
    text.appendChild(textValue);
    widget.appendChild(text);
  }

  // Text edits:
  //-------------
  QList<QTextEdit *> allTextEdits = parentWidget->findChildren<QTextEdit *>();
  
  for(int i = 0; i < allTextEdits.size(); i++) {
    QTextEdit *te = allTextEdits.at(i);

    if(!te)
      continue;

    QString teObjectName = te->objectName();
    QString teValue = te->toPlainText();

    if(teObjectName.isEmpty())
      continue;

    QDomElement widget = projectDoc->createElement("widget");
    widget.setAttribute("type", "TextEdit");
    item->appendChild(widget);
    
    QDomElement objectName = projectDoc->createElement("objectName");
    QDomText objectNameValue = projectDoc->createTextNode(teObjectName);
    objectName.appendChild(objectNameValue);
    widget.appendChild(objectName);

    QDomElement text = projectDoc->createElement("text");
    QDomText textValue = projectDoc->createTextNode(teValue);
    text.appendChild(textValue);
    widget.appendChild(text);
  }

  // Combo boxes:
  //--------------
  QList<QComboBox *> allComboBoxes = parentWidget->findChildren<QComboBox *>(); 
  
  for(int i = 0; i < allComboBoxes.size(); i++) {
    QComboBox *cx = allComboBoxes.at(i);

    if(!cx)
      continue;

    QString cxObjectName = cx->objectName();
    QString cxValue = QString::number(cx->currentIndex());

    if(cxObjectName.isEmpty())
      continue;

    QDomElement widget = projectDoc->createElement("widget");
    widget.setAttribute("type", "ComboBox");
    item->appendChild(widget);
    
    QDomElement objectName = projectDoc->createElement("objectName");
    QDomText objectNameValue = projectDoc->createTextNode(cxObjectName);
    objectName.appendChild(objectNameValue);
    widget.appendChild(objectName);

    QDomElement currentIndex = projectDoc->createElement("currentIndex");
    QDomText currentIndexValue = projectDoc->createTextNode(cxValue);
    currentIndex.appendChild(currentIndexValue);
    widget.appendChild(currentIndex);
  }
}

void ProjectIO::readFromProject(QDomDocument *projectDoc, QDomElement *item)
{
  // Radio buttons:
  //----------------
#if WITH_QT5 
  QList<QRadioButton *> allRadioButtons = parentWidget->findChildren<QRadioButton *>(QString()); 
#else
  QList<QRadioButton *> allRadioButtons = parentWidget->findChildren<QRadioButton *>(); 
#endif

  QList<QString> rbObjectNames;
  for(int i = 0; i < allRadioButtons.size(); i++) 
    rbObjectNames.append(allRadioButtons.at(i)->objectName());

  QDomElement widget = item->firstChildElement("widget");
  for( ; !widget.isNull(); widget = widget.nextSiblingElement()) {

    QString type = widget.attribute("type").trimmed();

    if(type != "RadioButton")
      continue;

    QString objectName = widget.firstChildElement("objectName").text().trimmed();
    bool isChecked = (widget.firstChildElement("isChecked").text().toInt() > 0);

    if(objectName.isEmpty())
      continue;

    int index = rbObjectNames.indexOf(objectName);
    
    if(index < 0) {
      cout << "Load project: RadioButton: mismatch with object name" << endl;
#if WITH_QT5
      cout << "*** " << string(objectName.toLatin1()) << " ***" << endl;
#else
      cout << "*** " << string(objectName.toAscii()) << " ***" << endl;
#endif
      return;
    }

    QRadioButton *rb = allRadioButtons.at(index);
    rb->setChecked(isChecked);
  }

  // Check boxes:
  //--------------
  QList<QCheckBox *> allCheckBoxes = parentWidget->findChildren<QCheckBox *>(); 

  QList<QString> cbObjectNames;
  for(int i = 0; i < allCheckBoxes.size(); i++)
    cbObjectNames.append(allCheckBoxes.at(i)->objectName());

  widget = item->firstChildElement("widget");
  for( ; !widget.isNull(); widget = widget.nextSiblingElement()) {

    QString type = widget.attribute("type").trimmed();

    if(type != "CheckBox")
      continue;

    QString objectName = widget.firstChildElement("objectName").text().trimmed();
    bool isChecked = (widget.firstChildElement("isChecked").text().toInt() > 0);

    if(objectName.isEmpty())
      continue;

    int index = cbObjectNames.indexOf(objectName);
    
    if(index < 0) {
      cout << "Load project: Check box: mismatch with object name" << endl;
#if WITH_QT5
      cout << "*** " << string(objectName.toLatin1()) << " ***" << endl;
#else
      cout << "*** " << string(objectName.toAscii()) << " ***" << endl;
#endif
      return;
    }

    QCheckBox *cb = allCheckBoxes.at(index);
    cb->setChecked(isChecked);
  }

  // Line edits:
  //-------------
  QList<QLineEdit *> allLineEdits = parentWidget->findChildren<QLineEdit *>(); 

  QList<QString> leObjectNames;
  for(int i = 0; i < allLineEdits.size(); i++)
    leObjectNames.append(allLineEdits.at(i)->objectName());

  widget = item->firstChildElement("widget");
  for( ; !widget.isNull(); widget = widget.nextSiblingElement()) {

    QString type = widget.attribute("type").trimmed();

    if(type != "LineEdit")
      continue;

    QString objectName = widget.firstChildElement("objectName").text().trimmed();
    QString text = widget.firstChildElement("text").text().trimmed();

    if(objectName.isEmpty())
      continue;

    int index = leObjectNames.indexOf(objectName);
    
    if(index < 0) {
      cout << "Load project: LineEdit: mismatch with object name" << endl;
#if WITH_QT5
      cout << "*** " << string(objectName.toLatin1()) << " ***" << endl;
#else
      cout << "*** " << string(objectName.toAscii()) << " ***" << endl;
#endif
      return;
    }

    QLineEdit *le = allLineEdits.at(index);
    le->setText(text);
  }

  // Text edits:
  //-------------
  QList<QTextEdit *> allTextEdits = parentWidget->findChildren<QTextEdit *>(); 

  QList<QString> teObjectNames;

  for(int i = 0; i < allTextEdits.size(); i++)
    teObjectNames.append(allTextEdits.at(i)->objectName());

  widget = item->firstChildElement("widget");
  for( ; !widget.isNull(); widget = widget.nextSiblingElement()) {

    QString type = widget.attribute("type").trimmed();

    if(type != "TextEdit")
      continue;

    QString objectName = widget.firstChildElement("objectName").text().trimmed();
    QString text = widget.firstChildElement("text").text().trimmed();

    if(objectName.isEmpty())
      continue;

    int index = teObjectNames.indexOf(objectName);
    
    if(index < 0) {
      cout << "Load project: TextEdit: mismatch with object name" << endl;
#if WITH_QT5
      cout << "*** " << string(objectName.toLatin1()) << " ***" << endl;
#else
      cout << "*** " << string(objectName.toAscii()) << " ***" << endl;
#endif
      return;
    }

    QTextEdit *te = allTextEdits.at(index);
    te->clear();
    te->append(text);
  }

  // Combo boxes:
  //--------------
#if WITH_QT5
  QList<QComboBox *> allComboBoxes = parentWidget->findChildren<QComboBox *>(QString()); 
#else
  QList<QComboBox *> allComboBoxes = parentWidget->findChildren<QComboBox *>(); 
#endif

  QList<QString> cxObjectNames;
  for(int i = 0; i < allComboBoxes.size(); i++)
    cxObjectNames.append(allComboBoxes.at(i)->objectName());

  widget = item->firstChildElement("widget");
  for( ; !widget.isNull(); widget = widget.nextSiblingElement()) {

    QString type = widget.attribute("type").trimmed();

    if(type != "ComboBox")
      continue;

    QString objectName = widget.firstChildElement("objectName").text().trimmed();
    int currentIndex = widget.firstChildElement("currentIndex").text().toInt();

    if(objectName.isEmpty())
      continue;

    int index = cxObjectNames.indexOf(objectName);
    
    if(index < 0) {
      cout << "Load project: Combo box: mismatch with object name" << endl;
#if WITH_QT5
      cout << "*** " << string(objectName.toLatin1()) << " ***" << endl;
#else
      cout << "*** " << string(objectName.toAscii()) << " ***" << endl;
#endif
      return;
    }

    QComboBox *cx = allComboBoxes.at(index);
    cx->setCurrentIndex(currentIndex);
  }
}
