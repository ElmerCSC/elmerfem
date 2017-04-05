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

#include <QtGui>
#include <iostream>
#include "dynamiceditor.h"

using namespace std;

DynLineEdit::DynLineEdit(QWidget *parent) : QWidget(parent)
{
   name = "";
   label = NULL;
   frame  = NULL;
   layout = NULL;
   textEdit = NULL;
   closeButton = NULL; 
   lineEdit = new QLineEdit;
}

DynLineEdit::~DynLineEdit()
{
}

DynamicEditor::DynamicEditor(QWidget *parent)
  : QWidget(parent)
{
  newIcon = QIcon(":/icons/document-new.png");
  addIcon = QIcon(":/icons/list-add.png");
  okIcon = QIcon(":/icons/dialog-ok-apply.png");
  removeIcon = QIcon(":/icons/list-remove.png");
  setWindowFlags(Qt::Window);

  menuAction = NULL;
  ID = -1;
  touched = false;

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

//----------------------------------------------------------------------------
DynamicEditor::~DynamicEditor()
{
}

//----------------------------------------------------------------------------
void DynamicEditor::setupTabs(QDomDocument *elmerDefs, const QString &Section, int ID)
{
  // Clear:
  //-------
  this->ID = ID;

  hash.clear();

  QLayout *layout = this->layout();
  if(layout != NULL) {
    QLayoutItem *item;
    while((item = layout->takeAt(0)) != 0)
      delete item;
    if(tabWidget != NULL) {
      tabWidget->clear();
      delete tabWidget;
    }
    delete layout;
  }

  // Get root element of elmerDefs:
  //-------------------------------
  root = elmerDefs->documentElement();

  tabWidget = new QTabWidget;
  //tabWidget->setTabShape(QTabWidget::Triangular);
  tabWidget->setUsesScrollButtons(true);
  tabWidget->setElideMode(Qt::ElideNone);
  all_stuff = root.firstChildElement("ALL");
  element = root.firstChildElement("PDE");

  tabs = 0;

  while(!element.isNull()) {

    name = element.firstChildElement("Name");

    QGridLayout *grid = new QGridLayout;

    int params = 0;

    for( int iter=0; iter<2; iter++ )
    {
      if ( iter==0 ) {
        if ( name.text().trimmed() == "General" ) continue;
        section = all_stuff.firstChildElement(Section);
      } else  {
        section = element.firstChildElement(Section);
      }

      param = section.firstChildElement("Parameter");
      
      // ML: Added argument "Parameter" for nextSiblingElement(), 5. August 2010:
      for( ; !param.isNull(); param=param.nextSiblingElement("Parameter"), params++ ) {

        // label
        QString widget_type = param.attribute("Widget","Edit");
        QString widget_enabled = param.attribute("Enabled","True");
        QString widget_visible = param.attribute("Visible","True");
        QString paramType = param.firstChildElement("Type").text().trimmed();
        QString labelName = param.firstChildElement("Name").text().trimmed();
        QString sifName   = param.firstChildElement("SifName").text().trimmed();
        if ( sifName == "" ) sifName = labelName;
        QString paramDefault = param.firstChildElement("DefaultValue").text().trimmed();
        QString whatis    = param.firstChildElement("Whatis").text().trimmed();
        QString statusTip = param.firstChildElement("StatusTip").text().trimmed();
        QString fullName  = "/"+name.text().trimmed()+"/"+Section+"/"+labelName+"/"+QString::number(ID);
        h.widget = NULL;

        if ( widget_type == "Edit" ) {
          DynLineEdit *edit = new DynLineEdit;
          h.widget = edit->lineEdit;
          edit->lineEdit->setText(paramDefault);
          edit->name = fullName;
          connect(edit->lineEdit, SIGNAL(returnPressed()),
		  edit, SLOT(editSlot()));
          connect(edit->lineEdit, SIGNAL(textChanged(QString)),
		  this, SLOT(textChangedSlot(QString)));

        } else if (widget_type == "TextEdit") {
	  QTextEdit *textEdit = new QTextEdit;
	  // set height to 5..8 lines of current font:
	  QFont currentFont = textEdit->currentFont();
	  QFontMetrics fontMetrics(currentFont);
	  int fontHeight = fontMetrics.height();
	  textEdit->setMinimumHeight(5*fontHeight);
	  textEdit->setMaximumHeight(8*fontHeight);
	  h.widget = textEdit;

	} else if ( widget_type == "Combo" ) {
          QComboBox *combo = new QComboBox;
          h.widget = combo;

          // combo->setObjectName(labelName);  // removed 30. sept. 2008, ML
          int count = 0, active=0;

          QDomElement item = param.firstChildElement("Item");
          for( ; !item.isNull(); item=item.nextSiblingElement("Item") ) {
            QString itemType = item.attribute( "Type", "" );
            if ( itemType == "Active" ) active=count;
            QDomElement itemName = item.firstChildElement("Name");
            combo->insertItem(count++,itemName.text().trimmed() );
          } 
          combo->setCurrentIndex(active);
          // connect(combo, SIGNAL(activated(QString)), this, SLOT(comboSlot(QString)));
	  connect(combo, SIGNAL(currentIndexChanged(QString)), this, SLOT(comboSlot(QString)));

        } else if ( widget_type == "CheckBox" ) {
          QCheckBox *l = new QCheckBox;
          h.widget = l;

          l->setText("");
          l->setChecked(false);
          if ( paramDefault == "True" ) l->setChecked(true);
          connect(l, SIGNAL(stateChanged(int)), this, SLOT(lSlot(int)));

        } else if ( widget_type == "Label" ) {
          QLabel *label = new QLabel;
          QFont font;
          font.setBold(true);
          font.setUnderline(true);
          label->setFont(font);
          label->setText(labelName);
          h.widget = label;
        }

        if ( h.widget ) {
          h.widget->setWhatsThis(whatis);
          h.widget->setStatusTip(statusTip);
          h.widget->setProperty( "dom address",fullName);
          h.elem = param;

          if ( widget_enabled == "False" ) h.widget->setEnabled(false);

	  if(widget_type != "TextEdit") h.widget->setFixedHeight(18);

	  if(widget_type == "TextEdit") {
            QLabel *textEditLabel = new QLabel;
            textEditLabel->setText(labelName);
            h.label = textEditLabel;
	    grid->addWidget(h.widget, params, 0, 1, 2);

            if ( widget_visible == "False" ) {
              h.label->hide();
              h.widget->hide();
            }

	  } else if ( widget_type != "Label" ) {
            QLabel *label = new QLabel;
            label->setText(labelName);
            h.label = label;
            grid->addWidget(h.label,  params, 0);
            grid->addWidget(h.widget, params, 1);

            if ( widget_visible == "False" ) {
              h.label->hide();
              h.widget->hide();
            }
          } else {
            h.label = NULL;
            grid->addWidget(h.widget, params, 0);
          }
          hash[fullName] = h;
        }
      }
    }

    // add a dummy widget in the left bottom corner:
    QWidget *dummyWidget = new QWidget;
    grid->addWidget(dummyWidget, params, 0);
    grid->setRowStretch(params, 1);
    
    // put the grid in a widget:
    QWidget *frmWidget = new QWidget;
    frmWidget->setLayout(grid);
    
    // set up the scroll area:
    QScrollArea *src = new QScrollArea;
    src->setWidget(frmWidget);
    src->setMinimumHeight(300);
    src->setWidgetResizable(true);
    
    // add the scroll area to the tab:
    if (params>0)
      tabWidget->addTab(src, name.text().trimmed());
    
    tabs++;
    element = element.nextSiblingElement("PDE"); // ML: Added "PDE" 5. August 2010
  }

  // Buttons:
  //----------
  QLabel *lbl = new QLabel;
  lbl->setText("Name:");

  nameEdit  = new QLineEdit;
  nameEdit->setText(Section + " " + QString::number(ID+1));

  applyButton = new QPushButton(tr("&Add"));
  applyButton->setIcon(addIcon);
  connect(applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClicked()));
  
  discardButton = new QPushButton(tr("&Remove"));
  discardButton->setIcon(removeIcon);
  connect(discardButton, SIGNAL(clicked()), this, SLOT(discardButtonClicked()));

  okButton = new QPushButton(tr("&OK"));
  okButton->setIcon(okIcon);
  connect(okButton, SIGNAL(clicked()), this, SLOT(okButtonClicked()));

  newButton = new QPushButton(tr("&New"));
  newButton->setIcon(newIcon);
  connect(newButton, SIGNAL(clicked()), this, SLOT(newButtonClicked()));

  QHBoxLayout *nameLayout = new QHBoxLayout;
  nameLayout->addWidget(lbl);
  nameLayout->addWidget(nameEdit);

  QHBoxLayout *buttonLayout = new QHBoxLayout;
  buttonLayout->addWidget(newButton);
  buttonLayout->addWidget(applyButton);
  buttonLayout->addWidget(okButton);
  buttonLayout->addWidget(discardButton);

  QHBoxLayout *spareButtonLayout = new QHBoxLayout;
  spareButton = new QPushButton(tr("SpareButton"));;
  spareButton->setVisible(false);
  spareButtonLayout->addWidget(spareButton);
  connect(spareButton, SIGNAL(clicked()), this, SLOT(spareButtonClicked()));  

  spareScroll = new QScrollArea;
  spareScroll->hide();

  // Main layout:
  //-------------
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(tabWidget);
  mainLayout->addWidget(spareScroll);
  mainLayout->addLayout(spareButtonLayout);
  mainLayout->addLayout(nameLayout);
  mainLayout->addLayout(buttonLayout);
  setLayout(mainLayout);

  // Window title:
  //---------------
  setWindowTitle(Section);
}

//----------------------------------------------------------------------------
void DynLineEdit::editSlot()
{
  QLineEdit *q =  lineEdit;
  QString s = q->text();
#if WITH_QT5
  cout << string(s.toLatin1()) << endl;
#else
  cout << string(s.toAscii()) << endl;
#endif

  if ( frame ) {
    frame->show();
    frame->raise();
    return;
  }

  textEdit = new QTextEdit;
  textEdit->setLineWrapMode(QTextEdit::NoWrap);

  s.replace( ';', '\n' );
  textEdit->append(s);

  closeButton = new QPushButton(tr("&Close"));
  connect(closeButton, SIGNAL(clicked()), this, SLOT(lineEditClose()));

  label = new QLabel;
  label->setText(name);
  
  layout = new QVBoxLayout;
  layout->addWidget(label);
  layout->addWidget(textEdit);
  layout->addWidget(closeButton);

  frame = new QFrame;
  frame->setLayout(layout);
  frame->show();
  frame->setWindowTitle(name);
}

//----------------------------------------------------------------------------
void DynLineEdit::lineEditClose()
{
  QString q = textEdit->toPlainText();
  q.replace( '\n', ';' );

  lineEdit->setText(q);

  frame->close();

  name = "";

  delete label;
  label = NULL;

  delete textEdit;
  textEdit = NULL;

  delete closeButton;
  closeButton = NULL;

  delete layout;
  layout = NULL;

  delete frame;
  frame = NULL;
}

//----------------------------------------------------------------------------
void DynamicEditor::lSlot(int state)
{
  QDomElement param;
  QString q = QObject::sender()->property("dom address").toString();

  int ind = q.lastIndexOf( '/', -1); 
  QString ID = q.mid(ind,-1);

  param = hash[q].elem.firstChildElement("Activate");
  for( ;!param.isNull(); param=param.nextSiblingElement("Activate") ) {
    q = param.text().trimmed() + ID;
    hash[q].widget->setEnabled(state);
    QString widget_visible = hash[q].elem.attribute("Visible","Unknown");
    if ( state == false && widget_visible != "Unknown" ) {
      hash[q].label->hide(); hash[q].widget->hide();
    } else {
      hash[q].label->show(); hash[q].widget->show();
    }
  }

  param = hash[q].elem.firstChildElement("Deactivate");
  for( ;!param.isNull(); param=param.nextSiblingElement("Deactivate") ) {
    q = param.text().trimmed() + ID;
    hash[q].widget->setEnabled(!state);
    QString widget_visible = hash[q].elem.attribute("Visible","Unknown");
    if ( state == true && widget_visible != "Unknown" ) {
      hash[q].label->hide(); hash[q].widget->hide();
    } else {
      hash[q].label->show(); hash[q].widget->show();
    }
  }
}

//----------------------------------------------------------------------------
void DynamicEditor::textChangedSlot(QString text)
{
  QDomElement param;
  QString q = QObject::sender()->property("dom address").toString();

  int ind = q.lastIndexOf( '/', -1); 
  QString ID = q.mid(ind,-1);

  param = hash[q].elem.firstChildElement("Activate");
  for( ;!param.isNull(); param=param.nextSiblingElement("Activate") ) {
    q = param.text().trimmed() + ID;
    QString widget_visible = hash[q].elem.attribute("Visible","Uknown");

    if ( text != "" ) {
      hash[q].widget->setEnabled(true);
      hash[q].widget->show();
      hash[q].label->show();
    } else {
      hash[q].widget->setEnabled(false);
      if ( widget_visible != "Unknown" ) {
        hash[q].label->hide();
        hash[q].widget->hide();
      }
    }
  }
}

//----------------------------------------------------------------------------
void DynamicEditor::comboSlot(QString select)
{
  QString q = QObject::sender()->property("dom address").toString();
  QDomElement item;

  int ind = q.lastIndexOf( '/', -1); 
  QString ID = q.mid(ind,-1);

  item = hash[q].elem.firstChildElement("Item");
  for( ;!item.isNull(); item=item.nextSiblingElement("Item") ) {
    QDomElement itemName = item.firstChildElement("Name");
    if ( itemName.text().trimmed() != select ) {
      QDomElement activ;

      activ = item.firstChildElement("Activate");
      for( ;!activ.isNull(); activ=activ.nextSiblingElement("Activate") ) {
        QString s=activ.text().trimmed() + ID;
        hash_entry_t h = hash[s];

        QString widget_enabled = h.elem.attribute("Enabled","True");
        QString widget_visible = h.elem.attribute("Visible","Unknown");

        h.widget->setEnabled(false);
        if ( widget_visible != "Unknown" ) {
          h.label->hide(); h.widget->hide();
        }
      }
    }
  }

  item = hash[q].elem.firstChildElement("Item");
  for( ;!item.isNull(); item=item.nextSiblingElement("Item") ) {
    QDomElement itemName = item.firstChildElement("Name");
    if ( itemName.text().trimmed() == select ) {
      QDomElement activ;
      activ = item.firstChildElement("Activate");
      for( ;!activ.isNull(); activ=activ.nextSiblingElement("Activate") ) {
        QString s=activ.text().trimmed() + ID;
        hash_entry_t h = hash[s];
        h.widget->setEnabled(true);
        h.label->show(); h.widget->show();
      }
    }
  }
  // this->show();  // Removed 30. sept. 2008, ML
}


//----------------------------------------------------------------------------
QSize DynamicEditor::minimumSizeHint() const
{
  return QSize(128, 128);
}

//----------------------------------------------------------------------------
QSize DynamicEditor::sizeHint() const
{
  return QSize(400, 500);
}

//----------------------------------------------------------------------------
void DynamicEditor::spareButtonClicked()
{
  emit(dynamicEditorSpareButtonClicked(tabWidget->currentIndex(), ID));
}

//----------------------------------------------------------------------------
void DynamicEditor::applyButtonClicked()
{

  cout << "Dynamic editor: apply-button clicked" << endl;
  cout.flush();

  touched = true;

  emit(dynamicEditorReady(MAT_APPLY, ID));
}


//----------------------------------------------------------------------------
void DynamicEditor::discardButtonClicked()
{

  cout << "Dynamic editor: Remove-button clicked" << endl;
  cout.flush();

  touched = false;

  emit(dynamicEditorReady(MAT_DELETE, ID));
}

//----------------------------------------------------------------------------
void DynamicEditor::okButtonClicked()
{

  // cout << "Dynamic editor: ok-button clicked" << endl;
  // cout.flush();

  touched = false;

  emit(dynamicEditorReady(MAT_OK, ID));
}

//----------------------------------------------------------------------------
void DynamicEditor::newButtonClicked()
{
  cout << "Dynamic editor: next-button clicked" << endl;
  cout.flush();

  touched = false;

  emit(dynamicEditorReady(MAT_NEW, ID));
}

//----------------------------------------------------------------------------
void DynamicEditor::dumpHash(QDomDocument *projectDoc, QDomElement *item)
{
  for(int j = 0; j < this->hash.count(); j++) {
    QString key = this->hash.keys().at(j);
    hash_entry_t value = this->hash.values().at(j);
    QDomElement elem = value.elem;
    QWidget *widget = value.widget;
    
    QDomElement itemWidget = projectDoc->createElement("widget");
    item->appendChild(itemWidget);
    
    QDomElement itemKey = projectDoc->createElement("key");
    QDomText itemKeyValue = projectDoc->createTextNode(key);
    itemKey.appendChild(itemKeyValue);
    itemWidget.appendChild(itemKey);
    
    if(elem.attribute("Widget") == "CheckBox") {
      QCheckBox *checkBox = (QCheckBox*)widget;
      QDomElement itemCheckBox = projectDoc->createElement("value");
      QDomText itemCheckBoxValue = projectDoc->createTextNode(QString::number(checkBox->isChecked()));
      itemCheckBox.appendChild(itemCheckBoxValue);
      itemWidget.appendChild(itemCheckBox);
      itemWidget.setAttribute("type", "CheckBox");
      
    } else if(elem.attribute("Widget") == "Edit") {
      QLineEdit *lineEdit = (QLineEdit*)widget;
      QDomElement itemLineEdit = projectDoc->createElement("value");
      QDomText itemLineEditValue = projectDoc->createTextNode(lineEdit->text().trimmed());
      itemLineEdit.appendChild(itemLineEditValue);
      itemWidget.appendChild(itemLineEdit);
      itemWidget.setAttribute("type", "Edit");

    } else if(elem.attribute("Widget") == "TextEdit") {
      QTextEdit *textEdit = (QTextEdit*)widget;
      QDomElement itemTextEdit = projectDoc->createElement("value");
      QDomText itemTextEditValue = projectDoc->createTextNode(textEdit->toPlainText());
      itemTextEdit.appendChild(itemTextEditValue);
      itemWidget.appendChild(itemTextEdit);
      itemWidget.setAttribute("type", "TextEdit");
      
    } else if(elem.attribute("Widget") == "Combo") {
      QComboBox *comboBox = (QComboBox*)widget;
      QDomElement itemComboBox = projectDoc->createElement("value");
      QDomText itemComboBoxValue = projectDoc->createTextNode(comboBox->currentText().trimmed());
      itemComboBox.appendChild(itemComboBoxValue);
      itemWidget.appendChild(itemComboBox);
      itemWidget.setAttribute("type", "Combo");
      
    } else if(elem.attribute("Widget") == "Label") {
      QLabel *label = (QLabel*)widget;
      QDomElement itemLabel = projectDoc->createElement("value");
      QDomText itemLabelValue = projectDoc->createTextNode(label->text().trimmed());
      itemLabel.appendChild(itemLabelValue);
      itemWidget.appendChild(itemLabel);
      itemWidget.setAttribute("type", "Label");
    }
  }
}

void DynamicEditor::populateHash(QDomElement *item)
{
  QDomElement widget = item->firstChildElement("widget");
  
  // ML: Added argument "widget" for nextSiblingElement(), 5. August 2010:
  for(; !widget.isNull(); widget = widget.nextSiblingElement("widget")) {
    QString type = widget.attribute("type").trimmed();
    QString key = widget.firstChildElement("key").text().trimmed();
    QString value = widget.firstChildElement("value").text().trimmed();
    
    if(value.isEmpty())
      continue;
    
    QStringList splittedKey = key.split("/");
    
    // Compare with current hash:
    //----------------------------
    bool match_found = false;
    for(int j = 0; j < this->hash.count(); j++) {
      QString hashkey = this->hash.keys().at(j);
      QStringList splittedHashKey = hashkey.split("/");
      hash_entry_t hashvalue = this->hash.values().at(j);
      QWidget *widget = hashvalue.widget;
      QDomElement elem = hashvalue.elem;
      
      if((splittedKey.at(1) == splittedHashKey.at(1)) &&
	 (splittedKey.at(2) == splittedHashKey.at(2)) &&
	 (splittedKey.at(3) == splittedHashKey.at(3))) {
	
	match_found = true;
	
	if(elem.attribute("Widget") == "CheckBox") {
	  if(type != "CheckBox")
	    cout << "Load project: type mismatch with checkBox" << endl;
	  QCheckBox *checkBox = (QCheckBox*)widget;
	  if(value.toInt() == 1)
	    checkBox->setChecked(true);
	  else
	    checkBox->setChecked(false);
	  
	} else if(elem.attribute("Widget") == "Edit") {
	  if(type != "Edit")
	    cout << "Load project: type mismatch with lineEdit" << endl;
	  QLineEdit *lineEdit = (QLineEdit*)widget;
	  lineEdit->setText(value);
	  
	} else if(elem.attribute("Widget") == "TextEdit") {
	  if(type != "TextEdit")
	    cout << "Load project: type mismatch with textEdit" << endl;
	  QTextEdit *textEdit = (QTextEdit*)widget;
	  textEdit->clear();
	  textEdit->append(value);
	  
	} else if(elem.attribute("Widget") == "Combo") {
	  if(type != "Combo")
	    cout << "Load project: type mismatch with comboBox" << endl;
	  QComboBox *comboBox = (QComboBox*)widget;
	  for(int k = 0; k < comboBox->count(); k++) {
	    QString current = comboBox->itemText(k).trimmed();
	    if(current == value.trimmed())
	      comboBox->setCurrentIndex(k);
	  }
	}
      }
    }
    
    if(!match_found) {
#if WITH_QT5
      cout << "Error: Unable to set menu entry: key: " << key.toLatin1().data() << endl;
#else
      cout << "Error: Unable to set menu entry: key: " << key.toAscii().data() << endl;
#endif
      cout.flush();
    }
  }
}
