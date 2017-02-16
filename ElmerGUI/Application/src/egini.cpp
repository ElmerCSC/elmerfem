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
 *  ElmerGUI egini                                                           *
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

#include <iostream>
#include "egini.h"

using namespace std;

EgIni::EgIni(QWidget *parent)
  : QDialog(parent)
{
  iniLoaded = false;

  // Determine ini-file location and name:
  //--------------------------------------
  QString elmerGuiHome;
  QString paraViewHome;
  
#ifdef __APPLE__DONTGOHERE_TODO
  //QString iniFileName = this->homePath +  "/edf/egini.xml";          
  QString iniFileName = QDir::homePath() +  "/edf/egini.xml";          
#else
  QString iniFileName = QCoreApplication::applicationDirPath() + "/../share/ElmerGUI/edf/egini.xml";  // @TODO: fix path to share/ElmerGUI/edf

  elmerGuiHome = QString(getenv("ELMERGUI_HOME"));

  if(!elmerGuiHome.isEmpty()) 
    iniFileName = elmerGuiHome + "/edf/egini.xml";  

  paraViewHome = QString(getenv("PARAVIEW_HOME"))+"/bin";
#endif
  
  // Load initialization file:
  //---------------------------
#if WITH_QT5
  cout << "Load " << string(iniFileName.toLatin1()) << "...";
#else
  cout << "Load " << string(iniFileName.toAscii()) << "...";
#endif
  cout.flush();
  
  QFile file(iniFileName);
  QString errStr;
  int errRow;
  int errCol;

  if(!file.exists()) {

    QMessageBox::information(window(), tr("Eg ini-file loader: ") + iniFileName,
			     tr("Initialization file does not exist"));
    return;
    
  } else {  
    
    if(!iniDoc.setContent(&file, true, &errStr, &errRow, &errCol)) {

      QMessageBox::information(window(), tr("Eg ini-file loader: ") + iniFileName,
			       tr("Parse error at line %1, col %2:\n%3")
			       .arg(errRow).arg(errCol).arg(errStr));
      file.close();
      return;
      
    } else {

      if(iniDoc.documentElement().tagName() != "egini") {
	QMessageBox::information(window(), tr("Eg ini-file loader: ") + iniFileName,
				 tr("This is not an eg initialization file"));
	file.close();	
	return;
      }
    }
  }
  
  cout << " done" << endl;
  file.close();
  iniLoaded = true;
}


EgIni::~EgIni()
{
}


bool EgIni::isPresent(QString tag)
{
  if(!iniLoaded)
    return false;

  root = iniDoc.documentElement();
  element = root.firstChildElement(tag);
  
  if(element.isNull())
    return false;

  return true;
}


bool EgIni::isSet(QString tag)
{
  if(!iniLoaded)
    return false;

  root = iniDoc.documentElement();
  element = root.firstChildElement(tag);
  
  if(element.isNull())
    return false;

  if(element.text().trimmed() != "0")
    return true;
  
  return false;
}


QString EgIni::value(QString tag)
{
  if(!iniLoaded)
    return "";

  root = iniDoc.documentElement();
  element = root.firstChildElement(tag);
  
  if(element.isNull())
    return "";

  return element.text().trimmed();
}
