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
 *  ElmerGUI readepfile                                                      *
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
#if WITH_QT5
  #include <QtWidgets>
#endif
#include <QtGui>
#include <iostream>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include "readepfile.h"

using namespace std;

ReadEpFile::ReadEpFile(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.browseButton, SIGNAL(clicked()), this, SLOT(browseButtonClickedSlot()));
  connect(ui.applyButton, SIGNAL(clicked()), this, SLOT(applyButtonClickedSlot()));
  connect(ui.cancelButton, SIGNAL(clicked()), this, SLOT(cancelButtonClickedSlot()));
  connect(ui.okButton, SIGNAL(clicked()), this, SLOT(okButtonClickedSlot()));
  connect(ui.allButton, SIGNAL(clicked()), this, SLOT(allButtonClickedSlot()));

  ui.nodesEdit->setEnabled(false);
  ui.elementsEdit->setEnabled(false);
  ui.timestepsEdit->setEnabled(false);
  ui.dofsEdit->setEnabled(false);

  setWindowTitle("Read input file");
  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

ReadEpFile::~ReadEpFile()
{
}

void ReadEpFile::browseButtonClickedSlot()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Select input file"), "", tr("Postprocessor files (*.vtu *.ep);;Paraview files (*.vtu);;ElmerPost files (*.ep)"));

  ui.fileName->setText(fileName.trimmed());

  readHeader();
}

void ReadEpFile::applyButtonClickedSlot()
{
  QString fileName = ui.fileName->text().trimmed();
  
  if(fileName.isEmpty()) return;

  int start = ui.start->value();
  int end = ui.end->value();
  int maxSteps = ui.timestepsEdit->text().toInt();

  if(end > maxSteps) {
    end = maxSteps;
    ui.end->setValue(maxSteps);
  }

  if(start > end) {
    start = end;
    ui.start->setValue(start);
  }
  
  repaint();

  emit(readPostFileSignal(fileName));
}

void ReadEpFile::cancelButtonClickedSlot()
{
  close();
}

void ReadEpFile::okButtonClickedSlot()
{
  applyButtonClickedSlot();
  cancelButtonClickedSlot();
}

void ReadEpFile::readHeader()
{ 
  QString fileName = ui.fileName->text().trimmed();

  QFile postFile(fileName);
  
  if(!postFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
    ui.fileName->setText("");
    return;
  }

  int nodes, elements, timesteps, components;

  if(ui.fileName->text().endsWith(".vtu", Qt::CaseInsensitive)){
  
    vtkXMLUnstructuredGridReader* reader =  vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(ui.fileName->text().toLatin1().data());
    reader->Update();

//potential = output->GetPointData().GetArray("Magnetization")
    //vtkInformation* outInfo = reader()->GetExecutive()->GetOutputInformation(0);
    
	nodes = reader->GetNumberOfPoints();
    elements = reader->GetNumberOfCells();
	components = 1;
    timesteps = reader->GetNumberOfTimeSteps();
	if(timesteps == 0) timesteps = 1;
	components = 0;
    vtkUnstructuredGrid *output = reader->GetOutput();
	vtkPointData *pointData = output->GetPointData();
	vtkCellData *cellData = output->GetCellData();
    
//    cout << "VTU: NumberOfPieces: " << reader->GetNumberOfPieces();
//    cout << "VTU: NumberOfPoints: " << reader->GetNumberOfPoints();
//    cout << "VTU: NumberOfCells: " << reader->GetNumberOfCells();
//    cout << "VTU: NumberOfTimeSteps: " << reader->GetNumberOfTimeSteps();
//    cout << "VTU: NumberOfPointArrays: " << reader->GetNumberOfPointArrays() << endl;
      for(int i = 0; i < reader->GetNumberOfPointArrays(); i++){
//        cout << "VTU: PointArray [" << i << "] " << reader->GetPointArrayName(i) << ", "<< pointData->GetArray(reader->GetPointArrayName(i))->GetNumberOfComponents()<< endl; 
		  components += pointData->GetArray(reader->GetPointArrayName(i))->GetNumberOfComponents();
      }
//      cout << "VTU: NumberOfCellArrays: " <<  reader->GetNumberOfCellArrays() << endl;
//      for(int i = 0; i < reader->GetNumberOfCellArrays(); i++){
//        cout << "VTU: CellArray [" << i << "] " << reader->GetCellArrayName(i) << ", "<< cellData->GetArray(reader->GetCellArrayName(i))->GetNumberOfComponents()<< endl; 
//      }        	
    
  
  }else if(ui.fileName->text().endsWith(".ep", Qt::CaseInsensitive)){

    QTextStream post(&postFile);

    QTextStream txtStream;
    QString tmpLine = post.readLine().trimmed();
    while(tmpLine.isEmpty() || (tmpLine.at(0) == '#'))
      tmpLine = post.readLine().trimmed();
    txtStream.setString(&tmpLine);

    txtStream >> nodes >> elements >> components >> timesteps;

    postFile.close();
  }

  ui.nodesEdit->setText(QString::number(nodes));
  ui.elementsEdit->setText(QString::number(elements));
  ui.timestepsEdit->setText(QString::number(timesteps));
  ui.dofsEdit->setText(QString::number(components));
}

void ReadEpFile::allButtonClickedSlot()
{
  ui.start->setValue(1);
  ui.end->setValue(ui.timestepsEdit->text().toInt());

  repaint();
}
