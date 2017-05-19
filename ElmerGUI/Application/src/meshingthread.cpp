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
 *  ElmerGUI meshingthread                                                   *
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
#include <stdio.h>
#include "meshingthread.h"

using namespace std;

MeshingThread::MeshingThread(QObject *parent)
  : QThread(parent)
{
  this->setTerminationEnabled(true);
}


MeshingThread::~MeshingThread()
{
}

nglib::Ng_Mesh* MeshingThread::getNgMesh()
{
  return this->ngmesh;
}

void MeshingThread::generate(int generatorType,
			     QString cs,
			     TetlibAPI *tetlibAPI,
			     nglib::Ng_Mesh *ngmesh,
			     nglib::Ng_STL_Geometry *nggeom,
			     nglib::Ng_Geometry_2D *nggeom2d,
			     int ngDim,
			     nglib::Ng_Meshing_Parameters *mp)
{
  this->generatorType = generatorType;
  
  this->tetgenControlString = cs;
  this->tetlibAPI = tetlibAPI;
  this->in = tetlibAPI->in;
  this->out = tetlibAPI->out;

  this->delegate_tetrahedralize = tetlibAPI->delegate_tetrahedralize;

  this->ngmesh = ngmesh;
  this->nggeom = nggeom;
  this->nggeom2d = nggeom2d;
  this->mp = mp;
  this->ngDim = ngDim;

  if (!isRunning()) {
    cout << "Starting meshing thread with low priority" << endl;
    cout.flush();
    start(LowPriority);
  } else {
    cout << "Meshing thread is already running" << endl;
    cout.flush();
    return;
  }
}

void MeshingThread::stopMeshing()
{
  cout << "Terminating meshing thread... ";
  cout.flush();
  
  terminate();
  // wait();

  cout << "done" << endl;
  cout.flush();
}

void MeshingThread::run()
{
  QString qs;
  char ss[1024];

  if(generatorType == GEN_TETLIB) {
    
    cout << "tetlib: control string: " 
#if WITH_QT5
	 << string(tetgenControlString.toLatin1()) << endl;
#else
	 << string(tetgenControlString.toAscii()) << endl;
#endif
    cout << "tetlib: input points: " << in->numberofpoints << endl;
    cout.flush();
    
    out->deinitialize();
    out->initialize();
    
#if WITH_QT5
    sprintf(ss, "%s", (const char*)(tetgenControlString.toLatin1()));
#else
    sprintf(ss, "%s", (const char*)(tetgenControlString.toAscii()));
#endif
    
    if(delegate_tetrahedralize) 
      delegate_tetrahedralize(1, NULL, ss, in, out, NULL, NULL);      
    
    cout << "tetlib: nodes: " << out->numberofpoints << endl;
    cout << "tetlib: elements: " << out->numberoftetrahedra << endl;
    cout << "tetlib: boundary elements: " << out->numberoftrifaces << endl;
    cout.flush();
    
  } else if(generatorType == GEN_NGLIB) {
    
    int rv = 0;
    
    if(ngDim == 3) {

      cout << "3D meshing..." << endl;

      rv = nglib::Ng_STL_GenerateSurfaceMesh(nggeom, ngmesh, mp);
      cout << "Generate Surface Mesh: Ng_result=" << rv << endl;
      
      rv = nglib::Ng_GenerateVolumeMesh(ngmesh, mp);
      cout << "Generate Volume Mesh: Ng_result=" << rv << endl;
      
      int np = nglib::Ng_GetNP(ngmesh);
      cout << "Meshing thread: nodes: " << np << endl;
      
      int ne = nglib::Ng_GetNE(ngmesh);
      cout << "Meshing thread: elements: " << ne << endl;
      
      int nse = nglib::Ng_GetNSE(ngmesh);
      cout << "Meshing thread: boundary elements: " << nse << endl;      
      cout.flush();

    } else if(ngDim == 2) {

      cout << "2D meshing..." << endl;

      rv = nglib::Ng_GenerateMesh_2D(nggeom2d, &ngmesh, mp);
      cout << "Generate 2D Mesh: Ng_result=" << rv << endl;

      int np = nglib::Ng_GetNP_2D(ngmesh);
      cout << "Meshing thread: nodes: " << np << endl;
      
      int ne = nglib::Ng_GetNE_2D(ngmesh);
      cout << "Meshing thread: elements: " << ne << endl;
      
      int nse = nglib::Ng_GetNSeg_2D(ngmesh);
      cout << "Meshing thread: boundary elements: " << nse << endl;      
      cout.flush();
      
    } else {

      cout << "Illegal spatial dimension: " << ngDim << endl;

    }
    
  } else {
    
    cout << "Meshgen: unknown generator type\n";
    cout.flush();
    
  }
}
