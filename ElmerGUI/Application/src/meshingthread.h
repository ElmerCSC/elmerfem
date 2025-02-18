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

#ifndef MESHINGTHREAD_H
#define MESHINGTHREAD_H

#include <QThread>

#ifdef WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include "plugins/tetlib_api.h"
#include "plugins/nglib_api.h"

namespace nglib {
#include "nglib.h"
}

class MeshingThread : public QThread
{
  Q_OBJECT

public:
  MeshingThread(QObject *parent = 0);
  ~MeshingThread();

  void generate(int generatorType, QString cs,
		TetlibAPI *tetlibAPI, 
		nglib::Ng_Mesh *ngmesh,
		nglib::Ng_STL_Geometry *nggeom, 
		nglib::Ng_Geometry_2D *nggeom2d,
		int ngDim, nglib::Ng_Meshing_Parameters *mp);

  void stopMeshing();

  nglib::Ng_Mesh *getNgMesh();

protected:
  void run();
  
private:
  int generatorType;

  // tetlib:
  QString tetgenControlString;
  TetlibAPI *tetlibAPI;
  tetgenio *in;
  tetgenio *out;
  delegate_tetrahedralize_t delegate_tetrahedralize;

  // nglib:
  NglibAPI *nglibAPI;
  nglib::Ng_Mesh *ngmesh;
  nglib::Ng_STL_Geometry *nggeom;
  nglib::Ng_Geometry_2D *nggeom2d;
  nglib::Ng_Meshing_Parameters *mp;
  int ngDim;
};

#endif // MESHINGTHREAD_H
