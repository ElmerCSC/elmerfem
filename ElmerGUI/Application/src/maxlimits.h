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
 *  ElmerGUI maxlimits                                                       *
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

#ifndef MAXLIMITS_H
#define MAXLIMITS_H

class Limit {
 public:
  Limit();
  ~Limit();

  int maxEquations();
  int maxMaterials();
  int maxBodyforces();
  int maxInitialconditions();
  int maxBcs();
  int maxBodies();
  int maxBoundaries();
  int maxSolvers();

  void setMaxEquations(int);
  void setMaxMaterials(int);
  void setMaxBodyforces(int);
  void setMaxInitialconditions(int);
  void setMaxBcs(int);
  void setMaxBodies(int);
  void setMaxBoundaries(int);
  void setMaxSolvers(int);

 private:
  int max_equations;
  int max_materials;
  int max_bodyforces;
  int max_initialconditions;
  int max_bcs;
  int max_bodies;
  int max_boundaries;
  int max_solvers;
};

#endif // MAXLIMITS_H
