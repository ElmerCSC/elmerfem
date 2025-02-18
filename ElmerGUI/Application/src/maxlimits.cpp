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

#include "maxlimits.h"

Limit::Limit()
{
  // Default limits:
  max_solvers = 10;
  max_equations = 10;
  max_materials = 10;
  max_bodyforces = 10;
  max_initialconditions = 10;
  max_bcs = 500;
  max_bodies = 100;
  max_boundaries = 500;
}

Limit::~Limit()
{
}

int Limit::maxEquations()
{
  return max_equations;
}

int Limit::maxMaterials()
{
  return max_materials;
}

int Limit::maxBodyforces()
{
  return max_bodyforces;
}

int Limit::maxInitialconditions()
{
  return max_initialconditions;
}

int Limit::maxBcs()
{
  return max_bcs;
}

int Limit::maxBodies()
{
  return max_bodies;
}

int Limit::maxBoundaries()
{
  return max_boundaries;
}

int Limit::maxSolvers()
{
  return max_solvers;
}

void Limit::setMaxEquations(int value)
{
  max_equations = value;
}

void Limit::setMaxMaterials(int value)
{
  max_materials = value;
}

void Limit::setMaxBodyforces(int value)
{
  max_bodyforces = value;
}

void Limit::setMaxInitialconditions(int value)
{
  max_initialconditions = value;
}

void Limit::setMaxBcs(int value)
{
  max_bcs = value;
}

void Limit::setMaxBodies(int value)
{
  max_bodies = value;
}

void Limit::setMaxBoundaries(int value)
{
  max_boundaries = value;
}

void Limit::setMaxSolvers(int value)
{
  max_solvers = value;
}
