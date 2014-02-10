/*  
   Elmer, A Finite Element Software for Multiphysical Problems
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/***********************************************************************
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#ifndef EIOSOLVERAGENT_H
#define EIOSOLVERAGENT_H

#include "EIOModelManager.h"

const int solverFiles = 4;

class EIOSolverAgent
{
public:
  EIOSolverAgent(EIOModelManager *mm);
  ~EIOSolverAgent();

  int createSolver();
  int openSolver();
  int closeSolver();

  int writeDescription(int& linsys, int& procs);
  int readDescription(int& linsys, int& procs);
  int writeSolverRecord(int& equation,
			int& main_type,
			int& sub_type,
			int& precond_type,
			int& stabilization,
			int& max_iter,
			double& stop_tol,
			double& steady_stop_tol,
			int& linearization,
			int& lin_max_iter,
			double& lin_stop_tol,
			int& lin_use_picard,
			int& lin_use_newton,
			int& newton_after_iter,
			double& newton_after_tol);
  int readSolverRecord(int& equation,
		       int& main_type,
		       int& sub_type,
		       int& precond_type,
		       int& stabilization,
		       int& max_iter,
		       double& stop_tol,
		       double& steady_stop_tol,
		       int& linearization,
		       int& lin_max_iter,
		       double& lin_stop_tol,
		       int& lin_use_picard,
		       int& lin_use_newton,
		       int& newton_after_iter,
		       double& newton_after_tol);
  int writeTimestepDescription(int& dependence,int& rlen);
  int readTimestepDescription(int& dependence,int& rlen);
  int writeTimestepRecord(int& type,
			       int *nof_timesteps,
			       double *timestep_sizes,
			       int *output_intervals,
			       int& steady_max_iter);
  int readTimestepRecord(int& type,
			       int *nof_timesteps,
			       double *timestep_sizes,
			       int *output_intervals,
			       int& steady_max_iter);
private:
  // We "use" ModelManager in every Agent.
  EIOModelManager *manager;

  // All streams
  fstream solverFileStream[solverFiles];
  // Sizes
  int len;
};

#endif /* EIOSOLVERAGENT_H */
