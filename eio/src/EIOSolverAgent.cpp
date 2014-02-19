/*  
   Elmer, A Finite Element Software for Multiphysical Problems
  
   Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library (in file ../../LGPL-2.1); if not, write 
   to the Free Software Foundation, Inc., 51 Franklin Street, 
   Fifth Floor, Boston, MA  02110-1301  USA
*/

/***********************************************************************
Program:    ELMER Data base interface (EIO)
Author(s):  Harri Hakula 10.03.98
************************************************************************/

#include "EIOSolverAgent.h"
#include <string.h>

extern
void make_filename(char *buf, const char *model, const char *suffix);

static const char *extension[] = {
  "solver.header",
  "solver.records",
  "timestep.header",
  "timestep.records"
};

enum { HEADER = 0, SOLVERS, TIMEHEAD, TIMESTEP};

EIOSolverAgent::EIOSolverAgent(EIOModelManager *mm)
{
  manager = mm;
}

EIOSolverAgent::~EIOSolverAgent()
{
}

int EIOSolverAgent::
createSolver()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < solverFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(solverFileStream[i], extension[i], std::ios::out);
    }

  return 0;
}

int EIOSolverAgent::
openSolver()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < solverFiles; ++i)
    {
      //      make_filename(filename, manager->name(), extension[i]);
      manager->openStream(solverFileStream[i], extension[i], std::ios::in);
    }

  return 0;
}

int EIOSolverAgent::
closeSolver()
{
  int i;
  char filename[PATH_MAX];

  for(i = 0; i < solverFiles; ++i)
    {
      manager->closeStream(solverFileStream[i]);
    } 
  return 0;
}

int EIOSolverAgent::
writeDescription(int& linsys, int& procs)
{
  fstream& str = solverFileStream[HEADER];
  str << linsys << ' ' << procs << '\n';
  return 0;
}

int EIOSolverAgent::
readDescription(int& linsys, int& procs)
{
  fstream& str = solverFileStream[HEADER];
  str >> linsys >> procs;
  return 0;
}

int EIOSolverAgent::
writeSolverRecord(int& equation,
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
		  double& newton_after_tol)
{
  fstream& str = solverFileStream[SOLVERS];
  str << equation << '\n'
    << main_type << '\n' 
    << sub_type << '\n' 
    << precond_type << '\n' 
    << stabilization << '\n'
    << max_iter << '\n' 
    << stop_tol << '\n' 
    << steady_stop_tol << '\n'
    << linearization << '\n' 
    << lin_max_iter << '\n' 
    << lin_stop_tol << '\n' 
    << lin_use_picard << '\n'
    << lin_use_newton << '\n'
    << newton_after_iter << '\n'
    << newton_after_tol << '\n';
  return 0;
}

int EIOSolverAgent::
readSolverRecord(int& equation,
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
		 double& newton_after_tol)
{
  fstream& str = solverFileStream[SOLVERS];
  str >> equation 
    >> main_type  
    >> sub_type  
    >> precond_type  
    >> stabilization 
    >> max_iter  
    >> stop_tol  
    >> steady_stop_tol 
    >> linearization  
    >> lin_max_iter  
    >> lin_stop_tol  
    >> lin_use_picard 
    >> lin_use_newton 
    >> newton_after_iter 
    >> newton_after_tol;
  return 0;
}


int EIOSolverAgent::
writeTimestepDescription(int& dependence, int& reclen)
{
  fstream& str = solverFileStream[TIMEHEAD];
  len = reclen;
  str << dependence << ' ' << len << '\n';
  return 0;
}


int EIOSolverAgent::
readTimestepDescription(int& dependence, int& reclen)
{
  fstream& str = solverFileStream[TIMEHEAD];
  str >> dependence >> len;
  reclen = len;
  return 0;
}


int EIOSolverAgent::
writeTimestepRecord(int& type,
		    int *nof_timesteps,
		    double *timestep_sizes,
		    int *output_intervals,
		    int& steady_max_iter)
{
  int i;
  fstream& str = solverFileStream[TIMESTEP];
  str << type << '\n';
  for(i = 0; i < len; ++i)
    {
      str << nof_timesteps[i] << ' ';
    }
  str << '\n';
  for(i = 0; i < len; ++i)
    {
      str << timestep_sizes[i] << ' ';
    }
  str << '\n';
  for(i = 0; i < len; ++i)
    {
      str << output_intervals[i] << ' ';
    }
  str << '\n';
  str << steady_max_iter;
  str << '\n';
  return 0;
}


int EIOSolverAgent::
readTimestepRecord(int& type,
		   int *nof_timesteps,
		   double *timestep_sizes,
		   int *output_intervals,
		   int& steady_max_iter)
{
  int i;
  fstream& str = solverFileStream[TIMESTEP];

  str >> type;
  for(i = 0; i < len; ++i)
    {
      str >> nof_timesteps[i];
    }
  for(i = 0; i < len; ++i)
    {
      str >> timestep_sizes[i];
    }
  for(i = 0; i < len; ++i)
    {
      str >> output_intervals[i];
    }
  str >> steady_max_iter;

  return 0;
}
