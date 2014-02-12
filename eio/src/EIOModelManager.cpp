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
            Martti Verho 08.10.98  (Win32 related changes)
************************************************************************/

#include "EIOModelManager.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
#include <iostream>
#include <fstream>

#if defined(MINGW32) 
  # include <direct.h> 
  # include <io.h>     
#else                
  #include <unistd.h>
#endif

/*
  The ModelManager presides over the filesystem access:
  - the correctness of paths supplied by the user is verified
  - the location in the filesystem is maintained
  - streams are managed (opened/created/closed) here
 */

/*
  eio_mkdir
  eio_chdir
  eio_checkmodel

  are cover functions for the system level calls.
  Note, that in the current version the error value is handled locally.
  THIS WILL BE CHANGED IN THE FUTURE!

  Note, that there are other system calls in the default constructor.
  They will get their own cover functions later.
 */
int eio_mkdir(const char *dir)
{
  int rc;
  /* #ifndef LINUX_PC 
  extern int errno;
  #endif */

#if defined(MINGW32)
  rc = _mkdir(dir);
#else
  rc = mkdir(dir, S_IRWXU|S_IRWXG);
#endif

  if(rc == -1)
  {
    switch(errno)
	  {
	  case EEXIST:
	    return 1;
	    break;
	  default:
	    std::cerr << "Unexpected error at mkdir" << std::endl;
	    break;
	  }

    return 0;
  }

  return 1;
}

int eio_chdir(const char *dir)
{
  int rc;
  /* #ifndef LINUX_PC
  extern int errno;
  #endif */

#if defined(MINGW32)
  rc = _chdir(dir);
#else
  rc = chdir(dir);
#endif


  if(rc == -1)
    {
      switch(errno)
	{
	case EACCES:
	  std::cerr << "Check permissions: dir " << std::endl;
	  break;
	case EIO:
	  std::cerr << "I/O error: dir " << std::endl;
	  break;	  
	case ENOENT:
	  std::cerr << "No such dir" << std::endl;
	  break;
	case ENOTDIR:
	  std::cerr << "Check path: dir" << std::endl;
	  break;
	default:
	  std::cerr << "Unexpected error at chdir" << std::endl;
	  break;
	}
      return 0;
    }
  return 1;
}

/*
  eio_checkmodel

  has two parts:
   (i) verify the integrity of the path
  (ii) verify the permissions
 */
int eio_checkmodel(const char *model)
{
  int rc;
  /* #ifndef LINUX_PC
  extern int errno;
  #endif */

#if defined(MINGW32)
  struct _stat buf;
  rc = _stat(model, &buf);
#else
  struct stat buf;
  rc = stat(model, &buf);
#endif

  if(rc == -1)
    {
      switch(errno)
	{
	case EACCES:
	  std::cerr << "Check permissions: model " << std::endl;
	  break;

	case EIO:
	  std::cerr << "I/O error: model " << std::endl;
	  break;	  

	case ENOENT:
	  std::cerr << "No such model" << std::endl;
	  break;

	case ENOTDIR:
	  std::cerr << "Check path: model" << std::endl;
	  break;
	  
	default:
	  std::cerr << "Unexpected error at stat" << std::endl;
	  break;
	}
      return 0;
    }

  /*
    Is model a directory?
    */
  int rc_access;

#if defined(MINGW32)
  rc = buf.st_mode & _S_IFDIR;
  if (rc)
    rc_access = _access(model, 06);
#else
  rc = S_ISDIR(buf.st_mode);
    rc_access = access(model, R_OK | W_OK | X_OK);
#endif

  if(rc)
    {
      /*
	We need read/write/exec permissions, however, since we could stat,
	we can search.
	*/
      if(rc_access == -1)
	{
	  std::cerr << "No permission to operate: model" << std::endl;
	  return 0;
	}
    }
  else
    {
      std::cerr << model << " is not a directory" << std::endl;
      return 0;
    }
  return 1;
}



EIOModelManager::EIOModelManager()
{
  /*
    We must remember the current directory so that we can safely return
    to there after the database has been closed.

    We should also get the mask and use it in opening the streams.
    TO BE FIXED SOON.
   */
#if defined(MINGW32)
  _getcwd(rundir, PATH_MAX);
/*  _umask(_S_IWRITE | _S_IREAD);*/
  _umask(0);
#else
  getcwd(rundir, PATH_MAX);
  umask(S_IRWXO);
#endif

}


EIOModelManager::~EIOModelManager()
{
  eio_chdir(rundir);
}

int EIOModelManager::createModel(const char *dir)
{
  strcpy(modeldir, dir);
  //  strcpy(modelname, model);

  if(!eio_chdir(modeldir))
    {
      return -1;
    }
  if(!eio_mkdir(modeldir))
    {
      return -1;
    }
  if(!eio_chdir(modeldir))
    {
      return -1;
    }
  return 0;
}

int EIOModelManager::openModel(const char *dir)
{
  strcpy(modeldir, dir);
  //  strcpy(modelname, dir);

  if(!eio_chdir(modeldir))
    {
      return -1;
    }
  if(!eio_checkmodel(modeldir))
    {
      return -1;
    }
  if(!eio_chdir(modeldir))
    {
      return -1;
    }
  return 0;
}

int EIOModelManager::
closeModel()
{
  return 0;
}

int EIOModelManager::openStream(fstream& fstr, const char *name, int mode)
{
  fstr.open(name, (std::ios::openmode) mode);
//  if(!fstr)
  if(fstr.fail())
    {
      std::cerr << "Could not open " << name << std::endl;
      return 0;
    }
  return 1;
}

int EIOModelManager::closeStream(fstream& fstr)
{
  fstr.close();
  return 1;
}

int EIOModelManager::makeDirectory(const char *dir)
{
  return eio_mkdir(dir);
}

