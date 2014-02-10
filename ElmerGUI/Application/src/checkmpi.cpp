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
 *  ElmerGUI checkmpi                                                        *
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
#include "checkmpi.h"

using namespace std;

CheckMpi::CheckMpi()
{
}

CheckMpi::~CheckMpi()
{
}

int CheckMpi::findSmpd()
{
  // Check if there is an instance of smpd running:

#ifdef WIN32

  unsigned int i;
  TCHAR szProcessName[50] = TEXT("");
  DWORD ProcessesIDs[MAX_PROCIDS], cbNeeded, cProcesses;
  bool found = false;

  cout << "Checking whether smpd is running... ";

  if(!EnumProcesses(ProcessesIDs, sizeof(ProcessesIDs), &cbNeeded)) {
    cout << "unable to enumerate processes - disabling parallel features" << endl;
    return -1;
  }
  
  cProcesses = cbNeeded / sizeof(DWORD);
  
  for(i = 0; i < cProcesses; i++) {
    HANDLE hProcess = OpenProcess( PROCESS_QUERY_INFORMATION | 
				   PROCESS_VM_READ, FALSE, 
				   ProcessesIDs[i] );
    
    if(hProcess != NULL)
      GetModuleBaseName(hProcess, NULL, szProcessName, 
			sizeof(szProcessName)/sizeof(TCHAR));
    
    if(!wcscmp(szProcessName, TEXT("smpd.exe"))) {
      found = true;
      cout << "yes (PID " << ProcessesIDs[i] << ")" << endl;
    }
    
    CloseHandle(hProcess);
  }
  
  if(!found) {
    cout << "no - disabling parallel features" << endl;
    return -1;
  }
  
#endif
  
  return 0;
}
