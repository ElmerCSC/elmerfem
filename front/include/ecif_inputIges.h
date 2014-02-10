/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_inputIges.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   A base class for input of IGES-format CAD input files. 

************************************************************************/

#ifndef _ECIF_INPUT_IGES 
#define _ECIF_IGES_IGES


#include "ecif_input.h" 


struct IgesStatusField {
public:
  IgesStatusField();
  void setValues(char* data_str);
  int blankFlag;
  int subordFlag;
  int useFlag;
  int hrchFlag;
};


struct IgesDirectoryEntry {
public:
  IgesDirectoryEntry();
  int id;                 // Fld-10
  int entNbr;             // Fld-01
  int startLine;          // Fld-02
  int nofLines;           // Fld-14
  colorIndices colorNbr;  // Fld-13
  int formNbr;            // Fld-15
  int transfId;           // Fld-07
  IgesStatusField status; // Fld-09
  bool canBeBody;         // Flag: can this entry define a body (liked a closed circular arc in 2D)
  int referingId;         // The directory entry id which is using this entry
  Body* body;             // The body this entry defines or belongs to
};


class BodyElelment;
class Model;
class GcPoint;


//*****
class InputIges : public Input
{
public:
  InputIges(enum ecif_modelDimension m_dim,
            ifstream& in_file, char* in_filename);
  ~InputIges() {};
protected:
  bool createNewBody;
  IgesDirectory* directory;
  int paramSecLine;
  int paramSecStart;
  int addToDirectory(IgesDirectoryEntry* dir_entry);
  bool checkEntryReferences();
  bool createBodies();
  bool dataLineStrmIsEmpty(istrstream& data_line);
  enum ecif_modelDimension findCadModelDimension();
  void getDataField(istrstream*& data_line, char* field_buffer);
  IgesDirectoryEntry* getDirectoryEntry(int entry_id);
  void getDirectoryField(int fld_nbr, char* line_buffer, char* fld_buffer);
  void locateParamEntry(IgesDirectoryEntry* de);
  bool read_100(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_102(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_106(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_110(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_126(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_128(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_141(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_142(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_143(IgesDirectoryEntry* de, bool check_only_status = false);
  bool read_144(IgesDirectoryEntry* de, bool check_only_status = false);
  bool readCadGeometry(); 
  bool readCadHeader();
  bool readBodyElements(); 
  int readDirectory();
  IgesDirectoryEntry* readDirectoryEntry(char* first_line);
  void readDataLine(char* line_buf, char* data_buf);
  int readDoubleFields(istrstream*& data_line, int nof_fields, double* buffer);
  int readIntFields(istrstream*& data_line, int nof_fields, int* buffer);
  bool readLine(Body* body); 
  bool readNurbs(Body* body); 
  bool readPoint(istrstream*& data_line, Point3& p);
  BodyElement* readVertex(istrstream*& data_line);
} ; 



#endif
