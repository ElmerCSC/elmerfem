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
Program:    ELMER model file (emf) reader defs
Module:     front_CPPdefs.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   Header file for the C++ interface for reading emf-files.

************************************************************************/

#ifndef _FRONT_DEFS_
#define _FRONT_DEFS_

extern "C" {
#include "front_Cdefs.h"
}



// =========================
// LibFront namespace starts
//==========================

namespace LibFront {

// Data value types 
// ================
enum dataValueType {
  EMF_INTEGER, EMF_REAL, EMF_LOGICAL, EMF_STRING, EMF_FILE
};


// Get numeric value data array
extern void setNumericData(bool*& target);
extern void setNumericData(int*& target);
extern void setNumericData(long*& target);
extern void setNumericData(float*& target);
extern void setNumericData(double*& target);
extern void setStringData(char*& buffer);

// Get numeric value from data values array by index
extern void setNumericData(bool& target, int value_index);
extern void setNumericData(int& target, int value_index);
extern void setNumericData(long& target, int value_index);
extern void setNumericData(float& target, int value_index);
extern void setNumericData(double& target, int value_index);
extern void setStringData(char* buffer, int value_index);
extern void setStringData(int nof_strings, char**& data_strings);

// Get numeric value from variable values array by index
extern void setNumericVariable(bool& target, int value_index);
extern void setNumericVariable(int& target, int value_index);
extern void setNumericVariable(long& target, int value_index);
extern void setNumericVariable(float& target, int value_index);
extern void setNumericVariable(double& target, int value_index);
extern void setStringVariable(char* buffer, int value_index);


// Misc utilities
// ==============
// A shorthand for isName(..)
bool in(const char* field_name, const char* test_name);

// Utilities for indented output
extern char* indent(short indent_size, short indent_level);
extern ostream& indent(ostream& out, short indent_size, short indent_level);

// String utilities
extern bool is_number(const char* str);
extern bool nocase_equal(char* str1, char* str2);
extern bool nocase_equal_partial(char* str1, char* str2);
extern void trim(char* buffer);
extern void trim_left(char* buffer);
extern void trim_right(char* buffer);



// Data output functions
// =====================
// Simple utility for debugging 
extern ostream&
output_data(ostream& out, char* data_type, void* data);


// Output field name properly indented
// -----------------------------------
extern ostream&
output_string(ostream& out, short indent_size, short indent_level,
              const char* str, bool output_eol = true);


// Output named scalar data
// ------------------------
extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const bool data);


extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const short data);


extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const int data);


extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const long data);


extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const float data);

extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const double data);


extern ostream&
output_scalar(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              const char* data, bool quoted = true);


// Output named vector data
// ------------------------
extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const bool* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const short* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const int* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const long* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const float* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const double* data, bool output_size = true);

extern ostream&
output_vector(ostream& out, short indent_size, short indent_level,
              const char* field_name, const char* field_type,
              int dim, const char** data, bool output_size = true, bool quoted = true);


// Output named table data (dim1*dim2)
// -----------------------------------
extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const bool** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const short** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const int** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const long** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const float** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const double** data);

extern ostream&
output_table(ostream& out, short indent_size, short indent_level,
             const char* field_name, const char* field_type,
             int dim1, int dim2, const char*** data, bool quoted = true);


}
//=====================
//End namespace LibFront


#endif
