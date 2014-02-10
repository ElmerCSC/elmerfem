/*  
   ElmerFront - A graphical user interface of Elmer software
   Copyright (C) 1995- , CSC - IT Center for Science Ltd.

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
Program:    ELMER model file (emf) reader
Module:     frontlib..cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   Implementation: Parser for ELMER model file (.emf files).

************************************************************************/

using namespace std;
#include <ctype.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>

#include "frontlib.h"


// Elmer MATC do-math call
extern "C" char* mtc_domath(char* s);


// *************************************************************************
// *************************************************************************
//
// START: Library globals
//
// *************************************************************************
// *************************************************************************
//
// Global interface variables
//
EMF_DECL int emf_ReadMode = 1;
EMF_DECL struct emf_ObjectData_X*  emf_ObjectData = 0;
EMF_DECL char* emf_UnknownData = new char[1 + emf_MAX_STRING_LEN];


// CallBack functions
//
extern "C" int (*emf_readDataCB) (void** user_data) = 0;
extern "C" int (*emf_readMessageCB) (char* msg_buffer) = 0;


// Global function definitions

// This function creates the reader and starts it
// NOTE: It is called by the user!
//
EMF_DECL int emf_readData(char* in_file_name, void** user_data,
                          int msg_buffer_size, char* msg_buffer,
                          int evaluate_only_matc)
{
  bool eval_only_matc = ( 0 != evaluate_only_matc);

  LibFrontReader* reader = new LibFrontReader(in_file_name, user_data,
                                    msg_buffer_size, msg_buffer,
                                    eval_only_matc);

  int rc = 0;  // Ok

  if ( eval_only_matc ) return rc;

  if ( reader == NULL || !reader->isOk() || !reader->start() ) {
    rc = 1;  // Error
  }

  delete reader;

  return rc;
}


// This function is called from the reader.
// It calls the CallBack function given by the user
// and thus informs about new available data.
//
int sendData(void** user_data)
{
  return emf_readDataCB(user_data);
}
// =========================================================================
// END: Library globals
// =========================================================================




// *************************************************************************
// *************************************************************************
// START Library locals
// *************************************************************************
// *************************************************************************

// Data handling functions
// -----------------------

//-Data arries
template <class T> void
setNumericData_impl(T*& target)
{
  if (LibFront::isLogicalData)
    target = (T*) LibFront::logicalData;

  else if (LibFront::isIntegerData)
    target = (T*) LibFront::integerData;

  else if (LibFront::isRealData)
    target = (T*) LibFront::realData;

}

//-Data values by index
template <class T> void
setNumericData_impl(T& target, int value_index)
{
  if (LibFront::isLogicalData)
    target = (T) LibFront::logicalData[value_index];

  else if (LibFront::isIntegerData)
    target = (T) LibFront::integerData[value_index];

  else if (LibFront::isRealData)
    target = (T) LibFront::realData[value_index];

}

//-Variable values
template <class T> void
setNumericVariable_impl(T& target, int value_index)
{
  if (LibFront::isRealVariable)
    target = (T) LibFront::realVariable[value_index];

}

// Data output functions
// ---------------------

// Scalar data
//
template <class T> ostream&
output_scalar_impl(ostream& out, short indent_size, short indent_level,
                   const char* field_name, const char* field_type,
                   const T data, bool quoted)
{
  char QM = '\"';

  // Field name
  if ( field_name != NULL ) {
    LibFront::output_string(out, indent_size, indent_level, field_name, false);
  }


  // Field type
  if ( field_type != NULL ) {
    out << endl;
    LibFront::output_string(out, indent_size, 1 + indent_level, field_type, false);
  }

  // If no name/type was output
  if ( field_name == NULL && field_type == NULL ) {
    LibFront::indent(out, indent_size, indent_level);

  // Data right after name/type string
  } else {
    out << " ";
  }

  //---Possible starting quote
  if (quoted)
    out << QM;

  // Null strings cannot be printed!
  if ( !(LibFront::ncEqual((char*)field_type, "string")) ||
       data != 0
     ) {

    //---Logical value as True/False
    if (LibFront::ncEqual((char*)field_type, "logical")) {
      if (data)
        out << "True";
      else
        out << "False";
    }
    //---Other data types
    else {
      out << data;
    }
  }

  //---Possible ending quote
  if (quoted)
    out << QM;

  out << endl;

  return out;
}


// Vector data
//
template <class T> ostream&
output_vector_impl(ostream& out, short indent_size, short indent_level,
                   const char* field_name, const char* field_type,
                   int dim, const T* data, bool output_size, bool quoted)
{
  char QM = '\"';

  // Field name
  if ( field_name != NULL ) {
    LibFront::output_string(out, indent_size, indent_level, field_name, false);
  }

  // Field size
  if (output_size) {
    out << endl;
    LibFront::indent(out, indent_size, 1 + indent_level);
    out << "Size " << dim << endl;
  }

  // Field type
  if (field_type != NULL ) {

    // If Size is not printed, we need a newline
    if (!output_size)
      out << endl;

    LibFront::output_string(out, indent_size, 2 + indent_level, field_type, false);
  }


  if (data == NULL) {
    return out;
  }

  // If no name/type was output
  if ( field_name == NULL && field_type == NULL ) {
    LibFront::indent(out, indent_size, indent_level);

  // Data right after name/type string
  } else {
    out << " ";
  }

  // Field data after type
  for (short i = 0; i < dim; i++) {

    //---Possible starting quote
    if (quoted)
      out << QM;

    // Null strings cannot be printed!
    if ( !(LibFront::ncEqual((char*)field_type, "string")) ||
         data[i] != 0
       ) {

      //---Logical value as True/Flase
      if (LibFront::ncEqual((char*)field_type, "logical")) {
        if ( data[i] )
          out << "True";
        else
          out << "False";
      }

      //---Other data types
      else  {
         out << data[i];
      }
    }

    //---Possible ending quote
    if (quoted)
      out << QM;

    out <<  ' ';
  }

  out << endl;

  return out;
}


// Table data
template <class T> ostream&
output_table_impl(ostream& out, short indent_size, short indent_level,
                  const char* field_name, const char* field_type,
                  int dim1, int dim2, const T** data, bool quoted)
{
  char QM = '\"';

  // Field name
  LibFront::output_string(out, indent_size, indent_level, field_name, true);

  // Field size
  LibFront::indent(out, indent_size, 1 + indent_level);
  out << "Size " << dim1 << " " << dim2 << endl;

  // Field type
  if (field_type != NULL )
    LibFront::output_string(out, indent_size, 2 + indent_level, field_type, true);

  if (data == NULL) {
    return out;
  }

  // Data lines after type line
  for (int i = 0; i < dim1; i++) {

    LibFront::indent(out, indent_size, 3 + indent_level);

    for (short j = 0; j < dim2; j++) {
      out << " ";

      // Possible starting quote
      if (quoted)
        out << QM;

      // Null strings cannot be printed!
      if ( !(LibFront::ncEqual((char*)field_type, "string")) ||
           data[i][j] != 0

         ) {

        //---Logical value as True/Flase
        if (LibFront::ncEqual((char*)field_type, "logical")) {
          if ( data[i][j] )
            out << "True";
          else
            out << "False";
        }
        //---Other data types
        else {
          out << data[i][j];
        }

        // Possible ending quote
        if (quoted)
          out << QM;

      } // If data_ij != NULL

    } // For dim2

    out << endl;

  } // For dim1

  return out;
}

ostream& output_data(ostream& out, char* data_type, void* data)
{
  if (data == NULL) {
    out << "NULL" << endl;
    return out;
  }

  // Associate to the data according to the data type
  if ( LibFront::ncEqual(data_type, "logical") ) {
    out << (int*)data;
  }
  else if ( LibFront::ncEqual(data_type, "integer") ) {
    out << (int*)data;
  }
  else if ( LibFront::ncEqual(data_type, "real") ) {
    out << (double*)data;
  }
  else if ( LibFront::ncEqual(data_type, "string") ) {
    out << (char*)data;
  }
  else if ( LibFront::ncEqual(data_type, "file") ) {
    out << (char*)data;
  }
  else {
    out << "Unknown datatype: " << data_type;
  }

  out << endl;

  return out;
}
// =========================================================================
// END: Library locals
// =========================================================================



// *************************************************************************
// *************************************************************************
//
// START: LibFront namespace implementation
//
// *************************************************************************
// *************************************************************************
//
// Data values related flags set in the library by the Reader and
// LibFront global data used to transfer data
//
bool LibFront::isFileData = false;
bool LibFront::isIntegerData = false;
bool LibFront::isLogicalData = false;
bool LibFront::isRealData = false;
bool LibFront::isStringData = false;
bool LibFront::isVoidData = false;
int* LibFront::integerData = NULL;
int* LibFront::logicalData = NULL;
double* LibFront::realData = NULL;
char* LibFront::stringData = NULL;
void* LibFront::voidData = NULL;
int LibFront::dimension1 = NULL;
int LibFront::dimension2 = NULL;
int* LibFront::dataLengths = NULL;

// Variables related flags set in the library by the Reader and
// LibFront global data used to transfer variables
//
bool LibFront::isRealVariable = false;
bool LibFront::isStringVariable = false;
bool LibFront::isVoidVariable = false;
double* LibFront::realVariable = NULL;
char* LibFront::stringVariable = NULL;
void* LibFront::voidVariable = NULL;
int* LibFront::variableValuesLengths = NULL;

char LibFront::stringSeparator = '\"';
char LibFront::matcSeparator = '$';


void
LibFront::setStringSeparator(char sep)
{
  stringSeparator = sep;
}

void
LibFront::setMatcSeparator(char sep)
{
  matcSeparator = sep;
}

// Read file into given buffer
//
// NOTE: Buffer is allocated in the function!
// Skip empty line, comment lines and eol comments
// Return: nof chars read or -1 if file not found
//
int
LibFront::readFile(char* filename, char*& buffer, bool skip_comments)
{
  // Initial size
  int bufferLen = 10000;
  int bufferPos = 0;

  buffer = new char[bufferLen];

  ifstream inputFile;

  inputFile.open(filename, ios::in );

  if ( inputFile.fail() ) return -1;

  bool is_string = false;
  bool skipping_comment = false;
  bool connecting_line = false;
  bool is_matc = false;
  int matc_start = -1;

	int fcount = 0;
  char c[1];


  // Read file until end-of-file
  //
  while ( !inputFile.eof() ) {

    c[0] = '\0';

    inputFile.read(c, 1);

    if ( c[0] == '\0' ) break;

    // Continuation marker (slash) irrelevant in a string and in comments
    if ( c[0] == '/' && !is_string && !skipping_comment) {
      connecting_line = true;
    }

    // String marker irrelevant in comments and after continuation marker
    if ( c[0] == stringSeparator && !(skipping_comment || connecting_line) ) {
      is_string = !is_string;
    }

    // Comment symbols irrelevant in strings
    if ( (c[0] == '!' || c[0] == '#') && skip_comments && !is_string) {
      skipping_comment = true;
    }

    // What newline means
    if ( c[0] == '\n' && skip_comments) {
      skipping_comment = false;

      // If we are connecting strings or we have an open string,
      // do not store the newline character
      if (connecting_line || is_string) {
        connecting_line = false;
        continue;
      }
    }

    // Skip carrage returns
    // NOTE: This way wwe can read dos-files (CR-LF) also in Unix
    // NOTE: Keep this after newline-handling, so that Dos files
    // are not ruined under Win32!
    if ( c[0] == '\r' && skip_comments) {
      continue;
    }

    // Comments and characters after continuation marker
    // in the same line can be skipped!
    if ( skipping_comment || connecting_line ) {
      continue;
    }

    // Allocate larget buffer is needed
    if ( fcount >= bufferLen - 1 ) {
      char* tmp = new char[bufferLen + 5000];
      for (int i = 0; i < bufferLen; i++)
        tmp[i] = buffer[i];
      bufferLen += 5000;
      delete[] buffer;
      buffer = tmp;
    }

    // Copy storable character to the buffer
    buffer[fcount++] = c[0];

  } // While !eof


  inputFile.close();

  // Check that data ends with eol!
  if ( buffer[fcount-1] != '\n' )  {
    buffer[fcount++] = '\n';
  }

  // Actual buffer size
  buffer[fcount] = '\0';

  return fcount;
}


// Read a Matc-expression from the source, starting from source_pos
// Return source-position after reading the expression
// NOTE: Source is NOT moved to this position actually, but the client
// can use this return value if its wants to advandec the source!
//
int
LibFront::readMatcExpression(char* source, int source_pos, int source_end, char* result, int max_result_len, bool skip_comments)
{
  int i,j;
  static char matc_peek_buffer[10001];

  bool is_string = false;
  bool skipping_comment = false;
  bool connecting_line = false;
  bool is_matc_end = false;
  bool is_matc_def = false;
  bool is_matc_func = false;

  bool check_matc_func = false;
  bool check_matc_end = false;

  int read_count = 0;
  int result_pos = 0;
  int matc_peek_pos = 0;
  int matc_lb_count = 0;
  int matc_rb_count = 0;
  int matc_lp_count = 0;
  int matc_rp_count = 0;

  char c_1 = '\0'; // Previous character
  char c = '\0';   // Current character

  // Read matc-expression in the result-buffer
  // advance the source-buffer to the end of the expression
  //
  while (true) {

    if ( source_pos >= source_end ) break;

    if ( result_pos >= max_result_len ) break;

    c = source[source_pos++];

    read_count++;

    // NOTE: We want to check if first non-blank text after $-sign
    // 'function'
    if ( matc_peek_pos < 10000 && c != ' ' ) {
      matc_peek_buffer[matc_peek_pos++] = c;
    }

    // Matc marker ($) is irrelevant in comments
    if ( c == matcSeparator && !skipping_comment) {
      is_matc_end = true;
      break;
    }

    if ( c == '\0' ) {
      check_matc_end = true;
    }

    // Continuation marker (slash) irrelevant in a string and in comments
    if ( c == '/' && !is_string && !skipping_comment) {
      connecting_line = true;
    }

    // Matc separators
    //
    // Left brace
    if ( c == '{' ) {
      matc_lb_count += 1;
      if ( matc_lb_count == 1 ) {
        check_matc_func = true;
      }

    // Right brace
    } else if ( c == '}' ) {
      matc_rb_count += 1;
      check_matc_end = true;

      // Left parenthesis
    } else if ( c == '(' ) {
      matc_lp_count += 1;
      if ( matc_lp_count == 1 ) {
        check_matc_func = true;
      }

    // Right parenthesis
    } else if ( c == ')' ) {
      matc_rp_count += 1;

    // Sign '=' is a definition, if lp==rp
    } else if ( c == '=' && matc_lp_count == matc_rp_count) {
      is_matc_def = true;

    // If not a definition, could be at end when lp==rp or when something
    // non-ws which is not =-sign
    //
    } else if ( !is_matc_def &&
                (c != '=' && c != ' ' && c != '\t') &&
                (c_1 == ' ' || c_1 == '\t' || matc_lp_count > 0) &&
                matc_lp_count == matc_rp_count) {
      check_matc_end = true;
    }

    // If we have the 'function'
    if ( check_matc_func ) {
      matc_peek_buffer[8] = '\0';
      if ( LibFront::in("function", matc_peek_buffer) ) {
        is_matc_func = true;
      }
    }

    // String marker (toggle)
    if ( c == stringSeparator && !(skipping_comment || connecting_line) ) {
      is_string = !is_string;
    }

    // Comment symbols irrelevant in strings
    if ( (c == '!' || c == '#') && skip_comments && !is_string) {
      skipping_comment = true;
    }

    // What newline means
    if ( c == '\n' ) {

      check_matc_end = true;
      skipping_comment = false;

      // If we are connecting strings or we have an open string,
      // do not store the newline character
      if (connecting_line || is_string) {
        connecting_line = false;
        continue;
      }
    }

    // Skip carrage returns
    // NOTE: This way we can read dos-files (CR-LF) also in Unix
    // NOTE: Keep this test after newline-handling, so that Dos files
    // are not ruined under Win32!
    if ( c == '\r' ) {
      continue;
    }

    // Comments and characters after continuation marker
    // in the same line can be skipped!
    if ( skipping_comment || connecting_line ) {
      continue;
    }

    // If the matc-expression is possibly at end
    if ( check_matc_end ) {

      // If we have the 'function'
      matc_peek_buffer[8] = '\0';
      if ( LibFront::in("function", matc_peek_buffer) ) {
        is_matc_func = true;
      }

      if ( !is_string &&
           ( matc_lb_count == matc_rb_count &&
             matc_lp_count == matc_rp_count &&
             ( (is_matc_func && matc_lb_count > 0) ||
               !is_matc_func
             )
           )
         ) {
        is_matc_end = true;
      }

      // If comments are not skipped, keep also the separating linefeeds!
      if ( !skip_comments && read_count == 1 && c == '\n' ) {
        result[result_pos++] = c;
        break;
      }
    }  // if check matc-end

    // When to store the charcter into result
    if ( !is_matc_end || c == ')' || c == '}' ) {
      result[result_pos++] = c;
    }

    if ( is_matc_end ) {

      // Move source backwards if we have been peeking too much to see
      // the end of the expression!!!
      //
      if ( c != ')' && c != '}' && c_1 == ' ' || c_1 == '\t' ) {
        while ( source[source_pos] != ' ' && source[source_pos] != '\t' ) {
          source_pos--;
        }
      }

      break;
    }

    // Store previous character
    c_1 = c;
  }

  result[result_pos] = '\0';

  // Return the source position just after the expression
  //
  return source_pos;
}


// Init Matc output format to 'plain' (ie. no Columns through) mode
// Precision set to 8
//
void LibFront::initMatcFrmt() {
  char matc_frmt[] = "format(8,\"rowform\")";
  mtc_domath( matc_frmt );
}


// Check that a Matc error message is Tcl-printable
void LibFront::formatMatcError(char* err_msg)
{
  // NOTE Tcl does not like [] !!!
  for(int i = 0; i < strlen(err_msg); i++) {
    if ( err_msg[i] == '[' ) err_msg[i] = ' ';
    if ( err_msg[i] == ']' ) err_msg[i] = ' ';
  }

}


// Check if the string is a matc function or variable definition
//
// NOTE: Whitespace must have been trimmed from 'str' before calling
// this function
//
bool
LibFront::isMatcDefinition(const char* str)
{
  static char buffer[9];
  int i;

  // Empty string or a comment or lf
  if ( str == NULL    ||
       str[0] == '\0' ||
       str[0] == '!'  ||
       str[0] == '#'  ||
       str[0] == '\n'
     ) return false;

  char* tmp = (char*)str;

  // Skip possible Matc-separator ($)
  if ( tmp[0] == matcSeparator ) tmp++;

  // Pick first 8 charcaters
  for (i = 0; i < 8; i++) {
    buffer[i] = tmp[i];
  }
  buffer[8] = '\0';
  toLower(buffer);

  // Function definition
  if ( 0 == strncmp(buffer, "function", 8) ) return true;

  // Variable definition
  for (i = 0; i < strlen(tmp); i++) {
    if ( tmp[i] == '=' ) return true;
  }

  return false;
}


// Checks if the string has matc variables
//
// NOTE: Whitespace must have been trimmed from 'str' before calling
// this function
//
bool
LibFront::hasMatcVariables(const char* str)
{
  if ( LibFront::isMatcDefinition(str) ) return false;

  int i;

  // Empty string or a comment or lf
  if ( str == NULL    ||
       str[0] == '\0' ||
       str[0] == '!'  ||
       str[0] == '#'  ||
       str[0] == '\n'
     ) return false;

  // Contains a variable name
  for (i = 0; i < strlen(str); i++) {
    if ( isalpha(str[i]) ) return true;
  }

  return false;
}


// Get a function or variable name from a Matc definition
//
// NOTE: Whitespace must have been trimmed from 'str' before calling
// this function
//
bool
LibFront::getMatcName(const char* str, char* name_buffer, int max_len, bool& is_var, bool& is_func)
{
  static char buffer[9];
  int i, j;

  is_var = false;
  is_func = false;
  name_buffer[0] = '\0';

  if ( !LibFront::isMatcDefinition(str) ) return false;

  char* tmp = (char*)str;

  // Skip possible Matc-separator ($)
  if ( tmp[0] == matcSeparator ) tmp++;

  // Pick first 8 characters to check if this
  // is a function definition
  for (i = 0; i < 8; i++) {
    buffer[i] = tmp[i];
  }
  buffer[8] = '\0';
  toLower(buffer);

  //--Function definition
  if ( 0 == strncmp(buffer, "function", 8) ) {
    tmp += 8;
    int start_pos = 0;
    LibFront::trimLeft(tmp, strlen(tmp), start_pos, true);
    tmp += start_pos;

    // Find function name end
    int len = 0;
    while (tmp[len] != '('  &&
           tmp[len] != ' '  &&
           tmp[len] != '\t' &&
           tmp[len] != '\n'
          ) len++;

    // Copy name
    for (j = 0; j < len && j < max_len; j++) {
      name_buffer[j] = tmp[j];
    }
    name_buffer[j] = '\0';
    is_func = true;
    return true;
  }

  //--Variable definition
  for (i = 0; i < strlen(tmp); i++) {
    if ( tmp[i] == '=' ) {

      // Find function name end
      int len = 0;
      while (tmp[len] != '='  &&
             tmp[len] != ' '  &&
             tmp[len] != '\t' &&
             tmp[len] != '\n'
            ) len++;
      // Copy name
      for (j = 0; j < len && j < max_len; j++) {
        name_buffer[j] = tmp[j];
      }
      name_buffer[j] = '\0';
      is_var = true;
      return true;
    }
  }

  return false;
}


// Check if result is a Matc error
//
// NOTE: Whitespace must have been trimmed from 'str' before calling
// this function
//
bool
LibFront::isMatcError(const char* str)
{
  static char buffer[11];
  int i;

  if ( str == NULL || str[0] == '\0' ) return false;

  int len = 10;

  // Pick first 10 charcaters
  for (i = 0; i < len; i++) {
    buffer[i] = str[i];
  }
  buffer[len] = '\0';
  toUpper(buffer);


  if ( 0 == strncmp("MATC ERROR", buffer, len) ) {
    return true;
  } else {
    return false;
  }
}


// Evaluate a string containing Matc-expressions
//
//
char*
LibFront::evalMatcString(const char* str)
{
  if ( str == NULL || str[0] == '\0' ) return NULL;

  int len = strlen(str);

  char* buffer = new char[1+len];

  strcpy(buffer, str);

  for (int i = 0; i < len; i++) {
    if ( buffer[i] == matcSeparator ) {
      buffer[i] = ' ';
    }
  }

  char* res = mtc_domath(buffer);

  delete[] buffer;

  if ( isMatcError(res) ) return NULL;

  return res;
}



//--Data reading functions

void
LibFront::setNumericData(bool*& target)
{
  setNumericData_impl(target);
}

void
LibFront::setNumericData(int*& target)
{
  setNumericData_impl(target);
}

void
LibFront::setNumericData(long*& target)
{
  setNumericData_impl(target);
}

void
LibFront::setNumericData(float*& target)
{
  setNumericData_impl(target);
}

void
LibFront::setNumericData(double*& target)
{
  setNumericData_impl(target);
}

void
LibFront::setStringData(char*& buffer)
{
  buffer = stringData;
}


void
LibFront::setNumericData(bool& target, int value_index)
{
  setNumericData_impl(target, value_index);
}

void
LibFront::setNumericData(int& target, int value_index)
{
  setNumericData_impl(target, value_index);
}

void
LibFront::setNumericData(long& target, int value_index)
{
  setNumericData_impl(target, value_index);
}

void
LibFront::setNumericData(float& target, int value_index)
{
  setNumericData_impl(target, value_index);
}

void
LibFront::setNumericData(double& target, int value_index)
{
  setNumericData_impl(target, value_index);
}


// NOTE: Buffer size must be large enough for the result!
//
void
LibFront::setStringData(char* buffer, int value_index)
{
  // Starting position in variable values string
  int start_pos = 0;
  for (int i = 0; i < value_index; i++)
    start_pos += dataLengths[i];

  // Copy string referred by the value_index
  int size = dataLengths[value_index];
  char* source = stringData + start_pos;
  strncpy(buffer, source, size);

  // End mark
  buffer[size] = '\0';
}


// NOTE: 'strings' argument is allocated here!
//
void
LibFront::setStringData(int nof_strings, char**& strings)
{
  strings = NULL;

  if ( nof_strings == 0 )
    return;

  strings = new char*[nof_strings];

  for (int i = 0; i < nof_strings; i++) {
    int len = dataLengths[i];
    strings[i] = new char[1 + len];
    setStringData(strings[i], i);
  }
}


void
LibFront::setNumericVariable(bool& target, int value_index)
{
  setNumericVariable_impl(target, value_index);
}

void
LibFront::setNumericVariable(int& target, int value_index)
{
  setNumericVariable_impl(target, value_index);
}

void
LibFront::setNumericVariable(long& target, int value_index)
{
  setNumericVariable_impl(target, value_index);
}

void
LibFront::setNumericVariable(float& target, int value_index)
{
  setNumericVariable_impl(target, value_index);
}

void
LibFront::setNumericVariable(double& target, int value_index)
{
  setNumericVariable_impl(target, value_index);
}

void
LibFront::setStringVariable(char* buffer, int value_index)
{
  // Starting position in varaible values string
  int start_pos = 0;
  for (int i = 0; i < value_index; i++)
    start_pos += variableValuesLengths[i];

  // Copy string referred by the value_index
  int size = variableValuesLengths[value_index];
  char* source = stringVariable + start_pos;
  strncpy(buffer, source, size);

  // End mark
  buffer[size] = '\0';
}



//--Field output functions

// Shorthand for "is_name" function
bool
LibFront::in(const char* field_name, const char* test_name)
{
  return LibFront::ncEqual((char*)field_name, (char*)test_name);
}


// Simple strin (field name etc.)
ostream&
LibFront::output_string(ostream& out, short indent_size, short indent_level,
                      const char* field_name, bool add_eol)
{
  LibFront::indent(out, indent_size, indent_level);

  out << field_name;

  if ( add_eol ) {
    out << endl;
  }

  return out;
}



// Callable cover functions
// ========================

ostream&
LibFront::output_matcDef(ostream& out, short indent_size, short indent_level,
                         const char* field_name, const char* field_type,
                         const char* def, bool add_dsign)
{
  strstream strm;

  if ( add_dsign ) {
    strm << '$';
  }
  strm << def << ends;


  return LibFront::output_scalar(out, indent_size, indent_level,
                                 field_name, field_type,
                                 strm.str(), false);
}


// Scalar data
ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const bool data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const short data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const int data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const long data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const float data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const double data)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, false);
}

ostream&
LibFront::output_scalar(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        const char* data, bool quoted)
{
  return output_scalar_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            data, quoted);
}


// Vector data
ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const bool* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}

ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const short* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}


ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const int* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}

ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const long* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}

ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const float* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}

ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const double* data, bool output_size)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, false);
}

ostream&
LibFront::output_vector(ostream& out, short indent_size, short indent_level,
                        const char* field_name, const char* field_type,
                        int dim, const char** data, bool output_size, bool quoted)
{
  return output_vector_impl(out, indent_size, indent_level,
                            field_name, field_type,
                            dim, data, output_size, quoted);
}


// Table data
ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const bool** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}


ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const short** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}


ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const int** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}


ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const long** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}

ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const float** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}

ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const double** data)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, data, false);
}

ostream&
LibFront::output_table(ostream& out, short indent_size, short indent_level,
                       const char* field_name, const char* field_type,
                       int dim1, int dim2, const char*** data, bool quoted)
{
  return output_table_impl(out, indent_size, indent_level,
                           field_name, field_type,
                           dim1, dim2, (const char**)data, quoted);
}


//****************
// Other utilities
//****************

// Check if string is a number using atof
bool
LibFront::isNumber(const char* str)
{
  if ( str == NULL || str[0] == '\0' )
    return false;

  int len = strlen(str);
  int index = 0;

  char* tmp = (char*)str;

  // Skip +/- sign
  if (tmp[0] == '-' || tmp[0] == '+' ) {
    tmp++;
    index++;
  }

  // Digits before decimal point are missing!
  if ( index == len || !isdigit(tmp[0]) )
    return false;

  // Skip digits before deciaml point
  while ( index < len && isdigit(tmp[0]) ) {
    tmp++;
    index++;
  }

  // OK: [+-]{digits}
  if ( index == len )
    return true;

  // Skip decimal point
  if ( tmp[0] == '.' ) {
    tmp++;
    index++;
  }


  // Skip digits after decimal point
  while ( index < len && isdigit(tmp[0]) ) {
    tmp++;
    index++;
  }

  // OK: [+-]{digits}.{digits}
  if ( index == len )
    return true;

  // Skip exponent indicator
  if ( tmp[0] == 'E' ||
       tmp[0] == 'e' ||
       tmp[0] == 'D' ||
       tmp[0] == 'd'
     ) {
    tmp++;
    index++;

  // Anything else is non-numeric!
  } else {
    return false;
  }

  // Digits after exponent are missing!
  if ( index == len )
    return false;

  // Skip +/- exponent's sign
  if (tmp[0] == '-' || tmp[0] == '+' ) {
    tmp++;
    index++;
  }

  // Digits after exponent's sign are missing!
  if ( index == len || !isdigit(tmp[0]) )
    return false;

  // Skip exponent's digits
  while ( index < len && isdigit(tmp[0]) ) {
    tmp++;
    index++;
  }

  // Now we should be at the end!

  // Not --> Error
  if ( index < len )
    return false;
  // Yes --> OK
  else
    return true;
}


// Compare two string (non-case sensitive)
// If string are same returns: true
// otherwise: false
bool
LibFront::ncEqual(char* str1, char* str2)
{
  if ( str1 == NULL || str2 == NULL ) {
    return false;
  }

  int i;
  int len1 = strlen(str1);
  int len2 = strlen(str2);

  if (len1 != len2)
    return false;

  for (i = 0; i < len1; i++)
    if ( toupper(str1[i]) != toupper(str2[i]) )
      return false;

  return true;
}


// Compare partially (according to the lengt of str2)
// two strings (non-case sensitive).
// If string are same returns: true
// otherwise: false
bool
LibFront::ncEqualPartial(char* str1, char* str2)
{
  int i;
  int len1 = strlen(str1);
  int len2 = strlen(str2);

  if (len1 < len2)
    return false;

  for (i = 0; i < len2; i++)
    if ( toupper(str1[i]) != toupper(str2[i]) )
      return false;

  return true;
}


// Indent stream by putting spaces into it
ostream&
LibFront::indent(ostream& out, short indent_size, short indent_level)
{
  for (int i = 0; i < indent_size * indent_level; i++)
    out << ' ';

  return out;
}

// Make a "indented" string by putting enough spaces into it
char*
LibFront::indent(short indent_size, short indent_level)
{
  const int MAX_INDENT = 80;
  static char buffer[1 + MAX_INDENT];

  int indent = indent_size * indent_level;

  if (indent > MAX_INDENT)
    indent = MAX_INDENT;

  for (int i = 0; i < indent; i++)
    buffer[i] = ' ';

  buffer[indent] = '\0';

  return buffer;
}


char*
LibFront::trim(char* buffer, bool trim_nl)
{
  LibFront::trimLeft(buffer, trim_nl);
  LibFront::trimRight(buffer, trim_nl);

  return buffer;
}


// Remove leading whitespace
char*
LibFront::trimLeft(char* buffer, bool trim_nl)
{
  int len = strlen(buffer);
  if ( len == 0)
    return buffer;

  int pos = 0;

  if ( trim_nl ) {
    while (buffer[pos] == '\t' || buffer[pos] == ' ' || buffer[pos] == '\n' ) {
      if (++pos == len)
        break;
    }

  } else {
    while (buffer[pos] == '\t' || buffer[pos] == ' ' ) {
      if (++pos == len)
        break;
    }
  }

  if (pos == 0)
    return buffer;

  // Start copying "back" from the first non-space position
  strcpy(buffer, buffer + pos);

  buffer[len-pos] = '\0';

  return buffer;
}


// Remove trailing whitespace
char*
LibFront::trimRight(char* buffer, bool trim_nl)
{
  int len = strlen(buffer);

  if ( len == 0)
    return buffer;

  // Put NULL after the last non-spce character in the buffer
  int pos = len - 1;

  if ( trim_nl ) {
    while (buffer[pos] == ' ' || buffer[pos] == '\t' || buffer[pos] == '\n' ) {
      buffer[pos] = '\0';
      if (--pos < 0)
        break;
    }

  } else {
    while (buffer[pos] == '\t' || buffer[pos] == ' ' ) {
      buffer[pos] = '\0';
      if (--pos < 0)
        break;
    }
  }

  return buffer;
}


// Find first non-ws position in the string
//
void
LibFront::trimLeft(const char* str, int max_pos, int& start_pos, bool trim_nl)
{
  // Pass by white-space at the end of the string
  while (1) {

    if (start_pos >= max_pos)
      break;

    char c = str[start_pos];

    if ( trim_nl ) {
      if (c != ' ' && c != '\t'  && c != '\n' )
        break;

    } else {
      if (c != ' ' && c != '\t' )
        break;
    }

    start_pos++;
  }
}


// Find last non-ws position in the string
//
void
LibFront::trimRight(const char* str, int min_pos, int& end_pos, bool trim_nl)
{
  // Pass by white-space at the end of the string
  while (1) {

    if (end_pos <= min_pos)
      break;

    char c = str[end_pos];

    if ( trim_nl ) {
      if (c != ' ' && c != '\t'  && c != '\n' )
        break;

    } else {
      if (c != ' ' && c != '\t' )
        break;
    }

    end_pos--;
  }
}


char*
LibFront::toLower(char* str)
{
  int i, len;
  len = strlen(str);
  for (i =0; i < len; i++) {
    str[i] = tolower(str[i]);
  }
  return str;
}


char*
LibFront::toLower(const char* source, char* target)
{
  int i, len;
  len = strlen(source);
  for (i =0; i <= len; i++) {
    target[i] = tolower(source[i]);
  }
  return target;
}


char*
LibFront::toUpper(char* str)
{
  int i, len;
  len = strlen(str);
  for (i =0; i < len; i++) {
    str[i] = toupper(str[i]);
  }
  return str;
}


char*
LibFront::toUpper(const char* source, char* target)
{
  int i, len;
  len = strlen(source);
  for (i =0; i <= len; i++) {
    target[i] = toupper(source[i]);
  }
  return target;
}


#if 0
// FILL GLOBAL FIELD NAME TABLE
// NOTE: This is strictly internal function
//
bool deleteFieldInfoTable()
{
  if (FieldInfoTable == NULL) {
    FieldInfoTableSize = 0;
    return true;
  }

  for (int i = 0; i < FieldInfoTableSize; i++)
    delete FieldInfoTable[i];

  delete[] FieldInfoTable;

  FieldInfoTableSize = 0;

  return true;
}


bool LibFront::createFieldInfoTable(int table_size)
{
  // Delete possible old table
  deleteFieldInfoTable();

  if (table_size <= 0)
    return false;

  // Allocate new table
  FieldInfoTable = new FieldInfo*[table_size];
  FieldInfoTableSize = table_size;

  for (int i = 0; i < table_size; i++)
    FieldInfoTable[i] = NULL;

  return true;
}
#endif
// =========================================================================
// END: LibFront namespace implementation
// =========================================================================



// *************************************************************************
// *************************************************************************
//
// START: LibFrontReaderFlags class implementation
//
// *************************************************************************
// *************************************************************************
//
LibFrontReaderFlags::LibFrontReaderFlags()
{
  init();
}


LibFrontReaderFlags::LibFrontReaderFlags(LibFrontReaderFlags& flags)
{
  copy(flags);
}


void
LibFrontReaderFlags::copy(LibFrontReaderFlags& flags)
{
  nonSizedDataAtEnd = flags.nonSizedDataAtEnd;
  canSkipSpace = flags.canSkipSpace;
  canSkipNewline = flags.canSkipNewline;
  endsField = flags.endsField;
  endsObject = flags.endsObject;
  equalSignPeeked = flags.equalSignPeeked;
  hasEvaluatedBuffer = flags.hasEvaluatedBuffer;
  matcPeeked = flags.matcPeeked;
  newlinePeeked = flags.newlinePeeked;
  hangingFieldNameNumber = flags.hangingFieldNameNumber;
  readingEndMarkedData = flags.readingEndMarkedData;
  readingFieldData = flags.readingFieldData;
  readingFileName = flags.readingFileName;
  readingGroupNames = flags.readingGroupNames;
  readingNonSizedData = flags.readingNonSizedData;
  readingProcedureName = flags.readingProcedureName;
  startsNewField = flags.startsNewField;
  startsNewObject = flags.startsNewObject;
}


void
LibFrontReaderFlags::init()
{
  nonSizedDataAtEnd = false;
  canSkipSpace = false;
  canSkipNewline = false;
  endsField = false;
  endsObject = false;
  equalSignPeeked = false;
  hasEvaluatedBuffer = false;
  matcPeeked = false;
  newlinePeeked = false;
  hangingFieldNameNumber = false;
  readingEndMarkedData = false;
  readingFieldData = false;
  readingFileName = false;
  readingGroupNames = false;
  readingNonSizedData = false;
  readingProcedureName = false;
  startsNewField = false;
  startsNewObject = false;
}
// =========================================================================
// END: LibFrontReaderFlags class implementation
// =========================================================================




// *************************************************************************
// *************************************************************************
// START: LibFrontReader class implementation
// *************************************************************************
// *************************************************************************
//
LibFrontReader::LibFrontReader(char* filename, void** user_data,
                     int msg_buffer_size, char* msg_buffer,
                     bool is_matc_file)
{
  readerOk = true;
  matcDefsWarned = false;

  inputFileName = new char[1 + strlen(filename)];
  strcpy(inputFileName, filename);

  // Attach stream to the file
  inputFile.open(inputFileName);
  inputFile.close();

  msgBufferSize = msg_buffer_size;
  msgBuffer = msg_buffer;

  stringSeparator = '\"';
  matcSeparator = '$';

  // Read file data into the readBuffer array (comments cleaned)
  readBufferPos = 0;
  readBufferLen = LibFront::readFile(inputFileName, readBuffer);

  // Read file data into the readBuffer array (comments cleaned)
  if ( readBufferLen <= 0 ) {
    readerOk = false;
  }

  if ( is_matc_file ) {
    if ( readerOk ) {
      evalMatcBuffer(readBuffer);
    }
    return;
  }

  dataBuffer = NULL;
  dataBufferSize = 0;
  dataLengths = NULL;

  dataAsStringBuffer[0] = '\0';
  dataAsStringPos = 0;

  dimension1 = 1;
  dimension2 = 1;

  evaluatedBufferEnd = -1;

  nameBuffer = new char[1 + emf_MAX_NAME_LEN];
  nameBufferMaxLen = emf_MAX_NAME_LEN;
  //nameBufferLen = 0;

  nofEntries = 0;
  nofVariables = 0;

  //tokenBuffer = new char[1 + emf_MAX_NAME_LEN];
  tokenBufferMaxLen = emf_MAX_NAME_LEN;
  tokenBufferLen = 0;

  userData = user_data;

  variableBuffer = NULL;
  variableBufferSize = 0;
  variableLengths = NULL;

  currentFlagsIndex = 0;
  flags = &flagsStack[currentFlagsIndex];
}


LibFrontReader::~LibFrontReader()
{
  inputFile.close();
  delete[] inputFileName;
  delete[] nameBuffer;
  delete[] readBuffer;

  delete_data_buffer();
  delete_variable_buffer();
}


bool
LibFrontReader::allocate_data_buffer(int current_size, int new_size)
{
  if (new_size == 0) {
    delete_data_buffer();
    return true;
  }

  void* new_buffer = NULL;

  int i;
  switch (dataBufferType) {
  case LibFront::EMF_INTEGER:
  case LibFront::EMF_LOGICAL:
    {
    int* new_data = new int[new_size];
    for (i = 0; i < current_size; i++)
        new_data[i] = ((int*)dataBuffer)[i];
    new_buffer = new_data;
    }
    break;
  case LibFront::EMF_REAL:
    {
    double* new_data = new double[new_size];
    for (i = 0; i < current_size; i++)
        new_data[i] = ((double*)dataBuffer)[i];
    new_buffer = new_data;
    }
    break;
  case LibFront::EMF_STRING:
  case LibFront::EMF_FILE:
    {
    char* new_data = new char[1 + new_size];
    new_data[new_size] = '\0';

    for (i = 0; i < current_size; i++)
        new_data[i] = ((char*)dataBuffer)[i];

    new_buffer = new_data;
    }
    break;
  }

  // Delete old buffer and set the new buffer
  delete_data_buffer();
  dataBufferSize = new_size;
  dataBuffer = new_buffer;

  return true;
}


// NOTE: Variable values are either Reals or Strings
bool
LibFrontReader::allocate_variable_buffer(int current_size, int new_size)
{
  if (new_size == 0) {
    delete_variable_buffer();
    return true;
  }

  void* new_buffer = NULL;

  int i;
  switch (variableBufferType) {
  case LibFront::EMF_REAL:
    {
    double* new_data = new double[new_size];
    for (i = 0; i < current_size; i++)
      new_data[i] = ((double*)variableBuffer)[i];
    new_buffer = new_data;
    }
    break;
  case LibFront::EMF_STRING:
    {
    char* new_data = new char[1 + new_size];
    new_data[new_size] = '\0';

    for (i = 0; i < current_size; i++)
      new_data[i] = ((char*)variableBuffer)[i];

    new_buffer = new_data;
    }
    break;
  }

  // Delete old buffer and set new
  delete_variable_buffer();
  variableBufferSize = new_size;
  variableBuffer = new_buffer;

  return true;
}


bool
LibFrontReader::appendReadBuffer(char c, int& bufferPos, bool init)
{
  int i,j;
  char* matc_result;
  static char matc_buffer[10001];

  static bool is_string = false;
  static bool skipping_comment = false;
  static bool connecting_line = false;
  static bool is_matc = false;
  static bool is_matc_def = false;
  static bool is_matc_func = false;
  static int matc_start = -1;
  static int matc_peek_pos = 0;
  static int matc_lb_count = 0;
  static int matc_rb_count = 0;
  static int matc_lp_count = 0;
  static int matc_rp_count = 0;

  if ( init ) {
    is_string = false;
    skipping_comment = false;
    connecting_line = false;
    is_matc = false;
    is_matc_func = false;
    is_matc_def = false;
    matc_start = -1;
    matc_peek_pos = 0;
    matc_lb_count = 0;
    matc_rb_count = 0;
    matc_lp_count = 0;
    matc_rp_count = 0;
  }

  bool check_matc_func = false;
  bool check_matc_end = false;
  bool do_matc = false;

  // NOTE: We want to check if first non-blank text after $-sign
  // 'function'
  if ( is_matc && matc_peek_pos < 10000 && c != ' ' ) {
    matc_buffer[matc_peek_pos++] = c;
  }

  // Matc marker ($) is irrelevant in a string and in comments
  //if ( c == '$' && !is_string && !skipping_comment) {
  if ( c == '$' && !skipping_comment) {

    if ( is_matc ) {
      do_matc = true;
    }
    is_matc = true;
    is_matc_def = false;
    is_matc_func = false;
    matc_start = bufferPos;
    matc_lb_count = 0;
    matc_rb_count = 0;
    matc_lp_count = 0;
    matc_rp_count = 0;
    matc_peek_pos = 0;
    return true;
  }

  if ( c == '\0' ) {
    check_matc_end = true;
    if ( is_matc ) {
      do_matc = true;
    }
  }

  // Continuation marker (slash) irrelevant in a string and in comments
  if ( c == '/' && !is_string && !skipping_comment) {
    connecting_line = true;
  }

  // Matc separators
  //
  if ( is_matc ) {

    // Left brace
    if ( c == '{' ) {
      matc_lb_count += 1;
      if ( matc_lb_count == 1 ) {
        check_matc_func = true;
      }

    // Right brace
    } else if ( c == '}' ) {
      matc_rb_count += 1;
      check_matc_end = true;

      // Left parenthesis
    } else if ( c == '(' ) {
      matc_lp_count += 1;
      if ( matc_lp_count == 1 ) {
        check_matc_func = true;
      }

    // Right parenthesis
    } else if ( c == ')' ) {
      matc_rp_count += 1;

    // Sign '=' is a definition, if lp==rp
    } else if ( c == '=' && matc_lp_count == matc_rp_count) {
      is_matc_def = true;

    // If a definition, could be at end when lp==rp
    } else if ( !is_matc_def &&
                (c != '=' && c != ' ' && c != '\t') &&
                matc_lp_count > 0 &&
                matc_lp_count == matc_rp_count) {
      check_matc_end = true;
    }

  }

  // If we have the 'function' after $-sign
  if ( check_matc_func ) {
    matc_buffer[8] = '\0';
    if ( LibFront::in("function", matc_buffer) ) {
      is_matc_func = true;
    }
  }

  // String marker irrelevant in comments and after continuation marker
  if ( c == stringSeparator && !(skipping_comment || connecting_line) ) {

    if ( is_string && is_matc ) {
      do_matc = true;
      is_string = false;

    } else {
      is_string = !is_string;
    }
  }

  // Comment symbols irrelevant in strings
  if ( (c == '!' || c == '#') && !is_string) {
    skipping_comment = true;
  }

  // What newline means
  if ( c == '\n' ) {
    check_matc_end = true;
    skipping_comment = false;

    // If we are connecting strings or we have an open string,
    // do not store the newline character
    if (connecting_line || is_string) {
      connecting_line = false;
      return true;
    }
  }

  // Skip carrage returns
  // NOTE: This way we can read dos-files (CR-LF) also in Unix
  // NOTE: Keep this test after newline-handling, so that Dos files
  // are not ruined under Win32!
  if ( c == '\r' ) {
    return true;
  }

  // Comments and characters after continuation marker
  // in the same line can be skipped!
  if ( skipping_comment || connecting_line ) {
    return true;
  }

  // If the matc-expression is possibly at end
  if ( is_matc && check_matc_end ) {

    // If we have the 'function' first after $-sign
    matc_buffer[8] = '\0';
    if ( LibFront::in("function", matc_buffer) ) {
      is_matc_func = true;
    }

    if ( !is_string &&
         ( matc_lb_count == matc_rb_count &&
           matc_lp_count == matc_rp_count &&
           ( (is_matc_func && matc_lb_count > 0) ||
             !is_matc_func
           )
         )
       ) {
      do_matc = true;
    }
  }

  // Allocate larger buffer is needed
  if ( bufferPos >= readBufferLen - 1 ) {
    char* tmp = new char[readBufferLen + 5000];
    for (int i = 0; i < readBufferLen; i++)
      tmp[i] = readBuffer[i];
    readBufferLen += 5000;
    delete[] readBuffer;
    readBuffer = tmp;
  }

  // Copy storable character to the buffer
  readBuffer[bufferPos++] = c;

  // Evaluate a Matc-expression
  // ------------------------
  if ( do_matc ) {

    // Copy all matc-stuff to the matc-buffer
    for (i = matc_start; i < bufferPos; i++) {
      matc_buffer[i-matc_start] = readBuffer[i];
    }

    matc_buffer[bufferPos - matc_start] = '\0';

    // EVALUATE Matc
    //
    matc_result = mtc_domath(matc_buffer);

    // If Matc-ERROR
    //
    if ( matc_result != NULL && 0 == strncmp(matc_result, "MATC ERROR", 10 ) ) {

      // NOTE Tcl does not like [] !!!
      for(j = 0; j < strlen(matc_result); j++) {
        if ( matc_result[j] == '[' ) matc_result[j] = ' ';
        if ( matc_result[j] == ']' ) matc_result[j] = ' ';
      }

      if ( matc_result != NULL ) {
        sendMessage(matc_result);
        sendMessage(matc_buffer);
      }

      is_matc = false;
      matc_start = -1;
      return false;
    }

    bufferPos = matc_start;

    is_matc = false;
    matc_start = -1;

    // If a Matc-expression result, insert result into read buffer
    // to be read 'normally'
    //
    for (i = 0; matc_result != NULL && i < strlen(matc_result); i++) {
      appendReadBuffer(matc_result[i], bufferPos, false);
    }

    return true;
  }

  return true;
}


bool
LibFrontReader::bufferAtEnd()
{
  if ( readBufferPos >= readBufferLen ||
       readBuffer[readBufferPos] == '\0'
     ) {
    return true;

  } else {
    return false;
  }
}




// NOTE: Equal sign (=) is also considered as white space!
//
void
LibFrontReader::bufferEatWs()
{
  // Skip spaces, tabs, newlines and equal sign
  //
  // NOTE: checkReadBuffer() handles also newlines
  while ( checkReadBuffer() &&
          ( readBuffer[readBufferPos] == ' '  ||
            readBuffer[readBufferPos] == '\t' ||
            readBuffer[readBufferPos] == '\n' ||
            readBuffer[readBufferPos] == '='
          )
          ) {

    readBufferPos++;
  }
}


// NOTE: Equal sign (=) is also considered as white space!
// NOTE: This version checks if an equal sign (token at end) or
// a newline was eaten (ie. data ended with a newline)
//
void
LibFrontReader::bufferEatWs(bool& eq_eaten, bool& nl_eaten)
{
  eq_eaten = false;
  nl_eaten = false;

  // Skip spaces, tabs, newlines and equal sign
  //
  // NOTE: checkReadBuffer() handles also newlines
  while ( checkReadBuffer() &&
          ( readBuffer[readBufferPos] == ' '  ||
            readBuffer[readBufferPos] == '\t' ||
            readBuffer[readBufferPos] == '\n' ||
            readBuffer[readBufferPos] == '='
          )
          ) {

    if ( readBuffer[readBufferPos] == '=' ) {
      eq_eaten = true;
    }

    if ( readBuffer[readBufferPos] == '\n' ) {
      nl_eaten = true;
    }

    readBufferPos++;
  }
}


bool
LibFrontReader::case_compare(char* str1, char* str2)
{
  int len1 = strlen(str1);
  int len2 = strlen(str2);
  if (len1 =! len2)
    return false;

  for (int i = 0; i < len1; i++) {
    if ( str1[i] != str2[i] )
      return false;
  }

  return true;
}


// Check if read buffer is at end
//
bool
LibFrontReader::checkReadBuffer()
{
  if ( bufferAtEnd() ) {
    return false;
  }

  bool ok = true;

  char c = readBuffer[readBufferPos];

  if ( c == '\n' ) {
    ok = handleEol();
  }

  if (bufferAtEnd()) {
    return false;
  }

  return ok;
}


void
LibFrontReader::constructStringFromTokens(char* buffer, int nof_string_tokens)
{
  char* tmp = buffer;

  for (int i = 0; i < nof_string_tokens; i++) {

    int old_len = strlen(tmp);
    getNextToken(tmp + old_len);
    int new_len = strlen(tmp);

    // Add space  separator between tokens
    if (i < nof_string_tokens - 1) {

      tmp[new_len] = ' ';
      tmp[++new_len] = '\0';
      tmp += new_len - old_len;
    }
  }
}


bool
LibFrontReader::data_at_end()
{
  bool at_end = false;

  //-these need two strings
  if ( flags->readingFileName || flags->readingProcedureName) {
    if (nofEntries == 2)
      at_end = true;

  //-if data end is marked with "End"
  } else if ( flags->readingEndMarkedData ) {
    //if ( next_is_text("end") )
    if ( next_is_end_keyword() )
      at_end = true;

  //-if there is not following a string starting with quote
  } else if ( !flags->readingEndMarkedData &&
              flags->readingNonSizedData
            ) {
    if ( !next_is_text("\"") )
      at_end = true;

  //-otherwise dimensions define data size
  // and it is always read to the end
  } else {
    at_end = true;
  }

  return at_end;
}


void
LibFrontReader::delete_data_buffer()
{
  if (dataBuffer == NULL)
    return;

  switch (dataBufferType) {
  case LibFront::EMF_INTEGER:
  case LibFront::EMF_LOGICAL:
    delete[] (int*)dataBuffer;
    break;
  case LibFront::EMF_REAL:
    delete[] (double*)dataBuffer;
    break;
  case LibFront::EMF_STRING:
  case LibFront::EMF_FILE:
    delete[] (char*)dataBuffer;
    break;
  }

  dataBuffer = NULL;
  dataBufferSize = 0;

}


void
LibFrontReader::delete_variable_buffer()
{
  if (variableBuffer == NULL)
    return;

  switch (dataBufferType) {
  case LibFront::EMF_REAL:
    delete[] (double*)variableBuffer;
    break;
  case LibFront::EMF_STRING:
    delete[] (char*)variableBuffer;
    break;
  }

  variableBuffer = NULL;
  variableBufferSize = 0;

}


// Evalute a source buffer containing possible Matc-expressions.
// Copy source buffer and all evaluated results into result-buffer
//
// NOTE: No $-sign accepted
//
bool
LibFrontReader::evalMatcBuffer(char* source)
{
  char matc_buffer[10001];

  int source_pos = 0;
  int source_len = strlen(source);

  // Read source buffer till end
  //
  while (source_pos < source_len ) {

    // Read expression and set new source-pos just after the expression
    //
    source_pos = LibFront::readMatcExpression(source, source_pos, source_len,
                                              matc_buffer, 10000);

    LibFront::trim(matc_buffer);

    if ( matc_buffer[0] == '\0' ) continue;

    // Evaluate the Matc-expression
    //
    char* matc_result = mtc_domath(matc_buffer);

    // An empty result
    if ( matc_result == NULL || matc_result[0] == '\0' ) continue;

    // Matc Error
    if ( LibFront::isMatcError(matc_result) ) {
      strstream strm, strm2;
      strm << "***WARNING: Matc evaluation error (" << matc_buffer << ")!" << ends;
      sendMessage(strm.str());
      strm2 << matc_result<< ends;
      sendMessage(strm2.str());
      continue;
    }
  }

  return true;
}


int
LibFrontReader::evalAndInsertNextMatc()
{
  readBufferPos = LibFront::readMatcExpression(readBuffer, readBufferPos, readBufferLen,
                                               matcBuffer, 10000);

  // WARNING: Matc definitions are not stored in the egf file, so a warning is in place!
  //
  if ( LibFront::isMatcDefinition(matcBuffer) && !matcDefsWarned ) {
    strstream strm;
    strm << "***WARNING: Matc definitions are not recommended in the egf-file! (" << matcBuffer << ")" << ends;

    if ( 1 == sendMessage(strm.str()) ) {
      matcDefsWarned = true;
    }
  }

  char* matc_result = mtc_domath(matcBuffer);

  if ( matc_result == NULL || matc_result[0] == '\0' ) return 0;

  // Matc error
  if ( LibFront::isMatcError(matc_result) ) {

    // NOTE Tcl does not like [] !!!
    for(int i = 0; i < strlen(matc_result); i++) {
      if ( matc_result[i] == '[' ) matc_result[i] = ' ';
      if ( matc_result[i] == ']' ) matc_result[i] = ' ';
    }

    if ( matc_result != NULL ) {
      sendMessage(matc_result);
      sendMessage(matcBuffer);
    }

    return -1;
  }

  LibFront::trim(matc_result);

  // Insert matc-expression
  int eof_insert = insertReadBuffer(matc_result);

  // If something was inserted, update eval-end-position
  if ( eof_insert > 0 ) {
    evaluatedBufferEnd = eof_insert;
  }

  storeDataAsString(matcBuffer, true);

  if ( LibFront::hasMatcVariables(matcBuffer) ) {
    dataHasMatcVars = true;
  }

  return strlen(matcBuffer);
}


void
LibFrontReader::evalAndSkipNextMatc()
{
  readBufferPos = LibFront::readMatcExpression(readBuffer, readBufferPos, readBufferLen,
                                               matcBuffer, 10000);

  // WARNING: Matc definitions are not stored nowhere if in the body of an egf file!
  if ( LibFront::isMatcDefinition(matcBuffer) && !matcDefsWarned ) {
    strstream strm;
    strm << "***WARNING: Matc definitions are not recommended in the egf-file! (" << matcBuffer << ")" << ends;

    if ( 1 == sendMessage(strm.str()) ) {
      matcDefsWarned = true;
    }
  }

  char* matc_result = mtc_domath(matcBuffer);
}


bool
LibFrontReader::extractSeparator(char sep)
{
  bufferEatWs();

  if (bufferAtEnd())
    return false;

  char c = readBuffer[readBufferPos];

  if ( c != sep ) {
    return false;
  }
  else {
    readBufferPos++;
    return true;
  }
}


bool
LibFrontReader::getNextNumber(int& number)
{
  return getNextNumber_impl(number);
}

bool
LibFrontReader::getNextNumber(double& number)
{
  return getNextNumber_impl(number);
}


template <class T> bool
LibFrontReader::getNextNumber_impl(T& number)
{
  if (bufferAtEnd()) {
    return false;
  }

  getNextToken(tokenBuffer);

  strstream strm;
  strm << tokenBuffer;

  return ( strm >> number );
}


// Basic version without any flag arguments
//
int
LibFrontReader::getNextToken(char* buffer)
{
  if (bufferAtEnd()) {
    return 0;
  }

  bool store_as_string = true;

  if ( !flags->readingFieldData ||
       ( flags->hasEvaluatedBuffer &&
         evaluatedBufferEnd >= readBufferPos
       )
       ) {
    store_as_string = false;
  }

  // A Matc-expression in the input
  //
  if ( readBuffer[readBufferPos] == matcSeparator ) {

    readBufferPos++;

    // Evaluate and insert the possible result
    //
    if ( flags->readingFieldData ) {

      int ec = evalAndInsertNextMatc();

      // Something as a result
      if ( ec > 0 ) {
        flags->hasEvaluatedBuffer = true;

      // Error, not a definition, but however no result value!
      } else if ( ec < 0 ) {
        return -1;
      }

    // Evaluate, but do not insert the possible result
    //
    } else {
      evalAndSkipNextMatc();
    }

    bufferEatWs();
    return getNextToken(buffer);
  }

  if ( readBuffer[readBufferPos] == stringSeparator ) {
    readBufferPos++;
    return getToSeparator(buffer, stringSeparator );
  }

  flags->nonSizedDataAtEnd = false;

  buffer[0] = '\0';

  int counter = 0;

  // Find next token separator
  //
  while ( checkReadBuffer() ) {

    char c = readBuffer[readBufferPos];

    // These stop reading (ws + nl + equal sign + $)
    if ( c == ' '  ||
         c == '\t' ||
         c == '\n' ||
         c == '='  ||
         c == matcSeparator
       ) {
      break;
    }

    readBufferPos++;
    buffer[counter++] = c;
  }

  // Clean ws from readBuffer
  bufferEatWs();

  buffer[counter] = '\0';

  if ( store_as_string ) {
    storeDataAsString(buffer);
  }

  return strlen(buffer);
}



// NOTE: Checks also if a newline or an equal sign was
// after the token
//
// NOTE: This version inserts only non-data related field items (names etc.) into
// readBuffer, field data related expression are inserted by the 'plain' getNextToken-function
//
// It also evaluates 'loose' Matc-expressions between fields!
//
// It is meant mainly for peeking tokens in the input!
//
int
LibFrontReader::getNextToken(char* buffer, bool& eq_found, bool& nl_found, bool& matc_found)
{
  static char matc_buffer[10001];

  if (bufferAtEnd()) {
    eq_found = false;
    nl_found = true;
    matc_found = false;
    return 0;
  }

  if ( readBuffer[readBufferPos] == matcSeparator ) {
    readBufferPos++;
    readBufferPos = LibFront::readMatcExpression(readBuffer, readBufferPos, readBufferLen, matc_buffer, 10000);

    char* matc_result = mtc_domath(matc_buffer);

    if ( matc_result != NULL && matc_result[0] != '\0' ) {
      strstream strm;
      strm << matc_result << ends;
      strm >> buffer;
      matc_found = true;

      return strlen(buffer);

    } else {
      bufferEatWs();
      return getNextToken(buffer, eq_found, nl_found, matc_found);
    }
  }

  if ( readBuffer[readBufferPos] == stringSeparator ) {
    readBufferPos++;
    return getToSeparator(buffer, stringSeparator );
  }

  flags->nonSizedDataAtEnd = false;

  buffer[0] = '\0';

  int counter = 0;

  // Find next token separator
  //
  while ( checkReadBuffer() ) {

    char c = readBuffer[readBufferPos];

    // These stop reading (ws + nl + equal sign + $)
    if ( c == ' '  ||
         c == '\t' ||
         c == '\n' ||
         c == '='  ||
         c == matcSeparator
       ) {
      matc_found = ( c == matcSeparator );
      break;
    }

    readBufferPos++;
    buffer[counter++] = c;
  }

  // Clean ws from readBuffer
  bufferEatWs(eq_found, nl_found);

  buffer[counter] = '\0';

  return strlen(buffer);
}


int
LibFrontReader::getToSeparator(char* buffer, char sep)
{
  buffer[0] = '\0';

  if (bufferAtEnd())
    return 0;

  bool store_as_string = true;

  if ( !flags->readingFieldData ||
       ( flags->hasEvaluatedBuffer &&
         evaluatedBufferEnd >= readBufferPos
       )
       ) {
    store_as_string = false;
  }

  flags->nonSizedDataAtEnd = false;

  int counter = 0;

  while ( checkReadBuffer() ) {

    char c = readBuffer[readBufferPos++];

    if ( c == sep )
      break;

    buffer[counter++] = c;
  }

  bufferEatWs();

  buffer[counter] = '\0';

  if ( store_as_string ) {
    storeDataAsString(buffer);
  }

  return strlen(buffer);
}


bool
LibFrontReader::handleEol()
{
  flags->nonSizedDataAtEnd = true;

  //readBufferPos++;

  return true;
}


void
LibFrontReader::init_transfer_info()
{
  LibFront::isLogicalData = false;
  LibFront::isIntegerData = false;
  LibFront::isRealData = false;
  LibFront::isStringData = false;
  LibFront::isFileData = false;
  LibFront::isVoidData = false;

  LibFront::logicalData = NULL;
  LibFront::integerData = NULL;
  LibFront::realData = NULL;
  LibFront::stringData = NULL;
  LibFront::voidData = NULL;

  LibFront::dimension1 = 0;
  LibFront::dimension2 = 0;
  LibFront::dataLengths = NULL;

  LibFront::isRealVariable = false;
  LibFront::isStringVariable = false;
  LibFront::isVoidVariable = false;

  LibFront::realVariable = NULL;
  LibFront::stringVariable = NULL;

  LibFront::variableValuesLengths = NULL;

}

// Insert new data into readBuffer
// This is typically used to insert an evaluated Matc-epxerssion
// into input buffer
//
// Return: end-of-insert position or -1 if nothing was inserted
//
int
LibFrontReader::insertReadBuffer(char* data)
{
  int eof_insert = -1;

  if ( data == NULL || data[0] == '\0' ) return eof_insert;

  int data_len = strlen(data);

  eof_insert = readBufferPos + data_len;

  // Allocate more space for the new read buffer
  char* tmp = new char[readBufferLen + data_len];

  int i;
  //-Copy old leading data
  for (i = 0; i < readBufferPos; i++) {
    tmp[i] = readBuffer[i];
  }

  //-Insert new data
  for (i = 0; i < data_len; i++) {
    tmp[readBufferPos + i] = data[i];
  }

  //-Copy old trailing data
  for (i = readBufferPos; i < readBufferLen; i++) {
    tmp[data_len + i] = readBuffer[i];
  }

  // Update readBuffer
  delete[] readBuffer;
  readBuffer = tmp;

  readBufferLen += data_len;

  return eof_insert;
}


bool
LibFrontReader::is_end_keyword(char* str)
{
  LibFront::toLower(str);
  if ( 0 == strcmp(str, "end") )
    return true;
  else
    return false;
}


bool
LibFrontReader::is_file_keyword(char* str)
{
  if ( 0 == strcmp(str, "file") )
    return true;
  else
    return false;
}


bool
LibFrontReader::is_keyword(char* str)
{
  LibFront::toLower(str);

  if ( is_size_keyword(str)     ||
       is_type_keyword(str)     ||
       is_variable_keyword(str)
     )
    return true;
  else
    return false;
}



bool
LibFrontReader::is_number(char* str)
{
  if ( str == NULL || str[0] == '\0' ) return false;

  return LibFront::isNumber(str);
}


bool
LibFrontReader::is_procedure_keyword(char* str)
{
  if ( 0 == strcmp(str, "procedure") )
    return true;
  else
    return false;
}


bool
LibFrontReader::is_separated_string()
{
  char c = inputFile.peek();
  if ( c == '\"' ) return true;
  return false;
}


bool
LibFrontReader::is_size_keyword(char* str)
{
  LibFront::toLower(str);
  if ( 0 == strcmp(str, "size") )
    return true;
  else
    return false;
}


bool
LibFrontReader::is_type_keyword(char* str)
{
  LibFront::toLower(str);
  if ( 0 == strcmp(str, "integer")  ||
       0 == strcmp(str, "real")     ||
       0 == strcmp(str, "logical")  ||
       0 == strcmp(str, "string")   ||
       0 == strcmp(str, "file")     ||
       0 == strcmp(str, "procedure")
     )
    return true;
  else
    return false;
}


bool
LibFrontReader::is_variable_keyword(char* str)
{
  LibFront::toLower(str);
  if ( 0 == strcmp(str, "variable") )
    return true;
  else
    return false;
}


// NOTE: End keyword must be the last token in the line!!!
//
bool
LibFrontReader::next_is_end_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(peek_buffer) ||
       !(flags->newlinePeeked && !flags->equalSignPeeked)
       ) {
    return false;
  }

  LibFront::toLower(peek_buffer);

  if ( is_end_keyword(peek_buffer) )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_group_names()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(peek_buffer) ) {
    return false;
  }

  if ( peek_buffer[0] == '[' ) {
    return true;
  } else {
    return false;
  }
}


bool
LibFrontReader::next_is_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( next_is_end_keyword() ) return true;

  if ( !peek_token(peek_buffer) ) {
    return false;
  }

  LibFront::toLower(peek_buffer);

  if ( is_keyword(peek_buffer) ) {
    return true;
  } else {
    return false;
  }
}


bool
LibFrontReader::next_is_matc()
{
  if (bufferAtEnd()) return false;

  push_flags();
  int old_pos = readBufferPos;

  bufferEatWs();

  char c = readBuffer[readBufferPos];

  readBufferPos = old_pos;
  pop_flags();

  if ( c == matcSeparator )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_number()
{
  static char buffer[emf_MAX_STRING_LEN];
  if ( !peek_token(buffer) ) return false;
  return is_number(buffer);
}


bool
LibFrontReader::next_is_separated_string()
{
  if (bufferAtEnd()) return false;

  push_flags();
  int old_pos = readBufferPos;

  bufferEatWs();

  char c = readBuffer[readBufferPos];

  readBufferPos = old_pos;
  pop_flags();

  if ( c == stringSeparator )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_size_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(peek_buffer) )
    return false;

  LibFront::toLower(peek_buffer);

  if ( is_size_keyword(peek_buffer) )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_text(char* text)
{
  static char buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(buffer) )
    return false;

  if ( LibFront::ncEqual(buffer, text) )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_type_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(peek_buffer) )
    return false;

  LibFront::toLower(peek_buffer);

  if ( is_type_keyword(peek_buffer) )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_is_variable_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_token(peek_buffer) )
    return false;

  LibFront::toLower(peek_buffer);

  if ( is_variable_keyword(peek_buffer) )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_non_number_is_data_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( !peek_next_string(peek_buffer) )
    return false;

  LibFront::toLower(peek_buffer);

  if ( is_size_keyword(peek_buffer)     ||
       is_type_keyword(peek_buffer)     ||
       is_variable_keyword(peek_buffer)
     )
    return true;
  else
    return false;
}


bool
LibFrontReader::next_non_number_is_keyword()
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  if ( next_is_end_keyword() ) return true;

  if ( !peek_next_string(peek_buffer) )
    return false;

  LibFront::toLower(peek_buffer);

  if ( is_keyword(peek_buffer) )
    return true;
  else
    return false;
}


bool
LibFrontReader::object_at_end()
{
  if (bufferAtEnd())
    return false;

  if ( next_is_end_keyword() ) {
    return true;
  } else {
    return false;
  }
}


// Method checks the number of string tokens before possible TWO
// consecutive numeric tokens.
// NOTE: the total number of numeric tokens is not counted, it is
// enough to know if there is 0, 1 or >= 2 of them!!!
// NOTE: Keywords are accepted as name tokens! ==> a LINEFEED must
//       be after the object/field name before a possible keyword!
bool
LibFrontReader::peek_field_name_tokens(int& nof_tokens)
{
  nof_tokens = 0;

  bool nl_found = false;
  bool eq_found = false;
  bool mc_found = false;

  push_flags();
  int old_pos = readBufferPos;

  flags->nonSizedDataAtEnd = false;

  int nof_picked_tokens = 0;
  int nof_name_tokens = 0;

  bool next_found = false;

  bool tmp = flags->readingFieldData;
  flags->readingFieldData = false;

  // Read all tokens in the line
  while( !bufferAtEnd() && !flags->nonSizedDataAtEnd ) {

    if ( next_is_separated_string() ) {
      flags->readingNonSizedData = true;
      break;
    }

    // Try to read next token
    next_found = ( 0 != getNextToken(tokenBuffer, eq_found, nl_found, mc_found));

    // If data at end or next starst with matc-$ sign, field name must end
    if ( !next_found || mc_found ) {
      break;
    }

    if ( next_found ) {
      nof_picked_tokens++;
    }

    // If next still was text, add to the field name
    //
    if ( next_found && !is_number(tokenBuffer) ) {
      nof_name_tokens = nof_picked_tokens;
    }

    // If an equal sign was found, field name must end
    if ( eq_found ) {
      break;
    }

  } // while

  flags->readingFieldData = tmp;

  // If next is keyword, all picked tokens belong to field name!
  if ( next_is_keyword() && !next_is_end_keyword() ) {
    nof_tokens = nof_picked_tokens;
  } else {
    nof_tokens = nof_name_tokens;
  }

  readBufferPos = old_pos;
  pop_flags();

  return true;
}


// Peeks first "real" string token in readStream, so all numbers are first skipped
bool
LibFrontReader::peek_next_string(char* buffer)
{
  buffer[0] = '\0';

  if (bufferAtEnd()) {
    return false;
  }

  push_flags();
  int old_pos = readBufferPos;

  bool tmp = flags->readingFieldData;
  flags->readingFieldData = false;

  // Skip numbers
  double test_number;
  while (!bufferAtEnd()) {

    getNextToken(tokenBuffer);

    strstream strm;
    strm << tokenBuffer;

    if ( !(strm >> test_number) ) {
      break;
    }
  }

  flags->readingFieldData = tmp;

  strcpy(buffer, tokenBuffer);

  readBufferPos = old_pos;
  pop_flags();

  if ( strlen(buffer) == 0 ) {
    return false;
  } else {
    return true;
  }
}


// Method checks the number of string tokens before possible TWO
// consecutive numeric tokens.
// NOTE: the total number of numeric tokens is not counted, it is
// enough to know if there is 0, 1 or >= 2 of them!!!
// NOTE: Keywords are accepted as name tokens! ==> a LINEFEED must
//       be after the object/field name before a possible keyword!
bool
LibFrontReader::peek_object_name_tokens(int& nof_tokens)
{
  nof_tokens = 0;

  push_flags();
  int old_pos = readBufferPos;

  flags->nonSizedDataAtEnd = false;

  int nof_picked_tokens = 0;
  int nof_consecutive_nbr_tokens = 0;

  bool tmp = flags->readingFieldData;
  flags->readingFieldData = false;

  while( !bufferAtEnd() && !flags->nonSizedDataAtEnd ) {

    // No separated strings accepted in object name
    if ( next_is_separated_string() ) {
      return false;
    }

    getNextToken(tokenBuffer);

    nof_picked_tokens++;

    // Include keyword is a bit exceptional!!!
    //
    if ( nof_picked_tokens == 1 && LibFront::in(tokenBuffer, "Include") ) {
      break;
    }

    if ( is_number(tokenBuffer) ) {
      nof_consecutive_nbr_tokens++;
    } else {
      nof_consecutive_nbr_tokens = 0;
    }

    // Only one number accepted after object name (object id)
    if ( nof_consecutive_nbr_tokens > 1 ) {
      return false;
    }

  } // while

  flags->readingFieldData = tmp;

  nof_tokens = nof_picked_tokens - nof_consecutive_nbr_tokens;

  readBufferPos = old_pos;
  pop_flags();

  return true;
}


// Peeks next token in the readBuffer
// NOTE: Checks also if an equal sign or a newline was ahed
//
bool
LibFrontReader::peek_token(char* buffer)
{
  bool eq_found = false;
  bool nl_found = false;
  bool mc_found = false;

  flags->equalSignPeeked = eq_found;
  flags->newlinePeeked = nl_found;
  flags->matcPeeked = mc_found;

  buffer[0] = '\0';

  if (bufferAtEnd()) {
    return false;
  }

  push_flags();
  int old_pos = readBufferPos;

  bool tmp = flags->readingFieldData;
  flags->readingFieldData = false;

  getNextToken(buffer, eq_found, nl_found, mc_found);

  flags->readingFieldData = tmp;

  readBufferPos = old_pos;
  pop_flags();

  flags->equalSignPeeked = eq_found;
  flags->newlinePeeked = nl_found;
  flags->matcPeeked = mc_found;

  if ( strlen(buffer) == 0 ) {
    return false;
  } else {
    return true;
  }
}


void
LibFrontReader::pop_flags()
{
  currentFlagsIndex--;

  flags = &flagsStack[currentFlagsIndex];
}


void
LibFrontReader::push_flags(bool copy_current)
{
  currentFlagsIndex++;

  if (copy_current)
    flagsStack[currentFlagsIndex].copy(flagsStack[currentFlagsIndex - 1]);
  else
    flagsStack[currentFlagsIndex].init();

  flags = &flagsStack[currentFlagsIndex];
}


// Read data in the entry
bool
LibFrontReader::read_entry_data(struct emf_ObjectData_X* od,
                           int data_alloc_count, int data_size,
                           int& read_count)
{
  // Buffer sizes
  int old_data_buffer_size, new_data_buffer_size;

  // Extra size allocated for this entry
  int data_alloc_size;

  data_alloc_size = data_alloc_count * data_size;
  old_data_buffer_size = dataBufferSize;
  new_data_buffer_size = dataBufferSize + data_alloc_size;

  allocate_data_buffer(old_data_buffer_size, new_data_buffer_size);

  currentBuffer = dataBuffer;
  currentBufferSize = dataBufferSize;
  currentBufferType = dataBufferType;

  int* read_sizes = new int[data_alloc_count];

  // Read values
  if ( !read_entry_data_values(old_data_buffer_size, new_data_buffer_size,
                               data_alloc_count, read_sizes, read_count)
     ) {

    // If no succes: ERROR MSG!
    sendMessage("\n");
    sendMessage("*** ERROR ***: When reading data for:");
    //-- Object name FieldName
    {
      strstream strm;
      strm << "OBJECT: " << od->object_name;

      if (od->object_id != -1) {
        strm << " " << od->object_id;
      }
      strm << "    FIELD:  " << od->field_name;
      strm << ends;
      sendMessage(strm.str());
    }
    //--Size info
    {
      strstream strm;
      strm << "Expected size: " << data_alloc_size;
      strm << "      Actual size: " << read_count;
      strm << ends;
      sendMessage(strm.str());
    }
    sendMessage("\n");
    return false;
  }

#if 0
// Debugging
    cerr << endl;
    cerr << "Object: " << od->object_name;
    if (od->object_id != -1)
      cerr << " " << od->object_id;
    cerr << endl;
    cerr << "Field:  " << od->field_name << endl;
    cerr << "Datatype: " << od->data_type << endl;
    cerr << "Data: ";
    output_data(cerr, od->data_type, od->data);
    cerr << "Buffer: " << readBuffer << endl;
    cerr << "Is sep string: " << isSeparatedStrings << endl;
#endif

  nofDataValues += read_count;

  int data_read_size = 0;

  // Actual total size of the read data
  for (int i = 0; i < read_count; i++)
    data_read_size += read_sizes[i];

  // Remove possible extra DATA entry space (relevant for string data!)
  if ( data_read_size < data_alloc_size ) {
    new_data_buffer_size += data_read_size - data_alloc_size;
    allocate_data_buffer(new_data_buffer_size, new_data_buffer_size);
  }

  update_data_lengths(read_count, read_sizes);

  delete[] read_sizes;

  if (bufferAtEnd()) {
    return false;
  }

  return true;
}


bool
LibFrontReader::read_entry_data_values(int old_data_buffer_size, int new_data_buffer_size, int data_count,
                                  int* data_sizes, int& read_count)
{
  if (data_count == 0) {
    read_count = 0;
    return true;
  }

  // Default data value sizes (for logical and numeric data)
  for (int i = 0; i < data_count; i++)
    data_sizes[i] = 1;

  // Character data
  if ( flags->readingProcedureName             ||
       currentBufferType == LibFront::EMF_STRING  ||
       currentBufferType == LibFront::EMF_FILE
     ) {


    if ( !read_string_data(old_data_buffer_size, new_data_buffer_size, data_count,
                           data_sizes, read_count)
       ) {
      return false;
    }

  // Logical data
  } else if ( currentBufferType == LibFront::EMF_LOGICAL ) {

    if ( !read_logical_data(old_data_buffer_size, new_data_buffer_size, read_count) ) {
      return false;
    }

  // Numeric data
  } else {

    if ( !read_numeric_data(old_data_buffer_size, new_data_buffer_size, read_count) ) {
      return false;
    }
  }

 return true;
}


// Read variables in the entry
bool
LibFrontReader::read_entry_variables(struct emf_ObjectData_X* od,
                                int var_alloc_count, int var_size)
{
  // Buffer sizes
  int old_var_buffer_size, new_var_buffer_size;

  // Extra size allocated for this entry
  int var_alloc_size;

  // Actual nof of read values
  int var_read_count;

  var_alloc_size = var_alloc_count * var_size;
  old_var_buffer_size = variableBufferSize;
  new_var_buffer_size = variableBufferSize + var_alloc_size;

  allocate_variable_buffer(old_var_buffer_size, new_var_buffer_size);

  currentBuffer = variableBuffer;
  currentBufferSize = variableBufferSize;
  currentBufferType = variableBufferType;

  // Sizes of the actually read values (normally 1, but string length for string data)
  int* read_sizes = new int[var_alloc_count];

  if ( !read_entry_variable_values(old_var_buffer_size, new_var_buffer_size,
                                   var_alloc_count, read_sizes, var_read_count)
     ) {

    // If no succes: ERROR MSG!
    sendMessage("\n");
    sendMessage("*** ERROR ***: When reading variable data for:");
    //--Object name, Field name
    {
      strstream strm;
      strm << "OBJECT: " << od->object_name;

      if (od->object_id != -1) {
        strm << " " << od->object_id;
      }
      strm << "   FIELD:  " << od->field_name;
      strm << ends;
      sendMessage(strm.str());
    }
    // Size info
    {
      strstream strm;
      strm << "Expected size: " << var_alloc_size;
      strm << "      Actual size: " << var_read_count;
      strm << ends;
      sendMessage(strm.str());
    }
    sendMessage("\n");

    return false;

  } else {
#if 0
// Debugging
    cerr << endl;
    cerr << "Object: " << od->object_name;
    if (od->object_id != -1)
      cerr << " " << od->object_id;
    cerr << endl;
    cerr << "Field:  " << od->field_name << endl;
    cerr << "Datatype: " << od->data_type << endl;
    cerr << "Data: ";
    output_data(cerr, od->data_type, od->data);
    cerr << "Buffer: " << readBuffer << endl;
    cerr << "Hanging sep string: " << hangingSeparatedString << endl;
    cerr << "Contains sep string: " << containsSeparatedStrings << endl;
#endif
  }

  nofVariableValues += var_read_count;

  int var_read_size = 0;
  // Actual total size of the read variable values
  for (int i = 0; i < var_read_count; i++)
    var_read_size += read_sizes[i];

  // Remove possible extra VARIABLE value space (relevant for name variables!)
  if ( var_read_size < var_alloc_size ) {
    new_var_buffer_size += var_read_size - var_alloc_size;
    allocate_variable_buffer(new_var_buffer_size, new_var_buffer_size);
  }
  update_variable_lengths(var_read_count, read_sizes);


  return true;
}


bool
LibFrontReader::read_entry_variable_values(int old_var_buffer_size, int new_var_buffer_size, int var_count,
                                      int* var_sizes, int& read_count)
{
  if (var_count == 0) {
    read_count = 0;
    return true;
  }

  // Default variable value sizes (for numeric data)
  for (int i = 0; i < var_count; i++)
    var_sizes[i] = 1;

  // Character variable values (a name string!)
  if ( currentBufferType == LibFront::EMF_STRING ) {

    if ( !read_string_data(old_var_buffer_size, new_var_buffer_size, var_count,
                           var_sizes, read_count) ) {
      return false;
    }

  // Numeric data
  } else {

    if ( !read_numeric_data(old_var_buffer_size, new_var_buffer_size, read_count) ) {
      return false;
    }
  }

 return true;
}


// Read field name (optional) and data entries
//
bool
LibFrontReader::read_field(struct emf_ObjectData_X*  od)
{
  static char peek_buffer[emf_MAX_STRING_LEN];

  flags->endsField = false;

  flags->readingNonSizedData = false;

  bool size_given, type_given, variable_given;

  if ( flags->startsNewField ) {

    dimension1 = 1;
    dimension2 = 1;
    nofEntries = 0;
    nofVariables = 0;
    od->field_name[0] = '\0';
    od->is_procedure = false;

    //---Read field name. NOTE: In some cases name is not needed!
    if (! read_field_name(od, od->field_name, emf_MAX_NAME_LEN) )
      return false;

    od->field_name_length = strlen(od->field_name);
    flags->startsNewField = false;

    // If there were group names (like [A B C] ) after field name
    if ( flags->readingGroupNames )
      read_group_names(od);
  }

  // All field data have been read
  if ( flags->endsField ) {
    flags->endsField = false;
    return true;
  }

  //---Check if Variable-keyword is given
  if (bufferAtEnd())
    return false;

  variable_given = false;

  if ( next_is_variable_keyword() ) {
    getNextToken(tokenBuffer);
    variable_given = true;
  }

  //---If Variable given, read variable names
  //   NOTE: Read until the Size-keyword or a Type keyword comes!
  if (variable_given) {
    while (1) {

      if (bufferAtEnd())
        return false;

      if (next_is_type_keyword() || next_is_size_keyword() )
        break;

      getNextToken(od->variable_names[nofVariables]);
      od->variable_name_lengths[nofVariables] = strlen(od->variable_names[nofVariables]);
      nofVariables++;

      if (nofVariables > emf_MAX_NOF_VARIABLES)
        return false;
    }
  }

  if ( variable_given && nofVariables == 0 )
    return false;

  //---Check if Size-keyword is given
  if (bufferAtEnd())
    return false;

  size_given = false;
  if ( next_is_size_keyword() ) {
    getNextToken(tokenBuffer);
    size_given = true;
  }

  if (bufferAtEnd())
    return false;

  //---If Size given, read dimensions
  if (size_given) {
    getNextNumber(dimension1);

    // If second (column) dimension is given
    // read it from the stream
    int test_nbr;
    if ( try_to_read(test_nbr) ) {
      getNextNumber(dimension2);
    }
  }

  // NOTE: When variable is given, data rows are ended
  // by the "end" mark and Size-parameters tell the size  of
  // the data entry rows!!!
  if ( nofVariables > 0 ) flags->readingEndMarkedData = true;

  //---Check data type keyword
  if (bufferAtEnd()) return false;

  type_given = false;
  if ( next_is_type_keyword() ) type_given = true;

  if ( size_given  &&
       !type_given
     ) {
    // ERROR MSG
    strstream strm;
    strm <<  "UNKNOWN DATA TYPE: " << peek_buffer << ends;
    sendMessage(strm.str());
    sendMessage("\n");

    return false;
  }

  //---If field-name is not yet resolved
  if ( flags->hangingFieldNameNumber ) {
    flags->hangingFieldNameNumber = false;
    // Hanging number belongs to the field name
    if (type_given) {
      strstream strm;
      strm << od->field_name << ' ' << (int)hangingNumber << ends;
      strm >> od->field_name;
      od->field_name_length = strlen(od->field_name);
    // We had one single data value
    } else {
     return storeNumber(od, LibFront::EMF_REAL, &hangingNumber);
    }
  }

  //---Read possible data type
  if (type_given) {
    getNextToken(od->data_type);
    LibFront::toLower(od->data_type);
    od->data_type_length = strlen(od->data_type);

    if ( 0 == strcmp(od->data_type, "procedure") ) {
      od->is_procedure = true;
    }
  }

  //---Check if we have a procedure definition (for int,real
  //   data types
  bool read_procedure = false;
  if (bufferAtEnd()) return false;
  if ( peek_token(peek_buffer) ) {
    LibFront::toLower(peek_buffer);
    read_procedure = is_procedure_keyword(peek_buffer);
  }

  //---Set possible data type
  if (type_given) {
    // integer
    if ( 0 == strcmp(od->data_type, "integer") ) {
      dataBufferType = LibFront::EMF_INTEGER;
      flags->readingProcedureName = read_procedure;

    // real
    } else if ( 0 == strcmp(od->data_type, "real") ) {
      dataBufferType = LibFront::EMF_REAL;
      flags->readingProcedureName = read_procedure;

    // logical
    } else if ( 0 == strcmp(od->data_type, "logical") ) {
      dataBufferType = LibFront::EMF_LOGICAL;
      flags->readingProcedureName = read_procedure;

    // string
    } else if ( 0 == strcmp(od->data_type, "string") ) {
      dataBufferType = LibFront::EMF_STRING;
      flags->readingProcedureName = read_procedure;

    // file name
    } else if ( 0 == strcmp(od->data_type, "file") ) {
      dataBufferType = LibFront::EMF_FILE;
      flags->readingFileName = true;

    // unknown --> string
    } else {
      dataBufferType = LibFront::EMF_STRING;
      strcpy(od->data_type, "string");
    }
  }

  // Reading field data
  // ==================

  flags->readingFieldData = true;
  dataAsStringBuffer[0] = '\0';
  dataAsStringPos = 0;

  //---Non-sized data
  if (!size_given) {
    return read_field_data_non_sized(od, type_given);
  }

  //---Sized data
  od->dimension1 = dimension1;
  od->dimension2 = dimension2;
  od->nof_variables = nofVariables;

  // Default sizes for variable value and data parts
  // in the entry row
  int var_count = nofVariables;
  int var_size = 1;
  int data_count = (dimension1 * dimension2);
  int data_size = 1;

  // Set size info for these data types
  if ( 0 == strcmp(od->data_type, "string")  ||
       0 == strcmp(od->data_type, "file")
     ) {
    data_size = readBufferLen; // we do not know the string size in advance!
  }

  int read_count = 0;

  return read_field_data_sized(od, var_count, var_size, data_count, data_size, read_count);
}


// Read LOOP, read variables and data in the entries
// *_count arguments tell the number of varaible values and data values
// *_size arguments tell the size to allocated for theses values
// NOTE: size is normally 1, but for strings it is the text buffer length
//       and this size must be updated when the value is read, because we
//       do not know the size in advance!
bool
LibFrontReader::read_field_data_sized(struct emf_ObjectData_X* od,
                                 int var_alloc_count, int var_size,
                                 int data_alloc_count, int data_size,
                                 int& read_count)
{
  nofDataValues = 0;
  nofEntries = 0;
  nofVariableValues = 0;

  // Check variable value types (NOTE: Real or String currently!)
  if ( nofVariables > 0 ) {

    if ( next_is_number() ) {
      variableBufferType = LibFront::EMF_REAL;
      strcpy(od->variable_type, "real");
      od->variable_type_length = strlen(od->variable_type);

    } else {
      variableBufferType = LibFront::EMF_STRING;
      strcpy(od->variable_type, "string");
      od->variable_type_length = strlen(od->variable_type);
      var_size = readBufferLen;
    }
  }

  delete_data_buffer();
  delete_variable_buffer();

  // Read LOOP
  while (1) {
    nofEntries++;

    // Read possible variable values
    if ( nofVariables > 0) {

      if (bufferAtEnd()) {
        return false;
      }

      read_entry_variables(od, var_alloc_count, var_size);
    }

    if (bufferAtEnd()) {
      return false;
    }

    // Read data values
    read_entry_data(od, data_alloc_count, data_size, read_count);


    if (bufferAtEnd()) {
      return false;
    }

    if (data_at_end()) {

      flags->readingFieldData = false;

      // Consume "End" mark for the data
      if ( flags->readingEndMarkedData ) {
        getNextToken(tokenBuffer);
        flags->readingEndMarkedData = false;
      }
      break;
    }

  }

  // Update transfer data
  od->data = dataBuffer;
  od->data_length = dataBufferSize;
  od->data_lengths = dataLengths;

  if ( currentBufferType == LibFront::EMF_STRING ) {
    od->nof_entries = read_count;
  } else {
    od->nof_entries = nofEntries;
  }

  od->variable_values = variableBuffer;
  od->variable_values_length = variableBufferSize;
  od->variable_values_lengths = variableLengths;

  return true;
}


// Read "simple" data entries when Size,Type is not given
// data type is either Real (for numbers) or String
bool
LibFrontReader::read_field_data_non_sized(struct emf_ObjectData_X*  od, bool type_given)
{
  flags->readingNonSizedData = true;
  flags->nonSizedDataAtEnd = false;

  delete_data_buffer();

  int data_count = 1024;
  int data_size  = 1;

  od->dimension1 = 1;
  od->dimension2 = 1;
  od->nof_variables = 0;

  if (!type_given) {
    if ( next_is_number() ) {
      dataBufferType = LibFront::EMF_REAL;
      strcpy(od->data_type, "real");
    }
    else {
      dataBufferType = LibFront::EMF_STRING;
      strcpy(od->data_type, "string");
      data_count = 10; // Max 10 strings per entry!
      data_size = readBufferLen;
    }
  }

  od->data_type_length = strlen(od->data_type);

  int read_count = 0;
  bool rc = read_field_data_sized(od, 0, 0, data_count, data_size, read_count);

  od->dimension1 = od->data_length;

  return rc;
}


bool
LibFrontReader::read_field_name(struct emf_ObjectData_X*  od,
                           char* name_buffer, int buffer_len)
{
  if (bufferAtEnd())
    return false;

  //-Check if End-keyword is given, then we have
  // no field-name or data,  just object and its data!!!
  if ( next_is_end_keyword() ) {
    flags->endsField = true;
    return true;
  }

  //-Check if Type or Variable-keyword is given, then we have
  // no field-name just object and data!!!
  if ( next_is_type_keyword() || next_is_variable_keyword() ) {
    return true;
  }

  //-We have to check if there is any number in the buffer
  // to resolve the the "Fieldname 1" problem
  int nbr_of_tokens = 0;
  if ( !peek_field_name_tokens(nbr_of_tokens) )
    return false;

  constructStringFromTokens(name_buffer, nbr_of_tokens);

  if ( name_buffer == NULL || name_buffer[0] == '\0' ) {
    return false;
  } else {
    return  true;
  }

}


bool
LibFrontReader::read_group_names(struct emf_ObjectData_X*  od)
{
  if (bufferAtEnd())
    return false;

  char groups_buffer[1024];
  char group_buffer[1024];
  int len;

  // Pick away group starting character '['
  readBufferPos++;

  // Read until next ']'
  getToSeparator(groups_buffer, ']');

  LibFront::trim(groups_buffer);
  len = strlen(groups_buffer);

  // Make input stream from the whole string between []
  istrstream groups_strm(groups_buffer, len);

  short group_index = 0;
  // Split string by group "grouper" characeter ';'
  while ( !groups_strm.eof() ) {
    groups_strm.getline(group_buffer, 1024, ';');

    LibFront::trim(group_buffer);
    len = strlen(group_buffer);

    // Make input stream from one group string between ';' characters
    istrstream group_strm(group_buffer, len);

    short name_index = 0;
    while( !group_strm.eof() ) {

      group_strm >> od->group_names[group_index][name_index];

      od->group_name_lengths[group_index][name_index] =
          strlen(od->group_names[group_index][name_index]);

      od->nof_group_names[group_index] = ++name_index;
    }

    if (++group_index > emf_MAX_NOF_GROUPS)
      break;
  }

  return true;
}


bool
LibFrontReader::read_logical_data(int start_pos, int max_pos, int& read_count)
{
  static char buffer[emf_MAX_STRING_LEN];
  int logical_value;
  double numeric_value;

  read_count = 0;

  for (int i = start_pos; i < max_pos; i++) {

    if (bufferAtEnd()) {
      return false;
    }

    //---Read logical value
    //-text True
    if ( next_is_text("True") ) {
      logical_value = 1;
      getNextToken(buffer);

    //-text False
    } else if ( next_is_text("False") ) {
      logical_value = 0;
      getNextToken(buffer);

    //-numeric (0/1) value
    } else if ( next_is_number() ) {
      getNextNumber(numeric_value);
      if (numeric_value > 0)
        logical_value = 1;
      else
        logical_value = 0;

    //-ERROR
    } else {
        return false;
    }

    ((int*)currentBuffer)[i] = logical_value;
    read_count++;

  } // for

  return true;
}


bool
LibFrontReader::read_numeric_data(int start_pos, int max_pos, int& read_count)
{

  read_count = 0;

  for (int i = start_pos; i < max_pos; i++) {

    if (bufferAtEnd()) {
      return false;
    }

    if ( flags->readingNonSizedData ) {
      if ( flags->nonSizedDataAtEnd || !next_is_number() )
        return true;
    }

    switch (currentBufferType) {

    case LibFront::EMF_INTEGER:
    case LibFront::EMF_LOGICAL:
      int int_number;
      if ( !getNextNumber(int_number) ) {
        return false;
      }
      ((int*)currentBuffer)[i] = int_number;
      break;

    case LibFront::EMF_REAL:
      // NonSized is always read as doubles and next non-number
      // ends the data and we should not try to tread it!!!
      // Normal Sized data should always succeed
      double dbl_number;
      if ( !getNextNumber(dbl_number) ) {
        return false;
      }
      ((double*)currentBuffer)[i] = dbl_number;
      break;

    } // switch

    read_count++;

  } // for

  return true;
}


bool
LibFrontReader::read_object(struct emf_ObjectData_X* od)
{
  //---Read object (section) name and a possible (id) number
  //   or a separated string (data like "string data")

  od->object_name[0] = '\0';
  od->object_id = -1;

  if ( !read_object_name(od, od->object_name, emf_MAX_NAME_LEN) ) {
    return false;
  }

  od->object_name_length = strlen(od->object_name);

   //-If there were group names (like [A B C] ) after object name
  if ( flags->readingGroupNames ) {
    read_group_names(od);
  }

#if 0
  // For debbugging
  cerr << od->object_name;
  if (od->object_id != -1)
    cerr << " " << od->object_id;
  cerr << endl;
#endif

  return true;
}


bool
LibFrontReader::read_object_name(struct emf_ObjectData_X*  od,
                            char* name_buffer, int buffer_len)
{
  if (next_is_separated_string()) {
    return false;
  }

  int nbr_of_tokens;

  if ( !peek_object_name_tokens(nbr_of_tokens) ) {
    return false;
  }


  constructStringFromTokens(name_buffer, nbr_of_tokens);

  // Remove white-space from the end of the string
  LibFront::trimRight(name_buffer);

  // Read possible id
  if ( next_is_number() ) {
    getNextNumber(od->object_id);
  }

  return true;
}


bool
LibFrontReader::read_string_data(int start_pos, int max_pos, int string_count,
                            int* string_sizes, int& read_count)
{
  // Result will be put into currentBuffer
  char* result_string = (char*)currentBuffer;

  read_count = 0;

  int string_pos = start_pos;
  int string_pos_prev = start_pos;

  for (int i = 0; i < string_count; i++) {

    bufferEatWs();

    if ( !checkReadBuffer() ) {
      return false;
    }

    // New line (ie. lf without continuation mark)
    // has stopped non-sized data reading!
    if ( flags->readingNonSizedData &&
         flags->nonSizedDataAtEnd
       )
      return true;

    // Append next token to the result_string
    char* buffer = result_string + string_pos;

    if ( -1 == getNextToken(buffer) ) return false;

    // Find current end position of the result
    //
    // NOTE: We will start writing from the NULL character, so
    // string-pos must be on that character --> -1.
    // If null charcaters are left in the buffer, they will
    // truncate the data!
    // BE CAREFUL WITH THIS!!!***!!!
    //
    string_pos = strlen(result_string) - 1;

    // Remove white-space from the end of the string
    //
    // NOTE: This function just tells the position of the last
    // non-ws character, it does NOT manipulate the string itself!!!
    //
    LibFront::trimRight(result_string, start_pos, string_pos);

    // This many was actually put into array
    int string_size = 1 + string_pos - string_pos_prev;

    string_sizes[i] = string_size;

    string_pos++;
    string_pos_prev = string_pos;

    read_count++;
  }

  return true;
}


void
LibFrontReader::reset_flags()
{
  flagsStack[currentFlagsIndex].init();
}


int
LibFrontReader::sendMessage(char* msg)
{
  cerr << msg;

  if (emf_readMessageCB == NULL || msgBuffer == NULL) return false;

  strncpy(msgBuffer, msg, msgBufferSize);

  return emf_readMessageCB(msgBuffer);
}


void
LibFrontReader::set_transfer_info()
{
  struct emf_ObjectData_X*  od = emf_ObjectData;

  init_transfer_info();

  // Data values
  if ( LibFront::ncEqual(od->data_type, "logical") ) {
    LibFront::logicalData = (int*) od->data;
    LibFront::isLogicalData = true;
  }
  else if ( LibFront::ncEqual(od->data_type, "integer") ) {
    LibFront::integerData = (int*) od->data;
    LibFront::isIntegerData = true;
  }
  else if ( LibFront::ncEqual(od->data_type, "real") ) {
    LibFront::realData = (double*) od->data;
    LibFront::isRealData = true;
  }
  else if ( LibFront::ncEqual(od->data_type, "string") ) {
    LibFront::stringData = (char*) od->data;
    LibFront::isStringData = true;
  }
  else if ( LibFront::ncEqual(od->data_type, "file") ) {
    LibFront::stringData = (char*) od->data;
    LibFront::isFileData = true;
  }
  else {
    LibFront::voidData = (void*) od->data;
    LibFront::isVoidData = true;
  }

  LibFront::dimension1 = od->dimension1;
  LibFront::dimension2 = od->dimension2;
  LibFront::dataLengths = od->data_lengths;

  // Variable values
  if ( LibFront::ncEqual(od->variable_type, "real") ) {
    LibFront::realVariable = (double*) od->variable_values;
    LibFront::isRealVariable = true;
  }
  else if ( LibFront::ncEqual(od->variable_type, "string") ) {
    LibFront::stringVariable = (char*) od->variable_values;
    LibFront::isStringVariable = true;
  }
  else {
    LibFront::voidVariable = (void*) od->variable_values;
    LibFront::isVoidVariable = true;
  }

  LibFront::variableValuesLengths = od->variable_values_lengths;
}


bool
LibFrontReader::start()
{
  struct emf_ObjectData_X*  od = emf_ObjectData;
  od->dataAsString = NULL;
  od->dataHasMatcVars = 0;

  init_transfer_info();

  flags->startsNewObject = true;

  for (int i = 0; i < emf_MAX_NOF_GROUPS; i++) {
    od->nof_group_names[i] = 0;
  }

  bufferEatWs();

  if (bufferAtEnd()) {
    return false;
  }

  while ( !bufferAtEnd() ) {

    //---New object is started
    if ( flags->startsNewObject ) {

      od->is_object_start = 1;
      od->is_object_end = 0;

      if ( !read_object(od) ) {
        strstream strm;
        strm <<  "***ERROR When reading object data for object: " << od->object_name << " " << od->object_id;
        strm << ends;
        sendMessage(strm.str());
        sendMessage("\n");
        return false;
      }

      flags->startsNewObject = false;

    //---New field for the object
    } else {

      flags->startsNewField = true;
      dataHasMatcVars = false;

      if ( !read_field(od) ) {
        strstream strm;
        strm <<  "***ERROR When reading field data for object: " << od->object_name << " " << od->object_id;
        strm << ". Reading field: " << od->field_name << ends;
        sendMessage(strm.str());
        sendMessage("\n");
        return false;
      }

      dataAsStringBuffer[dataAsStringPos] = '\0';

      flags->readingFieldData = false;
      flags->hasEvaluatedBuffer = false;

      if ( object_at_end() || LibFront::in(tokenBuffer, "Include") ) {

        od->is_object_end = 1;

        // Consume the possible object end marker "End"
        if ( object_at_end() ) {
          getNextToken(tokenBuffer);
        }

        flags->startsNewObject = true;
      }

      set_transfer_info();

      // Send data after each field
      delete[] od->dataAsString;
      od->dataAsString = new char[1+dataAsStringPos];
      strcpy(od->dataAsString, dataAsStringBuffer);
      od->dataHasMatcVars = int(dataHasMatcVars);

      int rc = sendData(userData);

      // ERROR!!!
      if ( rc != emf_OK ) {
        return false;
      }

      od->is_object_start = 0;

      // NOTE: Group data is "reset" after each send ==>
      //       If you want to read object-level group
      //       info, it must be done when "object-start = 1"!!!
      for (int i = 0; i < emf_MAX_NOF_GROUPS; i++)
        od-> nof_group_names[i] = 0;

    }

    flags->nonSizedDataAtEnd = false;
    flags->hangingFieldNameNumber = false;
    flags->readingNonSizedData = false;
    flags->readingFileName = false;
    flags->readingFieldData = false;
    flags->readingGroupNames = false;
    flags->readingProcedureName = false;

  }

  return true;
}


// Store field input data as a string.
// This is used to store Matc-expressions in input data etc.
//
void
LibFrontReader::storeDataAsString(const char* buffer, bool is_matc)
{
  // Add a separator before new data
  //
  if ( dataAsStringPos > 0 ) {
    dataAsStringBuffer[dataAsStringPos++] = ' ';
  }

  char* tmp = dataAsStringBuffer;
  tmp += dataAsStringPos;

  // Add $-sign in front of an Matc-expression, so it can be
  // used later as an input data!!!
  //
  if ( is_matc ) {
    strncpy(tmp, &matcSeparator, 1);
    strcpy(tmp+1, buffer);

  } else {
    strcpy(tmp, buffer);
  }

  // Add new data
  //
  dataAsStringPos += strlen(tmp);
}


// Read one given number
bool
LibFrontReader::storeNumber(struct emf_ObjectData_X* od,
                       LibFront::dataValueType data_type,
                       void* data)
{
  // Delete old buffer
  delete_data_buffer();

  // Allocate new with given type
  dataBufferType = data_type;
  allocate_data_buffer(dataBufferSize ,1);

  nofEntries = 1;
  nofDataValues = 1;

  int data_sizes[1] = {1};

  update_data_lengths(1, data_sizes);

  switch (data_type) {
  case LibFront::EMF_LOGICAL:
    strcpy(od->data_type,"logical");
    ((int*)dataBuffer)[0] = *((int*)data);
    break;
  case LibFront::EMF_INTEGER:
    strcpy(od->data_type, "integer");
    ((int*)dataBuffer)[0] = *((int*)data);
    break;
  case LibFront::EMF_REAL:
    strcpy(od->data_type, "real");
    ((double*)dataBuffer)[0] = *((double*)data);
    break;
  }

  od->data = dataBuffer;
  od->data_length = dataBufferSize;
  od->data_lengths = dataLengths;
  od->nof_entries = nofEntries;

  return true;
}


template <class T> bool
LibFrontReader::try_to_read(T& data)
{
  if (bufferAtEnd()) {
    return false;
  }

  push_flags();
  int old_pos = readBufferPos;

  bool rc = getNextNumber(data);

  readBufferPos = old_pos;
  pop_flags();

  return rc;
}


bool
LibFrontReader::update_data_lengths(int data_count, int* data_sizes)
{
  int old_nof_values = nofDataValues - data_count;

  int* data_lengths = new int[nofDataValues];

  int i;
  for (i = 0; i < old_nof_values; i++)
    // Copy old values
    data_lengths[i] = dataLengths[i];
  for (i = 0; i < data_count; i++)
    // Add new values
    data_lengths[old_nof_values + i] = data_sizes[i];

  delete[] dataLengths;

  dataLengths = data_lengths;

  return true;
}


bool
LibFrontReader::update_variable_lengths(int var_count, int* var_sizes)
{
  if (nofVariableValues == 0)
    return true;

  int old_nof_values = nofVariableValues - var_count;

  int* var_lengths = new int[nofVariableValues];

  int i;
  // Copy old values
  for (i = 0; i < old_nof_values; i++)
    var_lengths[i] = variableLengths[i];
  // Add new values
  for (i = 0; i < var_count; i++)
    var_lengths[old_nof_values + i] = var_sizes[i];

  delete[] variableLengths;

  variableLengths = var_lengths;

  return true;
}
// =========================================================================
// END: LibFrontReader class implementation
// =========================================================================



