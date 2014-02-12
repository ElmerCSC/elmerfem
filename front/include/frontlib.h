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
Program:    ELMER model file (emf) reader 
Module:     frontlib.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  
 
Abstract:   Header file for the emf-parser

************************************************************************/
 
#ifndef _FRONTLIB_
#define _FRONTLIB_


// OK/ERROR codes for Reader
#define emf_OK 0
#define emf_ERROR 1


// *************************************************************************
// *************************************************************************
//
// START: LibFront namespace definitions
//
// *************************************************************************
// *************************************************************************
//
namespace LibFront {

#if defined(WIN32)
  #define NSP_EXTERN extern
#else
  #define NSP_EXTERN extern
#endif

// Data value types 
// ================
enum dataValueType {
  EMF_INTEGER, EMF_REAL, EMF_LOGICAL, EMF_STRING, EMF_FILE
};


// Data values related flags set in the library by the Reader
extern bool isFileData;
extern bool isIntegerData;
extern bool isLogicalData;
extern bool isRealData;
extern bool isStringData;
extern bool isVoidData;
extern int* integerData; 
extern int* logicalData;
extern double* realData;
extern char* stringData;
extern void* voidData;
extern int  dimension1;
extern int  dimension2;
extern int* dataLengths;

// Variables related flags set in the library by the Reader
extern bool isRealVariable;
extern bool isStringVariable;
extern bool isVoidVariable;
extern double* realVariable;
extern char* stringVariable;
extern void* voidVariable;
extern int* variableValuesLengths;

extern char stringSeparator; // Default '\"';
extern char matcSeparator;   // Default = '$';

void setStringSeparator(char sep);
void setMatcSeparator(char sep);

// Read file into given buffer
// NOTE: Buffer is allocated in the function!
// Skip empty line, comment lines and eol comments
// NOTE: Comments can start with ! or #
// Return: nof chars read (= actual buffer size)
//
extern int readFile(char* filename, char*& buffer, bool skip_comments = true);


// Matc utilities
// ==============
extern void initMatcFrmt();
extern void formatMatcError(char* err_msg);

// NOTE: Whitespace must have been trimmed from 'str' before calling
// these functions
extern bool isMatcDefinition(const char* str);
extern bool hasMatcVariables(const char* str);
extern bool getMatcName(const char* str, char* name_buffer, int max_len, bool& is_var, bool& is_func);
extern bool isMatcError(const char* str);
extern char* evalMatcString(const char* str);

// Read a Matc-expression from the source, starting from source_pos
// Return source-position after reading the expression
// NOTE: Source is NOT moved to this position actually, but the client
// can use this return value if its wants to advandec the source!
//
extern int readMatcExpression(char* source, int source_pos, int source_end,
                              char* result, int max_result_len,
                              bool skip_commnets = true);


// Data reading utilities for getting data from global
// LibFront arries set by the LibFrontReader
// ===================================================
//
// NOTE: These function use global data arries as the source!!!
//
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
extern bool isNumber(const char* str);
extern bool ncEqual(char* str1, char* str2);
extern bool ncEqualPartial(char* str1, char* str2);
extern char* trim(char* buffer, bool trim_nl = false);
extern char* trimLeft(char* buffer, bool trim_nl = false);
extern char* trimRight(char* buffer, bool trim_nl = false);
extern char* toLower(char* buffer);
extern char* toLower(const char* source, char* target);
extern char* toUpper(char* buffer);
extern char* toUpper(const char* source, char* target);

// Find start/end of the first/last non-ws character
// NOTE: Initial start/end position mus be given in the reference
// arguments start_pos/end_pos
// NOTE: String itself is NOT modified!!!
//
extern void trimLeft(const char* str, int max_pos, int& start_pos, bool trim_nl = false);
extern void trimRight(const char* str, int min_pos, int& end_pos, bool trim_nl = false);



// Data output functions
// =====================
// Simple utility for debugging 
extern ostream&
output_data(ostream& out, char* data_type, void* data);


// Output a Matc definition
// ------------------------
// NOTE: add_dsign <--> add dollar sign
//
extern ostream&
output_matcDef(ostream& out, short indent_size, short indent_level,
               const char* field_name, const char* field_type,
               const char* def, bool add_dsign = false);


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
// =========================================================================
// END: LibFront namespace definitions
// =========================================================================



// *************************************************************************
// *************************************************************************
//
// START: Definition for LibFrontReader external interface
//
// *************************************************************************
// *************************************************************************
extern "C" {

// Max length of one single condition when delivered as 'raw' string
#define emf_MAX_STRING_LEN 8191

// Max length of object/field names etc
#define emf_MAX_NAME_LEN 64

// Max number of groups (name groups separated by ';' like [A B; C] afet object/field names
#define emf_MAX_NOF_GROUPS 2

// Max number of names in one group
#define emf_MAX_NOF_GROUP_NAMES 24

// Max number of variables in data parameters (like time, temperature etc)
#define emf_MAX_NOF_VARIABLES 12


// This struct stores the data for one parameter field
// ---------------------------------------------------
//
struct emf_ObjectData_X
{
  // OBJECT LEVEL RELATED STUFF
  // --------------------------

  // Flag for a new object 
  int is_object_start;

  // Flag for the end of the object 
  int is_object_end;

  // Object or section name like: Header, Initial Condition 
  char object_name[emf_MAX_NAME_LEN];

  // Actual length of the object-name character array 
  int object_name_length;

  // Possible id for the object, -1 if no id is used 
  int object_id;


  // FIELD LEVEL RELATED STUFF
  // -------------------------

  // NOTE: Field name is optional, all the data can also be at object level!

  // Name for the parameter field: Equation, Velocity1 etc.
  char field_name[emf_MAX_NAME_LEN];

  // Actual length of the field_name character array 
  int field_name_length;


  // GROUP LEVEL RELATED STUFF
  // Nof group names (like [H; N O] after object/field name
  // NOTE: character ';' is groups "grouper", currently we accept two different groups!
  // NOTE: if you have groups, you have to give same number of data values for each group!
  //
  int nof_group_names[emf_MAX_NOF_GROUPS];

  // Group names [ H2O N; C] etc.
  char group_names[emf_MAX_NOF_GROUPS][emf_MAX_NOF_GROUP_NAMES][emf_MAX_NAME_LEN];

  // Actual lengths of the group names
  int group_name_lengths[emf_MAX_NOF_GROUPS][emf_MAX_NOF_GROUP_NAMES];


  // DATA VALUES RELATED STUFF
  // -------------------------

  // Name for the data type: Integer, Real, Logical, String, File etc.
  char data_type[emf_MAX_NAME_LEN];

  // Actual length of the data_type character array
  int data_type_length;

  // If definition is for a dll-procedure
  int is_procedure;

  // Nof entries (rows of variables and data)
  int nof_entries;

  // Dimension1 for a data entry (default 1)
  int dimension1;

  // Dimension2 for a data entry (default 1) 
  int dimension2;

  // Lengths for data values.
  // SIZE: nof_entries * dimension1 * dimension2
  // NOTE: Normally each length is 1 (for a nuemric scalar item)
  //       But if data is of string type we normally have variable
  //       length strings and then the length of
  //       each string value can be read from this array
  //
  //       Actually you have to read string data using these
  //       lengths because there are no sperators between the strings
  //       in the "data" array !!!
  
  int* data_lengths;

  // Length of the data array = Sum(data_lengths), for convenience!
  int data_length;

  // Data array of homogenous data type
  // NOTE: cast to proprer type before using!
  void* data;

  // Field data as a string
  char* dataAsString;

  // If field data contains matc-variables
  int dataHasMatcVars;


  // VARIABLES RELATED STUFF
  // -----------------------

  // Name for the varaible value type. NOTE: Real or String only!
  char variable_type[emf_MAX_NAME_LEN];

  // Length of the variable_type character array
  int variable_type_length;

  // Nof argument variables
  int nof_variables;

  // Argument variable names (Temperature, Times etc)
  char variable_names[emf_MAX_NOF_VARIABLES][emf_MAX_NAME_LEN];

  // Lengths of the variable-name arrays
  int variable_name_lengths[emf_MAX_NOF_VARIABLES];

  // Lengths for variable values
  // SIZE = nof_entries * nofVariables
  //
  // NOTE: Normally these values are all 1, because a typical
  //       variable value is a (real) scalar. But we can also have
  //       names as variables and then the lengths of these character arries
  //       must be taken from this array (ref. data_lengths!)
   
  int* variable_values_lengths;

  // Total length of the variable values array = Sum(variable_lengths), for convenience!
  int variable_values_length;

  // Variable values array of homogenous data type (double or char only!)
  // NOTE: cast to proprer type before using!
  void* variable_values;

}; // end struct emf_ObjectData_X


// Declaration for library interface

#ifdef WIN32
  #ifdef BUILD_DLL
    #define EMF_DECL __declspec(dllimport)
  #else
    #define EMF_DECL extern
  #endif
#else
  #define EMF_DECL extern
#endif


// Pointer to the object data transfer structure.
// ----------------------------------------------
//
// Library function(s) use this pointer when transferin the data!
// Connect this pointer to your the variable which is allocated
// in your own module:
// declaration:  emf_ObjectData_X my_object_data;
// allocation:   emf_ObjectData = &my_object_data;

EMF_DECL  struct emf_ObjectData_X*  emf_ObjectData;

// Buffer for unknown data lines
// NOTE: This is dynamically allocated in the library, size: 1 + emf_MAX_STRING_LEN
EMF_DECL  char* emf_UnknownData;


// Interface FUNCTION decalrations
// ===============================

// This function starts parsing model-file
// ---------------------------------------
//
// Call it first in your module!!!
//
// Arguments:
// model_data_file: full path for the model file
// user_data      : this data is delivered back in the emf_readDataCB
// msg_buffer_size: maximum writeable size for the message buffer
// msg_buffer     : buffer for getting messages in the emf_messageReadCB
//
// Remember to allocate msg_buffer!
// If it is set NULL, call-back function is not used!

// Return values (int):
//   0 --> read ok
//   1 --> parsing was NOT ok, ERROR!!!

EMF_DECL int emf_readData(char* model_data_file, void** user_data,
                          int msg_buffer_size, char* msg_buffer,
                          int expand_matc);


// This call-back function does the actual parsing
// -----------------------------------------------
//
// It returns data using data-structures defined above.
//
// Connect this function to the call-back function
// defined in your module.
//
// LIKE: int my_read_function (void) {...};
//       emf_readDataCB = my_read_function;
//
// Return values (int):
//   -1 --> stop reading
//  anything else --> continue

EMF_DECL  int (*emf_readDataCB) (void** user_data);


// This call-back function to send error messages etc. from the parser.
// -------------------------------------------------------------------
//
// Connect this function to the call-back function
// defined in your module.
//
// LIKE: void my_msg_function (char* message_buffer) {...};
//       emf_readMessageCB = my_msg_function;
//
// Return values (int)
//   1 <--> ok
//   0 <--> no success

EMF_DECL  int (*emf_readMessageCB) (char* msg_buffer);

} // end extern "C"

// =========================================================================
// END: Definition for LibFrontReader external interface
// =========================================================================




// *************************************************************************
// *************************************************************************
//
// START: File reader class (= parser for egf/emf type files)
//
// *************************************************************************
// *************************************************************************
//
struct LibFrontReaderFlags {
  LibFrontReaderFlags();
  LibFrontReaderFlags(LibFrontReaderFlags& flags);
  void copy(LibFrontReaderFlags& flags);
  void init();

  bool nonSizedDataAtEnd;
  bool canSkipSpace;
  bool canSkipNewline;
  bool endsField;
  bool endsObject;
  bool equalSignPeeked;
  bool hasEvaluatedBuffer;
  bool newlinePeeked;
  bool hangingFieldNameNumber;
  bool matcPeeked;
  bool readingEndMarkedData;
  bool readingFieldData;
  bool readingFileName;
  bool readingGroupNames;
  bool readingNonSizedData;
  bool readingProcedureName;
  bool startsNewField;
  bool startsNewObject;
};


// Max length for LibFrontReader work buffers
// (matcBuffer, dataAsStringBuffer etc.)
#define MAX_BUFFER_LEN 640000

// LibFrontReader class
// ============
class LibFrontReader
{
public:
  LibFrontReader(char* filename, void** user_data,
           int msg_buffer_size, char* msg_buffer,
           bool is_matc_file = false);
  ~LibFrontReader();
  bool isOk() { return readerOk; }
  bool start();
protected:
  bool allocate_data_buffer(int buffer_size, int item_count);
  bool allocate_variable_buffer(int buffer_size, int item_count);
  bool appendReadBuffer(char c, int& bufferPos, bool init);
  bool bufferAtEnd();
  void bufferEatWs();
  void bufferEatWs(bool& eq_eaten, bool& nl_eaten);
  bool case_compare(char* str1, char* str2);
  bool checkReadBuffer();
  void constructStringFromTokens(char* buffer, int nof_string_token);
  bool data_at_end();
  void delete_data_buffer();
  void delete_variable_buffer();
  bool extractSeparator(char sep);
  bool evalMatcBuffer(char* source);
  void evalAll();
  int evalAndInsertNextMatc();
  void evalAndSkipNextMatc();
  bool getNextNumber(int& number);
  bool getNextNumber(double& number);
  template <class T> bool getNextNumber_impl(T& number);
  int getNextToken(char* buffer);
  int getNextToken(char* buffer, bool& eq_found, bool& nl_found, bool& matc_found);
  int getToSeparator(char* buffer, char sep);
  bool handleEol();
  void init_transfer_info();
  int insertReadBuffer(char* data);
  bool is_file_keyword(char* str);
  bool is_end_keyword(char* str);
  bool is_keyword(char* str);
  bool is_number(char* str);
  bool is_procedure_keyword(char* str);
  bool is_separated_string();
  bool is_size_keyword(char* str);
  bool is_type_keyword(char* str);
  bool is_variable_keyword(char* str);
  bool next_is_end_keyword();
  bool next_is_group_names();
  bool next_is_keyword();
  bool next_is_matc();
  bool next_is_number();
  bool next_is_separated_string();
  bool next_is_size_keyword();
  bool next_is_text(char* text);
  bool next_is_type_keyword();
  bool next_is_variable_keyword();
  bool next_non_number_is_data_keyword();
  bool next_non_number_is_keyword();
  bool object_at_end();
  bool peek_next_string(char* buffer);
  bool peek_token(char* buffer);
  bool peek_field_name_tokens(int& nof_tokens);
  bool peek_object_name_tokens(int& nof_tokens);
  void pop_flags();
  void push_flags(bool copy_current = true);
  bool read_entry_data(struct emf_ObjectData_X* pf,
                       int data_count, int data_size, int& total_read_count);
  bool read_entry_data_values(int old_data_buffer_size, int max_data_buffer_size, int data_count,
                              int* data_sizes, int& total_read_count);
  bool read_entry_variables(struct emf_ObjectData_X* pf, int var_count, int var_size);
  bool read_entry_variable_values(int old_var_buffer_size, int max_var_buffer_size, int var_count,
                                 int* var_sizes, int& total_read_count);
  bool read_field(struct emf_ObjectData_X* pf);
  bool read_field_data_sized(struct emf_ObjectData_X* pf, int var_count, int var_size,
                             int data_count, int data_size, int& total_read_count);
  bool read_field_data_non_sized(struct emf_ObjectData_X* pf, bool type_given = false);
  bool read_field_name(struct emf_ObjectData_X*  od,
                       char* name_buffer, int buffer_len);
  bool read_group_names(struct emf_ObjectData_X* pf);
  bool read_logical_data(int start_pos, int max_pos, int& read_count);
  bool read_numeric_data(int start_pos, int max_pos, int& read_count);
  bool read_object(struct emf_ObjectData_X* pf);
  bool read_object_name(struct emf_ObjectData_X*  od,
                        char* name_buffer, int buffer_len);
  bool read_string_data(int start_pos, int max_pos, int string_count, 
                        int* string_sizes, int& read_count);
  void reset_flags();
  int sendMessage(char* msg);
  void set_transfer_info();
  void storeDataAsString(const char* buffer, bool is_matc = false);
  bool storeNumber(struct emf_ObjectData_X* od, LibFront::dataValueType data_type,
                   void* data);
  template <class T> bool try_to_read(T& data);
  bool update_data_lengths(int data_count, int* data_sizes);
  bool update_variable_lengths(int var_count, int* var_sizes);

  bool readerOk;
  void* currentBuffer;
  int currentFlagsIndex;
  enum LibFront::dataValueType currentBufferType;
  int currentBufferSize;
  char dataAsStringBuffer[MAX_BUFFER_LEN];
  int dataAsStringPos;
  void* dataBuffer;
  bool dataHasMatcVars;
  enum LibFront::dataValueType dataBufferType;
  int dataBufferSize;
  int* dataLengths;
  int dimension1;
  int dimension2;
  int evaluatedBufferEnd;
  LibFrontReaderFlags* flags;
  LibFrontReaderFlags flagsStack[32];
  double hangingNumber;
  ifstream inputFile;
  char* inputFileName;
  char matcBuffer[MAX_BUFFER_LEN];
  bool matcDefsWarned;
  char matcSeparator;
  int msgBufferSize;    
  char* msgBuffer;      // for readMessageCB function!
  char* nameBuffer;
  int nameBufferMaxLen;
  int nofDataValues;
  int nofEntries;
  int nofVariables;
  int nofVariableValues;
  char* readBuffer;
  int readBufferPos;
  int readBufferLen;
  char resultBuffer[MAX_BUFFER_LEN];
  char stringSeparator;
  char tokenBuffer[MAX_BUFFER_LEN];
  int tokenBufferLen;
  int tokenBufferMaxLen;
  void** userData;
  void* variableBuffer;
  int variableBufferSize;
  enum LibFront::dataValueType variableBufferType;
  int* variableLengths;
};
// =========================================================================
// END: File reader (parser( class
// =========================================================================

#endif
