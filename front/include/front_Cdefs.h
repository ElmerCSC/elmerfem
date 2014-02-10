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
Module:     front_Cdefs.h
Language:   C
Date:       15.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:

Abstract:   A C-style header file for the interface for reading emf-files.
  All necessary data structures, functions and global variables needed
  for the interface are defined/declared here.

************************************************************************/

#ifndef _EMF_C_DEFS_
#define _EMF_C_DEFS_


/* Max length of one single condition when delivered as 'raw' string */
#define emf_MAX_STRING_LEN 8191

/* Max length of object/field names etc */
#define emf_MAX_NAME_LEN 64

/* Max number of groups (name groups separated by ';' like [A B; C] afet object/field names */
#define emf_MAX_NOF_GROUPS 2

/* Max number of names in one group */
#define emf_MAX_NOF_GROUP_NAMES 24

/* Max number of variables in data parameters (like time, temperature etc) */
#define emf_MAX_NOF_VARIABLES 12


/* **************************************************** */
/* This struct stores stats and control data for Emf2Db */
/* **************************************************** */
struct emf_Emf2DbData_X {
  int nof_bodies;
  int nof_bodyForces;
  int nof_bodyEquations;
  int nof_materials;
  int nof_boundaryConditions;
  int nof_initialConditions;
  int nof_meshParameters;
  char modeldir[emf_MAX_STRING_LEN];
};

/* *************************************** */
/* Callback function defined in emf2Db.lib */
/* *************************************** */
int Emf2DbCB(void**);

/* *************************************************** */
/* This struct stores the data for one parameter field */
/* *************************************************** */
struct emf_ObjectData_X
{
  /******************************/
  /* OBJECT LEVEL RELATED STUFF */
  /******************************/

  /* Flag for a new object */
  int is_object_start;

  /* Flag for the end of the object */
  int is_object_end;

  /* Object or section name like: Header, Initial Condition */
  char object_name[emf_MAX_NAME_LEN];

  /* Actual length of the object-name character array */
  int object_name_length;

  /* Possible id for the object, -1 if no id is used */
  int object_id;


  /*****************************/
  /* FIELD LEVEL RELATED STUFF */
  /*****************************/

  /* NOTE: Field name is optional, all the data can also be at object level! */

  /* Name for the parameter field: Equation, Velocity1 etc. */
  char field_name[emf_MAX_NAME_LEN];

  /* Actual length of the field_name character array */
  int field_name_length;


  /*****************************/
  /* GROUP LEVEL RELATED STUFF */
  /*****************************/
  /* Nof group names (like [H; N O] after object/field name */
  /* NOTE: character ';' is groups "grouper", currently we accept two different groups! */
  /* NOTE: if you have groups, you have to give same number of data values for each group! */
  int nof_group_names[emf_MAX_NOF_GROUPS];

  /* Group names [ H2O N; C] etc. */
  char group_names[emf_MAX_NOF_GROUPS][emf_MAX_NOF_GROUP_NAMES][emf_MAX_NAME_LEN];

  /* Actual lengths of the group names */
  int group_name_lengths[emf_MAX_NOF_GROUPS][emf_MAX_NOF_GROUP_NAMES];


  /*****************************/
  /* DATA VALUES RELATED STUFF */
  /*****************************/

  /* Name for the data type: Integer, Real, Logical, String, File etc. */
  char data_type[emf_MAX_NAME_LEN];

  /* Actual length of the data_type character array */
  int data_type_length;

  /* If definition is for a procedure */
  int is_procedure;

  /* Nof entries (rows of variables and data) */
  int nof_entries;

  /* Dimension1 for a data entry (default 1) */
  int dimension1;

  /* Dimension2 for a data entry (default 1) */
  int dimension2;

  /* Lengths for data values.
   * SIZE: nof_entries * dimension1 * dimension2
   *
   * NOTE: Normally each length is 1 (for a nuemric scalar item)
   *       But if data is of string type we normally have variable
   *       length strings and then the length of
   *       each string value can be read from this array
   *
   *       Actually you have to read string data using these
   *       lengths because there are no sperators between the strings
   *       in the "data" array !!!
  */
  int* data_lengths;

  /* Length of the data array = Sum(data_lengths), for convenience! */
  int data_length;

  /* Data array of homogenous data type */
  /* NOTE: cast to proprer type before using! */
  void* data;


  /*********************************/
  /* VARIABLES RELATED STUFF */
  /*********************************/

  /* Name for the varaible value type. NOTE: Real or String only! */
  char variable_type[emf_MAX_NAME_LEN];

  /* Length of the variable_type character array */
  int variable_type_length;

  /* Nof argument variables */
  int nof_variables;

  /* Argument variable names (Temperature, Times etc) */
  char variable_names[emf_MAX_NOF_VARIABLES][emf_MAX_NAME_LEN];

  /* Lengths of the variable-name arrays */
  int variable_name_lengths[emf_MAX_NOF_VARIABLES];

  /* Lengths for variable values.
   * SIZE = nof_entries * nofVariables
   *
   * NOTE: Normally these values are all 1, because a typical
   *       variable value is a (real) scalar. But we can also have
   *       names as variables and then the lengths of these character arries
   *       must be taken from this array (ref. data_lengths!)
   */
  int* variable_values_lengths;

  /* Total length of the variable values array = Sum(variable_lengths), for convenience! */
  int variable_values_length;

  /* Variable values array of homogenous data type (double or char only!) */
  /* NOTE: cast to proprer type before using! */
  void* variable_values;

};


/*** DECLARATION FOR LIBRARY INTERFACE ***/

#ifdef WIN32
  #ifdef BUILD_DLL
    #define EMF_DECL __declspec(dllimport)
  #else
    #define EMF_DECL extern
  #endif
#else
  #define EMF_DECL extern
#endif


/*Pointer to the object data transfer structure.
* Library function(s) use this pointer when transferin the data!
* Connect this pointer to your the variable which is allocated
* in your own module:
* declaration:  emf_ObjectData_X my_object_data;
* allocation:   emf_ObjectData = &my_object_data;
*/
EMF_DECL  struct emf_ObjectData_X*  emf_ObjectData;

/* Buffer for unknown data lines */
/* NOTE: This is dynamically allocated in the library, size: 1 + emf_MAX_STRING_LEN */
EMF_DECL  char* emf_UnknownData;


/*----- Interface FUNCTION decalrations */

/*This function starts parsing model-file
*
* Call it first in your module!!!
*
* Arguments:
* model_data_file: full path for the model file
* user_data      : this data is delivered back in the emf_readDataCB
* msg_buffer_size: maximum writeable size for the message buffer
* msg_buffer     : buffer for getting messages in the emf_messageReadCB
*
* Remember to allocate msg_buffer!
* If it is set NULL, call-back function is not used!

* Return values (int):
*   0 --> read ok
*   1 --> parsing was NOT ok, ERROR!!!
*/
EMF_DECL int emf_readData(char* model_data_file, void** user_data,
                          int msg_buffer_size, char* msg_buffer);


/*This call-back function does the actual parsing and
* returns data using data-structures defined above.
* Connect this function to the call-back function
* defined in your module.
*
* LIKE: int my_read_function (void) {...};
*       emf_readDataCB = my_read_function;
*
* Return values (int):
*   -1 --> stop reading
*  anything else --> continue
*/
EMF_DECL  int (*emf_readDataCB) (void** user_data);


/*This call-back function sends error messages etc.
* from the parser.
*
* Connect this function to the call-back function
* defined in your module.
*
* LIKE: void my_msg_function (char* message_buffer) {...};
*       emf_readMessageCB = my_msg_function;
*
* No return value
*
*
*/
EMF_DECL  void (*emf_readMessageCB) (char* msg_buffer);


#endif
