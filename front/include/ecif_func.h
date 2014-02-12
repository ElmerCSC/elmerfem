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

/**********************************************************************
Program:    ELMER Front 
Module:     ecif_func.h
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   File includes declarations for common non-class C++-functions. 

************************************************************************/

#ifndef _ECIF_FUNC_
#define _ECIF_FUNC_

#include "ecif_def.h"

// ===============
// Stream routines
// ===============
extern ostream& indent(ostream& strm, int isize);
extern istrstream& reset(istrstream& strm);
extern ostrstream& reset(ostrstream& strm);
extern strstream& reset(strstream& strm);


// ====================
// File system routines
// ====================
bool front_chdir(const char* dir);
bool front_mkdir(const char* dir);
void front_getcwd(char* dir, int max_len);


// ================
// String functions
// ================
extern void create_dyna_string(char*& string_var, const char* new_value);
extern void update_dyna_string(char*& string_var, const char* new_value);
extern void formTimeString(double seconds, char* buffer);
extern void getCurrentTs(char* ts, int max_len);

int countNumbersInString(char* str);

// ====================
// Arithmetic functions
// ====================
int mod(int num, int denom);

// ====================
// Comparison functions
// ====================
int compare(double a, double b);
bool isEqual(double x, double y);
bool isEqualE(double x, double y, double epsilon = 0.0);
bool isGreater(double x, double y);
bool isLess(double x, double y);
bool isZero(double x);
bool isZeroE(double x, double epsilon = 0.0);

bool idLoopsAreMatching(int nof_ids, int nof_ids_to_check, 
                        const int* ids1, const int* ids2,
                        int& direction, int& start1, int& start2,
                        bool signed_match = false);

// ==================
// Geometry functions
// ==================
void add3(Point3 p1, Point3 p2, Point3& sum);
void copy3(Point3 source, Point3& target);
Point3* create3(const char* str);
void cross3(Point3 v1, Point3 v2, Point3& product);
void diff3(Point3 p1, Point3 p2, Point3& diff); // p1 - p2
double dir2_to_angle(Point3 p);
double dist3(Point3 p1, Point3 p2);
double dot3(Point3 p1, Point3 p2);
void initPoint3(Point3& p, double value = 0.0);
double length3(Point3 v);
void normalize(Point3& v);
void rotate2(Point3 vector, double angle, Point3& result);
void perpnorm(Point3 p1, Point3 p2, Point3 p3, Point3& norm); //Unit normal for tri(p1,p2,p3)
bool samepoint(Point3 p1, Point3 p2);
bool samepointE(Point3 p1, Point3 p2, double epsilon);
void scalarmult(double scalar, Point3 v, Point3& result);

// Center of a Ccw 2D-circle with radius and through two points (p1, p2)
bool circle_center(double radius, Point3 p1, Point3 p2, Point3& center);
// Center of a Ccw 2D-circle through three points (p1, p2, p3)
bool circle_center(Point3 p1, Point3 p2, Point3 p3, Point3& center);

bool isCcwTetra(double* a, double* b, double* c, double* d);
bool isCcwTriangle(double* a, double* b, double* c);
bool pointInsideRectangle(Point3 point, Point3** rec_points, Point3 rec_center, double rec_radius);
bool pointInsideTriangle(Point3 point, Point3** tri_points, Point3 tri_center, double tri_radius);


// =============
// DLL utilities
// =============
//typedef int (*Hfunc)();
bool loadDllFunction(const char* library_name, char* func_name, Hdll& hDLL, Hfunc& hFunc, char*& err_msg);
Hdll openDllLibrary(const char* library_name, char*& err_msg);
void closeDllLibrary(Hdll hDLL);


// ==============
// Misc utilities
// ==============
// Check if filename ends with extension-string *str*.
// Like "filename.unv" etc
bool checkFname(char *fname, char *str);

// Called automatically if new fails to allocate memory.
void free_store_exception();

// User 'enter' from command-line
void pressAnyKey();


#endif
