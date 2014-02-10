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
Program:    ELMER Front 
Module:     ecif_func.cpp
Language:   C++
Date:       01.10.98
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Implementation

************************************************************************/
 
#include "ecif_func.h"

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(WIN32)   
  # include <direct.h> 
  # include <io.h>     
#else                
  #include <unistd.h>
  #include <dlfcn.h>
#endif



// ===============
// Stream routines
// ===============

ostream& indent(ostream& strm, int isize)
{
  for (int i = 0; i < isize; i++)
    strm << ' ';
 
  return strm;
}

istrstream& reset(istrstream& strm)
{
  strm.clear();
  strm.seekg(0);
  return strm;
}


ostrstream& reset(ostrstream& strm)
{
  strm.clear();
  strm.seekp(0);
  return strm;
}


strstream& reset(strstream& strm)
{
  strm.clear();
  strm.seekg(0);
  strm.seekp(0);
  strm.rdbuf()->freeze(false);
  return strm;
}


// ====================
// File system routines
// ====================

bool front_mkdir(const char* dir)
{
  if ( dir == NULL ) {
    cerr << "mkdir: No valid diretory given!" << endl;
    return false;
  }

  int rc;

#if defined(WIN32)
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
    case ENOENT:
	  default:
	    cerr << "mkdir: cannot create directory: " << dir << endl;
	    break;
	  }

    return false;
  }

  return true;
}



bool front_chdir(const char* dir)
{
  int rc;

  if ( dir == NULL ) {
    cerr << "chdir: No valid diretory given!" << endl;
    return false;
  }

#if defined(WIN32)
  rc = _chdir(dir);
#else
  rc = chdir(dir);
#endif

  if(rc == -1) {
    switch(errno) {
	    case EACCES:
	      cerr << "chdir: No permissions to:" << dir << endl;
	      break;
	    case EIO:
	      cerr << "chdir: I/O error with: " << dir << endl;
	      break;	  
	    case ENOENT:
	      cerr << "chdir: No such directory as: " << dir << endl;
	      break;
	    case ENOTDIR:
	      cerr << "chdir: Not a directory:" << dir << endl;
	      break;
	    default:
	      cerr << "chdir: Unexpected error with dir: " << dir << endl;
	      break;
	  }

    return false;
  }

  return true;
}


void front_getcwd(char* dir, int max_len)
{
#if defined(WIN32)
  _getcwd(dir, max_len);
#else
  getcwd(dir, max_len);
#endif
}


// Check if filename ends with str.
bool 
checkFname(char* fname, char* str)
{
  // Pick extension in the file name
  // (including the dot)
  char* ext = strrchr(fname, '.');

  if ( ext == NULL )
    return false;
  
  int i;
  int len1 = strlen(ext);
  int len2 = strlen(str);

  // Definitely different extensions!
  if (len1 != len2)
    return false;

  // Compare strings
  for (i = 0; i < len1; i++) {

    // Non case sensitive comparison
    if ( toupper(ext[i]) != toupper(str[i]) ) {

    // Case sensitive comparison
    //if ( ext[i] != str[i] ) {

      return false;
    }
  }

  // Ok, same they are!
  return true;
}  


// ===============
// String routines
// ===============

// Function calculates number (integer or real) count
// in the argument string
// NOTE: First non-number stops counting!
//
int countNumbersInString(char* str)
{
  strstream strm;

  double dummy;
  int count = 0;
  strm << str;

  while (!strm.eof()) {

    if ( strm >> dummy ) {
      count++;
    } else {
      break;
    }
  }

  return count;
}


// Function creates a time string from 'seconds' variable, 
// result in in the format HH:MM:SS
void formTimeString(double seconds, char* buffer)
{
  int hrs = int( seconds / 3600.0 );
  int mns = int( (seconds - hrs * 3600.0) / 60.0 );
  int scs = int( (seconds - hrs * 3600.0 - mns * 60.0) );
    
  ostrstream ostrm;
  ostrm.setf(ios::right, ios::adjustfield);
  ostrm << setw(4) << hrs << ":";
  ostrm << setfill('0');
  ostrm << setw(2) << mns << ":" << setw(2) << scs;
  int len = ostrm.pcount();
  ostrm << ends;

  for (int i=0; i< len; i++)
    buffer[i] = (ostrm.str())[i];

  buffer[len] = '\0';
}


// Current time as a string
void getCurrentTs(char* ts, int max_len)
{
  time_t ltime;
  time( &ltime );

  struct tm* today = localtime( &ltime );

  strftime( ts, max_len, "%a %b %d %H:%M:%S %Y", today );
}


// Create dynamically (with new) allocated string
// Allocate new string_var variable based on the size of the new_value
// Copy new_value into string_val
// NOTE: string_var is never NULL!, so it can always be out-stream argument
void create_dyna_string(char*& string_var, const char* new_value)
{
  string_var = NULL;
  update_dyna_string(string_var, new_value);
}
  

// Update dynamically (with new) allocated string
// Delete old string_var variable
// Allocate new string_var variable based on the size of the new_value
// Copy new_value into string_val
//
// NOTE: string_var should never be NULL!, so that it can be always as an out-stream argument
//
void update_dyna_string(char*& string_var, const char* new_value)
{
  
  delete[] string_var;
  string_var = NULL;

  // Calc new length
  int len = 0;

  if (new_value != NULL) {
    len = strlen(new_value);
  }

  // Allocate (always != NULL!)
  string_var = new char[1 + len];
  string_var[len] = '\0';

  // Copy new-value if not null
  if (new_value != NULL) {
    strcpy(string_var, new_value);
  }
}


//==========================
//  STL comparison operators
//==========================

// For struct Ids1
bool
operator<(const Ids1& p1, const Ids1& p2) {
  if (p1.id1 < p2.id1) 
     return true;
  else
    return false;
}


// For struct Ids1
bool
operator==(const Ids1& p1, const Ids1& p2) {
  if (p1.id1 == p2.id1) 
     return true;
  else
    return false;
}


// For struct Ids2
bool
operator<(const Ids2& p1, const Ids2& p2) {
  if (p1.id1 < p2.id1 || 
     p1.id1 == p2.id1 && p1.id2 < p2.id2)
     return true;
  else
    return false;
}


// For struct Ids2
bool
operator==(const Ids2& p1, const Ids2& p2) {
  if (p1.id1 == p2.id1 && p1.id2 == p2.id2)
     return true;
  else
    return false;
}

// For struct Ids3
bool
operator<(const Ids3& p1, const Ids3& p2) {
  if (p1.id1 < p2.id1 ||
     p1.id1 == p2.id1 && p1.id2 < p2.id2 ||
     p1.id1 == p2.id1 && p1.id2 == p2.id2 && p1.id3 < p2.id3)
    return true;
  else
    return false;
}


bool
operator==(const Ids3& p1, const Ids3& p2) {
  if (p1.id1 == p2.id1 &&
      p1.id2 == p2.id2 &&
      p1.id3 == p2.id3)
    return true;
  else
    return false;
}



// =============================
// Miscellanous helper functions
// =============================

// Called automatically if new fails to allocate memory.
void
free_store_exception()
{
  cerr << ": free store exhausted!" << endl;
  exit(1);
}


// User 'enter'
void
pressAnyKey()
{
  cerr << "Press any key ...";
  getchar();
  cerr << endl;
}


// =====================
// Dynamic library stuff
// =====================
 
bool
loadDllFunction(const char* library_name, char* func_name, Hdll& hDLL, Hfunc& hFunc, char*& err_msg)
{
  strstream err_strm;

  hDLL = NULL;
  hFunc = NULL;
  err_msg = NULL;

  // Open library
  // ------------
  hDLL = openDllLibrary(library_name, err_msg);

  // ERROR: library not found!
  //
  if ( hDLL == NULL ) {
    return false;
  }

  // Get function handle
  // -------------------
//--Windows
#if defined(WIN32)
  // Try first decorated name (_name@n)
  strstream strm;
  strm << "_" << func_name << "@28" << ends;
  hFunc = (Hfunc)GetProcAddress(hDLL, strm.str());

  // Next try undecorated name
  if ( hFunc == NULL ) {
    hFunc = (Hfunc)GetProcAddress(hDLL, func_name);
  }

//--Unix
#else
 hFunc= dlsym(hDLL,func_name);
#endif

  // ERROR: function not found!
  //
  if ( hFunc == NULL ) {
    strstream strm;
    strm << "Cannot open the library function: " << func_name << ends;;
    update_dyna_string(err_msg, strm.str());
    return false;
  }

  return true;
}


// Open dynamic library
//
Hdll openDllLibrary(const char* library_name, char*& err_msg)
{
  strstream err_strm;
  Hdll hDLL = NULL;

//--Windows
#if defined(WIN32)
  hDLL = LoadLibrary(library_name);

  // Try also with .dll extension
  if (hDLL == NULL ) {
    strstream strm;
    strm << library_name << ".dll" << ends;
    hDLL = LoadLibrary(strm.str());
  }

  // Error, library not found
  if (hDLL == NULL) {
    err_strm << "Cannot open dll-library: " << library_name;
  }

//--Unix
#else
   hDLL = dlopen(library_name,RTLD_NOW);

  // Error, library not found
  if (hDLL == NULL) {
    err_strm << "Cannot open dynamic library: " << library_name << " because  " << dlerror();
  }
#endif

  if ( hDLL == NULL ) {
    err_strm << ends;
    update_dyna_string(err_msg, err_strm.str());
  }

  return hDLL;
}


// Close dynamic library
//
void closeDllLibrary(Hdll hDLL)
{
//--Windows
#if defined(WIN32)
  FreeLibrary(hDLL);

//--Unix
#else
  dlclose(hDLL);
#endif
}


// ====================

// Function compares two double-values within EPSILON-limit
// Should be used in qsort etc.
//
int 
compare(double a, double b)
{
  double diff = a - b;

  /* Values are equal */
  if ( isZero(diff) )
    return 0;

  /* a > b */
  else if (diff > 0)
    return 1;   

 /* b > a */
  else
    return -1;
} 


bool
isEqual(double x, double y)
{
  if ( isZero(x-y) ) 
    return true;
  else
    return false;
}


bool
isEqualE(double x, double y, double epsilon)
{
  if ( epsilon == 0.0 ) {
    epsilon = EPSILON;
  }

  if ( isZeroE(x-y, epsilon) ) 
    return true;
  else
    return false;
}

bool 
isGreater(double x, double y)
{
  if ( x > (y + EPSILON) ) 
    return true;
  else
    return false;
}


bool 
isLess(double x, double y)
{
  if ( x < (y - EPSILON) ) 
    return true;
  else
    return false;
}

bool
isZero(double x)
{
  if ( x < EPSILON && x > -EPSILON ) 
    return true;
  else
    return false;
}


bool
isZeroE(double x, double epsilon)
{
  if ( epsilon == 0.0 ) {
    epsilon = EPSILON ;
  }

  if ( x < epsilon && x > -epsilon ) 
    return true;
  else
    return false;
}


int
mod(int num, int denom)
{  
  int quot, rem;    

  quot = num / denom;
  rem = num - denom * quot;
  return rem;
}



// Check if two ids loops are matching and find starting position
// in loop one and relative direction
//
bool idLoopsAreMatching(int nof_ids, int nof_ids_to_check, 
                        const int* ids1, const int* ids2,
                        int& direction, int& start1, int& start2,
                        bool signed_match)
{
  int i, pos, pos1;
  bool next_found;

  start1 = -1;
  start2 = 0;
  direction = 1;

  // Find possible start position in ids1
  for (i = 0; i < nof_ids; i++) {
    if ( signed_match && (ids1[i] == ids2[0]) ||
         !signed_match && (abs(ids1[i]) == abs(ids2[0]))
        ) {
      start1 = i;
    }
  }

  if (start1 == -1 ) return false;

  // Start found and only one element to check, so OK!
  //
  if ( nof_ids_to_check == 1 ) return true;

  // Peek next positions and find direction
  pos = start1;
  next_found = false;

  //--Try positive direction
  if (!next_found) {

    pos1 = pos + 1;

    if ( pos1 == nof_ids ) {
      pos1 = 0;
    }

    if ( signed_match && (ids1[pos1] == ids2[1]) ||
         !signed_match && (abs(ids1[pos1]) == abs(ids2[1]))
       ) {
      direction = 1;
      next_found = true;
    }
  }

  //--Try negative direction
  if (!next_found) {

    pos1 = pos - 1;

    if ( pos1 < 0 ) {
      pos1 = nof_ids - 1;
    }

    if ( signed_match && (ids1[pos1] == ids2[1]) ||
         !signed_match && (abs(ids1[pos1]) == abs(ids2[1]))
       ) {
      direction = -1;
      next_found = true;
    }
  }

  //--No match in either direction
  if (!next_found) return false;

  // Continue checking the match
  pos = pos1 + direction;

  // Loop enough ids to check the match
  for (i = 2; i < nof_ids_to_check; i++) {
    
    //-Check also table overflows
    // "dropping" from left
    if (pos < 0) {
      pos = nof_ids - 1;

    // "dropping" from right
    } else if (pos == nof_ids) {
      pos = 0;
    }

    //-If ids are different, we can stop
    if ( signed_match && (ids1[pos] != ids2[i]) ||
         !signed_match && (abs(ids1[pos]) != abs(ids2[i]))
       ) {
      return false;
    }

    //-Next check position
    pos += direction;
  }

  return true;
}


// ==================
// Geometry functions
// ==================

Point3* create3(const char* str)
{
  Point3* p = new Point3[1];

  strstream strm;
  strm << str << ends;

  for (int i = 0; i < 3; i++) {
    (*p)[i] = 0.0;
    strm >> (*p)[i];
  }

  return p;
}


// Point stuff
// -----------

void add3(Point3 p1, Point3 p2, Point3& sum)
{
  sum[0] = p1[0] + p2[0];
  sum[1] = p1[1] + p2[1];
  sum[2] = p1[2] + p2[2];
}


void cross3(Point3 v1, Point3 v2, Point3& product)
{
  Point3 p;	/* if product is same as v1 or v2 */

  p[0] = v1[1]*v2[2] - v2[1]*v1[2];
  p[1] = v1[2]*v2[0] - v2[2]*v1[0];
  p[2] = v1[0]*v2[1] - v2[0]*v1[1];
  product[0] = p[0]; product[1] = p[1]; product[2] = p[2];
}


// p1 - p2
void diff3(Point3 p1, Point3 p2, Point3& diff)
{
  diff[0] = p1[0] - p2[0];
  diff[1] = p1[1] - p2[1];
  diff[2] = p1[2] - p2[2];
}


// Direction (normalized vector) as radians (0...TWO_PI)
double dir2_to_angle(Point3 p)
{
  double base;
  int sign;

  // QUAD-1
  if ( p[0] >= 0.0 && p[1] >= 0.0 ) {
    base = 0.0;
    sign = 1;;
  // QUAD-2
  } else if ( p[0] < 0.0 && p[1] >= 0.0 ) {
    base = PI;
    sign = -1;
  // QUAD-3
  } else if ( p[0] < 0.0 && p[1] < 0.0 ) {
    base = PI;
    sign = 1;
  // QUAD-4
  } else {
    base = TWO_PI;
    sign = -1;
  }

  double prad = acos(fabs(p[0]));
  return base + sign * prad;
}


double dist3(Point3 p1, Point3 p2)
{
  Point3 d;
  diff3(p1, p2, d);
  return length3(d);
}


double dot3(Point3 p1, Point3 p2)
{
  return p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
}


void initPoint3(Point3& p, double value)
{
  p[0] = p[1] = p[2] = value;
}


double length3(Point3 v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


void copy3(Point3 source, Point3& target)
{
  target[0] = source[0], target[1] = source[1], target[2] = source[2];
}


void normalize(Point3& v)
{
  double d;

  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d == 0.0) {
	  v[0] = d = 1.0;
  }
  d = 1/d;
  v[0] *= d; v[1] *= d; v[2] *= d;
}


// Unit normal for triangle(p1,p2,p3)
void perpnorm(Point3 p1, Point3 p2, Point3 p3, Point3& norm)
{
  Point3 d1, d2;
  diff3(p2, p1, d1);
  diff3(p2, p3, d2);
  cross3(d1, d2, norm);
  normalize(norm);
}


// 2D rotation
void rotate2(Point3 vector, double angle, Point3& result)
{
  double sin_a = sin(angle);
  double cos_a = cos(angle);

  result[0] = vector[0] * cos_a - vector[1] * sin_a;
  result[1] = vector[0] * sin_a + vector[1] * cos_a;
  result[2] = vector[2];
}


bool samepoint(Point3 p1, Point3 p2)
{
  if (isZero(p1[0] - p2[0]) && isZero(p1[1] - p2[1]) && isZero(p1[2] - p2[2]))
	  return true;

  return false;
}


bool samepointE(Point3 p1, Point3 p2, double epsilon)
{
  if ( isZeroE(p1[0] - p2[0], epsilon) &&
       isZeroE(p1[1] - p2[1], epsilon) &&
       isZeroE(p1[2] - p2[2], epsilon)
     )
	  return true;
  else
    return false;
}


void scalarmult(double scalar, Point3 v, Point3& result)
{
  result[0] = scalar * v[0];
  result[1] = scalar * v[1];
  result[2] = scalar * v[2];
}


// Circle stuff
// ------------

// Center of a Ccw 2D-circle with radius and through two points (p1, p2)
bool circle_center(double radius, Point3 p1, Point3 p2, Point3& center)
{
  Point3 c, p3 , normal, to_cc;

  // Segment center
  diff3(p2, p1, c);
  scalarmult(0.5, c, c);

  // Normal to segment
  normal[0] = -c[1];
  normal[1] = c[0];
  normal[2] = c[2];

  // From segment center towards circle center
  cross3(normal, c, to_cc);
  normalize(to_cc);

  double half_seg = length3(c);

  double len = sqrt(radius * radius - half_seg * half_seg);

  scalarmult(len, to_cc, to_cc);

  add3(c, to_cc, center);

  return true;
}

// Center of a Ccw 2D-circle through three points (p1, p2, p3)
bool circle_center(Point3 p1, Point3 p2, Point3 p3, Point3& center)
{
  return false;
}



// Triangle stuff
// --------------

bool
pointInsideTriangle(Point3 point, Point3** tri_points, Point3 tri_center, double tri_radius)
{
  Point3 tmp, tmp1, tmp2;

  Point3* tp0 = tri_points[0];
  Point3* tp1 = tri_points[1];
  Point3* tp2 = tri_points[2];

 /* -First trivial distance check
  * Compare point distance from the center of the triangle
  * with the "radius" of the triangle
  */
  diff3(tri_center, point, tmp);

  if ( isGreater( dot3(tmp, tmp), tri_radius) ) {
    return false;
  }

 /*-Check if it remains between edges looked from
  * vertices
  * First look from p0
  */
  diff3(point, *tp0, tmp);
  diff3(*tp1, *tp0, tmp1);
  diff3(*tp2, *tp0, tmp2);
  normalize(tmp);
  normalize(tmp1);
  normalize(tmp2);

  /* if too much "left" */
  if ( dot3(tmp1, tmp2) > dot3(tmp1, tmp) ) {
    return false;
  }

  /* if too much "right" */
  if ( dot3(tmp2, tmp1) > dot3(tmp2, tmp) ) {
    return false;
  }

  /* Then look from p2 */
  diff3(point, *tp2, tmp);
  diff3(*tp1, *tp2, tmp1);
  diff3(*tp0, *tp2, tmp2);
  normalize(tmp);
  normalize(tmp1);
  normalize(tmp2);

  /* if too much "left" */
  if ( dot3(tmp2, tmp1) > dot3(tmp2, tmp) ) {
    return false;
  }

  /* if too much "right" */
  if ( dot3(tmp1, tmp2) > dot3(tmp1, tmp) ) {
    return false;
  }

  return true;
}


/* Checks if triangle is ccw-ordered */
bool isCcwTriangle(double* a, double* b, double* c)
{
 /*
	* Calculate triangle area from determinant (only x,y in use!)
	* O'Rourke, Computational geometry, p. 26)
 	* |a0 a1 1|
	* |b0 b1 1|
	* |c0 c1 1|
  */

	double area = 0.0 
		+a[0]*b[1] -a[1]*b[0] 
		+a[1]*c[0] -a[0]*c[1] 
		+b[0]*c[1] -c[0]*b[1]; 

 /*
  * Final orientation is concluded from sign of the are. 
  */
	int is_ccw = (area > 0.0);

	return is_ccw;
}



// Rectangle stuff
// ---------------

bool
pointInsideRectangle(Point3 point, Point3** rec_points, Point3 rec_center, double rec_radius)
{
  Point3 tmp, tmp1, tmp2;

  Point3* rp0 = rec_points[0];
  Point3* rp1 = rec_points[1];
  Point3* rp2 = rec_points[2];
  Point3* rp3 = rec_points[3];

 /* -First trivial distance check 
  * Compare point distance from the center of the rectangle
  * with the "radius" of the triangle
  */
  diff3(rec_center, point, tmp);

  if ( isGreater( dot3(tmp, tmp), rec_radius) ) {
    return false;
  }

 /*-Check if it remains between edges looked from
  * vertices

  *---First look from p0
  */
  diff3(point, *rp0, tmp);
  diff3(*rp1, *rp0, tmp1);
  diff3(*rp3, *rp0, tmp2);
  normalize(tmp);
  normalize(tmp1);
  normalize(tmp2);

  /* if too much "left" */
  if ( dot3(tmp1, tmp2) > dot3(tmp1, tmp) ) {
    return false;
  }

  /* if too much "right" */
  if ( dot3(tmp2, tmp1) > dot3(tmp2, tmp) ) {
    return false;
  }

  /*---Then look from p1 */
  diff3(point, *rp1, tmp);
  diff3(*rp2, *rp1, tmp1);
  diff3(*rp0, *rp1, tmp2);
  normalize(tmp);
  normalize(tmp1);
  normalize(tmp2);

  /* if too much "left" */
  if ( dot3(tmp2, tmp1) > dot3(tmp2, tmp) ) {
    return false;
  }

  /* if too much "right" */
  if ( dot3(tmp1, tmp2) > dot3(tmp1, tmp) ) {
    return false;
  }

  /*---Then look from p2 */
  diff3(point, *rp2, tmp);
  diff3(*rp1, *rp2, tmp1);
  diff3(*rp3, *rp2, tmp2);
  normalize(tmp);
  normalize(tmp1);
  normalize(tmp2);

  /* if too much "left" */
  if ( dot3(tmp2, tmp1) > dot3(tmp2, tmp) ) {
    return false;
  }

  /* if too much "right" */
  if ( dot3(tmp1, tmp2) > dot3(tmp1, tmp) ) {
    return false;
  }

  return true;
}


// Tetra stuff
// -----------

/* Checks if tetra is ccw-ordered */
bool isCcwTetra(double* a, double* b, double* c, double* d)
{

  /* Calculate tetrahedron volume from determinant
	 * O'Rourke, Computational geometry, p. 26)
	 * |a0 a1 a2 1|
	 * |b0 b1 b2 1|
	 * |c0 c1 c2 1|
	 * |d0 d1 d2 1|
   */

	double vol = 0.0
		-b[0]*c[1]*d[2] +a[0]*c[1]*d[2] +b[1]*c[0]*d[2] -a[1]*c[0]*d[2] -a[0]*b[1]*d[2] 
		+a[1]*b[0]*d[2] +b[0]*c[2]*d[1] -a[0]*c[2]*d[1] -b[2]*c[0]*d[1] +a[2]*c[0]*d[1] 
		+a[0]*b[2]*d[1] -a[2]*b[0]*d[1] -b[1]*c[2]*d[0] +a[1]*c[2]*d[0] +b[2]*c[1]*d[0] 
		-a[2]*c[1]*d[0] -a[1]*b[2]*d[0] +a[2]*b[1]*d[0] +a[0]*b[1]*c[2] -a[1]*b[0]*c[2] 
		-a[0]*b[2]*c[1] +a[2]*b[0]*c[1] +a[1]*b[2]*c[0] -a[2]*b[1]*c[0];

 /*
  * Final orientation is concluded from sign of the volume. 
  */
	int is_ccw = (vol > 0.0);

	return is_ccw;
}

