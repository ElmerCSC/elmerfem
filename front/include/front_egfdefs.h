#ifndef _FRONT_EGF_DEFS_
#define _FRONT_EGF_DEFS_



// ==============================================
// Structure for user egf DLL-procedure arguments
//===============================================

// Typedefs for safer function calls
//
enum egf_FuncType {
  EGF_CIRCLE,
  EGF_POLYLINE
};


typedef double egf_Point3[3];

// Result circle-structure for a circle egf-function
//
struct egf_Circle {
  bool radiusGiven;
  bool centerGiven;
  bool startGiven;
  bool endGiven;
  double radius;
  egf_Point3 center;
  egf_Point3 start;
  egf_Point3 end;
  int vflags[2]; // vertex flags for start/end points

  egf_Circle() {
    init();
  }

  ~egf_Circle() {
  }

  void init() {
    radiusGiven = false;
    centerGiven = false;
    startGiven = false;
    endGiven = false;
    radius = 0.0;
    vflags[0] = -1;
    vflags[1] = -1;
  }
};

// Result curve-structure for a curve egf-function
//
struct egf_PolyLine {
  int nof_points;
  egf_Point3* points;
  int* vflags; // size = nof-points; 1 <--> a point should become a vertex

  egf_PolyLine() {
    init();
  }

  ~egf_PolyLine() {
    delete[] points;
    delete[] vflags;
  }

  void init() {
    nof_points = 0;
    points = 0;
    vflags = 0;
  }
};

/*

*--PolyLine DLL-proc call:

 f(int arg_vec_size, double* arg_vec,
   egf_Point3 start, egf_Point3 end,
   int& nof_polylines, egf_Curve*& polylines);

*--Circle DLL-proc call:

 f(int arg_vec_size, double* arg_vec,
   egf_Point3 start, egf_Point3 end,
   int& nof_circles, egf_Circle*& circles);

*/

#endif
