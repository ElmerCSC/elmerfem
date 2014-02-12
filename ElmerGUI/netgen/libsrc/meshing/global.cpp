#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{
  ostream * testout = &cout;

  ostream * mycout = &cout;
  ostream * myerr = &cerr;


  // Flags globflags; // not used anymoure
  Flags parameters;


  int silentflag = 0;
  int testmode = 0;

  MeshingParameters mparam;
  volatile multithreadt multithread;

  string ngdir = ".";

  ARRAY<int> tets_in_qualclass;


  multithreadt :: multithreadt()
  {
    pause =0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }

  DebugParameters debugparam;
  bool verbose = 0;

  int timestamp = 0;
  int GetTimeStamp() 
  { 
    return timestamp; 
  }

  int NextTimeStamp()
  {
    timestamp++;
    return timestamp;
  }
}
