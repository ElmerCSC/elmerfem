/**************************************************************************/
/* File:   profiler.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/


#include <myadt.hpp>

namespace netgen
{
  //using namespace netgen;

  long int NgProfiler::tottimes[SIZE];
  long int NgProfiler::starttimes[SIZE];
  long int NgProfiler::counts[SIZE];
  string NgProfiler::names[SIZE];
  int NgProfiler::usedcounter[SIZE];
  

  NgProfiler :: NgProfiler()
  {
    for (int i = 0; i < SIZE; i++)
      {
	tottimes[i] = 0;
	usedcounter[i] = 0;
      }

    total_timer = CreateTimer ("total CPU time");
    StartTimer (total_timer);
  }

  NgProfiler :: ~NgProfiler()
  {
    StopTimer (total_timer);

    //ofstream prof;
    //prof.open("ng.prof");

    // ofstream-constructor may be called after STL-stuff is destructed,
    // which leads to an "order of destruction"-problem,
    // thus we use the C-variant:

    char filename[100];
#ifdef PARALLEL
    sprintf (filename, "netgen.prof.%d", id);
#else
    sprintf (filename, "netgen.prof");
#endif

    FILE *prof = fopen(filename,"w");
    Print (prof);
    fclose(prof);
  }


//   void NgProfiler :: Print (ostream & prof)
//   {
//     for (int i = 0; i < SIZE; i++)
//       if (counts[i] != 0 || usedcounter[i] != 0)
// 	{
// 	  prof.setf (ios::fixed, ios::floatfield);
// 	  prof.setf (ios::showpoint);

// 	  prof // << "job " << setw(3) << i 
// 	    << "calls " << setw(8) << counts[i] 
// 	    << ", time " << setprecision(2) << setw(6) << double(tottimes[i]) / CLOCKS_PER_SEC << " sec";

// 	  if (usedcounter[i]) 
// 	    prof << " " << names[i];
// 	  else
// 	    prof << " " << i;
	    
// 	  prof << endl;
// 	}
//   }


  void NgProfiler :: Print (FILE * prof)
  {
    for (int i = 0; i < SIZE; i++)
      if (counts[i] != 0 || usedcounter[i] != 0)
	{
	  //fprintf(prof,"job %3i calls %8i, time %6.2f sec",i,counts[i],double(tottimes[i]) / CLOCKS_PER_SEC);
	  fprintf(prof,"calls %8i, time %6.2f sec",counts[i],double(tottimes[i]) / CLOCKS_PER_SEC);
	  if(usedcounter[i])
	    fprintf(prof," %s",names[i].c_str());
	  else
	    fprintf(prof," %i",i);
	  fprintf(prof,"\n");
	}
  }

  int NgProfiler :: CreateTimer (const string & name)
  {
    for (int i = SIZE-1; i > 0; i--)
      if(names[i] == name)
	return i;

    for (int i = SIZE-1; i > 0; i--)
      if (!usedcounter[i])
	{
	  usedcounter[i] = 1;
	  names[i] = name;
	  return i;
	}
    return -1;
  }


  NgProfiler prof;
}
