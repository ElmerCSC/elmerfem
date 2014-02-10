//
//  Read solution file
//


#include <mystdlib.h>


#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "nginterface.h"

namespace netgen
{
#include "writeuser.hpp"


void ImportSolution (const char * filename)
{
  ifstream inf (filename);
  char buf[100], name[1000];
  int i, size, comps, order;
  bool iscomplex;
  const char * type;
  Flags flags;

  while (1)
    {
      buf[0] = 0;
      inf >> buf;
      if (strcmp (buf, "solution") == 0)
	{
	  inf >> name;
	  
	  inf >> buf[0];
	  flags.DeleteFlags ();
	  while (buf[0] == '-')
	    {
	      inf >> buf[1];
	      inf.putback (buf[1]);
	      if (!isalpha (buf[1]))
		{
		  break;
		}
	      inf >> (buf+1);
	      flags.SetCommandLineFlag (buf);
	      buf[0] = 0;
	      inf >> buf[0];
	    }
	  inf.putback (buf[0]);

	  (*testout) << "Flags: " << endl;
	  flags.PrintFlags (*testout);
	  (*testout) << "done" << endl;

	  size = int(flags.GetNumFlag ("size", Ng_GetNP()));
	  comps = int(flags.GetNumFlag ("components", 1));
	  type = flags.GetStringFlag ("type", "nodal");
	  order = int(flags.GetNumFlag ("order", 1));
	  iscomplex = flags.GetDefineFlag ("complex");

	  double * sol = new double[size*comps];
	  
	  (*testout) << "import solution " << name << " size = " << size << " comps = " << comps << " order = " << order << endl;

	  for (i = 0; i < size*comps; i++)
	    {
	      inf >> sol[i];
	      //	      (*testout) << "sol: " << sol[i] << endl;
	    }
	  
	  Ng_SolutionData soldata;
	  Ng_InitSolutionData (&soldata);
	  soldata.name = name;
	  soldata.data = sol;
	  soldata.dist = comps;
	  soldata.components = comps;
	  soldata.order = order;
	  soldata.iscomplex = iscomplex;
	  soldata.soltype = NG_SOLUTION_NODAL;
          soldata.draw_surface = 1;
          soldata.draw_volume = 1;
	  if (strcmp (type, "element") == 0)
            {
              soldata.soltype = NG_SOLUTION_ELEMENT;
              soldata.draw_surface = 0;
            }
	  if (strcmp (type, "surfaceelement") == 0)
            {
              soldata.soltype = NG_SOLUTION_SURFACE_ELEMENT;
              soldata.draw_volume = 0;
            }
	  if (strcmp (type, "noncontinuous") == 0)
	    soldata.soltype = NG_SOLUTION_NONCONTINUOUS;
	  if (strcmp (type, "surfacenoncontinuous") == 0)
	    soldata.soltype = NG_SOLUTION_SURFACE_NONCONTINUOUS;

	  Ng_SetSolutionData (&soldata);
	  }
      else
	{
	  //	  cout << "kw = (" << buf << ")" << endl;
	  (*testout) << "kw = (" << buf << ")" << endl;
	  break;
	}
    }
  /*
  struct Ng_SolutionData
    {
      char * name;      // name of gridfunction
      double * data;    // solution values
      int components;   // used components in solution vector
      int dist;         // num of doubles per entry (alignment!)
      Ng_SolutionType soltype;  // type of solution function
  };

  // initialize solution data with default arguments
  void Ng_InitSolutionData (Ng_SolutionData * soldata);
  // set solution data
  void Ng_SetSolutionData (Ng_SolutionData * soldata);
  */
}



}
