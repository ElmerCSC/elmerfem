/**************************************************************************/
/* File:   flags.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

/* 
   Datatype Flags
*/

#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{
  //using namespace netgen;

  Flags :: Flags ()
  {
    ;
  }
  
  Flags :: ~Flags ()
  {
    DeleteFlags ();
  }
  
  void Flags :: DeleteFlags ()
  {
    for (int i = 0; i < strflags.Size(); i++)
      delete [] strflags[i];
    for (int i = 0; i < numlistflags.Size(); i++)
      delete numlistflags[i];
    strflags.DeleteAll();
    numflags.DeleteAll();
    defflags.DeleteAll();
    strlistflags.DeleteAll();
    numlistflags.DeleteAll();
  }
  
  void Flags :: SetFlag (const char * name, const char * val)
  {
    char * hval = new char[strlen (val) + 1];
    strcpy (hval, val);
    strflags.Set (name, hval);
  }
  
  void Flags :: SetFlag (const char * name, double val)
  {
    numflags.Set (name, val);
  }
  
  void Flags :: SetFlag (const char * name)
  {
    defflags.Set (name, 1);
  }


  void Flags :: SetFlag (const char * name, const ARRAY<char*> & val)
  {
    ARRAY<char*> * strarray = new ARRAY<char*>;
    for (int i = 1; i <= val.Size(); i++)
      {
	strarray->Append (new char[strlen(val.Get(i))+1]);
	strcpy (strarray->Last(), val.Get(i));
      }
    strlistflags.Set (name, strarray);
  }

  void Flags :: SetFlag (const char * name, const ARRAY<double> & val)
  {
    ARRAY<double> * numarray = new ARRAY<double>;
    for (int i = 1; i <= val.Size(); i++)
      numarray->Append (val.Get(i));
    numlistflags.Set (name, numarray);
  }




  
  const char * 
  Flags :: GetStringFlag (const char * name, const char * def) const
  {
    if (strflags.Used (name))
      return strflags.Get(name);
    else
      return def;
  }

  double Flags :: GetNumFlag (const char * name, double def) const
  {
    if (numflags.Used (name))
      return numflags.Get(name);
    else
      return def;
  }
  
  const double * Flags :: GetNumFlagPtr (const char * name) const
  {
    if (numflags.Used (name))
      return & ((SYMBOLTABLE<double>&)numflags).Elem(name);
    else
      return NULL;
  }
  
  double * Flags :: GetNumFlagPtr (const char * name) 
  {
    if (numflags.Used (name))
      return & ((SYMBOLTABLE<double>&)numflags).Elem(name);
    else
      return NULL;
  }
  
  bool Flags :: GetDefineFlag (const char * name) const
  {
    return defflags.Used (name);
  }


  const ARRAY<char*> & 
  Flags :: GetStringListFlag (const char * name) const
  {
    if (strlistflags.Used (name))
      return *strlistflags.Get(name);
    else
      {
	static ARRAY<char*> hstra(0);
	return hstra;
      }
  }

  const ARRAY<double> & 
  Flags ::GetNumListFlag (const char * name) const
  {
    if (numlistflags.Used (name))
      return *numlistflags.Get(name);
    else
      {
	static ARRAY<double> hnuma(0);
	return hnuma;
      }
  }


  bool Flags :: StringFlagDefined (const char * name) const
  {
    return strflags.Used (name);
  }

  bool Flags :: NumFlagDefined (const char * name) const
  {
    return numflags.Used (name);
  }

  bool Flags :: StringListFlagDefined (const char * name) const
  {
    return strlistflags.Used (name);
  }

  bool Flags :: NumListFlagDefined (const char * name) const
  {
    return numlistflags.Used (name);
  }


  void Flags :: SaveFlags (const char * filename) const 
  {
    int i;
    ofstream outfile (filename);
  
    for (i = 1; i <= strflags.Size(); i++)
      outfile << strflags.GetName(i) << " = " << strflags.Get(i) << endl;
    for (i = 1; i <= numflags.Size(); i++)
      outfile << numflags.GetName(i) << " = " << numflags.Get(i) << endl;
    for (i = 1; i <= defflags.Size(); i++)
      outfile << defflags.GetName(i) << endl;
  }
 


  void Flags :: PrintFlags (ostream & ost) const 
  {
    int i;
  
    for (i = 1; i <= strflags.Size(); i++)
      ost << strflags.GetName(i) << " = " << strflags.Get(i) << endl;
    for (i = 1; i <= numflags.Size(); i++)
      ost << numflags.GetName(i) << " = " << numflags.Get(i) << endl;
    for (i = 1; i <= defflags.Size(); i++)
      ost << defflags.GetName(i) << endl;
  }
 

  void Flags :: LoadFlags (const char * filename) 
  {
    char name[100], str[100];
    char ch;
    double val;
    ifstream infile(filename);

    //  (*logout) << "Load flags from " << filename << endl << endl;
    while (infile.good())
      {
	infile >> name;
	if (strlen (name) == 0) break;

	if (name[0] == '/' && name[1] == '/')
	  {
	    //	  (*logout) << "comment: ";
	    ch = 0;
	    while (ch != '\n' && infile.good())
	      {
		ch = infile.get();
		//	      (*logout) << ch;
	      }
	    continue;
	  }

	//      (*logout)  << name;
	ch = 0;
	infile >> ch;
	if (ch != '=')
	  {
	    //	  (*logout) << endl;
	    infile.putback (ch);
	    SetFlag (name);
	  }
	else
	  {
	    infile >> val;
	    if (!infile.good())
	      {
		infile.clear();
		infile >> str;
		SetFlag (name, str);
		//	      (*logout) << " = " << str << endl;
	      }
	    else
	      {
		SetFlag (name, val);
		//	      (*logout) << " = " << val << endl;
	      }
	  }
      }
    //  (*logout) << endl;
  }


  void Flags :: SetCommandLineFlag (const char * st)
  {
    //  cout << "clflag = " << st << endl;
    istringstream inst( (char *)st);
    // istrstream defined with char *  (not const char *  ?????)

    char name[100];
    double val;


    if (st[0] != '-')
      {
	cerr << "flag must start with '-'" << endl;
	return;
      }
  
    const char * pos = strchr (st, '=');
  
    if (!pos)
      {
	//      (cout) << "Add def flag: " << st+1 << endl;
	SetFlag (st+1);
      }
    else
      {
	//      cout << "pos = " << pos << endl;

	strncpy (name, st+1, (pos-st)-1);
	name[pos-st-1] = 0;

	//      cout << "name = " << name << endl;

	pos++;
	char * endptr = NULL;

	val = strtod (pos, &endptr);

	//      cout << "val = " << val << endl;

	if (endptr == pos)
	  {
	    //	  (cout) << "Add String Flag: " << name << " = " << pos << endl;
	    SetFlag (name, pos);
	  }
	else
	  {
	    //	  (cout) << "Add Num Flag: " << name << " = " << val << endl;
	    SetFlag (name, val);
	  }
      }


    /*
      inst >> name;
      (*mycout) << "name = " << name << endl;

      ch = 0;
      inst >> ch;
      if (ch != '=')
      {
      SetFlag (name);
      }
      else
      {
      inst >> val;
      if (!inst.good())
      {
      inst.clear();
      inst >> str;
      SetFlag (name, str);
      (*mycout) << "str = " << str << endl;
      }
      else
      {
      SetFlag (name, val);
      (*mycout) << "val = " << val << endl;
      }
      }
    */
  }
}
