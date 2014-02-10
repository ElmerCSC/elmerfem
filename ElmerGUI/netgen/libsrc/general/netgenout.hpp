#ifndef NETGEN_OUT_STREAM_HPP__
#define NETGEN_OUT_STREAM_HPP__

// #include <ostream>
// #include <mystdlib.h>
// #include <meshing.hpp>

#ifdef PARALLEL
extern int id;
extern int ntasks;
#endif
extern int printmessage_importance;



class Imp
{
  int importance;
public:
  Imp () : importance(0) { ; }

  Imp ( int aimportance ) : importance(aimportance) { ; }

  int GetImp () const { return importance; }
};


class Proc
{
  int proc;
public:
  Proc () : proc(0) { ; }

  Proc ( int aproc ) : proc(aproc) { ; }

  int GetProc () const { return proc; }
};

class Procs
{
  const netgen::FlatArray<int> procs;

public:

  Procs ( const netgen::FlatArray<int> & aprocs ) : procs (aprocs) { ; }

  const netgen::FlatArray<int> & GetProcs () const { return procs; }
};



class NetgenOutStream
{
  ostream * out;

  bool print;
  bool printheader;


public:
  NetgenOutStream() :
    out(&std::cout),
    print(1),
    printheader(1)
  {
    ;
  }  

  NetgenOutStream(ostream * aout, Imp imp ) :
    out(aout),
    printheader(1)
  { 
    if ( netgen::printmessage_importance >= imp.GetImp() )
      print = true;
    else
      print = false;
  }

  NetgenOutStream(ostream * aout, Proc proc ) :
    out(aout),
    printheader(1)
  { 
#ifdef PARALLEL
    if ( netgen::id == proc.GetProc() )
      print = true;
    else
      print = false;
#else
    if ( 0 == proc.GetProc() )
      print = true;
    else
      print = false;

#endif
  }

  NetgenOutStream(ostream * aout, Procs & procs ) :
    out(aout),
    printheader(1)
  { 
#ifdef PARALLEL
    if ( procs.GetProcs().Contains(netgen::id) )
      print = true;
    else
      print = false;
#else
    if ( procs.GetProcs().Contains(0) )
      print = true;
    else
      print = false;

#endif
  }

  ostream & OStream ()
  {
    return *out;
  }

  template <typename T>
  NetgenOutStream & operator<< (T & var)
  {
    if ( print )
      {
#ifdef PARALLEL
	if ( printheader )
	  {
	    *out << "proc " << netgen::id << ": ";
	    printheader = false;
	  }
#endif
	*out << var;
      }
    return (*this); 
  }

  NetgenOutStream& operator<< (ostream& ( *pf )(ostream&))
  {
    if ( print )
      *out << (*pf) ;

    return (*this);
  }

  NetgenOutStream& operator<< (ios& ( *pf )(ios&))
  {
    if ( print)
      *out << (*pf) ;

    printheader = 1;

    return (*this);
  }

  NetgenOutStream& operator<< (ios_base& ( *pf )(ios_base&))
  {
    if (print )
      *out << (*pf) ;
    return (*this);
  }


};


NetgenOutStream operator<< ( ostream & ost, Imp  imp );
NetgenOutStream operator<< ( ostream & ost, Proc proc );
NetgenOutStream operator<< ( ostream & ost, Procs & procs );
// {
//   return ( NetgenOutStream ( &ost, imp.GetImp() ) );
// }

// template <typename T>
// NetgenOutStream& operator<< (NetgenOutStream& out, T c )
// {
//   out.OStream() << c << endl;
//   return out;
// }





#endif
