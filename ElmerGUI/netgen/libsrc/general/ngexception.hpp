#ifndef FILE_NGEXCEPTION
#define FILE_NGEXCEPTION

/**************************************************************************/
/* File:   ngexception.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 2002                                                  */
/**************************************************************************/


/// Base class for all ng exceptions
class NgException 
{
  /// verbal description of exception
  string what;
public:
  ///
  NgException (const string & s);
  ///
  virtual ~NgException ();

  /// append string to description
  void Append (const string & s);
  //  void Append (const char * s);
  
  /// verbal description of exception
  const string & What() const { return what; }
};

#endif
