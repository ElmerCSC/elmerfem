#ifndef FILE_FLAGS
#define FILE_FLAGS


/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

/** 
   Flag - Table.
   A flag table maintains string variables, numerical 
   variables and boolean flags.
*/
class Flags 
{
  ///
  SYMBOLTABLE<char *> strflags;
  ///
  SYMBOLTABLE<double> numflags;
  ///
  SYMBOLTABLE<int> defflags;
  ///
  SYMBOLTABLE<ARRAY<char*>*> strlistflags;
  ///
  SYMBOLTABLE<ARRAY<double>*> numlistflags;
public:
  ///
  Flags ();
  ///
  ~Flags ();
  
  /// Deletes all flags
  void DeleteFlags ();
  /// Sets string flag, overwrite if exists
  void SetFlag (const char * name, const char * val);
  /// Sets numerical flag, overwrite if exists
  void SetFlag (const char * name, double val);
  /// Sets boolean flag
  void SetFlag (const char * name);
  /// Sets string arary falg
  void SetFlag (const char * name, const ARRAY<char*> & val);
  /// Sets double array flag
  void SetFlag (const char * name, const ARRAY<double> & val);
  
  /// Save flags to file
  void SaveFlags (const char * filename) const;
  /// write flags to stream
  void PrintFlags (ostream & ost) const;
  /// Load flags from file
  void LoadFlags (const char * filename);
  /// set flag of form -name=hello -val=0.5 -defined
  void SetCommandLineFlag (const char * st);

  /// Returns string flag, default value if not exists
  const char * GetStringFlag (const char * name, const char * def) const;
  /// Returns numerical flag, default value if not exists
  double GetNumFlag (const char * name, double def) const;
  /// Returns address of numerical flag, null if not exists
  const double * GetNumFlagPtr (const char * name) const;
  /// Returns address of numerical flag, null if not exists
  double * GetNumFlagPtr (const char * name);
  /// Returns boolean flag
  bool GetDefineFlag (const char * name) const;
  /// Returns string list flag, empty array if not exist
  const ARRAY<char*> & GetStringListFlag (const char * name) const;
  /// Returns num list flag, empty array if not exist
  const ARRAY<double> & GetNumListFlag (const char * name) const;


  /// Test, if string flag is defined
  bool StringFlagDefined (const char * name) const;
  /// Test, if num flag is defined
  bool NumFlagDefined (const char * name) const;
  /// Test, if string list flag is defined
  bool StringListFlagDefined (const char * name) const;
  /// Test, if num list flag is defined
  bool NumListFlagDefined (const char * name) const;
};
  
#endif

