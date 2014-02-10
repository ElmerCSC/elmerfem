#ifndef FILE_SYMBOLTA
#define FILE_SYMBOLTA


/**************************************************************************/
/* File:   symbolta.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/**
   Base class for the generic SYMBOLTABLE.
   An array of identifiers is maintained.
*/
class BASE_SYMBOLTABLE
{
protected:
  /// identifiers
  ARRAY <char*> names;
  
public:
  /// Constructor
  BASE_SYMBOLTABLE ();
  ///
  ~BASE_SYMBOLTABLE ();
  ///
  void DelNames ();
  /// Index of symbol name, returns 0 if not used.
  int Index (const char * name) const;
};


/** 
    Abstract data type Symbol Table.
   
    To a string an value of the generic type T is associated.
    The string is not copied into the symbol table class!
*/
template <class T>
class SYMBOLTABLE : public BASE_SYMBOLTABLE
{
private:
  /// Associated data
  ARRAY <T> data;
  
public:
  /// Creates a symboltable
  inline SYMBOLTABLE ();
  /// Returns size of symboltable
  inline INDEX Size() const;
  /// Returns reference to element, error if not used
  inline T & Elem (const char * name);
  /// Returns reference to i-th element
  inline T & Elem (int i) 
  { return data.Elem(i); }
  /// Returns element, error if not used
  inline const T & Get (const char * name) const;
  /// Returns i-th element
  inline const T & Get (int i) const;
  /// Returns name of i-th element
  inline const char* GetName (int i) const;
  /// Associates el to the string name, overrides if name is used
  inline void Set (const char * name, const T & el);
  /// Checks whether name is used
  inline bool Used (const char * name) const;
  /// Deletes symboltable
  inline void DeleteAll ();

  inline T & operator[] (int i) 
  { return data[i]; }
  inline const T & operator[] (int i) const
  { return data[i]; }

private:
  /// Prevents from copying symboltable by pointer assignment
  SYMBOLTABLE<T> & operator= (SYMBOLTABLE<T> &);
};




template <class T>
inline SYMBOLTABLE<T> :: SYMBOLTABLE () 
{ 
  ;
}


template <class T>
inline INDEX SYMBOLTABLE<T> :: Size() const
{
  return data.Size();
}

template <class T>
inline T & SYMBOLTABLE<T> :: Elem (const char * name)
{
  int i = Index (name);
  if (i) 
    return data.Elem (i);
  else 
    return data.Elem(1);
}

template <class T>
inline const T & SYMBOLTABLE<T> :: Get (const char * name) const
{
  int i;
  i = Index (name);
  if (i) 
    return data.Get(i);
  else 
    return data.Get(1);
}

template <class T>
inline const T & SYMBOLTABLE<T> :: Get (int i) const
{
  return data.Get(i);
}

template <class T>
inline const char* SYMBOLTABLE<T> :: GetName (int i) const
{
  return names.Get(i);
}

template <class T>
inline void SYMBOLTABLE<T> :: Set (const char * name, const T & el)
{
  int i;
  i = Index (name);
  if (i) 
    data.Set(i, el);
  else
    {
      data.Append (el);
      char * hname = new char [strlen (name) + 1];
      strcpy (hname, name);
      names.Append (hname);
    }
}

template <class T>
inline bool SYMBOLTABLE<T> :: Used (const char * name) const
{
  return (Index(name)) ? true : false;
}

template <class T>
inline void SYMBOLTABLE<T> :: DeleteAll () 
{	
  DelNames();
  data.DeleteAll();
}


#endif
