#ifndef FILE_HASHTABL
#define FILE_HASHTABL

/**************************************************************************/
/* File:   hashtabl.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/**
   Abstract data type HASHTABLE.
   Hash is done by one INDEX
*/
class BASE_INDEX_HASHTABLE
{
protected:
  /// keys are stored in this table
  TABLE<INDEX> hash;
  
public:
  ///
  BASE_INDEX_HASHTABLE (int size)
    : hash (size) { };
  
protected:
  ///
  int HashValue (const INDEX & ind) const
    {
      return ind % hash.Size() + 1;
    }

  ///
  int Position (int bnr, const INDEX & ind) const;
};

///
template <class T>
class INDEX_HASHTABLE : private BASE_INDEX_HASHTABLE
{
  ///
  TABLE<T> cont;
  
public:
  ///
  inline INDEX_HASHTABLE (int size);
  ///
  inline void Set (const INDEX & hash, const T & acont);
  ///
  inline const T & Get (const INDEX & ahash) const;
  ///
  inline bool Used (const INDEX & ahash) const;
  ///
  inline int GetNBags () const;
  ///
  inline int GetBagSize (int bnr) const;
  ///
  inline void GetData (int bnr, int colnr, INDEX & ahash, T & acont) const;

  ///
  inline void PrintMemInfo (ostream & ost) const;
};











///
class BASE_INDEX_2_HASHTABLE
{
protected:
  ///
  TABLE<INDEX_2> hash;
  
public:
  ///
  BASE_INDEX_2_HASHTABLE (int size)
    : hash (size) { };

  ///
  void PrintStat (ostream & ost) const;
  void BaseSetSize(int s) {hash.SetSize(s);}
protected:
  ///
  int HashValue (const INDEX_2 & ind) const
    {
      return (ind.I1() + ind.I2()) % hash.Size() + 1;
    }
  ///
  int Position (int bnr, const INDEX_2 & ind) const
  {
    int i;
    for (i = 1; i <= hash.EntrySize (bnr); i++)
      if (hash.Get(bnr, i) == ind)
	return i;
    return 0;
  }
};


///
template <class T>
class INDEX_2_HASHTABLE : public BASE_INDEX_2_HASHTABLE
{
  ///
  TABLE<T> cont;
  
public:
  ///
 INDEX_2_HASHTABLE (int size)
   : BASE_INDEX_2_HASHTABLE (size), cont(size)
  { ; }  

  ///
  void SetSize(int s) 
  { 
    cont.SetSize(s); 
    BaseSetSize(s);
  }

  ///
  void Set (const INDEX_2 & ahash, const T & acont)
  {
    int bnr = HashValue (ahash);
      int pos = Position (bnr, ahash);
      if (pos)
	cont.Set (bnr, pos, acont);
      else
	{
	  hash.Add1 (bnr, ahash);
	  cont.Add1 (bnr, acont);
	}    
  }
  
  ///
  const T & Get (const INDEX_2 & ahash) const
  {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    return cont.Get (bnr, pos);
  }
  
  ///
  bool Used (const INDEX_2 & ahash) const
  {
    return (Position (HashValue (ahash), ahash)) ? 1 : 0;
  }
 ///
  int GetNBags () const
  {
    return cont.Size();
  }
  
  ///
  int GetBagSize (int bnr) const
  {
    return cont.EntrySize (bnr);
  }
    
  ///
  void GetData (int bnr, int colnr, 
		INDEX_2 & ahash, T & acont) const
  {
    ahash = hash.Get(bnr, colnr);
    acont = cont.Get(bnr, colnr);
  }

  ///
  void SetData (int bnr, int colnr, 
		const INDEX_2 & ahash, const T & acont) 
  {
    hash.Set(bnr, colnr, ahash);
    cont.Set(bnr, colnr, acont);
  }
  
  ///
  void PrintMemInfo (ostream & ost) const
  {
    ost << "Hash: " << endl;
    hash.PrintMemInfo (ost);
    ost << "Cont: " << endl;
    cont.PrintMemInfo (ost);
  }


  void DeleteData ()
  {
    int n = hash.Size();
    hash.SetSize (n);
    cont.SetSize (n);
  }


  class Iterator
  {
    const INDEX_2_HASHTABLE & ht;    
    int bagnr, pos;
  public:
    Iterator (const INDEX_2_HASHTABLE & aht,
	      int abagnr, int apos)
      : ht(aht), bagnr(abagnr), pos(apos)
    { ; }

    int BagNr() const { return bagnr; }
    int Pos() const { return pos; }

    void operator++ (int)
    {
      // cout << "begin Operator ++: bagnr = " << bagnr << " -  pos = " << pos << endl;
      pos++;
      while (bagnr < ht.GetNBags() && 
	     pos == ht.GetBagSize(bagnr+1))
	{
	  pos = 0;
	  bagnr++;
	}
      // cout << "end Operator ++: bagnr = " << bagnr << " - pos = " << pos << endl;
    }

    bool operator != (int i) const
    {
      return bagnr != i;
    }

  };
  
  Iterator Begin () const
  {
    Iterator it(*this, 0, -1);
    it++;
    return it;
  }

  int End() const
  {
    return GetNBags();
  }

  void GetData (const Iterator & it,
		INDEX_2 & ahash, T & acont) const
  {
    ahash = hash[it.BagNr()][it.Pos()];
    acont = cont[it.BagNr()][it.Pos()];
  }

  const INDEX_2 & GetHash (const Iterator & it) const
  { return hash[it.BagNr()][it.Pos()]; }

  const T & GetData (const Iterator & it) const
  { return cont[it.BagNr()][it.Pos()]; }
};



template <typename T> 
inline ostream & operator<< (ostream & ost, const INDEX_2_HASHTABLE<T> & ht)
{
  for (typename INDEX_2_HASHTABLE<T>::Iterator it = ht.Begin();
       it != ht.End(); it++)
    {
      ost << ht.GetHash(it) << ": " << ht.GetData(it) << endl;
    }

  return ost;
}







///
class BASE_INDEX_3_HASHTABLE
{
protected:
  ///
  TABLE<INDEX_3> hash;

public:
  ///
  BASE_INDEX_3_HASHTABLE (int size)
    : hash (size) { };

protected:
  ///
  int HashValue (const INDEX_3 & ind) const
    {
      return (ind.I1() + ind.I2() + ind.I3()) % hash.Size() + 1;
    }

  ///
  int Position (int bnr, const INDEX_3 & ind) const
  {
    const INDEX_3 * pi = &hash.Get(bnr, 1);
    int n = hash.EntrySize(bnr);
    for (int i = 1; i <= n; ++i, ++pi)
      {
	if (*pi == ind)
	return i;
      }
    
    return 0;
  }


};


///
template <class T>
class INDEX_3_HASHTABLE : private BASE_INDEX_3_HASHTABLE
{
  ///
  TABLE<T> cont;

public:
  ///
  inline INDEX_3_HASHTABLE (int size);
  ///
  inline void Set (const INDEX_3 & ahash, const T & acont);
  ///
  inline const T & Get (const INDEX_3 & ahash) const;
  ///
  inline bool Used (const INDEX_3 & ahash) const;
  ///
  inline int GetNBags () const;
  ///
  inline int GetBagSize (int bnr) const;
  ///
  inline void SetData (int bnr, int colnr, const INDEX_3 & ahash, const T & acont);
  ///
  inline void GetData (int bnr, int colnr, INDEX_3 & ahash, T & acont) const;
  /// returns position, if not existing, will create (create == return 1)
  inline int PositionCreate (const INDEX_3 & ahash, int & bnr, int & colnr);
  ///
  inline void SetSize (int size);

  ///
  inline void PrepareSet (const INDEX_3 & ahash);
  ///
  inline void AllocateElements ();

  ///
  inline void PrintMemInfo (ostream & ost) const;
  ///
  inline void DeleteData ();









  class Iterator
  {
    const INDEX_3_HASHTABLE & ht;    
    int bagnr, pos;
  public:
    Iterator (const INDEX_3_HASHTABLE & aht,
	      int abagnr, int apos)
      : ht(aht), bagnr(abagnr), pos(apos)
    { ; }

    int BagNr() const { return bagnr; }
    int Pos() const { return pos; }

    void operator++ (int)
    {
      // cout << "begin Operator ++: bagnr = " << bagnr << " -  pos = " << pos << endl;
      pos++;
      while (bagnr < ht.GetNBags() && 
	     pos == ht.GetBagSize(bagnr+1))
	{
	  pos = 0;
	  bagnr++;
	}
      // cout << "end Operator ++: bagnr = " << bagnr << " - pos = " << pos << endl;
    }

    bool operator != (int i) const
    {
      return bagnr != i;
    }

  };
  
  Iterator Begin () const
  {
    Iterator it(*this, 0, -1);
    it++;
    return it;
  }

  int End() const
  {
    return GetNBags();
  }

  void GetData (const Iterator & it,
		INDEX_3 & ahash, T & acont) const
  {
    ahash = hash[it.BagNr()][it.Pos()];
    acont = cont[it.BagNr()][it.Pos()];
  }

  const INDEX_3 & GetHash (const Iterator & it) const
  { return hash[it.BagNr()][it.Pos()]; }

  const T & GetData (const Iterator & it) const
  { return cont[it.BagNr()][it.Pos()]; }



};






template <typename T> 
inline ostream & operator<< (ostream & ost, const INDEX_3_HASHTABLE<T> & ht)
{
  for (typename INDEX_3_HASHTABLE<T>::Iterator it = ht.Begin();
       it != ht.End(); it++)
    {
      ost << ht.GetHash(it) << ": " << ht.GetData(it) << endl;
    }

  return ost;
}
























/// Closed Hashing HT

class BASE_INDEX_CLOSED_HASHTABLE
{
protected:
  ///
  MoveableArray<INDEX> hash;
  ///
  int invalid;
public:
  ///
  BASE_INDEX_CLOSED_HASHTABLE (int size);

  int Size() const { return hash.Size(); }
  int UsedPos (int pos) const { return ! (hash.Get(pos) == invalid); }
  int UsedElements () const;

  ///
  int HashValue (const INDEX & ind) const
  {
    return ind % hash.Size() + 1;
  }


  int Position (const INDEX & ind) const
  {
    int i = HashValue(ind);
    while (1)
      {
	if (hash.Get(i) == ind) return i;
	if (hash.Get(i) == invalid) return 0;
	i++;
	if (i > hash.Size()) i = 1;
      }
  }

  // returns 1, if new postion is created
  int PositionCreate (const INDEX & ind, int & apos)
  {
    int i = HashValue (ind);
    if (hash.Get(i) == ind) 
      {
	apos = i;
	return 0;
      }
    if (hash.Get(i) == invalid)
      {
	hash.Elem(i) = ind; 
	apos = i;
	return 1;
      }
    return PositionCreate2 (ind, apos);    
  }

protected:
  int Position2 (const INDEX & ind) const;
  int PositionCreate2 (const INDEX & ind, int & apos);
  void BaseSetSize (int asize);
};


template <class T>
class INDEX_CLOSED_HASHTABLE : public BASE_INDEX_CLOSED_HASHTABLE
{
  ///
  MoveableArray<T> cont;

public:
  ///
  INDEX_CLOSED_HASHTABLE (int size)
    : BASE_INDEX_CLOSED_HASHTABLE(size), cont(size)
  {
    cont.SetName ("ind-hashtable, contents");
  }


  void Set (const INDEX & ahash, const T & acont)
  {
    int pos;
    PositionCreate (ahash, pos);
    hash.Elem(pos) = ahash;
    cont.Elem(pos) = acont;
  }

  const T & Get (const INDEX & ahash) const
  {
    int pos = Position (ahash);
    return cont.Get(pos);
  }

  ///
  bool Used (const INDEX & ahash) const
  {
    int pos = Position (ahash);
    return (pos != 0);
  }
  
  
  ///
  inline void SetData (int pos, const INDEX & ahash, const T & acont)
  {
    hash.Elem(pos) = ahash;
    cont.Elem(pos) = acont;
  }

  ///
  void GetData (int pos, INDEX & ahash, T & acont) const
  {
    ahash = hash.Get(pos);
    acont = cont.Get(pos);
  }
  
  ///
  inline void SetData (int pos, const T & acont)
  {
    cont.Elem(pos) = acont;
  }
  
  ///
  void GetData (int pos, T & acont) const
  {
    acont = cont.Get(pos);
  }
  
  ///
  const T & GetData (int pos) { return cont.Get(pos); }
  ///
  inline void SetSize (int size)
  {
    BaseSetSize(size);
  cont.SetSize(size);
  }
  
  ///
  inline void DeleteData ()
  { SetSize (cont.Size()); }

  void SetName (const char * aname)
  {
    cont.SetName(aname);
    hash.SetName(aname);
  }
};







/// Closed Hashing HT

class BASE_INDEX_2_CLOSED_HASHTABLE
{
protected:
  ///
  MoveableArray<INDEX_2> hash;
  ///
  int invalid;
public:
  ///
  BASE_INDEX_2_CLOSED_HASHTABLE (int size);

  int Size() const { return hash.Size(); }
  int UsedPos (int pos) const { return ! (hash.Get(pos).I1() == invalid); }
  int UsedElements () const;

  ///
  int HashValue (const INDEX_2 & ind) const
    {
      return (ind.I1() + 71 * ind.I2()) % hash.Size() + 1;
    }


  int Position (const INDEX_2 & ind) const
  {
    int i = HashValue(ind);
    while (1)
      {
	if (hash.Get(i) == ind) return i;
	if (hash.Get(i).I1() == invalid) return 0;
	i++;
	if (i > hash.Size()) i = 1;
      }
  }

  // returns 1, if new postion is created
  int PositionCreate (const INDEX_2 & ind, int & apos)
  {
    int i = HashValue (ind);
    if (hash.Get(i) == ind) 
      {
	apos = i;
	return 0;
      }
    if (hash.Get(i).I1() == invalid)
      {
	hash.Elem(i) = ind; 
	apos = i;
	return 1;
      }
    return PositionCreate2 (ind, apos);    
  }

protected:
  ///

  int Position2 (const INDEX_2 & ind) const;
  int PositionCreate2 (const INDEX_2 & ind, int & apos);
  void BaseSetSize (int asize);
};


template <class T>
class INDEX_2_CLOSED_HASHTABLE : public BASE_INDEX_2_CLOSED_HASHTABLE
{
  ///
  MoveableArray<T> cont;

public:
  ///
  inline INDEX_2_CLOSED_HASHTABLE (int size);
  ///
  inline void Set (const INDEX_2 & ahash, const T & acont);
  ///
  inline const T & Get (const INDEX_2 & ahash) const;
  ///
  inline bool Used (const INDEX_2 & ahash) const;
  ///
  inline void SetData (int pos, const INDEX_2 & ahash, const T & acont);
  ///
  inline void GetData (int pos, INDEX_2 & ahash, T & acont) const;
  ///
  inline void SetData (int pos, const T & acont);
  ///
  inline void GetData (int pos, T & acont) const;
  ///
  const T & GetData (int pos) { return cont.Get(pos); }
  ///
  inline void SetSize (int size);
  ///
  inline void PrintMemInfo (ostream & ost) const;
  ///
  inline void DeleteData ()
  { SetSize (cont.Size()); }

  void SetName (const char * aname)
  {
    cont.SetName(aname);
    hash.SetName(aname);
  }
};

class BASE_INDEX_3_CLOSED_HASHTABLE
{
protected:
  MoveableArray<INDEX_3> hash;
  int invalid;

protected: 
  // public:    //SZ
  BASE_INDEX_3_CLOSED_HASHTABLE (int size)
    : hash(size)
  {
    hash.SetName ("i3-hashtable, hash");
    invalid = -1;
    for (int i = 0; i < size; i++)
      hash[i].I1() = invalid;
  }

public:
  int Size() const 
  {
    return hash.Size(); 
  }

  bool UsedPos (int pos) const 
  { 
    return ! (hash[pos].I1() == invalid); 
  }

  int UsedElements () const
  {
    int n = hash.Size();
    int cnt = 0;
    for (int i = 0; i < n; i++)
      if (hash[i].I1() != invalid)
	cnt++;
    return cnt;
  }

  int HashValue (const INDEX_3 & ind) const
  {
    return (ind.I1() + 15 * ind.I2() + 41 * ind.I3()) % hash.Size();
  }
  
  int Position (const INDEX_3 & ind) const
  {
    int i = HashValue(ind);
    while (1)
      {
	if (hash[i] == ind) return i;
	if (hash[i].I1() == invalid) return -1;
        i = (i+1) % hash.Size();
      }
  }

  int Costs (const INDEX_3 & ind) const
  {
    int i = HashValue(ind);
    int c = 1;
    while (1)
      {
	if (hash[i] == ind) return c;
	if (hash[i].I1() == invalid) return c;
        i = (i+1) % hash.Size();
        c++;
      }
  }


  
  // returns true, if new postion is created
  bool PositionCreate (const INDEX_3 & ind, int & apos)
  {
    int i = HashValue (ind);
    if (hash[i] == ind) 
      {
	apos = i;
	return false;
      }
    if (hash[i].I1() == invalid)
      {
	hash[i] = ind; 
	apos = i;
	return true;
      }
    return PositionCreate2 (ind, apos);    
  }

protected:
  bool PositionCreate2 (const INDEX_3 & ind, int & apos);
  void BaseSetSize (int asize);
};



template <class T>
class INDEX_3_CLOSED_HASHTABLE : public BASE_INDEX_3_CLOSED_HASHTABLE
{
  MoveableArray<T,0> cont;

public:
  INDEX_3_CLOSED_HASHTABLE (int size)
    : BASE_INDEX_3_CLOSED_HASHTABLE(size), cont(size)
  {
    cont.SetName ("i3-hashtable, contents");
  }
  
  void Set (const INDEX_3 & ahash, const T & acont)
  {
    int pos;
    PositionCreate (ahash, pos);
    hash[pos] = ahash;
    cont[pos] = acont;
  }

  const T & Get (const INDEX_3 & ahash) const
  {
    return cont[Position (ahash)];
  }

  bool Used (const INDEX_3 & ahash) const
  {
    return (Position (ahash) != -1);
  }

  void SetData (int pos, const INDEX_3 & ahash, const T & acont)
  {
    hash[pos] = ahash;
    cont[pos] = acont;
  }

  void GetData (int pos, INDEX_3 & ahash, T & acont) const
  {
    ahash = hash[pos];
    acont = cont[pos];
  }

  void SetData (int pos, const T & acont)
  {
    cont[pos] = acont;
  }

  void GetData (int pos, T & acont) const
  {
    acont = cont[pos];
  }

  const T & GetData (int pos) const
  {
    return cont[pos];
  }

   void SetSize (int size)
  {
    BaseSetSize(size);
    cont.SetSize(size);
  }

  void PrintMemInfo (ostream & ost) const
  {
    cout << "Hashtable: " << Size() 
         << " entries of size " << sizeof(INDEX_3) << " + " << sizeof(T) 
         << " = " << Size() * (sizeof(INDEX_3) + sizeof(T)) << " bytes" << endl;
    
  }

  void DeleteData ()
  { 
    SetSize (cont.Size()); 
  }

  void SetName (const char * aname)
  {
    cont.SetName(aname);
    hash.SetName(aname);
  }
};

















template<class T>
inline INDEX_3_HASHTABLE<T> :: INDEX_3_HASHTABLE (int size)
  : BASE_INDEX_3_HASHTABLE (size), cont(size)
{
  ;
}

template<class T>	
inline int INDEX_3_HASHTABLE<T> :: PositionCreate (const INDEX_3 & ahash, int & bnr, int & colnr)
{
  bnr = HashValue (ahash);
  colnr = Position (bnr, ahash);
  if (!colnr)
    {
      hash.Add (bnr, ahash);
      cont.AddEmpty (bnr);
      colnr = cont.EntrySize (bnr);
      return 1;
    }
  return 0;
}


template<class T>
inline void INDEX_3_HASHTABLE<T> :: Set (const INDEX_3 & ahash, const T & acont)
{
  int bnr = HashValue (ahash);
  int pos = Position (bnr, ahash);
  if (pos)
    cont.Set (bnr, pos, acont);
  else
    {
      hash.Add1 (bnr, ahash);
      cont.Add1 (bnr, acont);
    }
}

template<class T>
inline const T & INDEX_3_HASHTABLE<T> :: Get (const INDEX_3 & ahash) const
{
  int bnr = HashValue (ahash);
  int pos = Position (bnr, ahash);
  return cont.Get (bnr, pos);
}

template<class T>
inline bool INDEX_3_HASHTABLE<T> :: Used (const INDEX_3 & ahash) const
{
  return (Position (HashValue (ahash), ahash)) ? 1 : 0;
}

template<class T>
inline int INDEX_3_HASHTABLE<T> :: GetNBags () const
{
  return cont.Size();
}

template<class T>
inline int INDEX_3_HASHTABLE<T> :: GetBagSize (int bnr) const
{
  return cont.EntrySize (bnr);
}

template<class T>
inline void INDEX_3_HASHTABLE<T> :: GetData (int bnr, int colnr, INDEX_3 & ahash, T & acont) const
{
  ahash = hash.Get(bnr, colnr);
  acont = cont.Get(bnr, colnr);
}    

template<class T>
inline void INDEX_3_HASHTABLE<T> :: SetData (int bnr, int colnr, const INDEX_3 & ahash, const T & acont)
{
  hash.Set(bnr, colnr, ahash);
  cont.Set(bnr, colnr, acont);
}    

template<class T>
inline void INDEX_3_HASHTABLE<T> :: SetSize (int size)
{
  hash.SetSize (size);
  cont.SetSize (size);
}

template<class T>
inline void INDEX_3_HASHTABLE<T> :: DeleteData ()
{
  int n = hash.Size();
  hash.SetSize (n);
  cont.SetSize (n);
}

template<class T>
inline void INDEX_3_HASHTABLE<T> :: PrepareSet (const INDEX_3 & ahash)
{
  int bnr = HashValue (ahash);
  hash.IncSizePrepare (bnr-1);
  cont.IncSizePrepare (bnr-1);
}


template<class T>
inline void INDEX_3_HASHTABLE<T> :: AllocateElements ()
{
  hash.AllocateElementsOneBlock();
  cont.AllocateElementsOneBlock();
}



template<class T>
inline void INDEX_3_HASHTABLE<T> :: PrintMemInfo (ostream & ost) const
  {
    ost << "Hash: " << endl;
    hash.PrintMemInfo (ost);
    ost << "Cont: " << endl;
    cont.PrintMemInfo (ost);
  }





template<class T>
inline INDEX_HASHTABLE<T> :: INDEX_HASHTABLE (int size)
  : BASE_INDEX_HASHTABLE (size), cont(size)
  {
    ;
  }
	
template<class T>
inline void INDEX_HASHTABLE<T> :: Set (const INDEX & ahash, const T & acont)
    {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    if (pos)
      cont.Set (bnr, pos, acont);
    else
      {
      hash.Add (bnr, ahash);
      cont.Add (bnr, acont);
      }
    }

template<class T>
inline const T & INDEX_HASHTABLE<T> :: Get (const INDEX & ahash) const
    {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    return cont.Get (bnr, pos);
    }

template<class T>
inline bool INDEX_HASHTABLE<T> :: Used (const INDEX & ahash) const
    {
    return (Position (HashValue (ahash), ahash)) ? 1 : 0;
    }

template<class T>
inline int INDEX_HASHTABLE<T> :: GetNBags () const
    {
    return hash.Size();
    }

template<class T>
inline int INDEX_HASHTABLE<T> :: GetBagSize (int bnr) const
    {
    return hash.EntrySize(bnr);
    }

template<class T>
inline void INDEX_HASHTABLE<T> :: GetData (int bnr, int colnr, INDEX & ahash, T & acont) const
    {
    ahash = hash.Get(bnr, colnr);   
    acont = cont.Get(bnr, colnr);
    }
    
template<class T>
inline void INDEX_HASHTABLE<T> :: PrintMemInfo (ostream & ost) const
  {
    ost << "Hash: " << endl;
    hash.PrintMemInfo (ost);
    ost << "Cont: " << endl;
    cont.PrintMemInfo (ost);
  }


    
    
    
    
    








/* *********** Closed Hashing ************************* */



template<class T>
inline INDEX_2_CLOSED_HASHTABLE<T> :: 
INDEX_2_CLOSED_HASHTABLE (int size)
  : BASE_INDEX_2_CLOSED_HASHTABLE(size), cont(size)
{
  cont.SetName ("i2-hashtable, contents");
}

template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
Set (const INDEX_2 & ahash, const T & acont)
{
  int pos;
  PositionCreate (ahash, pos);
  hash.Elem(pos) = ahash;
  cont.Elem(pos) = acont;
}

template<class T>
inline const T & INDEX_2_CLOSED_HASHTABLE<T> :: 
Get (const INDEX_2 & ahash) const
{
  int pos = Position (ahash);
  return cont.Get(pos);
}

template<class T>
inline bool INDEX_2_CLOSED_HASHTABLE<T> :: 
Used (const INDEX_2 & ahash) const
{
  int pos = Position (ahash);
  return (pos != 0);
}

template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
SetData (int pos, const INDEX_2 & ahash, const T & acont)
{
  hash.Elem(pos) = ahash;
  cont.Elem(pos) = acont;
}
  
template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
GetData (int pos, INDEX_2 & ahash, T & acont) const
{
  ahash = hash.Get(pos);
  acont = cont.Get(pos);
}

template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
SetData (int pos, const T & acont)
{
  cont.Elem(pos) = acont;
}
  
template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
GetData (int pos, T & acont) const
{
  acont = cont.Get(pos);
}


template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
SetSize (int size)
{
  BaseSetSize(size);
  cont.SetSize(size);
}


  
template<class T>
inline void INDEX_2_CLOSED_HASHTABLE<T> :: 
PrintMemInfo (ostream & ost) const
{
  cout << "Hashtable: " << Size() 
       << " entries of size " << sizeof(INDEX_2) << " + " << sizeof(T) 
       << " = " << Size() * (sizeof(INDEX_2) + sizeof(T)) << " bytes." 
       << " Used els: " << UsedElements() 
       << endl;
}
















/*
template<class T>
inline INDEX_3_CLOSED_HASHTABLE<T> :: 
INDEX_3_CLOSED_HASHTABLE (int size)
  : BASE_INDEX_3_CLOSED_HASHTABLE(size), cont(size)
{
  cont.SetName ("i3-hashtable, contents");
}

template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
Set (const INDEX_3 & ahash, const T & acont)
{
  int pos;
  PositionCreate (ahash, pos);
  hash.Elem(pos) = ahash;
  cont.Elem(pos) = acont;
}

template<class T>
inline const T & INDEX_3_CLOSED_HASHTABLE<T> :: 
Get (const INDEX_3 & ahash) const
{
  int pos = Position (ahash);
  return cont[pos];
}

template<class T>
inline bool INDEX_3_CLOSED_HASHTABLE<T> :: 
Used (const INDEX_3 & ahash) const
{
  int pos = Position (ahash);
  return (pos != 0);
}

template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
SetData (int pos, const INDEX_3 & ahash, const T & acont)
{
  hash.Elem(pos) = ahash;
  cont.Elem(pos) = acont;
}
  
template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
GetData (int pos, INDEX_3 & ahash, T & acont) const
{
  ahash = hash.Get(pos);
  acont = cont.Get(pos);
}

template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
SetData (int pos, const T & acont)
{
  cont.Elem(pos) = acont;
}
  
template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
GetData (int pos, T & acont) const
{
  acont = cont.Get(pos);
}

template<class T>
inline const T & INDEX_3_CLOSED_HASHTABLE<T> :: 
GetData (int pos) const
{
  return cont.Get(pos);
}


template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
SetSize (int size)
{
  BaseSetSize(size);
  cont.SetSize(size);
}
  
template<class T>
inline void INDEX_3_CLOSED_HASHTABLE<T> :: 
PrintMemInfo (ostream & ost) const
{
  cout << "Hashtable: " << Size() 
       << " entries of size " << sizeof(INDEX_3) << " + " << sizeof(T) 
       << " = " << Size() * (sizeof(INDEX_3) + sizeof(T)) << " bytes" << endl;
}
*/

















#endif
