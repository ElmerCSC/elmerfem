
//**************************************************************
//
// filename:             mystring.cpp
//
// project:              doctoral thesis
//
// autor:                Dipl.-Ing. Gerstmayr Johannes
//
// generated:            20.12.98
// last change:          20.12.98
// description:          implementation for strings
// remarks:
//
//**************************************************************

// string class
#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <gprim.hpp>

/*
#include <iostream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream.h>
#include "mystring.hh"

 */

namespace netgen
{


  void ReadEnclString(istream & in, string & str, const char encl)
  {
    char currchar;
    str = "";

    in.get(currchar);
    while(in && currchar == ' ' || currchar == '\t' || currchar == '\n')
      in.get(currchar);
	
    if(currchar == encl)
      {
	in.get(currchar);
	while(in && currchar != encl)
	  {
	    str += currchar;
	    in.get(currchar);
	  }
      }
    else
      {
	in.putback(currchar);
	in >> str;
      }
  }
	    
	
    
  


void DefaultStringErrHandler()
{
  cerr << "Error : string operation out of range\n" << flush;
}

void (*MyStr::ErrHandler)() = DefaultStringErrHandler;

  /*     
MyStr::MyStr()
{
  length = 0;
  str = shortstr;
  str[0] = 0;
}
  */

MyStr::MyStr(const char *s)
{
  length = unsigned(strlen(s));

  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, s);
}

/*
MyStr::MyStr(char s)
{
  length = 1;
  str = shortstr;
  str[0] = s;
  str[1] = (char)0;
}
*/

MyStr::MyStr(const MyStr& s)
{
  length = s.length;
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, s.str);
}

MyStr::MyStr(int i)
{
  char buffer[32];
  sprintf(buffer, "%d", i);
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}

MyStr::MyStr(void * p)
{
  char buffer[32];
  sprintf(buffer, "%p", p);
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}


MyStr::MyStr(long l)
{
  char buffer[32];
  sprintf(buffer, "%ld", l);
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}

MyStr::MyStr(double d)
{
  char buffer[32];
  //if (fabs(d) < 1E-100) {d = 0;}
  sprintf(buffer, "%g", d);
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}

MyStr::MyStr(const Point3d& p)
{
  char buffer[80];
  //if (fabs(d) < 1E-100) {d = 0;}
  sprintf(buffer, "[%g, %g, %g]", p.X(), p.Y(), p.Z());
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}

MyStr::MyStr(const Vec3d& p)
{
  char buffer[80];
  //if (fabs(d) < 1E-100) {d = 0;}
  sprintf(buffer, "[%g, %g, %g]", p.X(), p.Y(), p.Z());
  length = unsigned(strlen(buffer));
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, buffer);
}

MyStr::MyStr(unsigned n, int)
{
  length = n;
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  str[n] = 0;
}

MyStr::MyStr(const string & st)
{
  length = unsigned(st.length());
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy (str, st.c_str());
}



MyStr MyStr::Left(unsigned r)
{
  if(r > length)
    {
      MyStr::ErrHandler();
      MyStr s;
      return s;
    }
  else
    {
      MyStr tmp(r, 0);
      strncpy(tmp.str, str, r);
      return tmp;
    }
}

MyStr MyStr::Right(unsigned l)
{
  if(l > length)
    {
      MyStr::ErrHandler();
      MyStr s;
      return s;
    }
  else
    {
      MyStr tmp(l, 0);
      strncpy(tmp.str, str + length - l, l);
      return tmp;
    }
}

MyStr& MyStr::InsertAt(unsigned pos, const MyStr& s)
{
  if(pos > length)
    {
      MyStr::ErrHandler();
      return *this;
    }
  int newLength = length + s.length;
  char *tmp = new char[newLength + 1];
  strncpy(tmp, str, pos);
  strcpy(tmp + pos, s.str);
  strcpy(tmp + pos + s.length, str + pos);

  if (length > SHORTLEN) delete [] str;
  length = newLength;
  if (length > SHORTLEN)
    str = tmp;
  else
    {
      strcpy (shortstr, tmp);
      delete [] tmp;
      str = shortstr;
    }
  return *this;
}

MyStr &MyStr::WriteAt(unsigned pos, const MyStr& s)
{
  if(pos > length)
  {
    MyStr::ErrHandler();
    return *this;
  }
  unsigned n = length - pos;
  if(s.length < n)
    n = s.length;
  strncpy(str + pos, s.str, n);
  return *this;
}

void MyStr::ConvertTextToExcel()
{
  /*
  for (int i = 0; i < Length(); i++)
    {
      if ((*this)[i]==',') {(*this)[i] = ';';}
      else if ((*this)[i]=='.') {(*this)[i] = ',';}
    }
  */
}

void MyStr::ConvertExcelToText()
{
  /*
  for (int i = 0; i < Length(); i++)
    {
      if ((*this)[i]==',') {(*this)[i] = '.';}
      else if ((*this)[i]==';') {(*this)[i] = ',';}
    }
  */
}

MyStr& MyStr::operator = (const MyStr& s)
{
  if (length > SHORTLEN) delete [] str;
  length = s.length;
  if (length > SHORTLEN)
    str = new char[length + 1];
  else
    str = shortstr;
  strcpy(str, s.str);
  return *this;
}

MyStr operator + (const MyStr& s1, const MyStr& s2)
{
  MyStr tmp(s1.length + s2.length, 0);
  if (s1.length != 0) strcpy(tmp.str, s1.str);
  if (s2.length != 0) strcpy(tmp.str + s1.length, s2.str);
  return tmp;
}

void MyStr::operator += (const MyStr& s)
{
  if (length+s.length <= SHORTLEN)
    {
      if (s.length != 0) strcpy(shortstr + length, s.str);
    }
  else
    {
      char *tmp = new char[length + s.length + 1];
      if (length != 0) strcpy(tmp, str);
      if (s.length != 0) strcpy(tmp + length, s.str);
      if (length > SHORTLEN) delete [] str;
      length += s.length;
      str = tmp;
    }
}

char& MyStr::operator [] (unsigned n)
{
  static char dummy;
  if(n < length)
    return str[n];
  else
  {
    MyStr::ErrHandler();
    return dummy;
  }
}

char MyStr::operator [] (unsigned n) const
{
  static char dummy;
  if(n < length)
    return str[n];
  else
  {
    MyStr::ErrHandler();
    return dummy;
  }
}

MyStr MyStr::operator () (unsigned l, unsigned r)
{
  if((l > r) || (r > length))
  {
    MyStr::ErrHandler();
    MyStr s;
    return s;
  }
  else
  {
    int n = r - l + 1;
    MyStr tmp(n, 0);
    strncpy(tmp.str, str + 1, n);
    return tmp;
  }
}

string MyStr::cpp_string(void) const
{
  string aux(str,length);
  return aux;
}

/*
istream& operator >> (istream& is, MyStr& s)
{
  const int buflen = 1000;
  char buffer[buflen+1];

  int end = 0;
  s = "";
  MyStr str;

  while (!end)
  {
    is.get(buffer, buflen);
    str = MyStr(buffer);
    s += str;
    if (is.peek() == EOF) {end = 1;}
  }

  return is;
}
*/
/*
#ifdef __borland
::ifstream& operator >> (::ifstream& is, MyStr& s)       // wb
{                                                        // wb
  const int buflen = 1000;                               // wb
  char buffer[buflen+1];                                 // wb
                                                         // wb
  int end = 0;                                           // wb
  s = "";                                                // wb
  MyStr str;                                             // wb
                                                         // wb
  while (!end)                                           // wb
  {                                                      // wb
    is.get(buffer, buflen);                              // wb
    str = MyStr(buffer);                                 // wb
    s += str;                                            // wb
    if (is.peek() == EOF) {end = 1;}                     // wb
  }                                                      // wb
                                                         // wb
  return is;                                             // wb
}          

#endif
*/
}
