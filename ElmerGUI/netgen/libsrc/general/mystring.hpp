
//**************************************************************
//
// filename:             mystring.h
//
// project:              doctoral thesis, program smart
//
// autor:                Dipl.-Ing. Gerstmayr Johannes
//
// generated:            20.12.98
// last change:          20.12.98
// description:          base class for strings
// remarks:              string with n characters has
//                       0..n-1 characters and at pos n a 0
//
//**************************************************************


#ifndef MYSTRING__H
#define MYSTRING__H

class Point3d;
class Vec3d;


// extract string str which is enclosed by the given character encl from a given string in
void ReadEnclString(istream & in, string & str, const char encl);


class MyStr;

MyStr operator + (const MyStr &, const MyStr &);
int operator == (const MyStr &, const MyStr &);
int operator < (const MyStr &, const MyStr &);
int operator <= (const MyStr &, const MyStr &);
int operator > (const MyStr &, const MyStr &);
int operator >= (const MyStr &, const MyStr &);
int operator != (const MyStr &, const MyStr &);
ostream& operator << (ostream &, const MyStr &);
istream& operator >> (istream &, MyStr &);

class MyStr
{
public:
  MyStr();
  MyStr(const char *);
  MyStr(char);
  MyStr(const MyStr &);
  MyStr(int);
  MyStr(void *);
  MyStr(long);
  MyStr(double);
  MyStr(const Point3d& p);
  MyStr(const Vec3d& p);
  MyStr(const string & st);

  ~MyStr();
  MyStr Left(unsigned);
  MyStr Right(unsigned);
  MyStr& InsertAt(unsigned, const MyStr &);
  MyStr& WriteAt(unsigned, const MyStr &);
  unsigned Length() const;
  int Find(const char);
  int Find(const char *);
  int Find(const MyStr &);
  MyStr& operator = (const MyStr &);
  friend MyStr operator + (const MyStr &, const MyStr &);
  void operator += (const MyStr &);
  char* c_str();
  string cpp_string(void) const;

  //change every ',' -> ';', '.' -> ','
  void ConvertTextToExcel();
  //change every ','->'.', ';'->','
  void ConvertExcelToText();

  MyStr operator () (unsigned, unsigned);
  operator int();
  operator double();
  operator long();
  operator char *();
  char& operator [] (unsigned int);
  char operator [] (unsigned int) const;

  friend int operator == (const MyStr &, const MyStr &);
  friend int operator < (const MyStr &, const MyStr &);
  friend int operator <= (const MyStr &, const MyStr &);
  friend int operator > (const MyStr &, const MyStr &);
  friend int operator >= (const MyStr &, const MyStr &);
  friend int operator != (const MyStr &, const MyStr &);
  friend ostream& operator << (ostream &, const MyStr &);
  friend istream& operator >> (istream &, MyStr &);
  static void SetToErrHandler(void (*)());
private:
  MyStr(unsigned, int);
  char *str;
  unsigned length;
  enum { SHORTLEN = 24 };
  char shortstr[SHORTLEN+1];
  static void(*ErrHandler)();
};


inline MyStr::MyStr()
{
  length = 0;
  str = shortstr;
  str[0] = 0;
}

inline MyStr::MyStr(char s)
{
  length = 1;
  str = shortstr;
  str[0] = s;
  str[1] = (char)0;
}

inline MyStr::~MyStr()
{
  if (length > SHORTLEN)
    delete [] str;
}

inline unsigned MyStr::Length() const
{
  return length;
}

inline int MyStr::Find(const char c)
{
  char *pos = strchr(str, int(c));
  return pos ? int(pos - str) : -1;
}

inline int MyStr::Find(const MyStr &s)
{
  char *pos = strstr(str, s.str);
  return pos ? int(pos - str) : -1;
}

inline int MyStr::Find(const char *s)
{
  char *pos = strstr(str, s);
  return pos ? int(pos - str) : -1;
}

inline MyStr::operator int()
{
  return atoi(str);
}

inline MyStr::operator double()
{
  return atof(str);
}

inline MyStr::operator long()
{
  return atol(str);
}

inline MyStr::operator char *()
{
  return str;
}

inline char* MyStr::c_str()
{
  return str;
}


inline int operator == (const MyStr &s1, const MyStr& s2)
{
  return strcmp(s1.str, s2.str) == 0;
}

inline int operator < (const MyStr &s1, const MyStr& s2)
{
  return strcmp(s1.str, s2.str) < 0;
}

inline int operator <= (const MyStr &s1, const MyStr& s2)
{
  return strcmp(s1.str, s2.str) <= 0;
}

inline int operator > (const MyStr &s1, const MyStr& s2)
{
  return strcmp(s1.str, s2.str) > 0;
}

inline int operator >= (const MyStr &s1, const MyStr& s2)
{
  return strcmp(s1.str, s2.str) >= 0;
}

inline int operator != (const MyStr &s1, const MyStr& s2)
{
  return !(s1 == s2);
}

inline ostream& operator << (ostream& os, const MyStr& s)
{
  return os << s.str;
}

inline void MyStr::SetToErrHandler(void (*Handler)())
{
  ErrHandler = Handler;
};

#endif

   
