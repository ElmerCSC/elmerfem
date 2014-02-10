//File for handling warnings, errors, messages
#include <meshing.hpp>

namespace netgen
{

int printmessage_importance = 5;
int printwarnings = 1;
int printerrors = 1;
int printdots = 1;
int printfnstart = 0;

// extern void Ng_PrintDest(const MyStr& s);
extern void Ng_PrintDest(const char * s);

//the dots for progression of program
void PrintDot(char ch)
{
  if (printdots)
    {
      char st[2];
      st[0] = ch;
      st[1] = 0;
      Ng_PrintDest(st);
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+MyStr("\n"));
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+MyStr("\n"));
    }
}

void PrintMessage(int importance, 
		  const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
    }
}

void PrintMessageCR(int importance, 
		    const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		    const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (importance <= printmessage_importance)
    {
      Ng_PrintDest(MyStr(" ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\r"));
    }
}

void PrintFnStart(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printfnstart)
    Ng_PrintDest(MyStr(" Start Function: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintWarning(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		  const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printwarnings)
    Ng_PrintDest(MyStr(" WARNING: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintFileError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		    const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" FILE ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintUserError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  Ng_PrintDest(MyStr(" USER ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintSysError(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
		const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printerrors)
    Ng_PrintDest(MyStr(" SYSTEM ERROR: ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}

void PrintTime(const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
	       const MyStr& s5, const MyStr& s6, const MyStr& s7, const MyStr& s8)
{
  if (printmessage_importance >= 3)
    Ng_PrintDest(MyStr(" Time = ")+s1+s2+s3+s4+s5+s6+s7+s8+MyStr("\n"));
}


static ARRAY<MyStr*> msgstatus_stack(0);
static ARRAY<double> threadpercent_stack(0);
static MyStr msgstatus = "";




void ResetStatus()
{
  SetStatMsg("idle");

  for (int i = 0; i < msgstatus_stack.Size(); i++)
    delete msgstatus_stack[i];
  msgstatus_stack.SetSize(0);
  threadpercent_stack.SetSize(0);

  // multithread.task = "";
  multithread.percent = 100.;
}

void PushStatus(const MyStr& s)
{
  msgstatus_stack.Append(new MyStr (s));  
  SetStatMsg(s);
  threadpercent_stack.Append(0);
}

void PushStatusF(const MyStr& s)
{
  msgstatus_stack.Append(new MyStr (s));
  SetStatMsg(s);
  threadpercent_stack.Append(0);
  PrintFnStart(s);
}

void PopStatus()
{
  if (msgstatus_stack.Size())
    {
      if (msgstatus_stack.Size() > 1)
	SetStatMsg (*msgstatus_stack.Last());
      else
	SetStatMsg ("");
      delete msgstatus_stack.Last();
      msgstatus_stack.DeleteLast();
      threadpercent_stack.DeleteLast();
      if(threadpercent_stack.Size() > 0)
	multithread.percent = threadpercent_stack.Last();
      else
	multithread.percent = 100.;
    }
  else
    {
      PrintSysError("PopStatus failed");
    }
}



/*
void SetStatMsgF(const MyStr& s)
{
  PrintFnStart(s);
  SetStatMsg(s);
}
*/

void SetStatMsg(const MyStr& s)
{
  msgstatus = s;
  multithread.task = msgstatus.c_str();  
}

void SetThreadPercent(double percent)
{
  multithread.percent = percent;
  if(threadpercent_stack.Size() > 0)
    threadpercent_stack.Last() = percent;
}


void GetStatus(MyStr & s, double & percentage)
{
  if(threadpercent_stack.Size() > 0)
    percentage = threadpercent_stack.Last();
  else
    percentage = multithread.percent;
  
  if ( msgstatus_stack.Size() )
    s = *msgstatus_stack.Last();
  else
    s = "idle";     
}


#ifdef SMALLLIB
#define SMALLLIBORNOTCL
#endif
#ifdef NOTCL
#define SMALLLIBORNOTCL
#endif

#ifdef SMALLLIBORNOTCL
void Ng_PrintDest(const char * s){cout << s <<flush;}
double GetTime(){return 0;}
void MyError(const char * ch)
{
  cerr << ch << endl;
}
#endif

}
