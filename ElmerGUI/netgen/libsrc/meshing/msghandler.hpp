#ifndef FILE_MSGHANDLER
#define FILE_MSGHANDLER

/**************************************************************************/
/* File:   msghandler.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/


extern void PrintDot(char ch = '.');


//Message Pipeline:

//importance: importance of message: 1=very important, 3=middle, 5=low, 7=unimportant
extern void PrintMessage(int importance, 
			 const MyStr& s1, const MyStr& s2=MyStr());
extern void PrintMessage(int importance, 
			 const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4=MyStr());
extern void PrintMessage(int importance, 
			 const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
			 const MyStr& s5, const MyStr& s6=MyStr(), const MyStr& s7=MyStr(), const MyStr& s8=MyStr());

// CR without line-feed
extern void PrintMessageCR(int importance, 
			   const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			   const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintFnStart(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			 const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintWarning(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			 const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
		       const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintFileError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
		       const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintSysError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
		       const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintUserError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
		       const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void PrintTime(const MyStr& s1="", const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
		      const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
extern void SetStatMsg(const MyStr& s);

extern void PushStatus(const MyStr& s);
extern void PushStatusF(const MyStr& s);
extern void PopStatus();
extern void SetThreadPercent(double percent);
extern void GetStatus(MyStr & s, double & percentage);


#endif

