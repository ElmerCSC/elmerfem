#include "glpfile.h"
#include <stdio.h>

GLPfile	*output;
GLPtext       *texthead;

extern "C" void initglp( char *name, int FitToPage )
{
  if ( FitToPage ) {
    output = new GLPfile( name,1 );
  } else {
    output = new GLPfile( name,0 );
  }
  if ( output != NULL ) output->StartPage();
}

extern "C" void stopglp()
{
  if (output != NULL)
  {
    if (output->EndPage() == 0)
    {
      delete output;
      output = NULL;
    };
  };
}
