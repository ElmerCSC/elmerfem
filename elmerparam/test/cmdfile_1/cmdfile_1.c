#include <stdlib.h>
#include <elmerparam.h>

int main ()
{
  int xi[2] = { 1, 2 };
  double xr[2] = { 3.0, 4.0 };
  double y;

  y = elmer_param(2, xr, 2, xi,  "FOO");
  if (y != 22.0) 
    return 1;
  else
    return 0;
}
