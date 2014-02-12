#include <stdlib.h>
#include <elmerparam.h>

int main()
{
    double xr[10], o[10];
    int i;

    for (i = 0; i < 10; i++)
        xr[i] = (double) i;

    elmer_param_vec(10, o, 10, xr, 0, NULL, NULL);

    for (i = 0; i < 10; i++)
        if (o[i] != xr[i]) return 1;

    return 0;
}
