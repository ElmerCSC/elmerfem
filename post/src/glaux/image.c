#include "../tk/tk.h"
#include "glaux.h"

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif


AUX_RGBImageRec *auxRGBImageLoad(char *fileName)
{
    return (AUX_RGBImageRec *)tkRGBImageLoad(fileName);
}
