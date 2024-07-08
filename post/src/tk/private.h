#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif


typedef struct _WINDOW_REC {
    int x, y, w, h;
    GLenum type;
    Window wMain, wOverlay;
    XVisualInfo *vInfoMain, *vInfoOverlay;
    Colormap cMapMain, cMapOverlay;
    GLXContext cMain, cOverlay;
    int cMapAllocated;  /*** ad@lms.be: some cMapMain's better not freed ***/
} WINDOW_REC;


extern Display *xDisplay;
extern int xScreen; 
extern Window wRoot;
extern WINDOW_REC w;
extern Atom deleteWindowAtom;

extern void (*ExposeFunc)(int, int);
extern void (*ReshapeFunc)(int, int);
extern void (*DisplayFunc)(void);
extern GLenum (*KeyDownFunc)(int, GLenum);
extern GLenum (*MouseDownFunc)(int, int, GLenum);
extern GLenum (*MouseUpFunc)(int, int, GLenum);
extern GLenum (*MouseMoveFunc)(int, int, GLenum);
extern void (*IdleFunc)(void);

extern GLenum drawAllowFlag;

extern int cursorNum;

