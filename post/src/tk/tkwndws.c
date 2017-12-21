/*
 * (c) Copyright 1993, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */

// Mesa Tweaking by: Mark E. Peterson (markp@ic.mankato.mn.us)

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tk.h"
/* #include "gl\wmesa.h" */

#define static

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif

#if DBG
#define TKASSERT(x)                                     \
if ( !(x) ) {                                           \
    PrintMessage("%s(%d) Assertion failed %s\n",        \
        __FILE__, __LINE__, #x);                        \
}
#else
#define TKASSERT(x)
#endif  /* DBG */

/******************************************************************************/

static struct _WINDOWINFO {
    int x, y;
    int width, height;
    GLenum type;
    GLenum dmPolicy;
    int ipfd;
    BOOL bDefPos;
} windInfo = {
    0, 0, 100, 100, TK_INDEX | TK_SINGLE, TK_MINIMUM_CRITERIA, 0, TRUE
};


static HWND     tkhwnd     = NULL;
static HDC      tkhdc      = NULL;
static HPALETTE tkhpalette = NULL;
GLboolean tkPopupEnable = TRUE;

// Fixed palette support.

#define BLACK   PALETTERGB(0,0,0)
#define WHITE   PALETTERGB(255,255,255)
#define NUM_STATIC_COLORS   (COLOR_BTNHIGHLIGHT - COLOR_SCROLLBAR + 1)

static void (*ExposeFunc)(int, int)              = NULL;
static void (*ReshapeFunc)(GLsizei, GLsizei)     = NULL;
static void (*DisplayFunc)(void)                 = NULL;
static GLenum (*KeyDownFunc)(int, GLenum)        = NULL;
static GLenum (*MouseDownFunc)(int, int, GLenum) = NULL;
static GLenum (*MouseUpFunc)(int, int, GLenum)   = NULL;
static GLenum (*MouseMoveFunc)(int, int, GLenum) = NULL;
static void (*IdleFunc)(void)                    = NULL;

static char     *lpszClassName = "tkLibWClass";
static WCHAR    *lpszClassNameW = L"tkLibWClass";

long FAR PASCAL _export tkWndProc(HWND hWnd, UINT message, DWORD wParam, LONG lParam);
static unsigned char ComponentFromIndex(int i, int nbits, int shift );
static void PrintMessage( const char *Format, ... );
//static PALETTEENTRY *FillRgbPaletteEntries( PIXELFORMATDESCRIPTOR *Pfd, PALETTEENTRY *Entries, UINT Count );
static HPALETTE CreateCIPalette( HDC Dc );
static HPALETTE CreateRGBPalette( HDC hdc );
static void DestroyThisWindow( HWND Window );
static void CleanUp( void );
static void DelayPaletteRealization( void );
static long RealizePaletteNow( HDC Dc, HPALETTE Palette);
static void ForceRedraw( HWND Window );
static void *AllocateMemory( size_t Size );
static void *AllocateZeroedMemory( size_t Size );
static void FreeMemory( void *Chunk );

/*
 *  Prototypes for the debugging functions go here
 */

#define DBGFUNC 0
#if DBGFUNC

static void DbgPrintf( const char *Format, ... );
static void pwi( void );
static void pwr(RECT *pr);
//static void ShowPixelFormat(HDC hdc);

#endif
#define NCOLORS 17
float tkRGBMap[NCOLORS][3] = {
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {0,0,0},
    {1,0,0},
    {0,1,0},
    {1,1,0},
    {0,0,1},
    {1,0,1},
    {0,1,1},
    {1,1,1}
};


/***************************************************************
 *                                                             *
 *  Exported Functions go here                                 *
 *                                                             *
 ***************************************************************/

void tkErrorPopups(GLboolean bEnable)
{
    tkPopupEnable = bEnable;
}

void tkCloseWindow(void)
{
    DestroyThisWindow(tkhwnd);
	/* if (w.cMain) {
	   XMesaDestroyContext(w.cMain);
	}*/
	WMesaDestroyContext(NULL);
}


void tkExec(void)
{
    MSG Message;

    /*
     *  WM_SIZE gets delivered before we get here!
     */

    if (ReshapeFunc)
    {
        RECT ClientRect;

        GetClientRect(tkhwnd, &ClientRect);
        (*ReshapeFunc)(ClientRect.right, ClientRect.bottom);
    }

    while (GL_TRUE)
    {
        /*
         *  Process all pending messages
         */

        while (PeekMessage(&Message, NULL, 0, 0, PM_NOREMOVE) == TRUE)
        {
            if (GetMessage(&Message, NULL, 0, 0) )
            {
                TranslateMessage(&Message);
                DispatchMessage(&Message);
            }
            else
            {
                /*
                 *  Nothing else to do here, just return
                 */

                return;
            }
        }

        /*
         *  If an idle function was defined, call it
         */

        if (IdleFunc)
        {
            (*IdleFunc)();
        }
    }
}

void tkExposeFunc(void (*Func)(int, int))
{
    ExposeFunc = Func;
}

void tkReshapeFunc(void (*Func)(GLsizei, GLsizei))
{
    ReshapeFunc = Func;
}

void tkDisplayFunc(void (*Func)(void))
{
    DisplayFunc = Func;
}

void tkKeyDownFunc(GLenum (*Func)(int, GLenum))
{
    KeyDownFunc = Func;
}

void tkMouseDownFunc(GLenum (*Func)(int, int, GLenum))
{
    MouseDownFunc = Func;
}

void tkMouseUpFunc(GLenum (*Func)(int, int, GLenum))
{
    MouseUpFunc = Func;
}

void tkMouseMoveFunc(GLenum (*Func)(int, int, GLenum))
{
    MouseMoveFunc = Func;
}

void tkIdleFunc(void (*Func)(void))
{
    IdleFunc = Func;
}

void tkInitPosition(int x, int y, int width, int height)
{
    if (x == CW_USEDEFAULT)
    {
        x = 0;
        y = 0;
        windInfo.bDefPos = TRUE;
    }
    else
        windInfo.bDefPos = FALSE;

    windInfo.x = x + GetSystemMetrics(SM_CXFRAME);
    windInfo.y = y + GetSystemMetrics(SM_CYCAPTION)
                 - GetSystemMetrics(SM_CYBORDER)
                 + GetSystemMetrics(SM_CYFRAME);
    windInfo.width = width;
    windInfo.height = height;
}

void tkInitDisplayMode(GLenum type)
{
    windInfo.type = type;
}

void tkInitDisplayModePolicy(GLenum type)
{
    windInfo.dmPolicy = type;
}

GLenum tkInitDisplayModeID(GLint ipfd)
{
    windInfo.ipfd = ipfd;
    return GL_TRUE;
}

GLenum tkInitWindowAW(char *title, BOOL bUnicode)
{
    WMesaContext Cur;
    WNDCLASS wndclass;
    RECT     WinRect;
    HANDLE   hInstance;
    ATOM     aRegister;
    GLenum   Result = GL_FALSE,RGB_Flag=GL_TRUE,DB_Flag=GL_FALSE;

    hInstance = GetModuleHandle(NULL);

    // Must not define CS_CS_PARENTDC style.
    wndclass.style         = CS_HREDRAW | CS_VREDRAW;
    wndclass.lpfnWndProc   = (WNDPROC)tkWndProc;
    wndclass.cbClsExtra    = 0;
    wndclass.cbWndExtra    = 0;
    wndclass.hInstance     = hInstance;
    wndclass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wndclass.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wndclass.hbrBackground = GetStockObject(BLACK_BRUSH);
    wndclass.lpszMenuName  = NULL;

    if (bUnicode)
        wndclass.lpszClassName = (LPCSTR)lpszClassNameW;
    else
        wndclass.lpszClassName = (LPCSTR)lpszClassName;

    if (bUnicode)
    {
        aRegister = RegisterClassW((CONST WNDCLASSW *)&wndclass);
    }
    else
    {
        aRegister = RegisterClass(&wndclass);
    }


    /*
     *  If the window failed to register, then there's no
     *  need to continue further.
     */

    if(0 == aRegister)
    {
        PrintMessage("Failed to register window class\n");
        return(Result);
    }


    /*
     *  Make window large enough to hold a client area as large as windInfo
     */

    WinRect.left   = windInfo.x;
    WinRect.right  = windInfo.x + windInfo.width;
    WinRect.top    = windInfo.y;
    WinRect.bottom = windInfo.y + windInfo.height;

    AdjustWindowRect(&WinRect, WS_OVERLAPPEDWINDOW, FALSE);

    /*
     *  Must use WS_CLIPCHILDREN and WS_CLIPSIBLINGS styles.
     */

    if (bUnicode)
    {
        tkhwnd = CreateWindowW(
                    (LPCWSTR)lpszClassNameW,
                    (LPCWSTR)title,
                    WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
                    (windInfo.bDefPos) ? CW_USEDEFAULT : WinRect.left,
                    (windInfo.bDefPos) ? CW_USEDEFAULT : WinRect.top,
                    WinRect.right - WinRect.left,
                    WinRect.bottom - WinRect.top,
                    NULL,
                    NULL,
                    hInstance,
                    NULL);
    }
    else
    {
        tkhwnd = CreateWindow(
                    lpszClassName,
                    title,
                    WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
                    (windInfo.bDefPos) ? CW_USEDEFAULT : WinRect.left,
                    (windInfo.bDefPos) ? CW_USEDEFAULT : WinRect.top,
                    WinRect.right - WinRect.left,
                    WinRect.bottom - WinRect.top,
                    NULL,
                    NULL,
                    hInstance,
                    NULL);
    }

    if ( NULL != tkhwnd )
    {
        // If default window positioning used, find out window position and fix
        // up the windInfo position info.

        if (windInfo.bDefPos)
        {
            GetWindowRect(tkhwnd, &WinRect);
            windInfo.x = WinRect.left + GetSystemMetrics(SM_CXFRAME);
            windInfo.y = WinRect.top  + GetSystemMetrics(SM_CYCAPTION)
                         - GetSystemMetrics(SM_CYBORDER)
                         + GetSystemMetrics(SM_CYFRAME);
        }

        tkhdc = GetDC(tkhwnd);

        if ( NULL != tkhdc )
            ShowWindow(tkhwnd, SW_SHOWDEFAULT);
        else
            PrintMessage("Could not get an HDC for window 0x%08lX\n", tkhwnd );
    }
    else
        PrintMessage("create window failed\n");
    if (windInfo.type & TK_INDEX)
    {
      RGB_Flag=GL_FALSE;
      tkSetRGBMap(NCOLORS,(float *) tkRGBMap);
    }
    if (windInfo.type & TK_DOUBLE)
      DB_Flag=GL_TRUE;
    Cur=WMesaCreateContext(tkhwnd,tkhpalette,RGB_Flag,DB_Flag);
    WMesaMakeCurrent(Cur);
    return GL_TRUE;
}

// Initialize a window, create a rendering context for that window
GLenum tkInitWindow(char *title)
{
    TKASSERT( NULL==tkhwnd      );
    TKASSERT( NULL==tkhdc       );
    TKASSERT( NULL==tkhrc       );
    TKASSERT( NULL==tkhpalette  );

    return tkInitWindowAW(title, FALSE);
}


/******************************************************************************/

/*
 * You cannot just call DestroyWindow() here.  The programs do not expect
 * tkQuit() to return;  DestroyWindow() just sends a WM_DESTROY message
 */

void tkQuit(void)
{
    DestroyThisWindow(tkhwnd);
    ExitProcess(0);
}

/******************************************************************************/

void tkSetOneColor(int index, float r, float g, float b)
{
    PALETTEENTRY PalEntry;
    HPALETTE Palette;
    if ( NULL != (Palette = CreateCIPalette( tkhdc )) )
    {
        PalEntry.peRed   = (BYTE)(r*(float)255.0 + (float)0.5);
        PalEntry.peGreen = (BYTE)(g*(float)255.0 + (float)0.5);
        PalEntry.peBlue  = (BYTE)(b*(float)255.0 + (float)0.5);
        PalEntry.peFlags = 0;
        SetPaletteEntries( Palette, index, 1, &PalEntry);
        DelayPaletteRealization();
    }
}

void tkSetFogRamp(int density, int startIndex)
{
    HPALETTE CurrentPal;
    PALETTEENTRY *pPalEntry;
    UINT n, i, j, k, intensity, fogValues, colorValues;

    if ( NULL != (CurrentPal = CreateCIPalette(tkhdc)) )
    {
        n = GetPaletteEntries( CurrentPal, 0, 0, NULL );

        pPalEntry = AllocateMemory( n * sizeof(PALETTEENTRY) );

        if ( NULL != pPalEntry)
        {
            fogValues = 1 << density;
            colorValues = 1 << startIndex;
            for (i = 0; i < colorValues; i++) {
                for (j = 0; j < fogValues; j++) {
                    k = i * fogValues + j;

                    intensity = i * fogValues + j * colorValues;
                    //mf: not sure what they're trying to do here
                    //intensity = (intensity << 8) | intensity; ???

                // This is a workaround for a GDI palette "feature".  If any of
                // the static colors are repeated in the palette, those colors
                // will map to the first occurrence.  So, for our case where there
                // are only two static colors (black and white), if a white
                // color appears anywhere in the palette other than in the last
                // entry, the static white will remap to the first white.  This
                // destroys the nice one-to-one mapping we are trying to achieve.
                //
                // There are two ways to workaround this.  The first is to
                // simply not allow a pure white anywhere but in the last entry.
                // Such requests are replaced with an attenuated white of
                // (0xFE, 0xFE, 0xFE).
                //
                // The other way is to mark these extra whites with PC_RESERVED
                // which will cause GDI to skip these entries when mapping colors.
                // This way the app gets the actual colors requested, but can
                // have side effects on other apps.
                //
                // Both solutions are included below.  The PC_RESERVED solution is
                // the one currently enabled.  It may have side effects, but taking
                // over the static colors as we are is a really big side effect that
                // should swamp out the effects of using PC_RESERVED.
                if (intensity > 0xFF)
                  intensity = 0xFF;
                pPalEntry[k].peRed =pPalEntry[k].peGreen = pPalEntry[k].peBlue = (BYTE) intensity;
                pPalEntry[k].peFlags = 0;

                }
            }

            SetPaletteEntries(CurrentPal, 0, n, pPalEntry);
            FreeMemory( pPalEntry );

            DelayPaletteRealization();
        }
    }
}

void tkSetGreyRamp(void)
{
    HPALETTE CurrentPal;
    PALETTEENTRY *Entries;
    UINT Count, i;
    float intensity;

    if ( NULL != (CurrentPal = CreateCIPalette( tkhdc )) )
    {
        Count   = GetPaletteEntries( CurrentPal, 0, 0, NULL );
        Entries = AllocateMemory( Count * sizeof(PALETTEENTRY) );

        if ( NULL != Entries )
        {
            for (i = 0; i < Count; i++)
            {
                intensity = (float)(((double)i / (double)(Count-1)) * (double)255.0 + (double)0.5);
                Entries[i].peRed =
                Entries[i].peGreen =
                Entries[i].peBlue = (BYTE) intensity;
                Entries[i].peFlags = 0;
            }
            SetPaletteEntries( CurrentPal, 0, Count, Entries );
            FreeMemory( Entries );

            DelayPaletteRealization();
        }
    }
}

void tkSetRGBMap( int Size, float *Values )
{
    HPALETTE CurrentPal;
    int i;
    if ( NULL != (CurrentPal = CreateCIPalette( tkhdc )) )
    {
      for (i=0; i<Size; i++)
        tkSetOneColor(i,Values[i*3],Values[i*3+1],Values[i*3+2]);
    }
}

/******************************************************************************/

void tkSwapBuffers(void)
{
  WMesaSwapBuffers();
}

/******************************************************************************/

GLint tkGetColorMapSize(void)
{
    CreateCIPalette( tkhdc );

    if ( NULL == tkhpalette )
        return( 0 );

    return( GetPaletteEntries( tkhpalette, 0, 0, NULL ) );
}

void tkGetMouseLoc(int *x, int *y)
{
    POINT Point;

    *x = 0;
    *y = 0;

    GetCursorPos(&Point);

    /*
     *  GetCursorPos returns screen coordinates,
     *  we want window coordinates
     */

    *x = Point.x - windInfo.x;
    *y = Point.y - windInfo.y;
}

HWND tkGetHWND(void)
{
    return tkhwnd;
}

HDC tkGetHDC(void)
{
    return tkhdc;
}
GLenum tkGetDisplayModePolicy(void)
{
    return windInfo.dmPolicy;
}

GLint tkGetDisplayModeID(void)
{
    return windInfo.ipfd;
}

GLenum tkGetDisplayMode(void)
{
    return windInfo.type;
}


/***********************************************************************
 *                                                                     *
 *  The Following functions are for our own use only. (ie static)      *
 *                                                                     *
 ***********************************************************************/

long FAR PASCAL _export
tkWndProc(HWND hWnd, UINT message, DWORD wParam, LONG lParam)
{
    int key;
    PAINTSTRUCT paint;
    HDC hdc;

    switch (message)
    {

    case WM_USER:

        if ( RealizePaletteNow( tkhdc, tkhpalette) > 0 )
            ForceRedraw( hWnd );
        return(0);

    case WM_SIZE:
        windInfo.width  = LOWORD(lParam);
        windInfo.height = HIWORD(lParam);

        if (ReshapeFunc)
        {
            (*ReshapeFunc)(windInfo.width, windInfo.height);

            ForceRedraw( hWnd );
        }
        return (0);

    case WM_MOVE:
        windInfo.x = LOWORD(lParam);
        windInfo.y = HIWORD(lParam);
        return (0);

    case WM_PAINT:

        /*
         *  Validate the region even if there are no DisplayFunc.
         *  Otherwise, USER will not stop sending WM_PAINT messages.
         */

        hdc = BeginPaint(tkhwnd, &paint);

        if (DisplayFunc)
        {
            (*DisplayFunc)();
        }

        EndPaint(tkhwnd, &paint);
        return (0);

    case WM_PALETTECHANGED:
        if ( hWnd != (HWND) wParam )
          RealizePaletteNow(tkhdc,tkhpalette);
        return (0);
    case WM_QUERYNEWPALETTE:

    // In the foreground!  Let RealizePaletteNow do the work--
    // if management of the static system color usage is needed,
    // RealizePaletteNow will take care of it.

        if ( NULL != tkhpalette )
        {
            if ( RealizePaletteNow(tkhdc, tkhpalette) > 0 )
                ForceRedraw( hWnd );

            return (1);
        }

        return (0);

    case WM_ACTIVATE:

    // If the window is going inactive, the palette must be realized to
    // the background.  Cannot depend on WM_PALETTECHANGED to be sent since
    // the window that comes to the foreground may or may not be palette
    // managed.

        if ( LOWORD(wParam) == WA_INACTIVE )
        {
            if ( NULL != tkhpalette )
            {
            // Realize as a background palette.  Need to call
            // RealizePaletteNow rather than RealizePalette directly to
            // because it may be necessary to release usage of the static
            // system colors.

                if ( RealizePaletteNow( tkhdc, tkhpalette) > 0 )
                    ForceRedraw( hWnd );
            }
        }

    // Allow DefWindowProc() to finish the default processing (which includes
    // changing the keyboard focus).

        break;

    case WM_MOUSEMOVE:

        if (MouseMoveFunc)
        {
            GLenum mask;

            mask = 0;
            if (wParam & MK_LBUTTON) {
                mask |= TK_LEFTBUTTON;
            }
            if (wParam & MK_MBUTTON) {
                mask |= TK_MIDDLEBUTTON;
            }
            if (wParam & MK_RBUTTON) {
                mask |= TK_RIGHTBUTTON;
            }

            if ((*MouseMoveFunc)( LOWORD(lParam), HIWORD(lParam), mask ))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_LBUTTONDOWN:

        SetCapture(hWnd);

        if (MouseDownFunc)
        {
            if ( (*MouseDownFunc)(LOWORD(lParam), HIWORD(lParam),
                 TK_LEFTBUTTON) )
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_LBUTTONUP:

        ReleaseCapture();

        if (MouseUpFunc)
        {
            if ((*MouseUpFunc)(LOWORD(lParam), HIWORD(lParam), TK_LEFTBUTTON))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_MBUTTONDOWN:

        SetCapture(hWnd);

        if (MouseDownFunc)
        {
            if ((*MouseDownFunc)(LOWORD(lParam), HIWORD(lParam),
                    TK_MIDDLEBUTTON))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_MBUTTONUP:

        ReleaseCapture();

        if (MouseUpFunc)
        {
            if ((*MouseUpFunc)(LOWORD(lParam), HIWORD(lParam),
                TK_MIDDLEBUTTON))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_RBUTTONDOWN:

        SetCapture(hWnd);

        if (MouseDownFunc)
        {
            if ((*MouseDownFunc)(LOWORD(lParam), HIWORD(lParam),
                TK_RIGHTBUTTON))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_RBUTTONUP:

        ReleaseCapture();

        if (MouseUpFunc)
        {
            if ((*MouseUpFunc)(LOWORD(lParam), HIWORD(lParam),
                TK_RIGHTBUTTON))
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_KEYDOWN:
        switch (wParam) {
        case VK_SPACE:          key = TK_SPACE;         break;
        case VK_RETURN:         key = TK_RETURN;        break;
        case VK_ESCAPE:         key = TK_ESCAPE;        break;
        case VK_LEFT:           key = TK_LEFT;          break;
        case VK_UP:             key = TK_UP;            break;
        case VK_RIGHT:          key = TK_RIGHT;         break;
        case VK_DOWN:           key = TK_DOWN;          break;
        default:                key = GL_FALSE;         break;
        }

        if (key && KeyDownFunc)
        {
            GLenum mask;

            mask = 0;
            if (GetKeyState(VK_CONTROL)) {
                mask |= TK_CONTROL;
            }

            if (GetKeyState(VK_SHIFT)) {

                mask |= TK_SHIFT;
            }

            if ( (*KeyDownFunc)(key, mask) )
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_CHAR:
        if (('0' <= wParam && wParam <= '9') ||
            ('a' <= wParam && wParam <= 'z') ||
            ('A' <= wParam && wParam <= 'Z')) {

            key = wParam;
        } else {
            key = GL_FALSE;
        }

        if (key && KeyDownFunc) {
            GLenum mask;

            mask = 0;

            if (GetKeyState(VK_CONTROL)) {
                mask |= TK_CONTROL;
            }

            if (GetKeyState(VK_SHIFT)) {
                mask |= TK_SHIFT;
            }

            if ( (*KeyDownFunc)(key, mask) )
            {
                ForceRedraw( hWnd );
            }
        }
        return (0);

    case WM_CLOSE:
        DestroyWindow(tkhwnd);
        return(0);

    case WM_DESTROY:
        CleanUp();
        PostQuitMessage(TRUE);
        return 0;
    }
    return(DefWindowProc( hWnd, message, wParam, lParam));
}
static HPALETTE CreateCIPalette( HDC Dc )
{
    LOGPALETTE *LogicalPalette;
    HPALETTE StockPalette;
    UINT PaletteSize, StockPaletteSize, EntriesToCopy;

    if ( (Dc != NULL) && (NULL == tkhpalette) )
    {
                PaletteSize = 256; //(Pfd.cColorBits >= 8) ? 256 : (1 << Pfd.cColorBits);

                LogicalPalette = AllocateZeroedMemory( sizeof(LOGPALETTE) +
                                        (PaletteSize * sizeof(PALETTEENTRY)) );

                if ( NULL != LogicalPalette )
                {
                    LogicalPalette->palVersion    = 0x300;
                    LogicalPalette->palNumEntries = PaletteSize;

                    StockPalette     = GetStockObject(DEFAULT_PALETTE);
                    StockPaletteSize = GetPaletteEntries( StockPalette, 0, 0, NULL );

                    /*
                     *  start by copying default palette into new one
                     */

                    EntriesToCopy = StockPaletteSize < PaletteSize ?
                                        StockPaletteSize : PaletteSize;

                    GetPaletteEntries( StockPalette, 0, EntriesToCopy,
                                        LogicalPalette->palPalEntry );

                    /*
                     *  If we are taking possession of the system colors,
                     *  must guarantee that 0 and 255 are black and white
                     *  (respectively).
                     */

                    tkhpalette = CreatePalette(LogicalPalette);

                    FreeMemory(LogicalPalette);

                    RealizePaletteNow( Dc, tkhpalette);
                }
            }
    return( tkhpalette );
}
static void
PrintMessage( const char *Format, ... )
{
    va_list ArgList;
    char Buffer[256];

    va_start(ArgList, Format);
    vsprintf(Buffer, Format, ArgList);
    va_end(ArgList);

    MessageBox(GetFocus(), Buffer, "Error", MB_OK);
}

static void
DelayPaletteRealization( void )
{
    MSG Message;

    TKASSERT(NULL!=tkhwnd);

    /*
     *  Add a WM_USER message to the queue, if there isn't one there already.
     */

    if (!PeekMessage(&Message, tkhwnd, WM_USER, WM_USER, PM_NOREMOVE) )
    {
        PostMessage( tkhwnd, WM_USER, 0, 0);
    }
}

/******************************Public*Routine******************************\
* RealizePaletteNow
*
* Select the given palette in background or foreground mode (as specified
* by the bForceBackground flag), and realize the palette.
*
* If static system color usage is set, the system colors are replaced.
*
* History:
*  26-Apr-1994 -by- Gilman Wong [gilmanw]
* Wrote it.
\**************************************************************************/

static long RealizePaletteNow( HDC Dc, HPALETTE Palette)
{
    long Result = -1;
    TKASSERT( NULL!=Dc      );
    TKASSERT( NULL!=Palette );
    if ( NULL != SelectPalette( Dc, Palette, FALSE ) )
    {
      Result = RealizePalette( Dc );
      WMesaPaletteChange(Palette);
    }
    return( Result );
}

static void
ForceRedraw( HWND Window )
{
    MSG Message;

    if (!PeekMessage(&Message, Window, WM_PAINT, WM_PAINT, PM_NOREMOVE) )
    {
        InvalidateRect( Window, NULL, FALSE );
    }
}
static void
DestroyThisWindow( HWND Window )
{
    if ( NULL != Window )
    {
        DestroyWindow( Window );
    }
}

/*
 *  This Should be called in response to a WM_DESTROY message
 */

static void
CleanUp( void )
{
    HPALETTE hStock;

// Cleanup the palette.

    if ( NULL != tkhpalette )
    {
    // If static system color usage is set, restore the system colors.

        if ((hStock = GetStockObject( DEFAULT_PALETTE ))!=NULL)
          SelectPalette( tkhdc, hStock, FALSE );
        DeleteObject( tkhpalette );
    }

// Cleanup the DC.

    if ( NULL != tkhdc )
        ReleaseDC( tkhwnd, tkhdc );
// Be really nice and reset global values.
    tkhwnd        = NULL;
    tkhdc         = NULL;
    tkhpalette    = NULL;

    ExposeFunc    = NULL;
    ReshapeFunc   = NULL;
    IdleFunc      = NULL;
    DisplayFunc   = NULL;
    KeyDownFunc   = NULL;
    MouseDownFunc = NULL;
    MouseUpFunc   = NULL;
    MouseMoveFunc = NULL;
}

static void *
AllocateMemory( size_t Size )
{
    return( LocalAlloc( LMEM_FIXED, Size ) );
}

static void *
AllocateZeroedMemory( size_t Size )
{
    return( LocalAlloc( LMEM_FIXED | LMEM_ZEROINIT, Size ) );
}


static void
FreeMemory( void *Chunk )
{
    TKASSERT( NULL!=Chunk );

    LocalFree( Chunk );
}
