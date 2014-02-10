/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/***********************************************************************
Program:    ELMER Front 
Module:     ecif_oglWIN32renderer.hpp
Language:   C++
Date:       08.04.97
Version:    1.00
Author(s):  Martti Verho
Revisions:  

Abstract:   Include type file. WIN specific parts for OGLrenderer class. 

************************************************************************/

#define WM_INIT WM_USER

// WIN32 codes:
//MK_LBUTTON          0x0001
//MK_RBUTTON          0x0002
//MK_SHIFT            0x0004
//MK_CONTROL          0x0008
//MK_MBUTTON          0x0010

void setDCPixelFormat(HDC hdc);

 
void
Renderer_OGL::createGLWindow(Hinst app_instance, const char* pName,
                             int xPos, int yPos, int pWidth, int pHeight,
                             WindowInfo& winfo)
{
  winfo.display = NULL;

  WNDCLASS wndclass;
  wndclass.lpszClassName = oglWinClass;
  wndclass.hInstance     = app_instance;
  wndclass.lpfnWndProc   = windowProcedure;
  wndclass.hCursor       = LoadCursor( NULL, IDC_ARROW );
  wndclass.hIcon         = NULL;
  wndclass.lpszMenuName  = NULL;
  wndclass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
  wndclass.style         = CS_HREDRAW | CS_VREDRAW | 
                           CS_BYTEALIGNCLIENT | CS_BYTEALIGNWINDOW |
                           CS_DBLCLKS;
  wndclass.cbClsExtra    = 0;
  wndclass.cbWndExtra    = sizeof(LONG);

  RegisterClass(&wndclass );

  HWND hWnd = CreateWindow( oglWinClass,
      pName,
      WS_OVERLAPPEDWINDOW |
      WS_CLIPCHILDREN |
      WS_CLIPSIBLINGS,
      xPos,//CW_USEDEFAULT,
      yPos,//CW_USEDEFAULT,
      pWidth,
      pHeight,
      HWND_DESKTOP,
      NULL,
      app_instance,
      NULL);

  winfo.window = hWnd;

  // Create GL-rendering contex and attach window to it.
  HDC hdc = GetDC(hWnd);
  setDCPixelFormat(hdc);
  HGLRC hrc = wglCreateContext(hdc);
  wglMakeCurrent(hdc, hrc);

}


//  First a helper functions.
void
Renderer_OGL::destroyWindow(Hdisplay* display, Hwindow window)
{
  DestroyWindow(window);
  status = HAS_NO_WINDOW;
  visible = FALSE;
}


// *** This must be defined, but not really used in WIN32 ***
void
Renderer_OGL::dummyWindowProc()
{
}


void
Renderer::findRendererWindowSize(int& w, int& h)
{
  RECT  rect;
  //int result = GetWindowRect(rendererWindow, &rect);
  int result = GetClientRect(rendererWindow, &rect);
  w = rect.right - rect.left;
  h = rect.bottom - rect.top ;
}


void 
Renderer_OGL::paintRenderer()
{ 
  SendMessage(rendererWindow, WM_INIT, 0, 0L);
  ShowWindow(rendererWindow, SW_SHOWNORMAL);
  //UpdateWindow(rendererWindow);
  Renderer::visible = TRUE;
  refresh();
}


//  Function sets the pixel format for a device context.
//  Needed in order to create a OpenGL-rendering context.
void setDCPixelFormat(HDC hdc)
{
  HANDLE       hHeap;
   int          nColors, i;
   LPLOGPALETTE lpPalette;
   BYTE         byRedMask, byGreenMask, byBlueMask;

   static PIXELFORMATDESCRIPTOR pfd = 
   {
        sizeof(PIXELFORMATDESCRIPTOR),  // Size of this structure
      1,                              // Version number
      PFD_DRAW_TO_WINDOW |            // Flags
      PFD_SUPPORT_OPENGL |
      PFD_DOUBLEBUFFER,
      PFD_TYPE_RGBA,                  // Use RGBA pixel values
      24,                             // Try to use 24-bit color
      0, 0, 0, 0, 0, 0,               // Don't care about these
      0, 0,                           // No alpha buffer
      32, 0, 0, 0, 0,                 // 32-bit accumulation buffer
      32,                             // 32-bit depth buffer
      0,                              // No stencil buffer
      0,                              // No auxiliary buffers
      PFD_MAIN_PLANE,                 // Layer type
      0,                              // Reserved (must be 0)
      0, 0, 0                         // No layer masks
   };

   int nPixelFormat;
    
   nPixelFormat = ChoosePixelFormat(hdc, &pfd);
   SetPixelFormat(hdc, nPixelFormat, &pfd);
  
   DescribePixelFormat(hdc, nPixelFormat, 
                        sizeof(PIXELFORMATDESCRIPTOR),
                        &pfd);

   if(pfd.dwFlags & PFD_NEED_PALETTE) 
   {
    nColors = 1 << pfd.cColorBits;
      hHeap = GetProcessHeap();

      lpPalette = (LPLOGPALETTE) HeapAlloc(hHeap, 0,
                   sizeof(LOGPALETTE) + (nColors * 
                     sizeof(PALETTEENTRY)));

      lpPalette->palVersion = 0x300;
      lpPalette->palNumEntries = nColors;

      byRedMask = (1 << pfd.cRedBits) - 1;
      byGreenMask = (1 << pfd.cGreenBits) - 1;
      byBlueMask = (1 << pfd.cBlueBits) - 1;

      for(i = 0; i < nColors; i++) 
      {
        lpPalette->palPalEntry[i].peRed =
          (((i >> pfd.cRedShift) & byRedMask) * 255) / 
            byRedMask;
         lpPalette->palPalEntry[i].peGreen =
          (((i >> pfd.cGreenShift) & byGreenMask) * 255) /
            byGreenMask;
         lpPalette->palPalEntry[i].peBlue =
          (((i >> pfd.cBlueShift) & byBlueMask) * 255) / 
            byBlueMask;
         lpPalette->palPalEntry[i].peFlags = 0;
    }

    // Create the palette and free the allocated memory
      HPALETTE m_hPalette = CreatePalette(lpPalette);
      HeapFree(hHeap, 0, lpPalette);

    // Realize the color palette
      if(m_hPalette != NULL) 
      {
        SelectPalette(hdc, m_hPalette, FALSE);
         RealizePalette(hdc);
      }
  }
}


void
Renderer_OGL::setWindowTitle(char* title)
{
  // Call Win32 API function
  SetWindowText(rendererWindow, title);
}


void
Renderer_OGL::swapBuffers()
{
  SwapBuffers(wglGetCurrentDC());
}


// *** WIN32 WindowProc for the OGLRenderer ***
LresCallback 
Renderer_OGL::windowProcedure(Hwindow hWnd, Event wMsg,
                    Wparam wParam, Lparam lParam)
{
  //RECT  rect;
  ////bool result = GetWindowRect(rendererWindow, &rect);
  //bool result = GetClientRect(hWnd, &rect);
  //int x_size = rect.right - rect.left;
  //int y_size = rect.bottom - rect.top;

  Renderer* renderer = theControlCenter->getRenderer();

  if (renderer == NULL)
    return DefWindowProc(hWnd, wMsg, wParam, lParam);

  int x_size, y_size;
  renderer->getRendererWindowSize(x_size, y_size);

  HDC hdc;
  PAINTSTRUCT ps;
  static int xPos1 = -1;
  static int yPos1 = -1;
  int xPos2, yPos2;
  int keyCode = -1;

  static bool ctrlPressed = false;

  switch (wMsg) {

    case WM_INIT: 
      break;

    case WM_PAINT:
      hdc = BeginPaint(hWnd, &ps);
      display(renderer);
      EndPaint(hWnd, &ps);
      return 0;

    case WM_SIZE:
       renderer->reshape();
       break;

    case WM_MOVE:
    case WM_MOVING:
      //theUI->generateEvent();
      break;

    case WM_MOUSEMOVE:
      if (  wParam != 0 && xPos1 != -1 ) {
        xPos2 = LOWORD(lParam);
        yPos2 = y_size - HIWORD(lParam);
        mouse_move_action(int(wParam), xPos1, xPos2, yPos1, yPos2);
        xPos1 = xPos2;
        yPos1 = yPos2;
      } else if ( wParam != 0 ) {
        xPos1 = LOWORD(lParam);
        yPos1 = y_size - HIWORD(lParam);
      }
      break;

    case WM_LBUTTONDOWN:
    case WM_RBUTTONDOWN:
      xPos1 = LOWORD(lParam);
      yPos1 = y_size - HIWORD(lParam);
      if ( wParam & (MK_SHIFT | MK_CONTROL) ) {
        mouse_clck_action(int(wParam), xPos1, yPos1);
      }
      break;

    case WM_MBUTTONDOWN:
      xPos1 = LOWORD(lParam);
      yPos1 = y_size - HIWORD(lParam);
      break;

    case WM_LBUTTONUP: 
    case WM_RBUTTONUP: 
    case WM_MBUTTONUP: 
      xPos1 = -1;
      yPos1 = -1;
      break;
     
    case WM_LBUTTONDBLCLK:
      xPos1 = LOWORD(lParam);
      yPos1 = y_size - HIWORD(lParam);
      mouse_dblclck_action(int(wParam), xPos1, yPos1);
      break;

    case WM_KEYDOWN:

      keyCode = wParam;

      // Ctrl pressed
      if ( keyCode == VK_CONTROL ) {
        ctrlPressed = true;
      }

      // A letter pressed while Ctrl-key is down
      if (ctrlPressed) {

        switch (keyCode) {

        case 'b':
        case 'B': theControlCenter->handleKeyAction(KEY_CTRL_B); break;
        case 'h':
        case 'H': theControlCenter->handleKeyAction(KEY_CTRL_H); break;
        case 'l':
        case 'L': theControlCenter->handleKeyAction(KEY_CTRL_L); break;
        case 'r':
        case 'R': theControlCenter->handleKeyAction(KEY_CTRL_R); break;
        case 'x':
        case 'X': theControlCenter->handleKeyAction(KEY_CTRL_X); break;
        case 'y':
        case 'Y': theControlCenter->handleKeyAction(KEY_CTRL_Y); break;
        case 'z':
        case 'Z': theControlCenter->handleKeyAction(KEY_CTRL_Z); break;
        }
      }

      return 0;

    case WM_KEYUP:

      keyCode = wParam;

      // Ctrl released
      if ( keyCode == VK_CONTROL ) {
        ctrlPressed = false;
      }

      return 0;

    case WM_CLOSE:
      theControlCenter->rendererIsClosed();
      break;

    case WM_DESTROY: {
      // Possible GL-rendering context is released
      Hglrc hrc = wglGetCurrentContext();
      if (hrc) {
        Hdc hdc = wglGetCurrentDC();
        wglMakeCurrent(NULL, NULL);
        wglDeleteContext(hrc);
        ReleaseDC(hWnd, hdc);
      }
      //status  = HAS_NO_WINDOW;
    }
      break;

    case WM_QUIT:
    break;
    default:
      break;
 
  }

  return DefWindowProc(hWnd, wMsg, wParam, lParam);
}
