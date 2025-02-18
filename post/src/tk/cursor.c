#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tk.h"
#include "private.h"

/******************************************************************************/

#define MAX_CURSOR 32

typedef struct _cursorRec {
    GLint id;
    Cursor cursor;
} cursorRec;

int cursorNum = 0;
cursorRec cursors[MAX_CURSOR];

/******************************************************************************/

void tkNewCursor(GLint id, GLubyte *shapeBuf, GLubyte *maskBuf, GLenum fgColor,
		 GLenum bgColor, GLint hotX, GLint hotY)
{
    GLubyte buf[32];
    Pixmap shapeMap, maskMap;
    XColor c1, c2;
    int i;

    if (cursorNum == MAX_CURSOR-1) {
	return;
    }

    for (i = 0; i < 32; i += 2) {
	buf[i] = shapeBuf[i+1];
	buf[i+1] = shapeBuf[i];
    }
    shapeMap = XCreatePixmapFromBitmapData(xDisplay, wRoot, (char*) buf,
					   16, 16, 1, 0, 1);
    for (i = 0; i < 32; i += 2) {
	buf[i] = maskBuf[i+1];
	buf[i+1] = maskBuf[i];
    }
    maskMap = XCreatePixmapFromBitmapData(xDisplay, wRoot, (char*) buf,
					  16, 16, 1, 0, 1);
    c1.red = (unsigned short)(tkRGBMap[fgColor][0] * 65535.0 + 0.5);
    c1.green = (unsigned short)(tkRGBMap[fgColor][1] * 65535.0 + 0.5);
    c1.blue = (unsigned short)(tkRGBMap[fgColor][2] * 65535.0 + 0.5);
    c1.flags = DoRed | DoGreen | DoBlue;
    c2.red = (unsigned short)(tkRGBMap[bgColor][0] * 65535.0 + 0.5);
    c2.green = (unsigned short)(tkRGBMap[bgColor][1] * 65535.0 + 0.5);
    c2.blue = (unsigned short)(tkRGBMap[bgColor][2] * 65535.0 + 0.5);
    c2.flags = DoRed | DoGreen | DoBlue;

    cursors[cursorNum].id = id;
    cursors[cursorNum].cursor = XCreatePixmapCursor(xDisplay, shapeMap, maskMap,
					            &c1, &c2, hotX, hotY);
    cursorNum++;
}

/******************************************************************************/

void tkSetCursor(GLint id)
{
    int i;

    for (i = 0; i < cursorNum; i++) {
	if (cursors[i].id == id) {
	    XDefineCursor(xDisplay, w.wMain, cursors[i].cursor);
	}
    }
}

/******************************************************************************/
