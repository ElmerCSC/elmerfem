#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tk.h"
#include "private.h"


/*
 * Return the index of the first bit set.  I'm using this instead of the
 * ffs() function because I'm not sure it's really portable.
 */
static int my_ffs( unsigned int n )
{
   if (n==0) {
      return 0;
   }
   else {
      unsigned int m = 1;
      int index = 1;

      while ((n&m)==0) {
	 m <<= 1;
	 index++;
      }
      return index;
   }
}


/******************************************************************************/

int tkGetColorMapSize(void)
{

    if (!xDisplay) {
	return 0;
    } else {
	return w.vInfoMain->colormap_size;
    }
}

/******************************************************************************/

void tkGetMouseLoc(int *x, int *y)
{
    int junk;

    *x = 0;
    *y = 0;
    XQueryPointer(xDisplay, w.wMain, (Window *)&junk, (Window *)&junk,
		  &junk, &junk, x, y, (unsigned int *)&junk);
}

/******************************************************************************/

void tkGetSystem(TKenum type, void *ptr)
{

    switch (type) {
      case TK_X_DISPLAY:
	*(Display **)ptr = xDisplay;
	break;
      case TK_X_WINDOW:
	*(Window *)ptr = w.wMain;
	break;
    }
}

/******************************************************************************/

void tkSetFogRamp(int density, int startIndex)
{
    XColor c[4096];
    int rShift, gShift, bShift, intensity, fogValues, colorValues;
    int i, j, k;

    switch (w.vInfoMain->class) {
      case DirectColor:
	fogValues = 1 << density;
	colorValues = 1 << startIndex;
	for (i = 0; i < colorValues; i++) {
	    for (j = 0; j < fogValues; j++) {
		k = i * fogValues + j;
		intensity = i * fogValues + j * colorValues;
		if (intensity > w.vInfoMain->colormap_size) {
		    intensity = w.vInfoMain->colormap_size;
		}
		intensity = (intensity << 8) | intensity;
		rShift = my_ffs((unsigned int)w.vInfoMain->red_mask) - 1;
		gShift = my_ffs((unsigned int)w.vInfoMain->green_mask) - 1;
		bShift = my_ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
		c[k].pixel = ((k << rShift) & w.vInfoMain->red_mask) |
			     ((k << gShift) & w.vInfoMain->green_mask) |
			     ((k << bShift) & w.vInfoMain->blue_mask);
		c[k].red = (unsigned short)intensity;
		c[k].green = (unsigned short)intensity;
		c[k].blue = (unsigned short)intensity;
		c[k].flags = DoRed | DoGreen | DoBlue;
	    }
	}
	XStoreColors(xDisplay, w.cMapMain, c, w.vInfoMain->colormap_size);
	break;
      case GrayScale:
      case PseudoColor:
	fogValues = 1 << density;
	colorValues = 1 << startIndex;
	for (i = 0; i < colorValues; i++) {
	    for (j = 0; j < fogValues; j++) {
		k = i * fogValues + j;
		intensity = i * fogValues + j * colorValues;
		if (intensity > w.vInfoMain->colormap_size) {
		    intensity = w.vInfoMain->colormap_size;
		}
		intensity = (intensity << 8) | intensity;
		c[k].pixel = k;
		c[k].red = (unsigned short)intensity;
		c[k].green = (unsigned short)intensity;
		c[k].blue = (unsigned short)intensity;
		c[k].flags = DoRed | DoGreen | DoBlue;
	    }
	}
	XStoreColors(xDisplay, w.cMapMain, c, w.vInfoMain->colormap_size);
	break;
    }

    XSync(xDisplay, 0);
}

/******************************************************************************/

void tkSetGreyRamp(void)
{
    XColor c[4096];
    float intensity;
    int rShift, gShift, bShift, i;

    switch (w.vInfoMain->class) {
      case DirectColor:
	for (i = 0; i < w.vInfoMain->colormap_size; i++) {
	    intensity = (float)i / (float)w.vInfoMain->colormap_size *
			65535.0 + 0.5;
	    rShift = my_ffs((unsigned int)w.vInfoMain->red_mask) - 1;
	    gShift = my_ffs((unsigned int)w.vInfoMain->green_mask) - 1;
	    bShift = my_ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
	    c[i].pixel = ((i << rShift) & w.vInfoMain->red_mask) |
			 ((i << gShift) & w.vInfoMain->green_mask) |
			 ((i << bShift) & w.vInfoMain->blue_mask);
	    c[i].red = (unsigned short)intensity;
	    c[i].green = (unsigned short)intensity;
	    c[i].blue = (unsigned short)intensity;
	    c[i].flags = DoRed | DoGreen | DoBlue;
	}
	XStoreColors(xDisplay, w.cMapMain, c, w.vInfoMain->colormap_size);
	break;
      case GrayScale:
      case PseudoColor:
	for (i = 0; i < w.vInfoMain->colormap_size; i++) {
	    intensity = (float)i / (float)w.vInfoMain->colormap_size *
			65535.0 + 0.5;
	    c[i].pixel = i;
	    c[i].red = (unsigned short)intensity;
	    c[i].green = (unsigned short)intensity;
	    c[i].blue = (unsigned short)intensity;
	    c[i].flags = DoRed | DoGreen | DoBlue;
	}
	XStoreColors(xDisplay, w.cMapMain, c, w.vInfoMain->colormap_size);
	break;
    }

    XSync(xDisplay, 0);
}

/******************************************************************************/

void tkSetOneColor(int index, float r, float g, float b)
{
    XColor c;
    int rShift, gShift, bShift;

    switch (w.vInfoMain->class) {
      case DirectColor:
	rShift = my_ffs((unsigned int)w.vInfoMain->red_mask) - 1;
	gShift = my_ffs((unsigned int)w.vInfoMain->green_mask) - 1;
	bShift = my_ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
	c.pixel = ((index << rShift) & w.vInfoMain->red_mask) |
		  ((index << gShift) & w.vInfoMain->green_mask) |
		  ((index << bShift) & w.vInfoMain->blue_mask);
	c.red = (unsigned short)(r * 65535.0 + 0.5);
	c.green = (unsigned short)(g * 65535.0 + 0.5);
	c.blue = (unsigned short)(b * 65535.0 + 0.5);
	c.flags = DoRed | DoGreen | DoBlue;
	XStoreColor(xDisplay, w.cMapMain, &c);
	break;
      case GrayScale:
      case PseudoColor:
	if (index < w.vInfoMain->colormap_size) {
	    c.pixel = index;
	    c.red = (unsigned short)(r * 65535.0 + 0.5);
	    c.green = (unsigned short)(g * 65535.0 + 0.5);
	    c.blue = (unsigned short)(b * 65535.0 + 0.5);
	    c.flags = DoRed | DoGreen | DoBlue;
	    XStoreColor(xDisplay, w.cMapMain, &c);
	}
	break;
    }

    XSync(xDisplay, 0);
}

/******************************************************************************/

void tkSetOverlayMap(int size, float *rgb)
{
    XColor c;
    unsigned long *buf;
    int max, i;

    if (w.vInfoOverlay->class == PseudoColor) {
	max = (size > w.vInfoOverlay->colormap_size) ?
	      w.vInfoOverlay->colormap_size : size;
	buf = (unsigned long *)calloc(max, sizeof(unsigned long));
	XAllocColorCells(xDisplay, w.cMapOverlay, True, NULL, 0, buf, max);
	for (i = 1; i < max; i++) {
	    c.pixel = i;
	    c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
	    c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
	    c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
	    c.flags = DoRed | DoGreen | DoBlue;
	    XStoreColor(xDisplay, w.cMapOverlay, &c);
	}
	free(buf);
    }

    XSync(xDisplay, 0);
}

/******************************************************************************/

void tkSetRGBMap(int size, float *rgb)
{
    XColor c;
    int rShift, gShift, bShift, max, i;

    switch (w.vInfoMain->class) {
      case DirectColor:
	max = (size > w.vInfoMain->colormap_size) ? w.vInfoMain->colormap_size
						  : size;
	rShift = my_ffs((unsigned int)w.vInfoMain->red_mask) - 1;
	gShift = my_ffs((unsigned int)w.vInfoMain->green_mask) - 1;
	bShift = my_ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
	c.flags = DoRed | DoGreen | DoBlue;
	for (i = 0; i < max; i++) {
	    c.pixel = ((i << rShift) & w.vInfoMain->red_mask) |
		      ((i << gShift) & w.vInfoMain->green_mask) |
		      ((i << bShift) & w.vInfoMain->blue_mask);
	    c.red = c.green = c.blue = (unsigned short) (i * 65535 / (max-1));
	    XStoreColor(xDisplay, w.cMapMain, &c);
	}
	break;
      case GrayScale:
      case PseudoColor:
	max = (size > w.vInfoMain->colormap_size) ? w.vInfoMain->colormap_size
						  : size;
	for (i = 0; i < max; i++) {
	    c.pixel = i;
	    c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
	    c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
	    c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
	    c.flags = DoRed | DoGreen | DoBlue;
	    XStoreColor(xDisplay, w.cMapMain, &c);
	}
	break;
    }

    XSync(xDisplay, 0);
}

/******************************************************************************/

GLenum tkSetWindowLevel(GLenum level)
{

    switch (level) {
      case TK_OVERLAY:
	if (TK_HAS_OVERLAY(w.type)) {
	    if (!glXMakeCurrent(xDisplay, w.wOverlay, w.cOverlay)) {
		return GL_FALSE;
	    }
	} else {
	    return GL_FALSE;
	}
	break;
      case TK_RGB:
      case TK_INDEX:
	if (!glXMakeCurrent(xDisplay, w.wMain, w.cMain)) {
	    return GL_FALSE;
	}
	break;
    }
    return GL_TRUE;
}

/******************************************************************************/
