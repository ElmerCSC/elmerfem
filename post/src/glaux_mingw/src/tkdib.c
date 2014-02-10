/******************************Module*Header*******************************\
* Module Name: dib.c
*
* (Brief description)
*
* Created: 21-Feb-1994 23:12:58
* Author: Gilman Wong [gilmanw]
*
* Copyright (c) 1994 Microsoft Corporation
*
* (General description of its use)
*
* Dependencies:
*
*   (#defines)
*   (#includes)
*
\**************************************************************************/

#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include "tk.h"

#define static

#define BFT_BITMAP  0x4d42  // 'BM' -- indicates structure is BITMAPFILEHEADER

// struct BITMAPFILEHEADER {
//      WORD  bfType
//      DWORD bfSize
//      WORD  bfReserved1
//      WORD  bfReserved2
//      DWORD bfOffBits
// }
#define OFFSET_bfType       0
#define OFFSET_bfSize       2
#define OFFSET_bfReserved1  6
#define OFFSET_bfReserved2  8
#define OFFSET_bfOffBits    10
#define SIZEOF_BITMAPFILEHEADER 14

// Read a WORD-aligned DWORD.  Needed because BITMAPFILEHEADER has
// WORD-alignment.
#define READDWORD(pv)   ( (DWORD)((PWORD)(pv))[0]               \
                          | ((DWORD)((PWORD)(pv))[1] << 16) )   \

// Computes the number of BYTES needed to contain n number of bits.
#define BITS2BYTES(n)   ( ((n) + 7) >> 3 )

/****************************************************************************
 *                                                                          *
 *  FUNCTION   : DibNumColors(VOID FAR * pv)                                *
 *                                                                          *
 *  PURPOSE    : Determines the number of colors in the DIB by looking at   *
 *               the BitCount filed in the info block.                      *
 *                                                                          *
 *  RETURNS    : The number of colors in the DIB.                           *
 *                                                                          *
 * Stolen from SDK ShowDIB example.                                         *
 ****************************************************************************/

WORD DibNumColors(VOID FAR * pv)
{
    WORD                bits;
    BITMAPINFOHEADER UNALIGNED *lpbi;
    BITMAPCOREHEADER UNALIGNED *lpbc;

    lpbi = ((LPBITMAPINFOHEADER)pv);
    lpbc = ((LPBITMAPCOREHEADER)pv);

    /*  With the BITMAPINFO format headers, the size of the palette
     *  is in biClrUsed, whereas in the BITMAPCORE - style headers, it
     *  is dependent on the bits per pixel ( = 2 raised to the power of
     *  bits/pixel).
     *
     *  Because of the way we use this call, BITMAPINFOHEADER may be out
     *  of alignment if it follows a BITMAPFILEHEADER.  So use the macro
     *  to safely access DWORD fields.
     */
    if (READDWORD(&lpbi->biSize) != sizeof(BITMAPCOREHEADER)){
        if (READDWORD(&lpbi->biClrUsed) != 0)
        {
            return (WORD) READDWORD(&lpbi->biClrUsed);
        }
        bits = lpbi->biBitCount;
    }
    else
        bits = lpbc->bcBitCount;

    switch (bits){
        case 1:
            return 2;
        case 4:
            return 16;
        case 8:
            return 256;
        default:
            /* A 24 bitcount DIB has no color table */
            return 0;
    }
}

/******************************Public*Routine******************************\
* tkDIBImageLoad
*
* ANSI version stub.  Only here for orthogonality with tkRGBImageLoad.
*
* History:
*  22-Feb-1994 -by- Gilman Wong [gilmanw]
* Wrote it.
\**************************************************************************/

TK_RGBImageRec *tkDIBImageLoadAW(char *fileName, BOOL bUnicode);

TK_RGBImageRec *tkDIBImageLoad(char *fileName)
{
    return tkDIBImageLoadAW(fileName, FALSE);
}

/******************************Public*Routine******************************\
* tkDIBImageLoadAW
*
* Loads a DIB file (specified as either an ANSI or Unicode filename,
* depending on the bUnicode flag) and converts it into a TK image format.
*
* The technique used is based on CreateDIBSection and SetDIBits.
* CreateDIBSection is used to create a DIB with a format easily converted
* into the TK image format (packed 24BPP RGB).  The only conversion 
* required is swapping R and B in each RGB triplet (see history below)
* The resulting bitmap is selected into a memory DC.
*
* The DIB file is mapped into memory and SetDIBits called to initialize
* the memory DC bitmap.  It is during this step that GDI converts the
* arbitrary DIB file format to RGB format.
*
* Finally, the RGB data in the DIB section is read out and repacked
* as 24BPP 'BGR'.
*
* Returns:
*   Pointer to TK_RGBImageRec.  If an error occurs, a diagnostic error
*   message is put into the error stream and tkQuit() is called,
*   terminating the app.
*
* History:
*  - 22-Feb-1994 -by- Gilman Wong [gilmanw]
*    Wrote it.
*
*  - 01-May-1995 : [marcfo]
*    Don't quit if can't open DIB file - return NULL.
*
*  - 27-Jul-1995 : [marcfo]
*    Changed CreateDIBSection to create an RGB mapping, so this would work
*    on win95.  Swapping of R and B required, to generate BGR (R=low byte)
*    for GL.  If the GL 'BGR' (R=high byte) extension becomes part of api
*    in OpenGL 2.0, then we can avoid swap and do a memcpy to grab the bits.
*
\**************************************************************************/

TK_RGBImageRec *tkDIBImageLoadAW(char *fileName, BOOL bUnicode)
{
    TK_RGBImageRec *final = (TK_RGBImageRec *) NULL; // Ptr to TK image struct
                                                     // to return.  Non-NULL
                                                     // only for success.

    WORD             wNumColors;    // Number of colors in color table
    BITMAPFILEHEADER *pbmf;         // Ptr to file header
    BITMAPINFOHEADER UNALIGNED *pbmihFile;
    BITMAPCOREHEADER UNALIGNED *pbmchFile; // Ptr to file's core header (if it exists)
    PVOID            pvBitsFile;    // Ptr to bitmap bits in file
    PBYTE            pjBitsRGB;     // Ptr to 24BPP RGB image in DIB section
    PBYTE            pjTKBits = (PBYTE) NULL;   // Ptr to final TK image bits
    PBYTE            pjSrc;         // Ptr to image file used for conversion
    PBYTE            pjDst;         // Ptr to TK image used for conversion

    // These need to be cleaned up when we exit:
    HANDLE     hFile = INVALID_HANDLE_VALUE;        // File handle
    HANDLE     hMap = (HANDLE) NULL;                // Mapping object handle
    PVOID      pvFile = (PVOID) NULL;               // Ptr to mapped file
    HDC        hdcMem = (HDC) NULL;                 // 24BPP mem DC
    HBITMAP    hbmRGB = (HBITMAP) NULL;             // 24BPP RGB bitmap
    BITMAPINFO *pbmiSource = (BITMAPINFO *) NULL;   // Ptr to source BITMAPINFO
    BITMAPINFO *pbmiRGB = (BITMAPINFO *) NULL;      // Ptr to file's BITMAPINFO

    int i, j;
    int padBytes;

// Map the DIB file into memory.

    hFile = bUnicode ? 
            CreateFileW((LPWSTR) fileName, GENERIC_READ, FILE_SHARE_READ, NULL,
                        OPEN_EXISTING, 0, 0) :
            CreateFileA((LPSTR) fileName, GENERIC_READ, FILE_SHARE_READ, NULL, 
                        OPEN_EXISTING, 0, 0);

    if (hFile == INVALID_HANDLE_VALUE)
        goto tkDIBLoadImage_cleanup;

    hMap = CreateFileMappingA(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!hMap)
        goto tkDIBLoadImage_cleanup;

    pvFile = MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, 0);
    if (!pvFile)
        goto tkDIBLoadImage_cleanup;

// Check the file header.  If the BFT_BITMAP magic number is there,
// then the file format is a BITMAPFILEHEADER followed immediately
// by either a BITMAPINFOHEADER or a BITMAPCOREHEADER.  The bitmap
// bits, in this case, are located at the offset bfOffBits from the
// BITMAPFILEHEADER.
//
// Otherwise, this may be a raw BITMAPINFOHEADER or BITMAPCOREHEADER
// followed immediately with the color table and the bitmap bits.

    pbmf = (BITMAPFILEHEADER *) pvFile;

    if ( pbmf->bfType == BFT_BITMAP )
    {
        pbmihFile = (BITMAPINFOHEADER *) ((PBYTE) pbmf + SIZEOF_BITMAPFILEHEADER);

    // BITMAPFILEHEADER is WORD aligned, so use safe macro to read DWORD
    // bfOffBits field.

        pvBitsFile = (PVOID *) ((PBYTE) pbmf
                                + READDWORD((PBYTE) pbmf + OFFSET_bfOffBits));
    }
    else
    {
        pbmihFile = (BITMAPINFOHEADER *) pvFile;

    // Determination of where the bitmaps bits are needs to wait until we
    // know for sure whether we have a BITMAPINFOHEADER or a BITMAPCOREHEADER.
    }

// Determine the number of colors in the DIB palette.  This is non-zero
// only for 8BPP or less.

    wNumColors = DibNumColors(pbmihFile);

// Create a BITMAPINFO (with color table) for the DIB file.  Because the
// file may not have one (BITMAPCORE case) and potential alignment problems,
// we will create a new one in memory we allocate.
//
// We distinguish between BITMAPINFO and BITMAPCORE cases based upon
// BITMAPINFOHEADER.biSize.

    pbmiSource = (BITMAPINFO *)
        LocalAlloc(LMEM_FIXED, sizeof(BITMAPINFO)
                               + wNumColors * sizeof(RGBQUAD));
    if (!pbmiSource)
    {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

    // Note: need to use safe READDWORD macro because pbmihFile may
    // have only WORD alignment if it follows a BITMAPFILEHEADER.

    switch (READDWORD(&pbmihFile->biSize))
    {
    case sizeof(BITMAPINFOHEADER):

    // Convert WORD-aligned BITMAPINFOHEADER to aligned BITMAPINFO.

        pbmiSource->bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
        pbmiSource->bmiHeader.biWidth         = READDWORD(&pbmihFile->biWidth);
        pbmiSource->bmiHeader.biHeight        = READDWORD(&pbmihFile->biHeight);
        pbmiSource->bmiHeader.biPlanes        = pbmihFile->biPlanes;
        pbmiSource->bmiHeader.biBitCount      = pbmihFile->biBitCount;
        pbmiSource->bmiHeader.biCompression   = 
                                        READDWORD(&pbmihFile->biCompression);
        pbmiSource->bmiHeader.biSizeImage     = 
                                        READDWORD(&pbmihFile->biSizeImage);
        pbmiSource->bmiHeader.biXPelsPerMeter = 
                                        READDWORD(&pbmihFile->biXPelsPerMeter);
        pbmiSource->bmiHeader.biYPelsPerMeter = 
                                        READDWORD(&pbmihFile->biYPelsPerMeter);
        pbmiSource->bmiHeader.biClrUsed       = 
                                        READDWORD(&pbmihFile->biClrUsed);
        pbmiSource->bmiHeader.biClrImportant  = 
                                        READDWORD(&pbmihFile->biClrImportant);

    // Copy color table.  It immediately follows the BITMAPINFOHEADER.

        memcpy((PVOID) &pbmiSource->bmiColors[0], (PVOID) (pbmihFile + 1),
               wNumColors * sizeof(RGBQUAD));

    // If we haven't already determined the position of the image bits,
    // we may now assume that they immediately follow the color table.

        if (!pvBitsFile)
            pvBitsFile = (PVOID) ((PBYTE) (pbmihFile + 1)
                         + wNumColors * sizeof(RGBQUAD));
        break;

    case sizeof(BITMAPCOREHEADER):
        pbmchFile = (BITMAPCOREHEADER *) pbmihFile;

    // Convert BITMAPCOREHEADER to BITMAPINFOHEADER.

        pbmiSource->bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
        pbmiSource->bmiHeader.biWidth         = (DWORD) pbmchFile->bcWidth;
        pbmiSource->bmiHeader.biHeight        = (DWORD) pbmchFile->bcHeight;
        pbmiSource->bmiHeader.biPlanes        = pbmchFile->bcPlanes;
        pbmiSource->bmiHeader.biBitCount      = pbmchFile->bcBitCount;
        pbmiSource->bmiHeader.biCompression   = BI_RGB;
        pbmiSource->bmiHeader.biSizeImage     = 0;
        pbmiSource->bmiHeader.biXPelsPerMeter = 0;
        pbmiSource->bmiHeader.biYPelsPerMeter = 0;
        pbmiSource->bmiHeader.biClrUsed       = wNumColors;
        pbmiSource->bmiHeader.biClrImportant  = wNumColors;

    // Convert RGBTRIPLE color table into RGBQUAD color table.

        {
            RGBQUAD *rgb4 = pbmiSource->bmiColors;
            RGBTRIPLE *rgb3 = (RGBTRIPLE *) (pbmchFile + 1);

            for (i = 0; i < wNumColors; i++)
            {
                rgb4->rgbRed   = rgb3->rgbtRed  ;
                rgb4->rgbGreen = rgb3->rgbtGreen;
                rgb4->rgbBlue  = rgb3->rgbtBlue ;
                rgb4->rgbReserved = 0;

                rgb4++;
                rgb3++;
            }
        }

    // If we haven't already determined the position of the image bits,
    // we may now assume that they immediately follow the color table.

        if (!pvBitsFile)
            pvBitsFile = (PVOID) ((PBYTE) (pbmihFile + 1)
                         + wNumColors * sizeof(RGBTRIPLE));
        break;

    default:
        MESSAGEBOX(GetFocus(), "Unknown DIB file format.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

// Fill in default values (for fields that can have defaults).

    if (pbmiSource->bmiHeader.biSizeImage == 0)
        pbmiSource->bmiHeader.biSizeImage = 
            BITS2BYTES( (DWORD) pbmiSource->bmiHeader.biWidth * 
                                pbmiSource->bmiHeader.biBitCount ) * 
                                pbmiSource->bmiHeader.biHeight;
    if (pbmiSource->bmiHeader.biClrUsed == 0)
        pbmiSource->bmiHeader.biClrUsed = wNumColors;

// Create memory DC.

    hdcMem = CreateCompatibleDC(NULL);
    if (!hdcMem) {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

// Create a 24BPP RGB DIB section and select it into the memory DC.

    pbmiRGB = (BITMAPINFO *)
              LocalAlloc(LMEM_FIXED|LMEM_ZEROINIT, sizeof(BITMAPINFO) );
    if (!pbmiRGB)
    {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

    pbmiRGB->bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
    pbmiRGB->bmiHeader.biWidth         = pbmiSource->bmiHeader.biWidth;
    pbmiRGB->bmiHeader.biHeight        = pbmiSource->bmiHeader.biHeight;
    pbmiRGB->bmiHeader.biPlanes        = 1;
    pbmiRGB->bmiHeader.biBitCount      = 24;
    pbmiRGB->bmiHeader.biCompression   = BI_RGB;
    pbmiRGB->bmiHeader.biSizeImage     = pbmiRGB->bmiHeader.biWidth
                                         * abs(pbmiRGB->bmiHeader.biHeight) * 3;

    hbmRGB = CreateDIBSection(hdcMem, pbmiRGB, DIB_RGB_COLORS, 
                              (PVOID *) &pjBitsRGB, NULL, 0);

    if (!hbmRGB)
    {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK); 
        goto tkDIBLoadImage_cleanup;
    }

    if (!SelectObject(hdcMem, hbmRGB))
    {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

// Slam the DIB file image into the memory DC.  GDI will do the work of
// translating whatever format the DIB file has into RGB format.

    if (!SetDIBits(hdcMem, hbmRGB, 0, pbmiSource->bmiHeader.biHeight, 
                   pvBitsFile, pbmiSource, DIB_RGB_COLORS))
    {
        MESSAGEBOX(GetFocus(), "Image file conversion error.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }
    GdiFlush();     // make sure that SetDIBits executes

// Convert to TK image format (packed RGB format).
// Allocate with malloc to be consistent with tkRGBImageLoad (i.e., app
// can deallocate with free()).

    pjTKBits = (PBYTE) malloc(pbmiRGB->bmiHeader.biSizeImage);
    if (!pjTKBits)
    {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }

    pjSrc = pjBitsRGB;
    pjDst = pjTKBits;
    // src lines end on LONG boundary - so need to skip over any padding bytes
    padBytes = pbmiSource->bmiHeader.biWidth % sizeof(LONG);
    for (i = 0; i < pbmiSource->bmiHeader.biHeight; i++)
    {
        for (j = 0; j < pbmiSource->bmiHeader.biWidth; j++)
        {
            // swap R and B
            *pjDst++ = pjSrc[2];
            *pjDst++ = pjSrc[1];
            *pjDst++ = pjSrc[0];
            pjSrc += 3;
        }
        pjSrc += padBytes;
    }

// Allocate and initialize the TK_RGBImageRec.
// Allocate with malloc to be consistent with tkRGBImageLoad (i.e., app
// can deallocate with free()).

    final = (TK_RGBImageRec *)malloc(sizeof(TK_RGBImageRec));
    if (final == NULL) {
        MESSAGEBOX(GetFocus(), "Out of memory.", "Error", MB_OK);
        goto tkDIBLoadImage_cleanup;
    }
    // If we get to here, we have suceeded!
    final->sizeX = pbmiSource->bmiHeader.biWidth;
    final->sizeY = pbmiSource->bmiHeader.biHeight;
    final->data = pjTKBits;

// Cleanup objects.

tkDIBLoadImage_cleanup:
    {
        if (hdcMem)
            DeleteDC(hdcMem);

        if (hbmRGB)
            DeleteObject(hbmRGB);

        if (pbmiRGB)
            LocalFree(pbmiRGB);

        if (pbmiSource)
            LocalFree(pbmiSource);

        if (pvFile)
            UnmapViewOfFile(pvFile);

        if (hMap)
            CloseHandle(hMap);

        if (hFile != INVALID_HANDLE_VALUE)
            CloseHandle(hFile);
    }

// Check for error.

    if (!final)
    {
        if (pjTKBits)
            free(pjTKBits);

        if ( (hFile == INVALID_HANDLE_VALUE) || !hMap || !pvFile )
        {
            CHAR ach[256];

            bUnicode ? wsprintf(ach, "Failed to open DIB file %ws.\n", fileName) :
                       wsprintf(ach, "Failed to open DIB file %s.\n", fileName);

            MESSAGEBOX(GetFocus(), ach, "Error", MB_OK);
        }
    }

    return final;
}
