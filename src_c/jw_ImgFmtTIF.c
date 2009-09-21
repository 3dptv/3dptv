/*
 * g_ImgFmtTIF.c  --
 *
 * S. Nebiker, 10.8.95
 *
 * taken from tkImgFmtPPM.c --
 *
 *	A photo image file handler for TIF files.
 *      (Microsoft Windows Bitmap Version 3.x)
 *
 * Copyright (c) 1994 The Australian National University.
 * Copyright (c) 1994-1995 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * Author: Paul Mackerras (paulus@cs.anu.edu.au),
 *	   Department of Computer Science,
 *	   Australian National University.
 */
/*
static char sccsid[] = "@(#) g_ImgFmtTIF.c 1.7 95/06/14 22:49:55";
*/
#include "ptv.h"


/*
 * The maximum amount of memory to allocate for data read from the
 * file.  If we need more than this, we do it in pieces.
 */

#define MAX_MEMORY	20000		/* don't allocate > 10KB */

/*
 * Define PGM and PPM, i.e. gray images and color images.
 */

#define TIF 1
#define PGM 2 /* dummies, to be removed */
#define PPM 3 /* dummies, to be removed */

/*
 * The format record for the TIF file format:
 */

static int FileMatchTIF(Tcl_Channel chan,
		  CONST char *fileName,
		  Tcl_Obj *formatString,
		  int *widthPtr,
		  int *heightPtr,
		  Tcl_Interp *interp);

static int FileReadTIF(Tcl_Interp *interp,
		Tcl_Channel chan,
		CONST char *fileName,
		Tcl_Obj *format,
		Tk_PhotoHandle imageHandle,
		int destX, int destY,
		int width, int height,
		int srcX, int srcY);

static int FileWriteTIF(Tcl_Interp *interp,
    CONST char *fileName,
    Tcl_Obj *formatString,
    Tk_PhotoImageBlock *blockPtr);


Tk_PhotoImageFormat tkImgFmtTIF = {
    "TIF",                      /* name */
    FileMatchTIF,               /* fileMatchProc */
    NULL,                       /* stringMatchProc */
    FileReadTIF,                /* fileReadProc */
    NULL,                       /* stringReadProc */
    FileWriteTIF,               /* fileWriteProc */
    NULL,                       /* stringWriteProc */
};

/*
 * Prototypes for local procedures defined in this file:
 *
 * tkTIFFErrorHandler cannot be declared as STATIC as it is used by
 * other (external) procedures, too.
 */

/*
void                    tkTIFFErrorHandler _ANSI_ARGS_((const char* module,
                            const char* fmt, va_list ap));
*/

/*
 *----------------------------------------------------------------------
 *
 * tkTIFFErrorHandler --
 *
 *	This procedure creates a new TIFF error handler replacing the
 *      default TIFFErrorHandler
 *
 * Results:
 *
 *
 * Side effects:
 *
 *
 * Created / modified:
 *      S. Nebiker, 16.8.95
 *
 *----------------------------------------------------------------------
 */

void
tkTIFFErrorHandler(const char* module, const char* fmt, va_list ap)
{
        if (module != NULL)
                fprintf(stderr, "(tkTIFFErrorHandler:) %s: ", module);
        vfprintf(stderr, fmt, ap);
        fprintf(stderr, ".\n");
}

/*
 *----------------------------------------------------------------------
 *
 * FileMatchTIF --
 *
 *	This procedure is invoked by the photo image type to see if
 *	a file contains image data in TIF format.
 *
 * Results:
 *	The return value is 1 if the first characters in file "f" look
 *	like TIF data, and 0 otherwise.
 *
 * Side effects:
 *
 *----------------------------------------------------------------------
 */

int FileMatchTIF(Tcl_Channel chan,
		  CONST char *fileName,
		  Tcl_Obj *formatString,
		  int *widthPtr,
		  int *heightPtr,
		  Tcl_Interp *interp)
{
    TIFF *tif;
    TIFFErrorHandler prev_handler; /* Current TIFF error handler  */
    static TIFFErrorHandler temp_handler = tkTIFFErrorHandler;
                                   /* Temporary TIFF error handler
                                      (is used to avoid output of
                                       libtiff error messages to stderr) */

    prev_handler = TIFFSetErrorHandler(temp_handler);

    /* tif = TIFFFdOpen(f, fileName, "r"); */ /* this did not work */
    tif = TIFFOpen(fileName, "r");
    if (tif != NULL) {
        TIFFClose(tif);
        /* printf("FileMatchTIF: TIFF file found !"); */
        TIFFSetErrorHandler(prev_handler);
        return (TIF);
    } else {
        /* printf("FileMatchTIF: No TIF file found !"); */
        TIFFSetErrorHandler(prev_handler);
        return (0);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * FileReadTIF --
 *
 *	This procedure is called by the photo image type to read
 *	TIF format data from a file and write it into a given
 *	photo image.
 *
 * Results:
 *	A standard TCL completion code.  If TCL_ERROR is returned
 *	then an error message is left in interp->result.
 *
 * Side effects:
 *	The access position in file f is changed, and new data is
 *	added to the image given by imageHandle.
 *
 *----------------------------------------------------------------------
 */


int FileReadTIF(Tcl_Interp *interp,
		Tcl_Channel chan,
		CONST char *fileName,
		Tcl_Obj *format,
		Tk_PhotoHandle imageHandle,
		int destX, int destY,
		int width, int height,
		int srcX, int srcY)
{
    int nBytes, h;
    int nLines;                   /* actual number of scan lines read
                                   * for each image block pointed at
                                   * by pixelPtr */
    int nDefBlockLines;           /* default number of scan lines to be read
                                   * for each image block pointed at
                                   * by pixelPtr */
    unsigned char *pixelPtr;
    Tk_PhotoImageBlock block;
    unsigned long fileImgHeight, fileImgWidth;
    unsigned short bitspersample;
    unsigned short photometrictype;

    short n;
    unsigned long row;
    unsigned long startrow;
    unsigned long pixelPtrIndex;
    TIFF *tif;
    TIFFErrorHandler prev_handler; /* Current TIFF error handler  */
    static TIFFErrorHandler temp_handler = tkTIFFErrorHandler;
                                   /* Temporary TIFF error handler
                                      (is used to avoid output of
                                       libtiff error messages to stderr) */

    prev_handler = TIFFSetErrorHandler(temp_handler);

    tif = TIFFOpen(fileName, "r");
    if (tif == NULL) {
        printf("FileReadTIF: Problem reading TIF file !");
	Tcl_AppendResult(interp, "couldn't read TIF header from file \"",
		fileName, "\"", NULL);
        TIFFSetErrorHandler(prev_handler);
        return TCL_ERROR;
    }

    /* test entry !!! */
    /*
    printf("\nFileReadTIF: file name     = %s", fileName);
    printf("\nFileReadTIF: height        = %d", height);
    printf("\nFileReadTIF: width         = %d", width);
    */
    /* obtain dimensions of (first) image in TIFF file */
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &fileImgHeight);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &fileImgWidth);
    /* test entry !!! */
    /*
     printf("\nFileReadTIF: image height  = %d", fileImgHeight);
     printf("\nFileReadTIF: image width   = %d\n", fileImgWidth);
    */


    if ((fileImgWidth <= 0) || (fileImgHeight <= 0)) {
	Tcl_AppendResult(interp, "TIF image file \"", fileName,
		"\" has dimension(s) <= 0", (char *) NULL);
	return TCL_ERROR;
    }

    if ((srcX + width) > fileImgWidth) {
	width = fileImgWidth - srcX;
    }
    if ((srcY + height) > fileImgHeight) {
	height = fileImgHeight - srcY;
    }
    if ((width <= 0) || (height <= 0)
	|| (srcX >= fileImgWidth) || (srcY >= fileImgHeight)) {
	return TCL_OK;
    }


    /*
       determine type of TIFF image: PhotometricInterpretation
       WhiteIsZero       => 0 (bw or grayscale image)
       BlackIsZero       => 1 (bw or grayscale image)
       RGB Model         => 2 (RGB full color model)
       Color Palette     => 3
       Transparency Mask => 4
    */
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photometrictype);
    /* test entry !!! */
    if ((photometrictype == 0) || (photometrictype == 1) || (photometrictype == 3)) {
        /*
           bw or grayscale image
        */
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
        /* test entry !!! */
	/*
        printf("\nFileReadTIF: bits / sample = %d", bitspersample);
	*/
        block.pixelSize = 1;
        block.offset[0] = 0;
        block.offset[1] = 0;
        block.offset[2] = 0;
    } else if (photometrictype == 2) {
        /*
           RGB full color image

           note: - bitspersample is always = 8
                 - there is no color map
        */
        block.pixelSize = 3;
        block.offset[0] = 0;
        block.offset[1] = 1;
        block.offset[2] = 2;
    } else {
        Tcl_AppendResult(interp, "photometric interpretation >2 ",
            "(palette color and transparency mask) ",
            "not yet supported", NULL);
        TIFFSetErrorHandler(prev_handler);
        return TCL_ERROR;
    }

    /* block.width = width; */
    height = fileImgHeight;
    width = fileImgWidth;
    block.width = fileImgWidth;
    block.pitch = block.pixelSize * fileImgWidth;
    /* test entry !!! */
    /*
    printf("\nFileReadTIF: height        = %d", height);
    printf("\nFileReadTIF: width         = %d", width);
    printf("\nFileReadTIF: fileImgWidth         = %d", fileImgWidth);
    */
    /* the following two values should be identical */
    /*
    printf("\nFileReadTIF: block.pitch     = %d", block.pitch);
    printf("\nFileReadTIF: block.pixelSize = %d", block.pixelSize);
    printf("\nFileReadTIF: scanline size   = %d", TIFFScanlineSize(tif));
    */
    Tk_PhotoExpand(imageHandle, destX + width, destY + height);

    /*
    if (srcY > 0) {
	fseek(f, (long) (srcY * block.pitch), SEEK_CUR);
    }
    */

    nDefBlockLines = (MAX_MEMORY + block.pitch - 1) / block.pitch;
    nLines = nDefBlockLines;
    nLines = height;
    if (nLines > height) {
        nLines = height;
    }
    if (nLines <= 0) {
        nLines = 1;
    }

    /* determine the required buffer size by calling TIFFScanlineSize */
    nBytes = nLines * block.pitch;
    /*
    printf("\nFileReadTIF: nLines        = %d", nLines);
    printf("\nFileReadTIF: nBytes        = %d", nBytes);
    */
    pixelPtr = (unsigned char *) ckalloc((unsigned)nBytes);
    block.pixelPtr = pixelPtr + srcX * block.pixelSize;

    n = 0;
    for (h = height; h > 0; h -= nLines) {
        n++;
        startrow = (n-1) * nLines;
        if (nLines > h) {
            nLines = h;
            nBytes = nLines * block.pitch;
        }
        for (row = startrow; row < (startrow + nLines); row++) {
            /* read the data for the spezified row into the buffer */
            pixelPtrIndex = row*block.pitch - (n-1)*nDefBlockLines*block.pitch;
            /* test entry !!! */
            /* printf("\nFileReadTIF: pixelPtr_Indx = %d", pixelPtrIndex); */
            TIFFReadScanline(tif, &(pixelPtr[pixelPtrIndex]), row, 1);
        }
	block.height = nLines;

	Tk_PhotoPutBlock(imageHandle, &block, destX, destY, width, nLines, TK_PHOTO_COMPOSITE_SET);
	destY += nLines;
    }


    ckfree((char *) pixelPtr);

    TIFFClose(tif);
    /* test entry !!! */
    /*
    printf("\nFileReadTIF: TIF file closed !");
    */
    TIFFSetErrorHandler(prev_handler);

    return TCL_OK;



}

/*
 *----------------------------------------------------------------------
 *
 * FileWriteTIF --
 *
 *	This procedure is invoked to write image data to a file in TIF
 *	format.
 *
 * Results:
 *	A standard TCL completion code.  If TCL_ERROR is returned
 *	then an error message is left in interp->result.
 *
 * Side effects:
 *	Data is written to the file given by "fileName".
 *static
 *----------------------------------------------------------------------
 */

int FileWriteTIF(interp, fileName, formatString, blockPtr)
    Tcl_Interp *interp;
    CONST char *fileName;
    Tcl_Obj *formatString;
    Tk_PhotoImageBlock *blockPtr;
{
  /* not implemented yet */
    FILE *f;
    int w, h;
    int greenOffset, blueOffset, nBytes;
    unsigned char *pixelPtr, *pixLinePtr;

    if ((f = fopen(fileName, "wb")) == NULL) {
	Tcl_AppendResult(interp, fileName, ": ", Tcl_PosixError(interp),
		(char *)NULL);
	return TCL_ERROR;
    }

    fprintf(f, "P6\n%d %d\n255\n", blockPtr->width, blockPtr->height);

    pixLinePtr = blockPtr->pixelPtr + blockPtr->offset[0];
    greenOffset = blockPtr->offset[1] - blockPtr->offset[0];
    blueOffset = blockPtr->offset[2] - blockPtr->offset[0];

    if ((greenOffset == 1) && (blueOffset == 2) && (blockPtr->pixelSize == 3)
	    && (blockPtr->pitch == (blockPtr->width * 3))) {
	nBytes = blockPtr->height * blockPtr->pitch;
	if (fwrite(pixLinePtr, 1, (unsigned) nBytes, f) != nBytes) {
	    goto writeerror;
	}
    } else {
	for (h = blockPtr->height; h > 0; h--) {
	    pixelPtr = pixLinePtr;
	    for (w = blockPtr->width; w > 0; w--) {
		if ((putc(pixelPtr[0], f) == EOF)
			|| (putc(pixelPtr[greenOffset], f) == EOF)
			|| (putc(pixelPtr[blueOffset], f) == EOF)) {
		    goto writeerror;
		}
		pixelPtr += blockPtr->pixelSize;
	    }
	    pixLinePtr += blockPtr->pitch;
	}
    }

    if (fclose(f) == 0) {
	return TCL_OK;
    }
    f = NULL;

 writeerror:
    Tcl_AppendResult(interp, "error writing \"", fileName, "\": ",
	    Tcl_PosixError(interp), (char *) NULL);
    if (f != NULL) {
	fclose(f);
    }
    return TCL_ERROR;
}
