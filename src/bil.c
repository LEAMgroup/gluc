/* 
** Support the BIL format.  Like the MAP2 format BIL is a multifile
** format for multiband raster images extended to support geolocation.
** GLUC supports only single band data and doesn't provide any form
** of interpolation so all grids are expexted to be of identical
** extent and resolutions.

** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
**/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "leam.h"
#include "bil.h"


/* Byteswap values when necessary.  Handles 16, 32, and 64-bit quantities.
**
** Credits: code lifted directly from Python Numeric module.
*/
static void byteswap(void *p, int n, int size) {
    char *a, *b, c;
    int x = 1;

    /* check for little-endianess */
    if (*(char*)&x) return;

    switch(size) {
    case 2:
        for (a = (char*)p ; n > 0; n--, a += 1) {
            b = a + 1;
            c = *a; *a++ = *b; *b   = c;
        }
        break;
    case 4:
        for (a = (char*)p ; n > 0; n--, a += 2) {
            b = a + 3;
            c = *a; *a++ = *b; *b-- = c;
            c = *a; *a++ = *b; *b   = c;
        }
        break;
    case 8:
        for (a = (char*)p ; n > 0; n--, a += 4) {
            b = a + 7;
            c = *a; *a++ = *b; *b-- = c;
            c = *a; *a++ = *b; *b-- = c;
            c = *a; *a++ = *b; *b-- = c;
            c = *a; *a++ = *b; *b   = c;
        }
        break;
    default:
        break;
    }
}

char *get_bil_type(char *mpi_name)
{
    if (!strcmp(mpi_name, "MPI_UNSIGNED_CHAR"))
        return "UNSIGNEDINT";
    else if (!strcmp(mpi_name, "MPI_INTEGER"))
        return "SIGNEDINT";
    else if (!strcmp(mpi_name, "MPI_FLOAT"))
        return "FLOAT";
    else
        return "SIGNEDINT";
}

/*
** converts BIL filename to the associated header file name
*/
static char *get_hdr_name(char *bilname)
{
    char *hdrname, *ptr;
    int x;

    x = strlen(bilname);
    hdrname = (char *)getMem(x+8, "BIL Header string dup");
    //hdrname = (char *)getMem(strlen(bilname)+8, "BIL Header string dup");
    strcpy(hdrname, bilname);

    if ((ptr = strrchr(hdrname, '.')) != NULL)
        *ptr = '\0';

    strcat(hdrname, ".hdr");

    return hdrname;
}


/* Read the map2 header file and return the most important
** data including rows, columns, and data size.
**
*/
int BILreadHeader(char *bilname, int *rows, int *cols, int *size)
{
    char line[1024], *hdrname, *tag;
    FILE *f;

    hdrname = get_hdr_name(bilname);

    /* Open the header file */
    if ((f = fopen((const char *)hdrname, "r")) == NULL)  {
        sprintf(estring, "Unable to open BIL header file %s", hdrname);
        errorExit(estring);
        return 0;
    }

    /* search for tags */
    *rows = *cols = *size = -1;
    while (fgets(line, sizeof line, f) != NULL)  {
        if (!strncasecmp(line, "NROWS", 5))
            *rows = strtol(line+5, NULL, 10);
        else if (!strncasecmp(line, "NCOLS", 5))
            *cols = strtol(line+5, NULL, 10);
        else if (!strncasecmp(line, "NBITS", 5))
            *size = strtol(line+5, NULL, 10);
    }

    fclose(f);

    /* make sure all the dimensional tags were found */
    if (*rows == -1 || *cols == -1 || *size == -1)  {
        sprintf(estring, 
            "Failed to located required tags in BIL header file %s", hdrname);
        errorExit(estring);
    }

    freeMem(hdrname);
    return 1;
}


/* Read the bil header file and return the extent information
** including ULXMAP, ULYMAP, XDIM, and YDIM
**
*/
int BILreadExtents(char *bilname, float *ulx, float *uly, 
                   float *xdim, float *ydim)
{
    char line[1024], *hdrname;
    FILE *f;

    hdrname = get_hdr_name(bilname);

    /* Open the header file */
    if ((f = fopen((const char *)hdrname, "r")) == NULL)  {
        sprintf(estring, "Unable to open BIL header file %s", hdrname);
        errorExit(estring);
        return 0;
    }

    /* search for tags */
    while (fgets(line, sizeof line, f) != NULL)  {
        if (!strncasecmp(line, "ULXMAP", 6))
            *ulx = strtod(line+6, NULL);
        else if (!strncasecmp(line, "ULYMAP", 6))
            *uly = strtod(line+6, NULL);
        else if (!strncasecmp(line, "XDIM", 4))
            *xdim = strtod(line+4, NULL);
        else if (!strncasecmp(line, "YDIM", 4))
            *ydim = strtod(line+4, NULL);
    }

    fclose(f);

    freeMem(hdrname);
    return 1;
}
/*
** BILwriteHeader: writes the .hdr file for the bil format.
**
** Inputs:
**   fname -- filename string
**   rows, cols -- number of rows and cols
**   type -- MPI defined data type
**
** The header format
**
** BYTEORDER      I
** LAYOUT         BIL
** NROWS          3771
** NCOLS          5149
** NBANDS         1
** NBITS          32
** BANDROWBYTES   20596
** TOTALROWBYTES  20596
** PIXELTYPE      FLOAT
** ULXMAP         641829
** ULYMAP         4322074
** XDIM           30
** YDIM           30
*/
int BILwriteHeader(char *bilname, int rows, int cols, int size, char *type,
                   float ulx, float uly, float xdim, float ydim)

{
    char *hdrname, *typename;
    FILE *f;

    hdrname = get_hdr_name(bilname);
    typename = get_bil_type(type);

    if ((f = fopen(hdrname, "w")) == NULL)  {
        sprintf(estring, "Unable write BIL header %s", hdrname);
        errorExit(estring);
    }

    /* write the header */
    fprintf(f, "BYTEORDER\tI\n");
    fprintf(f, "LAYOUT\t\tBIL\n");
    fprintf(f, "NROWS\t\t%d\nNCOLS\t\t%d\n", rows, cols);
    fprintf(f, "NBANDS\t\t1\n");
    fprintf(f, "NBITS\t\t%d\n", size * 8);
    fprintf(f, "BANDROWBYTES\t%d\n", size * cols);
    fprintf(f, "TOTALROWBYTES\t%d\n", size * cols);
    fprintf(f, "PIXELTYPE\t%s\n", typename);
    fprintf(f, "ULXMAP\t%f\nULYMAP\t%f\n", ulx, uly);
    fprintf(f, "XDIM\t%f\nYDIM\t%f\n", xdim, ydim);

    fclose(f);
    freeMem(hdrname);
    return 1;
}

/* 
** Opens the BIL data file.
**
** Note: odd naming and pretty much enpty function to remain
** compatible with the old Map2 format.
*/
FILE *BILopenBinary(char *fname, char *mode)
{
    FILE *f;

    if ((f = fopen(fname, mode)) == NULL)  {
        sprintf(estring, "Unable to open %s", fname);
        errorExit(estring);
    }

    return f;
}


/* Wrapper for fclose.  Not needed but provides a matching function
** MAP2openBinary and we might as well be consistant.
*/
void BILclose(FILE *f)
{
    fclose(f);
}


/* Open the Map2 binary file and seek to the correct location
** within the file.  Read the necessary number of elements and
** byteswap if necessary.
**
** ??? Should offset by elements or bytes?  (opted for bytes)
** ??? Should size be from header or from program?
*/
int BILreadBuffer(char *fname, int offset, char *ptr, int count, int size)
{
    int i;
    FILE *f;

    if ((f = BILopenBinary(fname, "rb")) == NULL)  {
        sprintf(estring, "file %s not found.", fname);
        errorExit(estring);
        return 0;
    }

    if (fseek(f, offset, SEEK_SET))  {
            sprintf(estring, "seek to %d in %s failed.", offset, fname);
            errorExit(estring);
            return 0;
        }

    if ((i = fread(ptr, size, count, f))  != count)  {
        sprintf(estring, "unable to read %d bytes from %s, got %d.\n", 
                count, fname, i);
        errorExit(estring);
        return 0;
    }

    BILclose(f);

    if (size > 1) byteswap(ptr, count, size);

    return 1;
}


int BILwriteBuffer(char *fname, int offset, void *data, int count, int size)
{
    FILE *f;

    f = BILopenBinary(fname, "wb");
    if (fseek(f, offset, SEEK_SET))  {
        sprintf(estring, "seek to %d in %s failed.", offset, fname);
        errorExit(estring);
        return 0;
    }

    /* swap data IN PLACE */
    if (size > 1) byteswap(data, count, size);

    if (fwrite(data, size, count, f) != count)  {
        sprintf(estring, "could not write %d elements to %s.", count, fname);
        errorExit(estring);
        return 0;
    }

    /* swap data IN PLACE (putting it back to original endianess */
    if (size > 1) byteswap(data, count, size);
    BILclose(f);

    return 1;
}
