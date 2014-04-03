/*
** Utility routines for handling MapII formatted files.  MapII is
** ugly file format not fit for production use BUT it's compatible
** with SME and seems to be one of the file formats that SME 
** consistently gets right.  So...in the name of backwards compatibility
** we keep this around.
**
** The default byte order for MapII seems to be little-endian although
** the format was originally developed on a Mac suggesting it should
** be big-endian.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* GLUC is defined when this file is used as part of the LEAM GLUC model.
** Otherwise, it it can be used as library of routines for dealing
** with Map2 files.
*/
#ifdef GLUC
#  include <mpi.h>
#  include "leam.h"
#else
extern char estring[];
extern void errorExit(char *s)
#endif

#include "map2.h"


/* Byteswap values when necessary.  Handles 16, 32, and 64-bit quantities.
**
** Note: code lifted directly from Python Numeric module.
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


/* Read the map2 header file and return the most important
** data including location, rows, columns, and data size.
**
** Note: location must be a pre-defined buffer, memory is not
** allocated for the location information.
*/
int MAP2readHeader(char *fname, char *location, 
                   int *rows, int *cols, int *size)
{
    char line[1024], *ptr;
    FILE *f;

    /* Open the header file */
    if ((f = fopen((const char *)fname, "r")) == NULL)  {
        sprintf(estring, "Unable to open Map2 header file %s", fname);
        errorExit(estring);
        return 0;
    }

    /* search for tags */
    *rows = *cols = *size = -1;
    if (location != NULL) location[0] = '\0';
    while (fgets(line, sizeof line, f) != NULL)  {
        if (!strncasecmp(line, "ROWS", 4))
            *rows = strtol(strchr(line, '=')+1, NULL, 10);
        else if (!strncasecmp(line, "COLUMNS", 7))
            *cols = strtol(strchr(line, '=')+1, NULL, 10);
        else if (!strncasecmp(line, "SIZE", 4))
            *size = strtol(strchr(line, '=')+1, NULL, 10);
        else if (!strncasecmp(line, "LOCATION", 8))  {

            if (location == NULL) continue;

            if ((ptr = strchr(line, 0xc8)) == NULL)  {
                sprintf(estring, "LOCATION without 0xc8 char in %s", fname);
                errorExit(estring);
            }
            else *ptr = '\0';
            if ((ptr = strchr(line, 0xc7)) == NULL) {
                sprintf(estring, "LOCATION without 0xc7 char in %s", fname);
                errorExit(estring);
            }
            else
                ptr += 1;
            strcpy(location, ptr);
        }
    }

    fclose(f);

    /* make sure all the dimensional tags were found */
    if (*rows == -1 || *cols == -1 || *size == -1)  {
        sprintf(estring, 
            "Failed to located required tags in Map2 header file %s", fname);
        errorExit(estring);
    }

    /* make sure the location information was found */
    if (location != NULL && location[0] == '\0')  {
        sprintf(estring, 
            "Failed to located required tags in Map2 header file %s", fname);
        errorExit(estring);
    }

    return 1;
}

/* Write the map2 header using the following format:
**
** FILETYPE=INTERCHANGE
** ROWS=y
** COLUMNS=x
** CELLSIZE=30
** UNITS=M
** FORMAT=BIN
** SIZE=bytes
** LOCATION=\0xc7stlme_boundary.bin\0xc8
*/
int MAP2writeHeader(char *fname, int rows, int cols, int size)
{
    char binfile[1024], *ptr;
    FILE *f;


    if ((f = fopen(fname, "w")) == NULL)  {
        sprintf(estring, "Unable write Map2 header %s", fname);
        errorExit(estring);
    }

    /* extract the base of the header file name (no path or extension) */
    ptr = strrchr(fname, '/');
    strcpy(binfile, (ptr == NULL) ? fname : ptr + 1);
    ptr = strrchr(binfile, '.');
    if (ptr != NULL) *ptr = '\0';
    /* if no extension was found we'll simply append .bin */

    /* write the header */
    fprintf(f, "FILETYPE=INTERCHANGE\n");
    fprintf(f, "ROWS=%d\nCOLUMNS=%d\n", rows, cols);
    fprintf(f, "CELLSIZE=30\nUNITS=m\n");
    fprintf(f, "SIZE=%d\n", size);
    fprintf(f, "LOCATION=\307%s.bin\310\n", binfile);

    fclose(f);
    return 1;
}

/* Open the Map2 file header, extract location of binary file
** and return the file point of the binary file opened in the 
** appropriate mode.
*/
FILE *MAP2openBinary(char *fname, char *mode)
{
    int  r, c, s;
    char location[128];
    char binname[1024], *ptr;
    FILE *f;

    MAP2readHeader(fname, location, &r, &c, &s);

    /* prepend the path information */
    strcpy(binname, fname);
    if ((ptr = strrchr(binname, '/')) != NULL)
        strcpy(ptr+1, location);
    else
        strcpy(binname, location);
   
    /* open the binary file */
    if ((f = fopen(binname, mode)) == NULL)  {
        sprintf(estring, "Unable to open %s", fname);
        errorExit(estring);
        return NULL;                    /* We'll never get here */
    }

    return f;
}


/* Wrapper for fclose.  Not needed but provides a matching function
** MAP2openBinary and we might as well be consistant.
*/
void MAP2close(FILE *f)
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
int MAP2readBuffer(char *fname, int offset, char *ptr, int count, int size)
{
    int i;
    FILE *f;

    if ((f = MAP2openBinary(fname, "rb")) == NULL)  {
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

    MAP2close(f);

    if (size > 1) byteswap(ptr, count, size);

    return 1;
}

int MAP2writeBuffer(char *fname, int offset, void *data, int count, int size)
{
    FILE *f;

    f = MAP2openBinary(fname, "wb");
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
    MAP2close(f);

    return 1;
}
