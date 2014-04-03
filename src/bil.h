/* bil.c header file
**
** Support the BIL format.  Like the MAP2 format BIL is a multifile
** format for multiband raster images extended to support geolocation.
** GLUC supports only single band data and doesn't provide any form
** of interpolation so all grids are expexted to be of identical
** extent and resolutions.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/

#ifndef BIL_H
#define BIL_H

extern int BILreadHeader(char *, int *, int *, int *);
extern int BILwriteHeader(char *, int, int, int, char *,
                          float, float, float, float);
extern FILE *BILopenBinary(char *, char *);
extern void BILclose(FILE *);
extern int BILreadBuffer(char *, int, char *, int, int);
extern int BILwriteBuffer(char *, int, void *, int, int);

#endif
