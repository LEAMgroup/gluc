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

#ifndef MAP2_H
#define MAP2_H

extern int MAP2readHeader(char *, char *, int *, int *, int *);
extern int MAP2writeHeader(char *, int, int, int);
extern FILE *MAP2openBinary(char *, char *);
extern void MAP2close(FILE *);
extern int MAP2readBuffer(char *, int, char *, int, int);
extern int MAP2writeBuffer(char *, int, void *, int, int);

#endif
