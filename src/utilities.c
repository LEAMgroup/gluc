/*
** This module provide wrapper functions for malloc (calloc really).
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "leam.h"

char estring[1024];


char *getMem(int bytes, char *description)
{
   char *ptr;

   if ((ptr=calloc(bytes, 1)) == NULL)  {
       sprintf(estring, "Insufficient Memory for %s, %d bytes requested",
               description, bytes);
       errorExit(estring);
       return ptr;          /* can't get here, but stops warnings */
   }
   else 
       return ptr;
}

/* simple wrapper, just to be consistant with getMem
*/
void freeMem(void *ptr)
{

    free(ptr);
}


/*
** Handle cases where we need to bail but we want
** to make sure that MPI_Finalize is called.
*/
void errorExit(char *str)
{
    fprintf(stderr, "P%d: Error: %s\n", myrank, str);
    MPI_Finalize();
    exit(-1);
}
