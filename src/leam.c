/*
** LEAM is a stocastic cellular automata for predicting land use change.
**
** This files provides the command line parsing and general setup.  The
** internals of the model are located in luc.c most other source files 
** provide utility functions.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.  See
** the LICENCE file for details.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>

#include "leam.h"
#include "bil.h"

static char *ID = "$Id: leam.c 376 2011-01-10 16:43:30Z jefft $";
static char *SVN = "$HeadURL: http://plone.leamgroup.com/svn/desktop/gluc/branches/nocompete/leam.c $";

int debug = 0;
int timing = 0;
int distribMethod = EQROWS;
int myrank, nproc;
char *runName = NULL;

float ulx, uly;
float xdim, ydim;
char *projection;

/* Distributed Grid by Rows - read the boundary map and distributed
** the rows an evenly as possible across all the processors.  
*/
static void distributeGridRows(char *bfile, int *gridrows, 
                int *rows, int *cols)
{
    int i, size;

    BILreadHeader(bfile, rows, cols, &size);
    BILreadExtents(bfile, &ulx, &uly, &xdim, &ydim);
    projection = strdup(bfile);

    /* handle single processor case */
    if (nproc == 1)  {
        gridrows[0] = 0;
        gridrows[1] = *rows;
        return;
    }

    gridrows[0] = 0;
    for (i=1; i<nproc; i+=1)  {
        if (i <= *rows % nproc)
            gridrows[i] = gridrows[i-1] + (*rows / nproc) + 1;
        else
            gridrows[i] = gridrows[i-1] + (*rows / nproc);
    }
    gridrows[nproc]=*rows;
}

/* Distribute Grid by Cells - read the boundary map and count the
** number of active cells on each row and split the rows so each
** processor gets roughly equal number of active cells.
*/
static void distributeGridCells(char *bfile, int *gridrows, 
                 int *rows, int *cols)
{
    int i, j, size;
    int partial = 0, total = 0;
    int currow, proc;
    unsigned char *data;
    FILE *f;

    BILreadHeader(bfile, rows, cols, &size);
    BILreadExtents(bfile, &ulx, &uly, &xdim, &ydim);
    projection = strdup(bfile);

    /* if single processor then return all rows */
    if (nproc == 1) {
       gridrows[0] = 0;
       gridrows[1] = *rows;
       return;
    }

    /* count the total number of active cells */
    data = (unsigned char *)getMem(*cols, "distributeGridCells");
    f = BILopenBinary(bfile, "rb");
    for (j=0; j<*rows; j+=1)  {
        if (fread(data, 1, *cols, f) != *cols)  {
            sprintf(estring, "failed reading from %s", bfile);
            errorExit(estring);
        }
        for (i=0; i<*cols; i+=1)
            total += (data[i] != 0) ? 1 : 0;
    }

    /* reset to the beginning of the file */
    fseek(f, 0, SEEK_SET); 

    /* read each row until the partial count of active cells > than
    ** 1/nproc.  Continue this stategy until each processor has roughly
    ** equal number of active cells.
    */
    currow = -1; proc = 1;
    while (proc < nproc)  {
        do {
            currow += 1;
            if (fread(data, 1, *cols, f) != *cols)  {
                sprintf(estring, "failed reading from %s", bfile);
                errorExit(estring);
            }
            for (i=0; i<*cols; i+=1)
                partial += (data[i] != 0) ? 1 : 0;

        } while (partial <= total/(nproc-proc+1));

        gridrows[proc++] = currow;
        total -= partial; partial = 0;
    }
    gridrows[nproc] = *rows;

    BILclose(f);
    freeMem(data);
}

void GLUCusage(char *prog)
{
   printf("Usage: %s [SME options] [options]\n\n", prog);

   printf("SME options:\n");
   printf(" -ppath <path> : sets the project path\n");
   printf(" -p <project>  : sets the project name (directory ");
   printf("under project path)\n");
   printf(" -ci <config>  : specify the configuration file\n");
   printf("\n");

   printf("Options:\n");
   printf(" --newk         : turns on improved K factor (default)\n");
   printf(" --oldk         : turns on the old version of K factor\n");
   printf(" --resoff       : turns off residential development\n");
   printf(" --comoff       : turns off commercial development\n");
   printf(" --osoff        : turns off openspace development\n");
   printf(" --eqrows       : distribute rows equally across processors\n");
   printf(" --eqcells      : distribute cells equally across processors\n");
   printf(" --version      : print version information and exit\n");
   printf(" --histogram    : print histograms as part of debugging info\n");
   printf(" -r || --random : randomly seeds random number generator\n");
   printf(" -f || --final  : generates only final landuse and change map\n");
   printf(" -d || --debug  : turns on debugging information\n");
   printf(" -t || --timing : turns on performance timing results\n");
   printf(" -h || --help   : display this information\n");

   printf("\nversion info = %s from %s\n", ID, SVN);
}


void parseOptions(int argc, char *argv[])
{
    int i;
    char *name, *value;
    struct timeval randseed;

    for (i=1; i<argc; i+=1)  {
        if (!strcmp(argv[i], "--newk"))
            setBetterK(1);
        else if (!strcmp(argv[i], "--oldk"))
            setBetterK(0);
        else if (!strcmp(argv[i], "--resoff"))
            setModel(RES_MODEL, 0);
        else if (!strcmp(argv[i], "--comoff"))
            setModel(COM_MODEL, 0);
        else if (!strcmp(argv[i], "--osoff"))
            setModel(OS_MODEL, 0);
        else if (!strcmp(argv[i], "--eqrows"))
            distribMethod = EQROWS;
        else if (!strcmp(argv[i], "--eqcells"))
            distribMethod = EQCELLS;
        else if (!strcmp(argv[i], "--histogram"))
            debug = 2;
        else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--random")) {
            gettimeofday(&randseed, NULL);
            srand48(randseed.tv_usec);
        }
        else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--graph"))
            GRAPHreadFile(argv[++i]);
        else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug"))
            debug = 1;
        else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--timing"))
            timing = 1;
        else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--quiet"))
            debug = 0;
        else if (!strncmp(argv[i], "--", 2))  {
            name = strtok(argv[i]+2, "=");
            value = strtok(NULL, "=");
            SMEaddVar(name, value);
        }
        
    }
}

/* quickArgs is checked before anything else.  If the program invocation
** only requires a usage statement then print that message and exit.
*/
void quickArgs(int argc, char *argv[])
{
    int i;

    if (argc < 2)  {
        fprintf(stderr, "USAGE: %s [options] %s\n", argv[0], 
                "-ppath <path> -p <project> -ci <config>");
        GLUCusage(argv[0]);
        exit(1);
    }

    for (i=0; i<argc; i+=1)  {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))  {
            GLUCusage(argv[0]);
            exit(0);
        }
        else if (!strcmp(argv[i], "--version"))  {
            printf("Repository = %s\n", SVN);
            printf("Revision = %s\n", ID);
            exit(0);
        }
    }
}

int main(int argc, char *argv[])
{
    int i;
    int rows, cols, *gridrows;
    struct timeval start, ioend, end;

    /* check for invokations only requiring usage and version info */
    quickArgs(argc, argv);

    /* begin timing */
    gettimeofday(&start, NULL);

    // initial MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // set the MPI type names to match BIL format
    MPI_Type_set_name(MPI_CHAR, "UNSIGNEDINT");
    MPI_Type_set_name(MPI_BYTE, "UNSIGNEDINT");
    MPI_Type_set_name(MPI_SHORT, "SIGNEDINT");
    MPI_Type_set_name(MPI_UNSIGNED_SHORT, "UNSIGNEDINT");
    MPI_Type_set_name(MPI_INT, "SIGNEDINT");
    MPI_Type_set_name(MPI_UNSIGNED, "UNSIGNEDINT");
    MPI_Type_set_name(MPI_LONG, "SIGNEDINT");
    MPI_Type_set_name(MPI_UNSIGNED_LONG, "UNSIGNEDINT");
    MPI_Type_set_name(MPI_FLOAT, "FLOAT");
    MPI_Type_set_name(MPI_DOUBLE, "FLOAT");

    // load default graphs
    GRAPHinit();

    SMEparseOptions(argc, argv);
    parseOptions(argc, argv);

    if (debug && myrank == 0) { 
        fprintf(stderr, "LEAM model running. Args = '");
        for (i=1; i<argc; i+=1)
           fprintf(stderr, "%s ", argv[i]);
        fprintf(stderr, "'\n");
    }

    /* Read the boundary file and spread the computation
    ** over the processors evenly.  Broadcast array defining
    ** which processor is responsible for which rows.
    */
    gridrows = (int *)getMem(sizeof (int) * (nproc + 1), "gridrows");
    if (myrank == 0)  {
        if (distribMethod == EQROWS)
            distributeGridRows(SMEgetFileName("BOUNDARY_MAP"),
                gridrows, &rows, &cols);
        else if (distribMethod == EQCELLS)
            distributeGridCells(SMEgetFileName("BOUNDARY_MAP"),
                gridrows, &rows, &cols);
    }

    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(gridrows, nproc+1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Process the SME config file and run the model.
    */
    LUCconfigGrids(rows, cols, gridrows);
    LUCinitGrids();
    gettimeofday(&ioend, NULL);
    
    /*
    ** if GA Engine is specified loop...indefinitely?
    */
    do {
        LUCresetGrids();
        LUCrun();
    }  while (SMEgetFileName("GA_ENGINE") != NULL);

        
    /* Clean up.  We'll make sure that everyone gets here before
    ** calling MPI_Finalize.  We have to report timing prior to 
    ** calling MPI_Finalize because the number of active processes
    ** following MPI_Finalize is undefined!
    */
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&end, NULL);
    if (timing && myrank == 0)  {
        printf("\n===== Timing Information ======\n");
        if (end.tv_usec - start.tv_usec < 0)  {
            printf("%d processor run complete, wallclock time = %ld.%ld\n", 
            nproc, (long)(end.tv_sec - start.tv_sec - 1), 
            end.tv_usec - start.tv_usec + 1000000);
        }
        else {
            printf("%d processor run complete, wallclock time = %ld.%ld\n",
            nproc, (long)(end.tv_sec-start.tv_sec), end.tv_usec-start.tv_usec);
        }
        if (end.tv_usec - ioend.tv_usec < 0)  {
            printf("computional component, wallclock time = %ld.%ld\n", 
            (long)(end.tv_sec - ioend.tv_sec - 1), 
            end.tv_usec - ioend.tv_usec + 1000000);
        }
        else {
            printf("computional components, wallclock time = %ld.%ld\n",
            (long)(end.tv_sec-ioend.tv_sec), end.tv_usec-ioend.tv_usec);
        }
    }
    MPI_Finalize();
}
