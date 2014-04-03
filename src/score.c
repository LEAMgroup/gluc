/*
** This module was used provided a mechanism to  score the model resultsi
** as part of an experiment using genetic algorithms to assign model weights.
** It should be removed at this point.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "leam.h"

/*
** Simple routine for reading in reference counts
** from another file.
*/
void readReferenceCounts(int *totals, int len, char *fname)
{
    int i, zone, count;
    char line[1024];
    FILE *fptr;

    if (fname == NULL)
        return;

    if (debug && myrank == 0)
        fprintf(stderr, "readReferenceCounts reading %s\n", fname);

    if ((fptr = fopen(fname, "r")) == NULL)  {
        sprintf(estring, "unable to open REFERENCE_COUNTS = %s\n", fname);
        errorExit(estring);
    }

    while (fgets(line, sizeof line, fptr) != NULL)  {
        sscanf(line, "%d\t%d", &zone, &count);
        if (zone < 0 || zone >= len)  {
            fprintf(stderr, "Reference count out of range: %d, %d\n", 
                    zone, count);
            continue;
        }
        if (zone >= 0 && zone < len)
            totals[zone] = count;
        else if (debug && myrank == 0)
            fprintf(stderr, "Warning: zone out of range, %s\n", line);
    }
    fclose(fptr);
    return;
}

readReferencePop(float *totals, int len, char *fname)
{
    int i;
    FILE *fptr;

    if ((fptr = fopen(fname, "r")) == NULL)  {
        sprintf(estring, "unable to open REFERENCE_COUNTS = %s\n", fname);
        errorExit(estring);
    }

    for (i=0; i<len; i+=1);

    fclose(fptr);
    return;
}

void printResultCounts(int *totals, int len)
{
    int i;
    char *cptr;
    FILE *fptr;

    cptr = SMEgetFileName("REFERENCE_RESULTS");
    if (myrank != 0 || cptr == NULL)
        return;

    if ((fptr = fopen(cptr, "w")) == NULL) {
        sprintf(estring, "failed opening REFERENCE_RESULTS = %s\n", cptr );
        errorExit(estring);
    }

    fprintf(fptr, "ZONE\tCOUNT\n");

    for (i=0; i<len; i+=1)
        if (totals[i] != 0)
            fprintf(fptr, "%d\t%d\n", i, totals[i]);

    if (fptr != stdout) fclose(fptr);
    return;
}

void printResultPop(float *totals, int len)
{
    int i;
    char *cptr;
    FILE *fptr;

    if ((cptr = SMEgetFileName("REFERENCE_RESULTS")) == NULL)
        fptr = stdout;

    else if ((fptr = fopen(cptr, "w+")) == NULL)  {
        sprintf(estring, "failed opening REFERENCE_RESULTS = %s\n", cptr);
        errorExit(estring);
    }

    fprintf(fptr, "ZONE\tCOUNT\n");

    for (i=0; i<len; i+=1)
            fprintf(fptr, "%d\t%f\n", i, totals[i]);

    if (fptr != stdout) fclose(fptr);
    return;
}

double scoreSumErrSquared(int *ref, int *totals, int *active, int len)
{
    int i;
    double sum = 0.0;

    for (i=0; i<len; i+=1)  {
        if (active[i] > 0)
            sum += (ref[i] - totals[i]) * (ref[i] - totals[i]);
    }

    /* return to GA engine */
    if (debug && myrank == 0)
        fprintf(stderr, "SumErrSquared: Score = %f\n", sum);

    if (myrank == 0 && SMEgetFileName("GA_ENGINE") != NULL)
        GAsendFit(sum);

    return sum;
}

/*
** scoreResults compares to the final land use change map
** to a reference map and reference counts.  The basic procedure
** counts the number of new residential cells within each
** region of a reference map.  This array is compared against
** the reference counts using a sum of square of the error.
**
** Optionally the population density map is used to compute
** the exact population change within each region.  Score value
** is computed in the same manner.
*/
double scoreResults(int *refcounts, int reflen, int *refmap, 
                  unsigned char *lu, int count)
{
    int *totals, *active;
    double score;

    /* if no reference map given return without scoring */
    if (reflen <= 0)
        return;

    if (debug && myrank == 0)
        fprintf(stderr, "Scoring results, reflen = %d\n", reflen);

    active = (int *)getMem(reflen * sizeof (int), "active zones");
    spatialHistogram(active, reflen, refmap, count);

    totals = (int *)getMem(reflen * sizeof (int), "score totals array");
    spatialCorrelatedCount(totals, reflen, refmap, lu, count, LU_LRES);

    printResultCounts(totals, reflen);

    /** Currently using Sum of the error squared
    ** but other methods are probably better.
    */
    score = scoreSumErrSquared(refcounts, totals, active, reflen);

    free(totals);
    free(active);
    return score;
}
