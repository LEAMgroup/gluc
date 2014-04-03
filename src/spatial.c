/*
** The spatial module provides basic utility function for performing
** raster operations in an MPI safe manner.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include "leam.h"

#define COMP_NW(p,val) ((*(p - cols - 1) == val) ? 1: 0)
#define COMP_N(p,val)  ((*(p - cols) == val) ? 1 : 0)
#define COMP_NE(p,val) ((*(p - cols + 1) == val) ? 1: 0)
#define COMP_E(p,val)  ((*(p + 1) == val) ? 1: 0)
#define COMP_W(p,val)  ((*(p - 1) == val) ? 1: 0)
#define COMP_SW(p,val) ((*(p + cols - 1) == val) ? 1: 0)
#define COMP_S(p,val)  ((*(p + cols) == val) ? 1: 0)
#define COMP_SE(p,val) ((*(p + cols + 1) == val) ? 1: 0)

#define GET_NW(p)  (*(p - cols - 1))
#define GET_N(p)   (*(p - cols))
#define GET_NE(p)  (*(p - cols + 1))
#define GET_W(p)   (*(p - 1))
#define GET_E(p)   (*(p + 1))
#define GET_SW(p)  (*(p + cols - 1))
#define GET_S(p)   (*(p + cols))
#define GET_SE(p)  (*(p + cols + 1))


/* Compute the number of nearest neighbors of a particular type,
** stores results in a grid.
*/
void nearestNeighbors(unsigned char *dst, unsigned char *tmp, int val,
                      unsigned char *src, int rows, int cols)
{
    int i, j, offset;
    int res=0, com=0;

    for (i=0; i<rows*cols; i+=1)  {
      switch (src[i])  {
      case LU_LRES: case LU_HRES:
        tmp[i] = (val & RES_FLAG) ? 1 : 0;
        break;
      case LU_COM:
        tmp[i] = (val & COM_FLAG) ? 1 : 0;
        break;
      case LU_ROAD:
        tmp[i] = (val & ROAD_FLAG) ? 1 : 0;
        break;
      case LU_OS:
        tmp[i] = (val & OS_FLAG) ? 1 : 0;
        break;
      case LU_WATER:
        tmp[i] = (val & WATER_FLAG) ? 1 : 0;
        break;
      default:
        tmp[i] = 0;
      }
    }

    for (j=0; j<rows; j+=1)  {

        /* row's western most cell */
        offset = j * cols;
        dst[offset] =  GET_N(tmp+offset) + GET_NE(tmp+offset)
                     +                     GET_E(tmp+offset)
                     + GET_S(tmp+offset) + GET_SE(tmp+offset)
                     ;

        /* row's middle cells */
        for (i=1; i<cols-1; i+=1)
            dst[offset+i] = 
                              GET_NW(tmp+offset+i)
                            + GET_N(tmp+offset+i)
                            + GET_NE(tmp+offset+i)
                            + GET_W(tmp+offset+i)
                            + GET_E(tmp+offset+i)
                            + GET_SW(tmp+offset+i)
                            + GET_S(tmp+offset+i)
                            + GET_SE(tmp+offset+i)
                            ;

        /* row's eastern most cell */
        offset = (j * cols) + cols - 1;
        dst[offset] = 
                        GET_NW(tmp+offset) + GET_N(tmp+offset)
                      + GET_W(tmp+offset)
                      + GET_SW(tmp+offset) + GET_S(tmp+offset)
                      ;
    }
}

/* Compute the weighted sum of all the cells that match a value,
** returns a scalar value.
*/
float SPATIALweightedSum(float *w, unsigned char *src, int count, int val)
{
    int i;
    float total = 0, gtotal;

    for (i=0; i<count; i+=1)
        total += (src[i] == val) ? w[i] : 0.0;

    MPI_Reduce(&total, &gtotal, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gtotal;
}


float SPATIALfalseDev(float w, float *probmap, float *density,
                      float *ranvals, int count)
{
    int i;
    float total = 0, gtotal;

    // check the density map to catch possible nodata values
    for (i=0; i<count; i+=1)  {
        if (ranvals[i] < powf(probmap[i], w))  {
            total += (density[i] > MIN_DENSITY) ? density[i] : 0.0;
        }
    }

    MPI_Reduce(&total, &gtotal, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gtotal;
}


float SPATIALfalseDev2(float w, float best, float *probmap, float *density,
                      float *ranvals, int count)
{
    int i, cells = 0;
    float p, total = 0, gtotal;

    // check the density map to catch possible nodata values
    for (i=0; i<count; i+=1)  {
        p = (probmap[i] * w > best) ? best : probmap[i] * w;
        if (ranvals[i] < p)  {
            cells += 1;
            total += (density[i] > MIN_DENSITY) ? density[i] : 0.0;
        }
    }

    MPI_Reduce(&total, &gtotal, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gtotal;
}


/* Compute the total number of cells that match a value, returns
** a scalar value.
*/
int spatialCount(unsigned char *src, int count, int val)
{
    int i, total = 0, gtotal;

    for (i=0; i<count; i+=1)
        total += (*(src+i) == val) ? 1 : 0;

    MPI_Reduce(&total, &gtotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return gtotal;
}

int spatialCountGreaterF(float *map, int count, float val)
{
   int i, total = 0, gtotal;

   for (i=0; i<count; i+=1)
       total += ((map[i] > val) ? 1 : 0);

    MPI_Reduce(&total, &gtotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return gtotal;
}



/* Compute the total number of cells that match a value, returns
** a array with the total number of cells associated with each
** class used within the map.
*/
void spatialCorrelatedCount(int *totals, int len, int *map,
                           unsigned char *src, int count, int val)
{
    int i, *gtotals = NULL;

    if (debug && myrank == 0)
        fprintf(stderr, "spatialCorrelatedCount called, len = %d\n", len);

    gtotals = (int *)getMem(len * sizeof (int), "CorrelatedCount array");

    for (i=0; i<count; i+=1)  {
        if (map[i] < 0 || map[i] >= len)  {
            continue;
        }
        if (src[i] == val) totals[map[i]] += 1;
    }

    MPI_Reduce(totals, gtotals, len, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    for (i=0; i<len; i+=1)  totals[i] = gtotals[i];
    MPI_Bcast(totals, len, MPI_INT, 0, MPI_COMM_WORLD);

    free(gtotals);
    return;
}

/* Compute the sum of all values that are associated with each class
** used within the map.
*/
void spatialCorrelatedSumF(float *totals, int len, int *map,
      float *src, int count)
{
      int i;
      float *sums, *gsums;

      sums = (float *)getMem(len * sizeof (float), 
                    "vector for spatialCorrelatedSumF");
      gsums = (float *)getMem(len * sizeof (float),
                    "vector for spatialCorrelatedSumF");

      for (i=0; i<count; i+=1)
           sums[map[i]] += src[i];

      MPI_Reduce(sums, gsums, len, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Bcast(gsums, len, MPI_FLOAT, 0, MPI_COMM_WORLD);

      if (myrank == 0)  {
          printf("GRID_ID,     Delta Pop\n");
          for (i=1; i<len; i+=1)
              printf("%6d, %6.3f\n", i, gsums[i]);
      }
        
      freeMem(sums);
      return;
}

/* Perform a summation over entire grid and return a scalar value.
*/
float spatialSumF(float *src, int count)
{
    int i;
    float total = 0.0, gtotal;

    for (i=0; i<count; i+=1)
        total += *(src+i);

    MPI_Reduce(&total, &gtotal, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gtotal;
}

/* Perform a summation over entire region where any value
** greater than 1 counts as 1.  Return a scalar value.
*/
float spatialSumOneF(float *src, int count)
{
    int i;
    float total = 0.0, gtotal;

    for (i=0; i<count; i+=1)
        total += (*(src+i) > 1.0) ? 1.0 : *(src+i);

    MPI_Reduce(&total, &gtotal, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gtotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gtotal;
}

/* Search for maximum value over entire region and return scalar value.
*/
float spatialMaxF(float *src, int count)
{
    int i;
    float max, gmax;

    max = src[0];
    for (i=1; i<count; i+=1)
        if (src[i] > max)
            max = src[i];

    MPI_Reduce(&max, &gmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gmax, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gmax;
}

float spatialMinF(float *src, int count, int gt_zero)
{
    int i;
    float min, gmin;

    // find min that's greater than 0.0
    if (gt_zero)  {
        for (i=0; i<count; i+=1)  // find starting min
            if (src[i] > 0.0)
                min = src[i];

        for (i; i<count; i+=1)    // find local min thats > 0.0
            if ((src[i] > 0.0) && (src[i] < min))
                min = src[i];
    }

    // find min any min
    else  {
        min = src[0];
        for (i=1; i<count; i+=1)
            if (src[i] < min)
                min = src[i];
    }

    MPI_Reduce(&min, &gmin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gmin, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return gmin;
}

/* Normalize grid */
void spatialNormalizeF(float *dst, float *src, int count)
{
    int i;
    float max;

    max = spatialMaxF(src, count);

    if (max != 0.0)
        for (i=0; i<count; i+=1)
             dst[i] = src[i] / max;

    return;
}


/* Search for maximum value over entire region and return scalar value.
*/
int spatialMax(int *src, int count)
{
    int i, max, gmax;

    max = src[0];
    for (i=1; i<count; i+=1)
        if (src[i] > max)
            max = src[i];

    MPI_Reduce(&max, &gmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return gmax;
}

/* Simple Histogram - histogram array is pre-allocated and should
** be zero'ed out prior to call.
*/
void spatialHistogram(int *hist, int len, int *src, int count)
{
    int i, *ghist = NULL;

    if (debug && myrank == 0)
        fprintf(stderr, "spatialHistogram called\n");

    ghist = (int *)getMem(len * sizeof (int), "temporary hist array");

    for (i=0; i<count; i+=1)  {
        if (src[i] < 0 || src[i] >= len)
            continue;
        else
            hist[src[i]] += 1;
    }

    MPI_Reduce(hist, ghist, len, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    for (i=0; i<len; i+=1)  hist[i] = ghist[i];
    MPI_Bcast(hist, len, MPI_INT, 0, MPI_COMM_WORLD);

    free(ghist);
    return;
}


// Logrithmic Histogram - performs a log function on floating point array
// and uses this to bin the cells.
// 
#define BINS 40
void SPATIALhistogramLog(float *src, int count)
{
  int i, bins[BINS], gbins[BINS], x, zeros=0, gt_one=0;
  float min, max;
  float ln2 = -1 * log(2.0);

  for (i=0; i<BINS; i+=1) bins[i] = 0;

  for (i=0; i<count; i+=1)  {
    if (src[i] <= 0.0)
      zeros += 1;
    else if (src[i] > 1.0)
      gt_one += 1;
    else {
      x = (int)(logf(src[i])/ln2);
      if (x >= BINS) x = BINS - 1;
      bins[x] += 1;
    }
  }

  MPI_Reduce(bins, gbins, BINS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  min = spatialMinF(src, count, 1);
  max = spatialMaxF(src, count);

  if (myrank == 0)  {
    fprintf(stderr, "Histogram: max = %g, min = %g ", max, min);
    fprintf(stderr, "with %d < 0.0 and %d > 1.0\n", zeros, gt_one);
    for (i=0; i<BINS; i+=1)  {
      if (gbins[i] > 0)
        fprintf(stderr, "\t%d:\t%d\n", i, gbins[i]);
    }
    fprintf(stderr, "Histogram done\n");
  }
}

/*
** Originally this was called the utility model but we're now calling
** it the diffision model with is more discriptive and accurate.
**
** Generate a utility model that diffuses energy from developed
** cells to nearby cells.  Energy continues to diffuse over time.
**
** UTILITY(t) = UTILITY(t - dt) + (U_DIFFUSE_IN - U_DIFFUSE_OUT) * dt
** INIT UTILITY = IF (UTILITY_MAP = 1) THEN U_LEVEL  ELSE 
** IF (UTILITY_MAP = 2) THEN U_LEVEL  ELSE 0
**
** INFLOWS:
** U_DIFFUSE_IN = IF UTILITY_MAP=0 THEN 
**                    (IF (UTILITY >= U_LEVEL) THEN
**                         U_DIFFUSE_OUT 
**                     ELSE 
**                         (IF (U_LEVEL-UTILITY) >=  TOTAL_U_IN THEN 
**                              TOTAL_U_IN 
**                          ELSE U_LEVEL-UTILITY)
**                    ) 
**                ELSE 
**                    U_LEVEL-UTILITY+U_DIFFUSE_OUT
** OUTFLOWS:
** U_DIFFUSE_OUT = UTILITY*U_DIFFUSE_RATE
**
** TOTAL_U_IN = U_E@W+U_N@S+U_NE@SW+U_NW@SE+U_S@N+U_SE@NW+U_SW@NE+U_W@E
** UTILITIES = RANDOM(0,UTILITY/U_LEVEL)
** UTILITY_MAP = RESIDENTIAL+2*COM_IND
**
** Ok, if this baffles you then you're not alone.  I think it's jumping
** through hoops to clamp utility at at U_LEVEL (which I'll hardcode
** at 1.0) because there always has to be an outflow from the stock?
**
** 'type' parameter allows diffusion from any set of land cover types.
** 
*/
int spatialDiffusion(float *src, float *tmp, float rate, int type,
                    unsigned char *luptr, int rows, int cols)
{
   int i, j, offset;

   /* Set utilities to max for specified cells.  This jumps newly 
   ** developed cells to the max level so they can begin diffusing 
   ** in earnest.
   */
   for (i=0; i<rows*cols; i+=1)  {
       tmp[i] = 0.0;
       switch (*(luptr+i))  {
       case LU_LRES: case LU_HRES:
           if (type & RES_FLAG) *(src+i) = 1.0;
           break;

       case LU_COM:
           if (type & COM_FLAG) *(src+i) = 1.0;
           break;

       case LU_ROAD:
           if (type & ROAD_FLAG) *(src+i) = 1.0;
           break;

       case LU_OS:
           if (type & OS_FLAG) *(src+i) = 1.0;
           break;
       }
   }

   for (j=0; j<rows; j+=1)  {

       /* row's left most cell */
       offset = j * cols;
       *(tmp+offset) = rate *
                      ( 
                      + GET_N(src+offset) + GET_NE(src+offset)
                      +                      GET_E(src+offset)
                      + GET_S(src+offset) + GET_SE(src+offset)
                      ) / 8.0;

       /* row's middle cells */
       for (i=1; i<cols-1; i+=1)
           *(tmp+offset+i) += rate * 
                      ( GET_NW(src+offset+i)
                      + GET_N(src+offset+i) 
                      + GET_NE(src+offset+i) 
                      + GET_W(src+offset+i) 
                      + GET_E(src+offset+i) 
                      + GET_SW(src+offset+i) 
                      + GET_S(src+offset+i) 
                      + GET_SE(src+offset+i) 
                      ) / 8.0;

       /* row's right most cell */
       offset = (j * cols) + cols - 1;
       *(tmp+offset) = rate *
                      ( 
                      + GET_NW(src+offset) + GET_N(src+offset)
                      + GET_W(src+offset)
                      + GET_SW(src+offset) + GET_S(src+offset)
                      ) / 8.0;
   }

   /* clamp utilities at 1 */
   for (i=0; i<rows * cols; i+=1)
#ifdef NOCLAMP
       src[i] += tmp[i];
#else
       src[i] = (src[i] + tmp[i] > 1.0) ? 1.0 : src[i] + tmp[i];
#endif

   shareGrid((char *)src, rows*cols, MPI_FLOAT);
}
