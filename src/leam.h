/*
**
** leam.h standard include file with constants that are generally useful
** for LEAM models and applications.
**
** leam.h also includes function prototypes and external variables
** definitions that should probably be moved elsewhere.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#ifndef _LEAM_
#define _LEAM_

#define EQROWS         1
#define EQCELLS        2

// prevents bad density values (typically NODATA) from
// getting into the demand calculation
#define MIN_DENSITY    -100.0

typedef int SUBMODELS;
#define RES_MODEL      ((SUBMODELS)1)
#define COM_MODEL      ((SUBMODELS)1 << 1)
#define OS_MODEL       ((SUBMODELS)1 << 2)

// exit codes
#define EXIT_INSUFFICIENT 1

/* Setup definitions of the various land use classifications.
** Use these instead of hardcoding land use classifications into
** the source code.
*/ 
typedef int LU_VALUE;
#define LU_WATER       ((LU_VALUE)11)
#define LU_LRES        ((LU_VALUE)21)
#define LU_HRES        ((LU_VALUE)22)
#define LU_COM         ((LU_VALUE)23)
#define LU_ROAD        ((LU_VALUE)24)
#define LU_OS          ((LU_VALUE)85)
#define LU_WET         ((LU_VALUE)91)
#define LU_HWET        ((LU_VALUE)92)

/* LU_FLAGs are used to pass groups of land use classifications
** without resorting to arrays or strings.
*/
typedef int LU_FLAG;
#define WATER_FLAG     ((LU_FLAG)1)              /* LU_VALUE 11 */
#define RES_FLAG       ((LU_FLAG)1 << 1)         /* LU_VALUE 21 and 22 */
#define COM_FLAG       ((LU_FLAG)1 << 2)         /* LU_VALUE 23 */
#define ROAD_FLAG      ((LU_FLAG)1 << 3)         /* LU_VALUE 24 */
#define OS_FLAG        ((LU_FLAG)1 << 4)         /* LU_VALUE 86 */
#define AG_FLAG        ((LU_FLAG)1 << 5)         /* LU_VALUE xx xx xx */
#define WETLAND_FLAG   ((LU_FLAG)1 << 6)         /* LU_VALUE 91 92 */
#define DEVELOPED_FLAG  (RES_FLAG | COM_FLAG | ROAD_FLAG)
#define NONDEVELOPABLE_FLAGS    \
        (WATER_FLAG | DEVELOPED_FLAG | WETLAND_FLAG )


/**** Function Prototypes *****/

/* leam.c */
extern int debug;
extern int myrank, nproc;
extern char *runName;
extern float ulx, uly;
extern float xdim, ydim;

/* utilities.c */
extern char estring[];
extern char *getMem(int, char *);
extern void freeMem(void *);
extern void errorExit(char *);

/* luc.c */
extern int gRows, gCols;
extern void selector(unsigned char *, float *, float *, float*);
extern void shareGrid(void *, int, MPI_Datatype);
extern char *initGridMaps(char *, int, int);
extern void readProbmap(float **, int, int, char *, int);
extern void LUCconfigGrids(int, int, int *);
extern void LUCinitGrids();
extern void LUCrun();

/* SME.c */
extern char *SMEgetBoundary(char*);
extern char *SMEgetFileName(char*);
extern char *SMEgetString(char*, char*);
extern float SMEgetData(char*, float);
extern double SMEgetFloat(char*, float);
extern int SMEgetInt(char*, int);
extern void SMEparseConfig(char*, char*, char*);
extern void SMEparseOptions(int, char**);

/* spatial.c */
extern float SPATIALfalseDev(float, float *, float *, float *, int);
extern float SPATIALfalseDev2(float, float, float *, float *, float *, int);
extern void SPATIALhistogramLog(float *, int);
extern float SPATIALweightedSum(float *, unsigned char *, int, int);
extern float spatialSumF(float *, int);
extern float spatialSumOneF(float *, int);
extern float spatialMaxF(float *, int);
extern void spatialNormalizeF(float *, float *, int);
extern int spatialCount(unsigned char *, int, int);
extern void spatialCorrelatedCount(int *, int , int *, unsigned char *, 
                                   int, int);
extern void spatialCorrelatedSumF(float *, int , int *, float *, int);
extern void nearestNeighbors(unsigned char *, unsigned char *, int, unsigned char *, int, int);
extern int spatialDiffusion(float *, float *, float, int, unsigned char *,
                            int, int);

/* score.c */
extern double scoreResults(int *, int , int *, unsigned char *, int );


/* graph.c */
extern void GRAPHinit();
extern void GRAPHreadFile(char *);
extern void GRAPHaddGraph(char *, int, float *);
extern void *GRAPHgetGraph(char *);
extern float *GRAPHgetGraphData(char *);
extern float GRAPHinterp(void *, float);
extern float GRAPHlookup(void *, float);

#endif /* _LEAM_ */
