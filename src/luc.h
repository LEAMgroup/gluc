/*
** Gereralized Land Use Change Model (gluc)
**
** This is primary code for the land use change model.  The model is
** primarily driven by a probability map and a land use map.  For each
** developed cell on the land use map a probability of development is
** calculated and a random number used to determine if the cell develops
** or not.  A feedback mechanism controls the rate of development increasing
** or decreasing the overall probability based on the desired demand.
**
** Note: he model has evolved significantly over the years and a lot cruf
** has accumulated.  It could use a good cleaning.
**
** Note: this file should be rewritten.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/

/* Pre-Processor Macros to reference varaibles */
#ifdef XXX
/* Scalar Variables */
int time = 0;

/* Spatial Model Variables */
unsigned char *boundary;
unsigned char *landuse_map;
unsigned char *landuse;
#endif

