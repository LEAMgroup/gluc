/*
** Early versions of the LEAM model handled all scoring of the drivers
** internally and this module provided many of utility function to
** perform these operations.  It's still used population and employment
** projections but not much more.  It could probably be simplified.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/

#define MAXGRAPHS    4096            // maximum number of graphs
#define GRAPHBUFSIZE 16 * 1024       // maximum size of a single graph


// default NearestNeighbors graph
static float NearestNeighbors[] = {         
  0, 1.000, 
  1, 1.000,
  2, 1.000,
  3, 1.000,
  4, 1.000,
  5, 0.100,
  6, 0.010,
  7, 0.001,
  8, 0.0001,
};

static float CellDemand[] = {
  1992, 0,
  2000, 10000,
  2050, 50000,
};

