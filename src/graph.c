/*
** Early versions of the LEAM model handled all scoring of the drivers 
** internally and this module provided many of utility function to
** perform these operations.  It's still used for population and employment
** projections but not much more.  It could probably be simplified.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "leam.h"
#include "graph.h"

// maximum size of graphs
#define DATABUFSIZE 16 * 1024

typedef struct {
  char name[256];
  int  len;
  float *data; 
} GRAPH_T;


static int gidx = 0;
static GRAPH_T graphs[4096];


#define G(x,y) *(data+2*(x)+(y))
static float graph(float *data, int count, float val)
{
   int i, rows = count / 2;

   if (val <= G(0,0))
       return G(0,1);
   else if (val >= G(rows-1,0)) 
       return G(rows-1,1);
   else 
       for (i=0; i<rows; i+=1)
           if (val >= G(i,0) && val < G(i+1,0))
               return (val - G(i,0)) * (G(i+1,1) - G(i,1)) /
                      (G(i+1,0) - G(i,0)) + G(i,1);

   errorExit("bad graph data?");
}



/*
** DEBUG ROUTINE -- dumps the current graph table
*/
void dumpGraphs(FILE *f, int verbose)
{
  int i, j;

  for (i=0; i<gidx; i+=1)  {
    fprintf(f, "GRAPH %s: len = %d\n", graphs[i].name, graphs[i].len);
    if (verbose)  {
      for (j=0; j<graphs[i].len; j+=2)
        fprintf(f, "%8.2f, %8.2f\n", graphs[i].data[j],
                graphs[i].data[j+1]);
      fprintf(f, "\n");
    }
  }
}


void *GRAPHgetGraph(char *name)
{
  int i;

  for (i=0; i<gidx; i+=1)  {
    if (!strcmp(graphs[i].name, name))
      return (void *)(graphs+i);
  }

  return NULL;
}

float *GRAPHgetGraphData(char *name)
{
  GRAPH_T *g;

  g = (GRAPH_T *)GRAPHgetGraph(name);
  if (g != NULL)
    return g->data;
  else
    return NULL;
}

int GRAPHdump(void *ptr)
{
  int i;
  GRAPH_T *g= (GRAPH_T *)ptr;

  fprintf(stderr, "GRAPHdump: ");
  for (i=0; i<g->len; i+=1)
    fprintf(stderr, "%f, ", (g->data)[i]);
  fprintf(stderr, "\n");

}

float GRAPHlookup(void *ptr, float v)
{
  GRAPH_T *g = (GRAPH_T *)ptr;

  return graph(g->data, g->len, v);
}

float GRAPHinterp(void *ptr, float v)
{
  GRAPH_T *g = (GRAPH_T *)ptr;

  return graph(g->data, g->len, v);
}

void GRAPHaddGraph(char *name, int count, float *data)
{
  GRAPH_T *g;

  if ((g = (GRAPH_T *)GRAPHgetGraph(name)) == NULL)  {
    strcpy(graphs[gidx].name, name);
    graphs[gidx].len = count;
    graphs[gidx].data = (float *)getMem(count * sizeof (float), "graph data");
    memcpy(graphs[gidx].data, data, count * sizeof (float));
    gidx += 1;
  }

  else  {
    freeMem(g->data);
    g->len = count;
    g->data = (float *)getMem(count * sizeof (float), "graph data");
    memcpy(g->data, data, count * sizeof (float));
  }
}


void GRAPHreadFile(char *fname)
{
  FILE *f;
  int count = 0, findname = 1;
  char name[256], l[256];
  float databuf[DATABUFSIZE];

  fprintf(stderr, "GRAPHreadFile: reading %s\n", fname);

  if ((f = fopen(fname, "r")) == NULL)  {
    sprintf(estring, "unable to open '%s'\n", fname);
    errorExit(estring);
  }

  while (fgets(l, sizeof l, f) != NULL)  {
    if (l[0] == '#') continue;

    // seeking the next graph name
    if (findname)  {
      if (sscanf(l, "%s", name) == 1)
        findname = 0;
    }
        
    // continue scanning of tuples
    else if (sscanf(l, "%f,%f", databuf+count, databuf+count+1) == 2) {
      count += 2;
      if (count > sizeof databuf/sizeof (float))  {
        sprintf(estring, "graph %s exceeds buffer size", name);
        errorExit(estring);
      }
    }

    // hit the end of tuples so store new graph
    else  {
      GRAPHaddGraph(name, count, databuf);
      count = 0;
      findname = 1;   // return to findname mode
    }
  }

  // write out the final graph if necessary
  if (count > 0)
    GRAPHaddGraph(name, count, databuf);

  fclose(f);
}

void GRAPHinit()
{
  int i;

  GRAPHaddGraph("CellDemandRes", sizeof CellDemand / sizeof (float), 
                CellDemand);
  GRAPHaddGraph("CellDemandCom", sizeof CellDemand / sizeof (float), 
                CellDemand);
  GRAPHaddGraph("CellDemandOS", sizeof CellDemand / sizeof (float), 
                CellDemand);
  GRAPHaddGraph("NearestNeighbors", sizeof NearestNeighbors / sizeof (float),
                NearestNeighbors);
}
