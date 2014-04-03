/*
** This code was an attempt to use genetic algorithms to determine
** code weights.  It should be removed and can safely be ignored.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <ctype.h>
#include <expat.h>
#include <mpi.h>

#include "leam.h"
#include "http.h"

/* the GA engine URL */
static char GA_URL[256];

/* the gene ID */
static char GA_ID[33];

/* some constants for constucting the url methods */
/* get_baby command looks like:
   GA_URL/cgi-bin/GAEngine.py?op=get_baby" */
/* mature_baby command looks like:
   GA_URL/cgi=bin/GAEngine.py?op=mature_baby&id=GA_ID&fit=ID_FIT */
#define MAX_URL_LEN 1024
#define ROOT_CMD "cgi-bin/GAEngine.py"
#define GET_BABY "op=get_baby"
#define MATURE_BABY "op=mature_baby"
#define FIT_STUB "fit="
#define ID_STUB "id="
#define MATURE_BABY_OK "OK\n"

/* structure for holding the weight/value pairs */
#define WEIGHT_LEN 256
typedef struct _weight_value_
{
  char weight[WEIGHT_LEN];
  double value;
} weight_value;
#define MAX_WEIGHTS 1024
static weight_value weights[MAX_WEIGHTS];
static unsigned nweights = 0;

/* string case normalizers...converts string to uppercase */
static char *upcase(char *str)
{
  char *cptr;
  cptr = str;
  while (cptr && *cptr)
    {
      *cptr = toupper(*cptr);
      cptr++;
    }
  return str;
}

/************************************************************
 * helper routines for accessing the array of weight_values *
 ************************************************************/

/* determine if a weight is in the array return the index if it is 
   or -1 if it is not */
static int isWeight(const char *str)
{
  int i;
  /* if weights is empty it isn't in it */
  if (nweights == 0) return -1;

  /* look for it in the array if it is there return the index */
  for (i = 0; i < nweights; i++)
    if (!strcasecmp(weights[i].weight, str)) 
      return i;

  /* if we get here it wasn't in the array */
  return -1;
}

/* insert a weight into the array overwriting value if it already exists */
/* return 1 on success, 0 on failure */
static int insertWeight(const char *str, double val)
{
  int i;
  char tmpstr[WEIGHT_LEN], *cptr;
  /* check to see if we already have filled it up */
  if (nweights >= MAX_WEIGHTS)
    {
      fprintf(stderr, "Error: InsertWeight - max values exceeded\n");
      return 0;
    }

  i = isWeight(str);
  /* if index is -1 then it is new weight else overwrite old value */
  if (i == -1)
    {
      /* normalize to upper case */
      strcpy(tmpstr, str);
      upcase(tmpstr);
      strcpy(weights[nweights].weight, tmpstr);
      weights[nweights].value = val;
      /* incr # of weights in array */
      nweights++;
    }
  else weights[i].value = val;
  return 1;
}

/***************************************
 * helper routines for parsing the XML *
 ***************************************/
/* global variable to know when we are in the variable section */
static unsigned char varFlag = 0;

/* event handler for beginning tags */
static void startHandler(void *data, const char *el, const char **attr)
{
  /* if we are the start of the variables section flag it */
  if (!strcasecmp(el, "VARIABLES")) varFlag = 1;
  /* if we are <ID> get value into GA_ID otherwise
     if we are in <variable> region insert the VALUES for the keys */
  else if (!strcasecmp(el, "ID"))
    { 
      if (!strcasecmp(attr[0], "VALUE")) strcpy(GA_ID,attr[1]);
    }
  else if (varFlag)
    if (!strcasecmp(attr[0], "VALUE"))
      insertWeight(el,atof(attr[1])); /* insert into array of weights*/
}
/* event handler for ending tags */
static void endHandler(void *data, const char *el)
{
  /* reset the varFlag when we close the </variables> section */
  if (!strcasecmp(el,"VARIABLES")) varFlag = 0;
}
/* actual string parser */
void XMLParseString(char *str, int len)
{
  XML_Parser p; /* the parser */

  /* create the XML SAX parser */
  if ((p = XML_ParserCreate(NULL)) == NULL)
    {
      sprintf(estring,"Error allocating memory for parser\n");
      errorExit(estring);
    }

  /* set the tag handlers */
  XML_SetElementHandler(p,startHandler,endHandler);

  /* now actually parse through the string */
  if (!XML_Parse(p, str, len, 1))
    {
      sprintf(estring, "Parse Error at line %d:\n%s\n",
	      XML_GetCurrentLineNumber(p),
	      XML_ErrorString(XML_GetErrorCode(p)));
      errorExit(estring);
    }
  XML_ParserFree(p);
}

/****************************************************************
 **                        Main API                            **
 ****************************************************************/
/* initialize the interface to the GA engine by contacting the engine
   at url and getting a set of weights in the form of an XML file. 
   Parse out the weights and the ID and put them into static global 
   variables */
void GAinit(const char *url)
{
  char url_cmd[MAX_URL_LEN];
  HTTP_Response hResponse;

  if (debug && myrank == 0)
      fprintf(stderr, "GAinit: url = %s\n", url);

  /* zero out the weight array and nweights */
  {
    unsigned i;
    for (i=0; i<MAX_WEIGHTS;i++)
      {  
	memset(weights[i].weight,0,WEIGHT_LEN);
	weights[i].value = 0.0;
      }
    nweights = 0;
  }

  /* copy url into GA_URL for future use */
  strcpy(GA_URL,url);
  if (GA_URL[strlen(GA_URL) - 1] == '/' ||
      GA_URL[strlen(GA_URL) - 1] == '\\' ) 
    GA_URL[strlen(GA_URL) - 1] = '\0';
  
  /* construct the get baby cmd url */
  sprintf(url_cmd, "%s/%s?%s", GA_URL, ROOT_CMD, GET_BABY);

  /* actually get the XML document from the server */
  hResponse = http_request(url_cmd, NULL, kHMethodGet, HFLAG_NONE);
  /* if we got the XML document so call the parser to fill in weights array */
  if (hResponse.lSize > 0)
      XMLParseString(hResponse.pData, hResponse.lSize);
  else /* otherwise there was an error */
    {
      sprintf(estring,"Error GAInit: Error retrieving baby\n");
      errorExit(estring);
    }

  /* we have to free up the pData from the response when done */
  if (hResponse.pData) free(hResponse.pData);
}

/* return the value for the weight given by str or dflt if it doesn't exist */
double GAgetData(const char *str, double dflt)
{
  double val = dflt;
  int i;
  /* get the index of the weight if it exists */
  i = isWeight(str);
  if (i != -1) /* if index not -1 we found weight set val to the value */
    val = weights[i].value;

  return val;
}

/* return the fitness to the GA engine */
void GAsendFit(double fit)
{
  char url_cmd[MAX_URL_LEN];
  HTTP_Response hResponse;

  if (debug)
     fprintf(stderr, "GAsendFit: fitness = %f\n", fit);

  /* construct the mature baby cmd url */
  sprintf(url_cmd, "%s/%s?%s&%s%s&%s%f",
	  GA_URL, ROOT_CMD, MATURE_BABY,
	  ID_STUB, GA_ID, FIT_STUB, fit);
  
  /* actually send the fitness for kid GA_ID to the server */
  hResponse = http_request(url_cmd, NULL, kHMethodGet, HFLAG_NONE);
  /* check to see if we got OK as a response */
  if (hResponse.lSize <= 0 || strcasecmp(hResponse.pData, MATURE_BABY_OK))
    {
      /* there was an error */
      sprintf(estring,"Error GAsendFit: Error sending fit\n");
      errorExit(estring);
    }
  /* free up the pData from the response */
  if (hResponse.pData) free(hResponse.pData);
}
