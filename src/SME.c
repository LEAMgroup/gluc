/*
** This module provides utility function to allow the model to emulate
** the original SME modeling approach.  SME was a tool that allowed Stella
** models to be automatically spatialized.  Because many configuration
** files had already been created the gluc model continued using the format.
**
** Note: the current configuration file is kludgy at best and should be
** replaced with something more modern.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "leam.h"

/* Start off badly by statically allocating the number of variables
** that can be read/processed from the SME configuration file.
*/
static struct {
    char  varname[80];
    char  *mapname;
    float fdata;
} vardata[4000]; 

static int idx = 0;
static char *datapath, *mappath;


const char *SMEgetDataPath()
{
    return datapath;
}

const char *SMEgetMapPath()
{
    return mappath;
}


void SMEaddVar(char *name, char *value)
{
    int i;

    /* check for existing instance of the variable and replace it */
    for (i=0; i<idx; i+=1)  {
        if (!strcmp(vardata[i].varname, name))  {
            free(vardata[i].mapname);
            vardata[i].mapname = strdup(value);
            if (vardata[i].mapname == NULL)  {
                sprintf(estring, "out of memory allocating %s\n", name);
                errorExit(estring);
            }
            return;
        }
    }

    /* new varible, so add it */
    strcpy(vardata[idx].varname, name);
    vardata[idx].mapname = strdup(value);
    if (vardata[idx].mapname == NULL)  {
        sprintf(estring, "out of memory allocating %s\n", name);
        errorExit(estring);
    }
    idx += 1;
}

/* Parse the SME configuration file and build variable database.
*/
void SMEparseConfig(char *path, char *proj, char *fname)
{
    char line[1024], fstring[1024], varname[80], *ptr;
    char *delimit = " \t\n(),";
    FILE *f;

    if (path == NULL || proj == NULL)
        strcpy(fstring, fname);
    else
        sprintf(fstring,"%s/%s/Config/%s", path, proj, fname);

    if ((f = fopen(fstring, "r")) == NULL)  {
        sprintf(estring, "unable to open SME config file %s", fstring);
        errorExit(estring);
    }

    while (fgets(line, sizeof line, f) != NULL)  {

        if ((ptr = strtok(line, delimit)) == NULL)
            continue;

        /* global parameters */
        if (!strcmp(ptr, "#") && !strcasecmp(strtok(NULL, delimit), 
                "global"))  {
            while ((ptr = strtok(NULL, delimit)) != NULL)  {
                if (!strcasecmp(ptr, "s"))
                    srand48(atol(strtok(NULL, delimit)));
                else if (!strcmp(ptr, "d"))
                    debug = atol(strtok(NULL, delimit)); 
                else if (!strcmp(ptr, "OT"))  {
                    SMEaddVar("TIMESTEP", strtok(NULL, delimit));
                    SMEaddVar("START_DATE", strtok(NULL, delimit));
                    SMEaddVar("END_DATE", strtok(NULL, delimit));
                }
            }
        }

        /* module level parameters */
        else if (!strcmp(ptr, "$"))  {
            while ((ptr = strtok(NULL, delimit)) != NULL )  {
                if (!strcasecmp(ptr, "g"))  {
                    strtok(NULL, delimit);  /* format */
                    sprintf(fstring, "%s/%s", datapath, strtok(NULL, delimit));
                    SMEaddVar("STUDY_AREA", fstring);
                }
            }
        }

        /* variable level parameters */
        else if (!strcmp(ptr, "*"))  {
            strcpy(varname, strtok(NULL, delimit));
            while ((ptr = strtok(NULL, delimit)) != NULL)  {
                if (!strcasecmp(ptr, "pm"))  {

                   // ugly hack to handle cases where
                   // config files where created with
                   // floating point values on _MAP variables
                   if (strstr(varname, "_MAP") != NULL)
                       *strstr(varname, "_MAP") = '\0';
                   SMEaddVar(varname, strtok(NULL, delimit));
                }
                else if (!strcasecmp(ptr, "c") || !strcasecmp(ptr, "d"))  {
                   strtok(NULL, delimit);  // format
                   sprintf(fstring, "%s/%s", datapath, strtok(NULL, delimit));
                   SMEaddVar(varname, fstring);
                }

                else if (!strcasecmp(ptr, "M"))  {
                   strtok(NULL, delimit);  // format
                   strtok(NULL, delimit);  // scale
                   sprintf(fstring, "%s/%s", mappath, strtok(NULL, delimit));
                   SMEaddVar(varname, fstring);
                }

                // non-standard variable for URL types
                else if (!strcasecmp(ptr, "URL"))  {
                   SMEaddVar(varname, strtok(NULL, delimit));
                }

                // non-standard variable for string types
                else if (!strcasecmp(ptr, "s"))  {
                   SMEaddVar(varname, strtok(NULL, delimit));
                }

                // parse date in "YYYY MMM DD" format, only year is kept
                // and it is kept in a string format.
                else if (!strcasecmp(ptr, "DATE")) {
                   SMEaddVar(varname, strtok(NULL, delimit));
                }

                else if (!strcasecmp(ptr, "GRAPH"))  {
                   sprintf(fstring, "%s/%s", datapath, strtok(NULL, delimit));
                   GRAPHreadFile(fstring);
                }
            }
        }
    }

    fclose(f);
}

char *SMEgetFileName(char *var)
{
    int i;

    for (i=0; i<idx; i+=1)
        if (!strcmp(var, vardata[i].varname))
            return vardata[i].mapname;

    return NULL;
}

char *SMEgetString(char *var, char *def)
{
    int i;

    for (i=0; i<idx; i+=1)  {
        if (!strcmp(var, vardata[i].varname))
            return strdup(vardata[i].mapname);
    }

    if (def != NULL)
        return strdup(def);

    else
        return NULL;
}

int SMEgetInt(char *var, int def)
{
    int i, x;
    char *endptr;

    for (i=0; i<idx; i+=1)  {
        if (!strcmp(var, vardata[i].varname))  {
            x = strtol(vardata[i].mapname, &endptr, 10);
            if (vardata[i].mapname == endptr)
                return def;
            else
                return x;
        }
    }

    return def;
}

double SMEgetFloat(char *var, float def)
{
    int i;
    double x;
    char *endptr;

    for (i=0; i<idx; i+=1)  {
        if (!strcmp(var, vardata[i].varname))  {
            x = strtod(vardata[i].mapname, &endptr);
            if (vardata[i].mapname == endptr)
                return def;
            else
                return x;
        }
    }

    return def;
}


float SMEgetData(char *var, float dval)
{
    int i;

    for (i=0; i<idx; i+=1)
        if (!strcmp(var, vardata[i].varname))  {
            return vardata[i].fdata;
        }

    return dval;
}

/* Parse out the options associated with SME libraries.
** Other options are ignored.
*/
void SMEparseOptions(int argc, char *argv[])
{
    char fstring[1024];
    char *ppath = NULL;
    char *project = NULL;
    char *model = NULL;
    char *scenario = NULL;
    char *cfile = NULL;
    char *cptr;
    int i;

    for (i=1; i<argc; i+=1)  {
        if (!strcmp(argv[i], "-ppath"))  {
            ppath = strdup(argv[++i]);
        }

        else if (!strcmp(argv[i], "-p"))  {
            project = strdup(argv[++i]);
        }

        else if (!strcmp(argv[i], "-m"))  {
            model = strdup(argv[++i]);
        }

        else if (!strcmp(argv[i], "-scen"))  {
            scenario = strdup(argv[++i]);
        }

        else if (!strcmp(argv[i], "-ci"))  {
            cfile = strdup(argv[++i]);
        }
    }

    if (ppath == NULL || project == NULL)  {
        datapath = mappath = strdup(".");
        fprintf(stderr, "P%d: Project and/or Project Path unspecified,", 
                myrank);
        fprintf(stderr, "working in local directory.\n");
    }

    else  {
        sprintf(fstring, "%s/%s/Data", ppath, project);
        datapath = strdup(fstring);
        sprintf(fstring, "%s/%s/DriverOutput/Maps", ppath, project);
        mappath = strdup(fstring);

        if (debug)  {
             fprintf(stderr, "P%d: Data path set to %s\n", myrank, datapath);
             fprintf(stderr, "P%d: Map path set to %s\n", myrank, mappath);
             fprintf(stderr, "P%d: model='%s', scenario='%s'\n", myrank, 
                     model, scenario);
        }
    }

    if (cfile == NULL)  {
        sprintf(estring, "no config file specified on command line");
        errorExit(estring);
    }
    else  {
        runName = strdup(cfile);
        if ((cptr = strrchr(runName, '.')) != NULL)
            *cptr = '\0';
        SMEparseConfig(ppath, project, cfile);
    }
}
