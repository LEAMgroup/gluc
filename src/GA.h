/*
** This code was an attempt to use genetic algorithms to determine
** code weights.  It should be removed and can safely be ignored.
**
** Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.
*/

#ifndef GA_H
#define GA_H
void GAinit(const char *url);
double GAgetData(const char *str, double dflt);
void GAsendFit(double fit);
#endif
