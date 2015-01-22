# Simple Makefile for building the LEAM land use change (gluc) model.
# The model is written using MPI for parallelization and has been compiled
# on a variety of platforms and libraries.  Some example flags used
# on various system are provided below.  

# Copyright (C) 2013 LEAMgroup, Inc. Released under GPL v2.


# Standard C-compiler
#CC = cc

# GNU C-compiler
#CC = gcc

CC = mpicc
#CC = /usr/lib64/lam/bin/mpicc
#MPIINC = -I/usr/local/mpich/include
#MPILIBS = -L/usr/local/mpich/lib -lmpich

# NCSA Platinum Cluster
#CC = gcc
#MPIINC = -I/usr/local/vmi/mpich/include
#MPILIB = -L/usr/local/vmi/mpich/lib/$(CC) -lmpich -lvmi -ldl -lpthread

# Alpha Cluster
#CC = cc
#MPIINC = -I.
#MPILIB = -lmpi -lelan

# SGI Origin
#CC = cc
#MPIINC = 
#MPILIB = -lmpi
#ARCH_CFLAGS = -64

# AIX C-compiler with MPI support
#CC = mpcc
#ARCH_CFLAGS = -Dpowf=pow

CFLAGS = $(ARCH_CFLAGS) $(MPIINC) -DGLUC -DMYGATHER
LIBS = $(MPILIB)  -lexpat -lm

SRCS = leam.c utilities.c luc.c SME.c bil.c spatial.c graph.c score.c GA.c http.c
OBJS = leam.o utilities.o luc.o SME.o bil.o spatial.o graph.o score.o GA.o http.o

gluc: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)
	@echo gluc `Built ./gluc --version`

graph.o: graph.c graph.h

clean:
	-rm gluc *.o

tags: $(SRCS)
	ctags $(SRCS)
