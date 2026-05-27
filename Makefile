#
# Makefile for MAIN
#
#SHELL = /bin/sh
#F77FLAGS =  -check_bounds

FFLAGS=-O -g
LDFLAGS=-g

OBJS = opm.o min.o rmsd.o solva.o readpdb.o read_small.o watface.o \
	 tilting.o locate.o profile.o hbcor.o deftm.o orient.o find_segm.o
EXEC = immers

FC=gfortran 

$(EXEC): $(OBJS)
	 $(FC)  $(LDFLAGS) -o $(EXEC) $(OBJS)

%.o: %.f
	 $(FC) $(FFLAGS) -c $<

clean:
	-rm $(EXEC) $(OBJS)
