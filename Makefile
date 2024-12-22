#
# Makefile for MAIN
#
SHELL = /bin/sh
#F77FLAGS =  -check_bounds
OBJS = opm.o min.o rmsd.o solva.o readpdb.o read_small.o watface.o \
	 tilting.o locate.o profile.o hbcor.o deftm.o orient.o find_segm.o
EXEC = immers

load: $(OBJS)
	 gfortran-9 -O -g -o $(EXEC) $(OBJS)
.f.o:
	 gfortran-9 $(F77FLAGS) -O -c $<
