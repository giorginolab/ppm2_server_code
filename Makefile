#
# Makefile for MAIN
#
#SHELL = /bin/sh
#F77FLAGS =  -check_bounds

FFLAGS=-O -g
LDFLAGS=-g
OMPFLAGS=-fopenmp

OBJS = opm.o min.o rmsd.o solva.o readpdb.o read_small.o watface.o \
	 tilting.o locate.o profile.o hbcor.o deftm.o orient.o find_segm.o
EXEC = immers

FC=gfortran 

$(EXEC): $(OBJS)
	 $(FC) $(OMPFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS)

locate.o solva.o orient.o: override FFLAGS += $(OMPFLAGS)

%.o: %.f
	 $(FC) $(FFLAGS) -c $<

.PHONY: test clean
test: $(EXEC)
	@set -eux; \
	refdir="$(CURDIR)/ref"; \
	workdir=$$(mktemp -d "$${TMPDIR:-/tmp}/ppm2-test.XXXXXX"); \
	trap 'rm -rf "$$workdir"' 0; \
	cp $(EXEC) res.lib test.inp 1rsy.pdb 1gzm.pdb "$$workdir/"; \
	cd "$$workdir"; \
	./$(EXEC) <test.inp >test.out; \
	for output in test.out 1rsyout.pdb 1gzmout.pdb; do \
		diff -u "$$refdir/$$output" "$$output"; \
	done; \
	set +x; \
	echo ok

clean:
	-rm $(EXEC) $(OBJS)
