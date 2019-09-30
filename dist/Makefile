# MAKEFILE FOR MARKOV CHAIN MONTE CARLO FOR A SPECIFIED DISTRIBUTION.

# Copyright (c) 1995-2004 by Radford M. Neal
#
# Permission is granted for anyone to copy, use, modify, or distribute this
# program and accompanying programs and documents for any purpose, provided 
# this copyright notice is retained and prominently displayed, along with
# a note saying that the original programs are available from Radford Neal's
# web page, and note is made of any changes made to the programs.  The
# programs and documents are distributed without any warranty, express or
# implied.  As the programs were written for research purposes only, they have
# not been tested to the degree that would be advisable in any important
# application.  All use of these programs is entirely at the user's own risk.

include ../make.include

programs:	dist-spec dist-mc dist-plt dist-tbl dist-hist dist-grad-test \
		dist-stepsizes dist-initial dist-display dist-est dist-gen \
		dist-genp dist-his
clean:
	rm -f *.o dist-spec dist-mc dist-plt dist-tbl dist-hist dist-grad-test \
	          dist-stepsizes dist-initial dist-display dist-est dist-gen \
		  dist-genp dist-his


include ../util/util.make
include ../mc/mc.make

dist-spec:	dist-spec.o	misc.o log.o formula.o digamma.o rand.o
		$(CC) dist-spec.o misc.o log.o formula.o digamma.o rand.o \
		  -lm -o dist-spec

dist-spec.o:	dist-spec.c	misc.h log.h mc.h dist.h formula.h

dist-grad-test:	mc-grad-test.o dist-mc.o dist-util.o dist-plt.o \
				dist-quantities.o formula.o digamma.o misc.o \
				log.o rand.o quantities.o dist-data.o numin.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o data-trans.o
		$(CC) $(LFLAGS) mc-grad-test.o dist-mc.o dist-util.o dist-plt.o\
		  dist-quantities.o formula.o misc.o log.o rand.o quantities.o \
		  dist-data.o numin.o mc-iter.o mc-heatbath.o mc-metropolis.o \
                  mc-hybrid.o mc-slice.o mc-traj.o mc-util.o mc-quantities.o \
		  data-trans.o digamma.o \
		  -lm -o dist-grad-test

dist-stepsizes:	dist-stepsizes.o dist-mc.o dist-util.o dist-plt.o digamma.o \
				dist-quantities.o formula.o misc.o log.o \
				rand.o quantities.o dist-data.o numin.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o data-trans.o 
		$(CC) $(LFLAGS) dist-stepsizes.o dist-mc.o dist-util.o \
		  dist-quantities.o formula.o misc.o log.o rand.o quantities.o \
		  dist-data.o numin.o mc-traj.o mc-iter.o mc-heatbath.o \
                  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o dist-plt.o\
                  mc-quantities.o data-trans.o digamma.o \
		  -lm -o dist-stepsizes

dist-stepsizes.o: dist-stepsizes.c misc.h log.h mc.h formula.h dist.h

dist-genp:	mc-genp.o       dist-mc.o dist-util.o dist-plt.o digamma.o \
				dist-quantities.o formula.o misc.o log.o \
				rand.o quantities.o dist-data.o numin.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o data-trans.o 
		$(CC) $(LFLAGS) mc-genp.o dist-mc.o dist-util.o \
		  dist-quantities.o formula.o misc.o log.o rand.o quantities.o \
		  dist-data.o numin.o mc-traj.o mc-iter.o mc-heatbath.o \
                  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o dist-plt.o\
                  mc-quantities.o data-trans.o digamma.o \
		  -lm -o dist-genp

dist-mc:	mc.o dist-mc.o	dist-plt.o dist-util.o dist-quantities.o \
				formula.o misc.o log.o rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o dist-data.o numin.o \
				data-trans.o digamma.o
		$(CC) $(LFLAGS) mc.o dist-mc.o dist-util.o dist-quantities.o \
		  formula.o dist-plt.o misc.o log.o rand.o \
		  quantities.o mc-traj.o mc-iter.o mc-heatbath.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o \
		  mc-quantities.o dist-data.o numin.o data-trans.o digamma.o \
		  -lm -o dist-mc

dist-his:	mc-his.o dist-mc.o dist-plt.o dist-util.o dist-quantities.o \
				formula.o misc.o log.o rand.o quantities.o \
				mc-heatbath.o mc-traj.o mc-util.o \
				mc-quantities.o dist-data.o numin.o \
				data-trans.o digamma.o
		$(CC) $(LFLAGS) mc-his.o dist-mc.o dist-util.o \
		  dist-quantities.o \
		  formula.o dist-plt.o misc.o log.o rand.o \
		  quantities.o mc-traj.o mc-heatbath.o mc-util.o \
		  mc-quantities.o dist-data.o numin.o data-trans.o digamma.o \
		  -lm -o dist-his

dist-mc.o:	dist-mc.c	misc.h log.h dist.h mc.h formula.h \
				data.h dist-data.h

dist-plt:	dist-plt.o	dist-mc.o dist-util.o dist-quantities.o \
				formula.o mc-quantities.o mc-util.o dist-data.o\
				numin.o misc.o plt.o quantities.o log.o rand.o \
				data-trans.o digamma.o
		$(CC) $(LFLAGS) dist-plt.o dist-mc.o dist-util.o \
		  dist-quantities.o formula.o mc-quantities.o mc-util.o \
		  dist-data.o numin.o misc.o plt.o quantities.o log.o rand.o \
                  data-trans.o digamma.o \
                  -lm -o dist-plt

dist-tbl:	dist-plt.o	dist-mc.o dist-util.o dist-quantities.o \
				formula.o mc-quantities.o mc-util.o dist-data.o\
				numin.o misc.o tbl.o quantities.o log.o rand.o \
				data-trans.o digamma.o
		$(CC) $(LFLAGS) dist-plt.o dist-mc.o dist-util.o \
		  dist-quantities.o formula.o mc-quantities.o mc-util.o \
		  dist-data.o numin.o misc.o tbl.o quantities.o log.o rand.o \
                  data-trans.o digamma.o \
                  -lm -o dist-tbl

dist-hist:	dist-plt.o	dist-mc.o dist-util.o dist-quantities.o \
				formula.o mc-quantities.o mc-util.o dist-data.o\
				numin.o misc.o hist.o quantities.o log.o rand.o\
				data-trans.o digamma.o
		$(CC) $(LFLAGS) dist-plt.o dist-mc.o dist-util.o \
		  dist-quantities.o formula.o mc-quantities.o mc-util.o \
		  dist-data.o numin.o misc.o hist.o quantities.o log.o rand.o \
		  data-trans.o digamma.o \
                  -lm -o dist-hist

dist-plt.o:	dist-plt.c	misc.h log.h quantities.h mc.h dist.h

dist-initial:	dist-initial.o	dist-util.o dist-data.o misc.o log.o formula.o \
				numin.o data-trans.o rand.o digamma.o
		$(CC) $(LFLAGS) dist-initial.o dist-util.o dist-data.o \
		  misc.o log.o formula.o numin.o data-trans.o rand.o digamma.o \
		  -lm -o dist-initial

dist-initial.o:	dist-initial.c	misc.h log.h mc.h dist.h formula.h rand.h

dist-display:	dist-display.o	dist-util.o dist-data.o misc.o log.o formula.o \
				numin.o data-trans.o digamma.o rand.o
		$(CC) $(LFLAGS) dist-display.o dist-util.o dist-data.o \
		  misc.o log.o formula.o numin.o data-trans.o digamma.o rand.o \
		  -lm -o dist-display

dist-display.o:	dist-display.c	misc.h log.h mc.h dist.h formula.h

dist-est:	dist-est.o	dist-util.o dist-data.o misc.o log.o formula.o \
				numin.o data-trans.o digamma.o rand.o
		$(CC) $(LFLAGS) dist-est.o dist-util.o dist-data.o \
		  misc.o log.o formula.o numin.o data-trans.o digamma.o rand.o \
		  -lm -o dist-est

dist-est.o:	dist-est.c	misc.h log.h mc.h dist.h formula.h

dist-gen:	dist-gen.o	dist-util.o dist-data.o misc.o log.o formula.o \
				numin.o data-trans.o rand.o digamma.o
		$(CC) $(LFLAGS) dist-gen.o dist-util.o dist-data.o \
		  misc.o log.o formula.o numin.o data-trans.o rand.o digamma.o \
		  -lm -o dist-gen

dist-gen.o:	dist-gen.c	misc.h log.h mc.h dist.h formula.h rand.h


dist-util.o:	dist-util.c	misc.h log.h mc.h formula.h dist.h \
				data.h dist-data.h
dist-quantities.o: dist-quantities.c misc.h log.h mc.h quantities.h formula.h \
                                     data.h dist.h dist-data.h
dist-data.o:	dist-data.c	misc.h log.h data.h numin.h dist-data.h
