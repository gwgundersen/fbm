# MAKEFILE FOR GAUSSIAN PROCESS PROGRAMS. 

# Copyright (c) 1996 by Radford M. Neal 
#
# Permission is granted for anyone to copy, use, or modify this program 
# for purposes of research or education, provided this copyright notice 
# is retained, and note is made of any changes that have been made. 
#
# This program is distributed without any warranty, express or implied.
# As this program was written for research purposes only, it has not been
# tested to the degree that would be advisable in any important application.
# All use of this program is entirely at the user's own risk.

include ../make.include

programs:	gp-spec gp-gen gp-display gp-eval gp-plt gp-hist gp-mc gp-pred \
		gp-stepsizes gp-grad-test

clean:
	rm -f *.o gp-spec gp-gen gp-display gp-eval gp-pred gp-mc \
	          gp-plt gp-hist gp-grad-test gp-stepsizes gp-grad-test


include ../util/util.make
include ../mc/mc.make


gp-spec:	gp-spec.o	gp-util.o prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) gp-spec.o gp-util.o \
		  prior.o misc.o log.o rand.o -lm -o gp-spec

gp-spec.o:	gp-spec.c	misc.h prior.h model.h log.h data.h gp.h

gp-gen:		gp-gen.o	gp-util.o gp-prior.o prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) gp-gen.o gp-util.o gp-prior.o \
		  prior.o misc.o log.o rand.o -lm -o gp-gen

gp-gen.o:	gp-gen.c	misc.h prior.h model.h log.h data.h gp.h

gp-display:	gp-display.o	gp-util.o prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) gp-display.o gp-util.o \
		  prior.o misc.o log.o rand.o -lm -o gp-display

gp-display.o:	gp-display.c	misc.h prior.h model.h log.h data.h gp.h

gp-eval:	gp-eval.o	gp-util.o gp-data.o gp-cov.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o matrix.o
		$(CC) $(LFLAGS) gp-eval.o gp-util.o gp-data.o gp-cov.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  matrix.o -lm -o gp-eval

gp-eval.o:	gp-eval.c	misc.h prior.h model.h log.h data.h matrix.h \
				gp.h gp-data.h

gp-plt:		gp-plt.o	gp-mc.o prior.o gp-data.o gp-cov.o model.o \
				gp-quantities.o gp-util.o gp-prior.o \
				mc-quantities.o mc-util.o plt.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) gp-plt.o gp-mc.o gp-cov.o prior.o \
		  gp-data.o model.o gp-quantities.o gp-util.o gp-prior.o \
		  mc-quantities.o mc-util.o plt.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o matrix.o ars.o \
		  -lm -o gp-plt

gp-hist:	gp-plt.o	gp-mc.o prior.o gp-data.o gp-cov.o model.o \
				gp-quantities.o gp-util.o gp-prior.o \
				mc-quantities.o mc-util.o hist.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) gp-plt.o gp-mc.o gp-cov.o prior.o \
		  gp-data.o model.o gp-quantities.o gp-util.o gp-prior.o \
		  mc-quantities.o mc-util.o hist.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o matrix.o ars.o \
		  -lm -o gp-hist

gp-plt.o:	gp-plt.c	misc.h log.h quantities.h mc.h

gp-pred:	gp-pred.o	gp-util.o gp-data.o model.o misc.o log.o \
				matrix.o rand.o numin.o data-trans.o prior.o \
				gp-cov.o 
		$(CC) $(LFLAGS) gp-pred.o gp-util.o gp-data.o gp-cov.o \
		  model.o misc.o log.o rand.o numin.o data-trans.o prior.o \
		  matrix.o -lm -o gp-pred

gp-pred.o:	gp-pred.c	misc.h prior.h model.h gp.h gp-data.h log.h \
				matrix.h data.h numin.h rand.h mc.h

gp-stepsizes:	mc-stepsizes.o	gp-mc.o gp-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o matrix.o ars.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-cov.o \
				model.o gp-util.o gp-quantities.o \
				mc-quantities.o quantities.o
		$(CC) $(LFLAGS) mc-stepsizes.o gp-mc.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-data.o model.o gp-util.o \
		  gp-plt.o gp-quantities.o mc-quantities.o quantities.o \
		  matrix.o ars.o gp-cov.o -lm -o gp-stepsizes

gp-grad-test:	mc-grad-test.o	gp-mc.o gp-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o matrix.o ars.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-cov.o \
				model.o gp-util.o gp-quantities.o \
				mc-quantities.o quantities.o
		$(CC) $(LFLAGS) mc-grad-test.o gp-mc.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-data.o model.o gp-util.o \
		  gp-plt.o gp-quantities.o mc-quantities.o quantities.o \
		  matrix.o ars.o gp-cov.o -lm -o gp-grad-test

gp-mc:		mc.o gp-mc.o	misc.o log.o rand.o numin.o data-trans.o \
				gp-plt.o mc-iter.o mc-traj.o mc-util.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-cov.o \
				model.o gp-util.o gp-quantities.o \
				mc-quantities.o quantities.o mc-heatbath.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) mc.o gp-mc.o misc.o log.o rand.o numin.o\
		  data-trans.o mc-iter.o mc-traj.o mc-util.o matrix.o ars.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-plt.o gp-quantities.o mc-quantities.o \
		  quantities.o gp-data.o model.o gp-util.o gp-cov.o \
		  -lm -o gp-mc

gp-mc.o:	gp-mc.c		misc.h rand.h log.h mc.h data.h prior.h \
				matrix.h model.h gp.h gp-data.h


gp-util.o:	gp-util.c	misc.h log.h prior.h model.h data.h gp.h rand.h
gp-cov.o:	gp-cov.c	misc.h log.h prior.h model.h data.h gp.h
gp-data.o:	gp-data.c	misc.h log.h data.h prior.h model.h numin.h \
				gp.h gp-data.h
gp-prior.o:	gp-prior.c	misc.h log.h prior.h model.h data.h gp.h rand.h
gp-quantities.o: gp-quantities.c misc.h log.h prior.h model.h gp.h data.h \
				gp-data.h quantities.h
