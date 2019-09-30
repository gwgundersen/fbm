# MAKEFILE FOR GAUSSIAN PROCESS PROGRAMS. 

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

programs:	gp-spec gp-gen gp-display gp-eval gp-plt gp-hist gp-mc \
		gp-pred gp-cov gp-eigen gp-stepsizes gp-grad-test gp-tbl \
		gp-genp

clean:
	rm -f *.o gp-spec gp-gen gp-display gp-eval gp-plt gp-hist gp-mc \
	          gp-pred gp-cov gp-eigen gp-stepsizes gp-grad-test gp-tbl \
		  gp-genp


include ../util/util.make
include ../mc/mc.make


gp-spec:	gp-spec.o	gp-util.o prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) gp-spec.o gp-util.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o gp-spec

gp-spec.o:	gp-spec.c	misc.h prior.h model.h log.h data.h gp.h

gp-gen:		gp-gen.o	gp-util.o gp-prior.o prior.o misc.o log.o \
				rand.o ars.o
		$(CC) $(LFLAGS) gp-gen.o gp-util.o gp-prior.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o gp-gen

gp-gen.o:	gp-gen.c	misc.h prior.h model.h log.h data.h gp.h

gp-display:	gp-display.o	gp-util.o prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) gp-display.o gp-util.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o gp-display

gp-display.o:	gp-display.c	misc.h prior.h model.h log.h data.h gp.h

gp-eval:	gp-eval.o	gp-util.o gp-data.o gp-covm.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o matrix.o ars.o
		$(CC) $(LFLAGS) gp-eval.o gp-util.o gp-data.o gp-covm.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  matrix.o ars.o -lm -o gp-eval

gp-eval.o:	gp-eval.c	misc.h prior.h model.h log.h data.h matrix.h \
				gp.h gp-data.h

gp-eigen:	gp-eigen.o	gp-util.o gp-data.o gp-covm.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o matrix.o ars.o
		$(CC) $(LFLAGS) gp-eigen.o gp-util.o gp-data.o gp-covm.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  matrix.o ars.o -lm -o gp-eigen

gp-eigen.o:	gp-eigen.c	misc.h prior.h model.h log.h data.h matrix.h \
				gp.h gp-data.h

gp-cov:		gp-cov.o	gp-util.o gp-data.o gp-covm.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o matrix.o ars.o
		$(CC) $(LFLAGS) gp-cov.o gp-util.o gp-data.o gp-covm.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  matrix.o ars.o -lm -o gp-cov

gp-cov.o:	gp-cov.c	misc.h prior.h model.h log.h data.h matrix.h \
				gp.h gp-data.h

gp-plt:		gp-plt.o	gp-mc.o prior.o gp-data.o gp-covm.o model.o \
				gp-quantities.o gp-util.o gp-prior.o gp-like.o \
				mc-quantities.o mc-util.o plt.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) gp-plt.o gp-mc.o gp-covm.o gp-like.o prior.o \
		  gp-data.o model.o gp-quantities.o gp-util.o gp-prior.o \
		  mc-quantities.o mc-util.o plt.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o matrix.o ars.o \
		  -lm -o gp-plt

gp-tbl:		gp-plt.o	gp-mc.o prior.o gp-data.o gp-covm.o model.o \
				gp-quantities.o gp-util.o gp-prior.o gp-like.o \
				mc-quantities.o mc-util.o tbl.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) gp-plt.o gp-mc.o gp-covm.o gp-like.o prior.o \
		  gp-data.o model.o gp-quantities.o gp-util.o gp-prior.o \
		  mc-quantities.o mc-util.o tbl.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o matrix.o ars.o \
		  -lm -o gp-tbl

gp-hist:	gp-plt.o	gp-mc.o prior.o gp-data.o gp-covm.o model.o \
				gp-quantities.o gp-util.o gp-prior.o gp-like.o \
				mc-quantities.o mc-util.o hist.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				matrix.o ars.o
		$(CC) $(LFLAGS) gp-plt.o gp-mc.o gp-covm.o gp-like.o prior.o \
		  gp-data.o model.o gp-quantities.o gp-util.o gp-prior.o \
		  mc-quantities.o mc-util.o hist.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o matrix.o ars.o \
		  -lm -o gp-hist

gp-plt.o:	gp-plt.c	misc.h log.h quantities.h mc.h

gp-pred:	pred.o		gp-pred.o gp-data.o gp-util.o gp-covm.o \
				model.o misc.o log.o ars.o \
				matrix.o rand.o numin.o data-trans.o prior.o 
		$(CC) $(LFLAGS) pred.o gp-pred.o gp-data.o gp-util.o gp-covm.o \
		  model.o misc.o log.o rand.o numin.o data-trans.o prior.o \
		  matrix.o ars.o -lm -o gp-pred

gp-pred.o:	gp-pred.c	misc.h prior.h model.h gp.h gp-data.h log.h \
				matrix.h data.h numin.h rand.h mc.h pred.h

gp-stepsizes:	mc-stepsizes.o	gp-mc.o gp-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o matrix.o ars.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-covm.o \
				model.o gp-util.o gp-quantities.o gp-like.o \
				mc-quantities.o quantities.o prior.o 
		$(CC) $(LFLAGS) mc-stepsizes.o gp-mc.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-data.o model.o gp-util.o gp-like.o \
		  gp-plt.o gp-quantities.o mc-quantities.o quantities.o \
		  matrix.o ars.o gp-covm.o -lm -o gp-stepsizes

gp-genp:	mc-genp.o	gp-mc.o gp-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o matrix.o ars.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-covm.o \
				model.o gp-util.o gp-quantities.o gp-like.o \
				mc-quantities.o quantities.o prior.o 
		$(CC) $(LFLAGS) mc-genp.o gp-mc.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-data.o model.o gp-util.o gp-like.o \
		  gp-plt.o gp-quantities.o mc-quantities.o quantities.o \
		  matrix.o ars.o gp-covm.o -lm -o gp-genp

gp-grad-test:	mc-grad-test.o	gp-mc.o gp-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o matrix.o ars.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-covm.o \
				model.o gp-util.o gp-quantities.o gp-like.o \
				mc-quantities.o quantities.o prior.o
		$(CC) $(LFLAGS) mc-grad-test.o gp-mc.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-data.o model.o gp-util.o gp-like.o \
		  gp-plt.o gp-quantities.o mc-quantities.o quantities.o \
		  matrix.o ars.o gp-covm.o -lm -o gp-grad-test

gp-mc:		mc.o gp-mc.o	misc.o log.o rand.o numin.o data-trans.o \
				gp-plt.o mc-iter.o mc-traj.o mc-util.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o gp-prior.o gp-data.o gp-covm.o \
				model.o gp-util.o gp-quantities.o gp-like.o \
				mc-quantities.o quantities.o mc-heatbath.o \
				matrix.o ars.o prior.o
		$(CC) $(LFLAGS) mc.o gp-mc.o misc.o log.o rand.o numin.o\
		  data-trans.o mc-iter.o mc-traj.o mc-util.o matrix.o ars.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o gp-prior.o gp-plt.o gp-quantities.o mc-quantities.o \
		  quantities.o gp-data.o model.o gp-util.o gp-covm.o gp-like.o \
		  -lm -o gp-mc

gp-mc.o:	gp-mc.c		misc.h rand.h log.h mc.h data.h prior.h ars.h \
				matrix.h model.h gp.h gp-data.h 


gp-util.o:	gp-util.c	misc.h log.h prior.h model.h data.h gp.h rand.h
gp-covm.o:	gp-covm.c	misc.h log.h prior.h model.h data.h gp.h
gp-data.o:	gp-data.c	misc.h log.h data.h prior.h model.h numin.h \
				gp.h gp-data.h
gp-prior.o:	gp-prior.c	misc.h log.h prior.h model.h data.h gp.h rand.h
gp-like.o:	gp-like.c	misc.h log.h prior.h model.h data.h gp.h rand.h
gp-quantities.o: gp-quantities.c misc.h log.h prior.h model.h gp.h data.h \
				gp-data.h quantities.h
