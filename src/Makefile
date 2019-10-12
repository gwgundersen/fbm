# MAKEFILE FOR SOURCE LOCATION MODULE

# Copyright (c) 2007 by Radford M. Neal
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

programs:	src-spec det-spec flow-spec src-initial src-gen src-dgen \
		src-display src-plt src-tbl src-hist src-stepsizes \
		src-mc src-pred src-intensity

clean:
	rm -f *.o *.exe \
		  src-spec det-spec flow-spec src-initial src-gen src-dgen \
		  src-display src-plt src-tbl src-hist src-stepsizes \
		  src-mc src-pred src-intensity

include ../util/util.make
include ../mc/mc.make

src-util.o:	src-util.c	log.h rand.h src.h
src-data.o:	src-data.c	misc.h numin.h src.h data.h src-data.h
src-prior.o:	src-prior.c	log.h rand.h src.h
src-quantities.o: src-quantities.c misc.h log.h quantities.h mc.h src.h \
				numin.h data.h src-data.h
src-flow.o:	src-flow.c	log.h src.h

src-spec:	src-spec.o	log.o src-util.o misc.o
		$(CC) $(LFLAGS) src-spec.o log.o src-util.o misc.o \
		  -lm -o src-spec

src-spec.o:	src-spec.c	misc.h log.h src.h

det-spec:	det-spec.o	log.o src-util.o
		$(CC) $(LFLAGS) det-spec.o log.o src-util.o misc.o \
		  -lm -o det-spec

det-spec.o:	det-spec.c	misc.h log.h src.h

flow-spec:	flow-spec.o	log.o src-util.o misc.o
		$(CC) $(LFLAGS) flow-spec.o log.o src-util.o misc.o \
		  -lm -o flow-spec

flow-spec.o:	flow-spec.c	misc.h log.h src.h

src-initial:	src-initial.o	misc.o rand.o log.o src-util.o src-prior.o
		$(CC) $(LFLAGS) src-initial.o misc.o rand.o log.o src-util.o \
		  src-prior.o -lm -o src-initial

src-initial.o:	src-initial.c	misc.h log.h src.h rand.h

src-gen:	src-gen.o	misc.o rand.o log.o src-util.o src-prior.o
		$(CC) $(LFLAGS) src-gen.o misc.o rand.o log.o src-util.o \
		  src-prior.o -lm -o src-gen

src-gen.o:	src-gen.c	misc.h log.h src.h rand.h

src-dgen:	src-dgen.o	log.o rand.o src-util.o src-data.o src-flow.o \
				numin.o misc.o data-trans.o phi.o
		$(CC) $(LFLAGS) src-dgen.o log.o rand.o src-util.o src-data.o \
		  src-flow.o numin.o misc.o data-trans.o phi.o \
		  -lm -o src-dgen

src-dgen.o:	src-dgen.c	misc.h log.h src.h numin.h data.h src-data.h \
				rand.h

src-display:	src-display.o	log.o src-util.o
		$(CC) $(LFLAGS) src-display.o log.o src-util.o \
		  -lm -o src-display

src-display.o:	src-display.c	log.h src.h

src-grad-test:	mc-grad-test.o src-mc.o src-plt.o src-util.o src-data.o numin.o\
				src-prior.o rand.o quantities.o misc.o log.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o src-quantities.o data-trans.o \
				src-flow.o phi.o
		$(CC) $(LFLAGS) mc-grad-test.o \
		  src-mc.o src-plt.o src-util.o src-data.o numin.o src-prior.o \
		  src-quantities.o rand.o quantities.o misc.o log.o phi.o \
		  mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o \
		  mc-slice.o mc-traj.o mc-util.o mc-quantities.o data-trans.o \
		  src-flow.o -lm -o src-grad-test

src-stepsizes:	mc-stepsizes.o src-mc.o src-plt.o src-util.o src-data.o numin.o\
				src-prior.o rand.o quantities.o misc.o log.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o src-quantities.o data-trans.o \
				src-flow.o phi.o
		$(CC) $(LFLAGS) mc-stepsizes.o \
		  src-mc.o src-plt.o src-util.o src-data.o numin.o src-prior.o \
		  rand.o quantities.o misc.o log.o data-trans.o phi.o \
		  mc-traj.o mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o\
		  mc-slice.o mc-util.o mc-quantities.o src-quantities.o \
		  src-flow.o -lm -o src-stepsizes

src-mc:		mc.o src-mc.o	src-plt.o src-util.o src-data.o numin.o phi.o \
				data-trans.o misc.o log.o rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o src-quantities.o src-prior.o \
				src-flow.o
		$(CC) $(LFLAGS) mc.o src-mc.o src-plt.o src-util.o src-data.o \
		  numin.o data-trans.o misc.o log.o rand.o phi.o \
		  quantities.o mc-traj.o mc-iter.o mc-heatbath.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o \
		  mc-quantities.o src-quantities.o src-prior.o src-flow.o \
		  -lm -o src-mc

src-mc.o:	src-mc.c	misc.h rand.h log.h numin.h mc.h src.h data.h \
				phi.h src-data.h

src-plt:	src-plt.o	src-mc.o src-data.o mc-quantities.o mc-util.o \
				misc.o plt.o quantities.o log.o rand.o \
				data-trans.o numin.o src-prior.o phi.o \
				src-quantities.o src-data.o src-util.o \
				src-flow.o
		$(CC) $(LFLAGS) src-plt.o src-mc.o mc-quantities.o mc-util.o \
		  misc.o plt.o quantities.o log.o rand.o numin.o data-trans.o \
		  src-quantities.o src-data.o src-util.o src-prior.o phi.o \
		  src-flow.o -lm -o src-plt

src-tbl:	src-plt.o	src-mc.o src-data.o mc-quantities.o mc-util.o \
				misc.o tbl.o quantities.o log.o rand.o \
				data-trans.o numin.o src-prior.o phi.o \
				src-quantities.o src-data.o src-util.o \
				src-flow.o
		$(CC) $(LFLAGS) src-plt.o src-mc.o mc-quantities.o mc-util.o \
		  misc.o tbl.o quantities.o log.o rand.o numin.o data-trans.o \
		  src-quantities.o src-data.o src-util.o src-prior.o phi.o \
		  src-flow.o -lm -o src-tbl

src-hist:	src-plt.o	src-mc.o src-data.o mc-quantities.o mc-util.o \
				misc.o hist.o quantities.o log.o rand.o \
				data-trans.o numin.o src-prior.o phi.o \
				src-quantities.o src-data.o src-util.o \
				src-flow.o
		$(CC) $(LFLAGS) src-plt.o src-mc.o mc-quantities.o mc-util.o \
		  misc.o hist.o quantities.o log.o rand.o numin.o data-trans.o \
		  src-quantities.o src-data.o src-util.o src-prior.o phi.o \
		  src-flow.o -lm -o src-hist

src-plt.o:	src-plt.c	misc.h log.h quantities.h mc.h src.h

src-pred:	pred.o		src-pred.o src-data.o src-util.o src-flow.o \
				misc.o phi.o log.o rand.o numin.o data-trans.o
		$(CC) $(LFLAGS)	pred.o src-pred.o src-data.o src-util.o \
		  src-flow.o \
		  misc.o phi.o log.o rand.o numin.o data-trans.o \
		  -lm -o src-pred

src-pred.o:	src-pred.c	misc.h rand.h log.h mc.h src.h numin.h \
				data.h src-data.h pred.h

src-intensity:	src-intensity.o	src-data.o src-util.o src-flow.o \
				misc.o phi.o log.o rand.o numin.o data-trans.o
		$(CC) $(LFLAGS)	src-intensity.o src-data.o src-util.o \
		  src-flow.o misc.o phi.o log.o rand.o \
		  numin.o data-trans.o -lm -o src-intensity

src-intensity.o: src-intensity.c misc.h rand.h log.h mc.h src.h
