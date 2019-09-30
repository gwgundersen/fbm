# MAKEFILE FOR TESTS ON BIVARIATE GAUSSIAN.

# Copyright (c) 1994 by Radford M. Neal

CFLAGS = -O
LFLAGS = 

programs:	bvg-spec bvg-mc bvg-plt bvg-hist bvg-grad-test bvg-stepsizes 

clean:
	rm -f *.o bvg-spec bvg-mc bvg-plt bvg-hist bvg-grad-test bvg-stepsizes 


include ../util/util.make
include ../mc/mc.make

bvg-spec:	bvg-spec.o	misc.o log.o
		$(CC) bvg-spec.o misc.o log.o -lm -o bvg-spec

bvg-spec.o:	bvg-spec.c	misc.h log.h mc.h bvg.h

bvg-grad-test:	mc-grad-test.o bvg-mc.o bvg-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-metropolis.o mc-hybrid.o \
				mc-slice.o mc-traj.o mc-util.o mc-quantities.o
		$(CC) $(LFLAGS) mc-grad-test.o bvg-mc.o bvg-plt.o misc.o log.o \
                  rand.o quantities.o \
		  mc-iter.o mc-metropolis.o mc-hybrid.o mc-slice.o \
		  mc-traj.o mc-util.o mc-quantities.o \
		  -lm -o bvg-grad-test

bvg-stepsizes:	mc-stepsizes.o bvg-mc.o bvg-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-metropolis.o mc-hybrid.o \
				mc-slice.o mc-traj.o mc-util.o mc-quantities.o
		$(CC) $(LFLAGS) mc-stepsizes.o bvg-mc.o bvg-plt.o misc.o log.o \
		  rand.o quantities.o \
		  mc-traj.o mc-iter.o mc-metropolis.o mc-hybrid.o \
		  mc-slice.o mc-util.o mc-quantities.o \
		  -lm -o bvg-stepsizes

bvg-mc:		mc.o bvg-mc.o	bvg-plt.o misc.o log.o rand.o quantities.o \
				mc-iter.o mc-metropolis.o mc-hybrid.o \
				mc-slice.o mc-traj.o mc-util.o mc-quantities.o
		$(CC) $(LFLAGS) mc.o bvg-mc.o bvg-plt.o misc.o log.o rand.o \
		  quantities.o mc-traj.o mc-iter.o mc-metropolis.o mc-hybrid.o \
		  mc-slice.o mc-util.o mc-quantities.o \
		  -lm -o bvg-mc

bvg-mc.o:	bvg-mc.c	misc.h log.h bvg.h mc.h

bvg-plt:	bvg-plt.o	bvg-mc.o mc-quantities.o mc-util.o \
				misc.o plt.o quantities.o log.o rand.o
		$(CC) $(LFLAGS) bvg-plt.o bvg-mc.o mc-quantities.o mc-util.o \
		  misc.o plt.o quantities.o log.o rand.o -lm -o bvg-plt

bvg-hist:	bvg-plt.o	bvg-mc.o mc-quantities.o mc-util.o \
				misc.o hist.o quantities.o log.o rand.o
		$(CC) $(LFLAGS) bvg-plt.o bvg-mc.o mc-quantities.o mc-util.o \
		  misc.o hist.o quantities.o log.o rand.o -lm -o bvg-hist

bvg-plt.o:	bvg-plt.c	misc.h log.h quantities.h mc.h bvg.h
