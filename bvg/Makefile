# MAKEFILE FOR TESTS ON BIVARIATE GAUSSIAN.

# Copyright (c) 1995, 1996 by Radford M. Neal
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

programs:	bvg-spec bvg-mc bvg-plt bvg-hist bvg-grad-test bvg-stepsizes \
		bvg-initial

clean:
	rm -f *.o bvg-spec bvg-mc bvg-plt bvg-hist bvg-grad-test bvg-stepsizes \
	          bvg-initial


include ../util/util.make
include ../mc/mc.make

bvg-spec:	bvg-spec.o	misc.o log.o
		$(CC) bvg-spec.o misc.o log.o -lm -o bvg-spec

bvg-spec.o:	bvg-spec.c	misc.h log.h mc.h bvg.h

bvg-grad-test:	mc-grad-test.o bvg-mc.o bvg-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o
		$(CC) $(LFLAGS) mc-grad-test.o bvg-mc.o bvg-plt.o misc.o log.o \
                  rand.o quantities.o \
		  mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o \
		  mc-slice.o mc-traj.o mc-util.o mc-quantities.o \
		  -lm -o bvg-grad-test

bvg-stepsizes:	mc-stepsizes.o bvg-mc.o bvg-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o
		$(CC) $(LFLAGS) mc-stepsizes.o bvg-mc.o bvg-plt.o misc.o log.o \
		  rand.o quantities.o \
		  mc-traj.o mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o\
		  mc-slice.o mc-util.o mc-quantities.o \
		  -lm -o bvg-stepsizes

bvg-mc:		mc.o bvg-mc.o	bvg-plt.o misc.o log.o rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o
		$(CC) $(LFLAGS) mc.o bvg-mc.o bvg-plt.o misc.o log.o rand.o \
		  quantities.o mc-traj.o mc-iter.o mc-heatbath.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o \
		  mc-quantities.o \
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

bvg-initial:	bvg-initial.o	misc.o log.o 
		$(CC) $(LFLAGS) bvg-initial.o misc.o log.o -lm -o bvg-initial

bvg-initial.o:	bvg-initial.c	misc.h log.h mc.h bvg.h
