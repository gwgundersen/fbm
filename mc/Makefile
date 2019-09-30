# MAKEFILE FOR MARKOV CHAIN MONTE CARLO MODULES.

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

programs:	mc-spec mc-temp-sched mc-temp-filter mc-ais

modules:	mc-traj.o mc-iter.o mc-metropolis.o mc-hybrid.o mc-slice.o \
		mc-util.o mc-quantities.o \
		mc-grad-test.o mc-stepsizes.o mc.o

clean:
	rm -f *.o mc-spec mc-temp-sched mc-temp-filter


include ../util/util.make
include mc.make

mc-spec:	mc-spec.o	misc.o log.o
		$(CC) $(LFLAGS) mc-spec.o misc.o log.o -lm -o mc-spec 

mc-spec.o:	mc-spec.c	misc.h log.h mc.h

mc-temp-sched:	mc-temp-sched.o	misc.o log.o
		$(CC) $(LFLAGS) mc-temp-sched.o misc.o log.o \
		      -lm -o mc-temp-sched 

mc-temp-sched.o: mc-temp-sched.c misc.h log.h mc.h

mc-temp-filter:	mc-temp-filter.o misc.o log.o
		$(CC) $(LFLAGS) mc-temp-filter.o misc.o log.o \
		      -lm -o mc-temp-filter 

mc-temp-filter.o: mc-temp-filter.c misc.h log.h mc.h

mc-ais:		mc-ais.o misc.o log.o
		$(CC) $(LFLAGS) mc-ais.o misc.o log.o \
		      -lm -o mc-ais 

mc-ais.o: 	mc-ais.c misc.h log.h mc.h
