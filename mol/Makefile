# MAKEFILE FOR MOLECULAR DYNAMICS MODULE.

# Copyright (c) 2005 by Radford M. Neal
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

programs:	mol-spec mol-mc mol-plt mol-hist mol-grad-test mol-stepsizes \
		mol-tbl mol-display mol-his mol-mis

clean:
	rm -f *.o mol-spec mol-mc mol-plt mol-hist mol-grad-test mol-stepsizes \
	          mol-tbl mol-display mol-his


include ../util/util.make
include ../mc/mc.make

mol-spec:	mol-spec.o	misc.o log.o
		$(CC) $(LFLAGS) mol-spec.o misc.o log.o -lm -o mol-spec

mol-spec.o:	mol-spec.c	misc.h log.h mc.h mol.h

mol-display:	mol-display.o	misc.o log.o
		$(CC) $(LFLAGS) mol-display.o misc.o log.o -lm -o mol-display

mol-display.o:	mol-display.c	misc.h log.h mc.h mol.h

mol-grad-test:	mc-grad-test.o mol-mc.o mol-util.o mol-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o mol-quantities.o
		$(CC) $(LFLAGS) mc-grad-test.o mol-mc.o mol-plt.o misc.o log.o \
                  mol-util.o rand.o quantities.o \
		  mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o \
		  mc-slice.o mc-traj.o mc-util.o mc-quantities.o \
		  mol-quantities.o -lm -o mol-grad-test

mol-stepsizes:	mc-stepsizes.o mol-mc.o mol-util.o mol-plt.o misc.o log.o \
				rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o mol-quantities.o
		$(CC) $(LFLAGS) mc-stepsizes.o mol-mc.o mol-plt.o misc.o log.o \
		  mol-util.o rand.o quantities.o \
		  mc-traj.o mc-iter.o mc-heatbath.o mc-metropolis.o mc-hybrid.o\
		  mc-slice.o mc-util.o mc-quantities.o mol-quantities.o \
		  -lm -o mol-stepsizes

mol-his:        mc-his.o mol-mc.o mol-util.o mol-plt.o misc.o \
				log.o rand.o quantities.o \
				mc-heatbath.o mc-traj.o mc-util.o \
				mc-quantities.o mol-quantities.o
		$(CC) $(LFLAGS) mc-his.o mol-mc.o mol-util.o mol-plt.o misc.o \
		  log.o rand.o quantities.o mc-traj.o mc-heatbath.o \
		  mc-util.o mc-quantities.o mol-quantities.o \
		  -lm -o mol-his

mol-mis:        mc-mis.o mol-mc.o mol-util.o mol-plt.o misc.o \
				log.o rand.o quantities.o \
				mc-util.o \
				mc-quantities.o mol-quantities.o
		$(CC) $(LFLAGS) mc-mis.o mol-mc.o mol-util.o mol-plt.o misc.o \
		  log.o rand.o quantities.o mc-util.o \
		  mc-quantities.o mol-quantities.o \
		  -lm -o mol-mis

mol-mc:		mc.o mol-mc.o	mol-plt.o mol-util.o misc.o \
				log.o rand.o quantities.o \
				mc-iter.o mc-heatbath.o mc-metropolis.o \
				mc-hybrid.o mc-slice.o mc-traj.o mc-util.o \
				mc-quantities.o mol-quantities.o
		$(CC) $(LFLAGS) mc.o mol-mc.o mol-util.o mol-plt.o misc.o \
		  log.o rand.o quantities.o mc-traj.o mc-iter.o mc-heatbath.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-util.o \
		  mc-quantities.o mol-quantities.o \
		  -lm -o mol-mc

mol-mc.o:	mol-mc.c	misc.h log.h mol.h mc.h

mol-plt:	mol-plt.o	mol-mc.o mol-util.o mc-quantities.o mc-util.o \
				misc.o plt.o quantities.o log.o rand.o \
				mol-quantities.o 
		$(CC) $(LFLAGS) mol-plt.o mol-mc.o mol-util.o mc-quantities.o \
		  mc-util.o \
		  misc.o plt.o quantities.o log.o rand.o mol-quantities.o \
		  -lm -o mol-plt

mol-tbl:	mol-plt.o	mol-mc.o mol-util.o mc-quantities.o mc-util.o \
				misc.o tbl.o quantities.o log.o rand.o \
				mol-quantities.o
		$(CC) $(LFLAGS) mol-plt.o mol-mc.o mol-util.o mc-quantities.o \
		  mc-util.o \
		  misc.o tbl.o quantities.o log.o rand.o mol-quantities.o \
		  -lm -o mol-tbl

mol-hist:	mol-plt.o	mol-mc.o mol-util.o mc-quantities.o mc-util.o \
				misc.o hist.o quantities.o log.o rand.o \
				mol-quantities.o
		$(CC) $(LFLAGS) mol-plt.o mol-mc.o mol-util.o mc-quantities.o \
		  mc-util.o \
		  misc.o hist.o quantities.o log.o rand.o mol-quantities.o \
		  -lm -o mol-hist

mol-plt.o:	mol-plt.c	misc.h log.h quantities.h mc.h mol.h
mol-util.o:	mol-util.c	misc.h rand.h log.h mc.h mol.h
mol-quantities.o: mol-quantities.c misc.h log.h quantities.h mc.h mol.h 
