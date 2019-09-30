# MAKEFILE FOR MARKOV CHAIN MONTE CARLO, FOR INCLUSION IN OTHER MAKEFILES

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

mc.h:
	ln -s ../mc/mc.h .
mc-traj.c:
	ln -s ../mc/mc-traj.c .
mc-iter.c:
	ln -s ../mc/mc-iter.c .
mc-heatbath.c:
	ln -s ../mc/mc-heatbath.c .
mc-metropolis.c:
	ln -s ../mc/mc-metropolis.c .
mc-hybrid.c:
	ln -s ../mc/mc-hybrid.c .
mc-slice.c:
	ln -s ../mc/mc-slice.c .
mc-util.c:
	ln -s ../mc/mc-util.c .
mc-quantities.c:
	ln -s ../mc/mc-quantities.c .
mc-grad-test.c:
	ln -s ../mc/mc-grad-test.c .
mc-stepsizes.c:
	ln -s ../mc/mc-stepsizes.c .
mc-genp.c:
	ln -s ../mc/mc-genp.c .
mc.c:
	ln -s ../mc/mc.c .
mc-his.c:
	ln -s ../mc/mc-his.c .

mc-traj.o:	mc-traj.c	misc.h rand.h log.h mc.h
mc-iter.o:	mc-iter.c	misc.h rand.h log.h mc.h quantities.h
mc-heatbath.o:	mc-heatbath.c	misc.h rand.h log.h mc.h
mc-metropolis.o:mc-metropolis.c	misc.h rand.h log.h mc.h 
mc-hybrid.o:	mc-hybrid.c	misc.h rand.h log.h mc.h
mc-slice.o:	mc-slice.c	misc.h rand.h log.h mc.h
mc-util.o:	mc-util.c	misc.h rand.h log.h mc.h
mc-quantities.o:mc-quantities.c misc.h log.h quantities.h mc.h
mc-grad-test.o:	mc-grad-test.c	misc.h rand.h log.h mc.h
mc-stepsizes.o:	mc-stepsizes.c	misc.h rand.h log.h mc.h
mc.o:		mc.c		misc.h rand.h log.h mc.h quantities.h
mc-his.o:	mc-his.c	misc.h rand.h log.h mc.h quantities.h
