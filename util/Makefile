# MAKEFILE FOR UTILITY PROGRAMS AND MODULES.

# Copyright (c) 1995 by Radford M. Neal 
#
# Permission is granted for anyone to copy, use, or modify this program 
# for purposes of research or education, provided this copyright notice 
# is retained, and note is made of any changes that have been made. 
#
# This program is distributed without any warranty, express or implied.
# As this program was written for research purposes only, it has not been
# tested to the degree that would be advisable in any important application.
# All use of this program is entirely at the user's own risk.

CFLAGS = -O

all:		programs test-progs modules

programs:	data-spec log-records log-copy log-last rand-seed \
		grid series extract

test-progs:	demo-plt demo-hist rand-test numin-test

modules:	misc.o log.o quantities.o plt.o hist.o rand.o numin.o \
		data-trans.o 

clean:
	rm -f *.o data-spec log-records log-copy log-last rand-seed \
	          grid series extract demo-plt demo-hist rand-test numin-test


include util.make


grid:		grid.o
		$(CC) grid.o -lm -o grid

grid.o:		grid.c


series:		series.o	misc.o
		$(CC) series.o misc.o -lm -o series

series.o:	series.c	misc.h


data-spec:	data-spec.o	log.o numin.o data-trans.o
		$(CC) data-spec.o log.o numin.o data-trans.o -lm -o data-spec

data-spec.o:	data-spec.c	data.h log.h


log-last:	log-last.o	log.o
		$(CC) log-last.o log.o -lm -o log-last

log-last.o:	log-last.c	log.h


log-records:	log-records.o	log.o
		$(CC) log-records.o log.o -lm -o log-records

log-records.o:	log-records.c	log.h


log-copy:	log-copy.o	log.o misc.o
		$(CC) log-copy.o log.o misc.o -lm -o log-copy

log-copy.o:	log-copy.c	misc.h log.h


rand-seed:	rand-seed.o	log.o rand.o
		$(CC) rand-seed.o log.o rand.o -lm -o rand-seed

rand-seed.o:	rand-seed.c	rand.h log.h


extract:	extract.o	rand.o
		$(CC) extract.o rand.o -lm -o extract

extract.o:	extract.c	rand.h


demo-plt:	demo-plt.o	demo-quantities.o plt.o quantities.o \
				misc.o log.o
		$(CC) demo-plt.o demo-quantities.o plt.o quantities.o \
		  misc.o log.o -lm -o demo-plt

demo-hist:	demo-plt.o	demo-quantities.o hist.o quantities.o \
				misc.o log.o
		$(CC) demo-plt.o demo-quantities.o hist.o quantities.o \
		  misc.o log.o -lm -o demo-hist

demo-plt.o:	demo-plt.c	log.h quantities.h 
demo-quantities.o: demo-quantities.c log.h quantities.h


rand-test:	rand-test.o	rand.o
		$(CC) rand-test.o rand.o -lm -o rand-test

rand-test.o:	rand-test.c	rand.h


numin-test:	numin-test.o	numin.o
		$(CC) numin-test.o numin.o -lm -o numin-test

numin-test.o:	numin-test.c	numin.h
