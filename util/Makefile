# MAKEFILE FOR UTILITY PROGRAMS AND MODULES.

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

all:		programs modules

programs:	data-spec model-spec log-records log-copy log-append log-last \
		log-equal rand-seed grid series extract combine find-min calc \
		case-indicators

modules:	misc.o matrix.o log.o quantities.o plt.o hist.o rand.o numin.o \
		data-trans.o model.o prior.o ars.o digamma.o 

clean:
	rm -f *.o data-spec model-spec log-records log-copy log-append log-last\
	          log-equal rand-seed grid series extract combine find-min calc\
		  case-indicators


include util.make


calc:		calc.o formula.o digamma.o rand.o
		$(CC) $(LFLAGS) calc.o formula.o digamma.o rand.o -lm -o calc

calc.o:		calc.c		formula.h


grid:		grid.o
		$(CC) $(LFLAGS) grid.o -lm -o grid

grid.o:		grid.c


case-indicators:	case-indicators.o
		$(CC) $(LFLAGS) case-indicators.o -lm -o case-indicators

case-indicators.o:	case-indicators.c


series:		series.o	misc.o
		$(CC) $(LFLAGS) series.o misc.o -lm -o series

series.o:	series.c	misc.h


data-spec:	data-spec.o	log.o numin.o data-trans.o misc.o
		$(CC) $(LFLAGS) data-spec.o log.o numin.o data-trans.o misc.o \
		-lm -o data-spec

data-spec.o:	data-spec.c	data.h log.h misc.h numin.h


model-spec:	model-spec.o	log.o rand.o prior.o matrix.o ars.o
		$(CC) $(LFLAGS) model-spec.o log.o rand.o prior.o matrix.o \
		ars.o -lm -o model-spec

model-spec.o:	model-spec.c	log.h prior.h model.h matrix.h


log-last:	log-last.o	log.o
		$(CC) $(LFLAGS) log-last.o log.o -lm -o log-last

log-last.o:	log-last.c	log.h


log-records:	log-records.o	log.o
		$(CC) $(LFLAGS) log-records.o log.o -lm -o log-records

log-records.o:	log-records.c	log.h


log-copy:	log-copy.o	log.o misc.o
		$(CC) $(LFLAGS) log-copy.o log.o misc.o -lm -o log-copy

log-copy.o:	log-copy.c	misc.h log.h


log-append:	log-append.o	log.o misc.o
		$(CC) $(LFLAGS) log-append.o log.o misc.o -lm -o log-append

log-append.o:	log-append.c	misc.h log.h


log-equal:	log-equal.o	log.o
		$(CC) $(LFLAGS) log-equal.o log.o -lm -o log-equal

log-equal.o:	log-equal.c	log.h


rand-seed:	rand-seed.o	log.o rand.o
		$(CC) $(LFLAGS) rand-seed.o log.o rand.o -lm -o rand-seed

rand-seed.o:	rand-seed.c	rand.h log.h


extract:	extract.o	rand.o
		$(CC) $(LFLAGS) extract.o rand.o -lm -o extract

extract.o:	extract.c	rand.h

combine:	combine.o
		$(CC) $(LFLAGS) combine.o -lm -o combine

combine.o:	combine.c

find-min:	find-min.o
		$(CC) $(LFLAGS) find-min.o -lm -o find-min

find-min.o:	find-min.c

formula.o:	formula.c	formula.h extfunc.h rand.h
