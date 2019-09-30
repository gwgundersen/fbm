# MAKEFILE FOR MIXTURE MODEL PROGRAMS. 

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

programs:	mix-spec mix-gen mix-display mix-cases mix-mc mix-pred \
		mix-plt mix-tbl mix-hist mix-genp 

clean:
	rm -f *.o mix-spec mix-gen mix-display mix-cases mix-mc mix-pred \
	          mix-plt mix-tbl mix-hist mix-genp


include ../util/util.make

Ponly=1		# Overrides setting of 0 in util.make

include ../mc/mc.make

mix-genp: 	# Just a stub, since standard mc methods aren't currently used
		echo >mix-genp; chmod 755 mix-genp

mix-spec:	mix-spec.o	mix-util.o prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) mix-spec.o mix-util.o ars.o \
		  prior.o misc.o log.o rand.o -lm -o mix-spec

mix-spec.o:	mix-spec.c	misc.h prior.h model.h log.h data.h mix.h

mix-gen:	mix-gen.o	mix-util.o prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) mix-gen.o mix-util.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o mix-gen

mix-gen.o:	mix-gen.c	misc.h prior.h model.h log.h data.h mix.h

mix-display:	mix-display.o	mix-util.o mix-print.o \
				prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) mix-display.o mix-util.o mix-print.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o mix-display

mix-display.o:	mix-display.c	misc.h prior.h model.h log.h data.h mix.h

mix-cases:	mix-cases.o	mix-util.o mix-data.o mix-print.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o ars.o
		$(CC) $(LFLAGS) mix-cases.o mix-util.o mix-data.o mix-print.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  ars.o -lm -o mix-cases

mix-cases.o:	mix-cases.c	misc.h prior.h model.h log.h data.h rand.h \
				mix.h mix-data.h

mix-tbl:	mix-plt.o	mix-mc.o mix-util.o prior.o mix-data.o model.o \
				mix-quantities.o \
				mc-quantities.o mc-util.o tbl.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				ars.o
		$(CC) $(LFLAGS) mix-plt.o mix-util.o mix-mc.o prior.o \
		  mix-data.o model.o mix-quantities.o \
		  mc-quantities.o mc-util.o tbl.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o ars.o \
		  -lm -o mix-tbl

mix-plt:	mix-plt.o	mix-mc.o mix-util.o prior.o mix-data.o model.o \
				mix-quantities.o \
				mc-quantities.o mc-util.o plt.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				ars.o
		$(CC) $(LFLAGS) mix-plt.o mix-util.o mix-mc.o prior.o \
		  mix-data.o model.o mix-quantities.o \
		  mc-quantities.o mc-util.o plt.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o ars.o \
		  -lm -o mix-plt

mix-hist:	mix-plt.o	mix-mc.o mix-util.o prior.o mix-data.o model.o \
				mix-quantities.o \
				mc-quantities.o mc-util.o hist.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				ars.o
		$(CC) $(LFLAGS) mix-plt.o mix-mc.o mix-util.o prior.o \
		  mix-data.o model.o mix-quantities.o \
		  mc-quantities.o mc-util.o hist.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o ars.o \
		  -lm -o mix-hist

mix-plt.o:	mix-plt.c	misc.h log.h quantities.h mc.h

mix-mc:		mc.o mix-mc.o	misc.o log.o rand.o numin.o data-trans.o \
				mix-plt.o mc-iter.o mc-traj.o mc-util.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o mix-util.o mix-data.o \
				model.o mix-quantities.o \
				mc-quantities.o quantities.o mc-heatbath.o \
				ars.o
		$(CC) $(LFLAGS) mc.o mix-mc.o misc.o log.o rand.o numin.o\
		  data-trans.o mc-iter.o mc-traj.o mc-util.o ars.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  prior.o mix-plt.o mix-quantities.o mc-quantities.o \
		  quantities.o mix-util.o mix-data.o model.o \
		  -lm -o mix-mc

mix-mc.o:	mix-mc.c	misc.h rand.h log.h mc.h data.h prior.h ars.h \
				model.h mix.h mix-data.h

mix-pred:	pred.o mix-pred.o misc.o log.o rand.o numin.o data-trans.o \
				mix-plt.o mc-util.o \
				prior.o mix-util.o mix-data.o \
				model.o ars.o
		$(CC) $(LFLAGS) pred.o mix-pred.o misc.o log.o rand.o numin.o \
		  data-trans.o prior.o mix-util.o mix-data.o model.o ars.o \
		  -lm -o mix-pred

mix-pred.o:	mix-pred.c	misc.h rand.h log.h model.h data.h prior.h \
				numin.h mix.h mix-data.h mc.h pred.h

mix-data.o:	mix-data.c	misc.h log.h data.h prior.h model.h numin.h \
				mix.h mix-data.h
mix-util.o:	mix-util.c	misc.h log.h data.h prior.h model.h mix.h 
mix-print.o:	mix-print.c	misc.h log.h data.h prior.h model.h mix.h 
mix-quantities.o: mix-quantities.c misc.h log.h prior.h model.h mix.h data.h \
				mix-data.h quantities.h
