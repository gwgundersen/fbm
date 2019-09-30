# MAKEFILE FOR MIXTURE MODEL PROGRAMS. 

# Copyright (c) 1997 by Radford M. Neal 
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

programs:	mix-spec mix-gen mix-display mix-cases mix-mc \
		mix-plt mix-hist 

clean:
	rm -f *.o mix-spec mix-gen mix-display mix-cases mix-mc \
	          mix-plt mix-hist 


include ../util/util.make
include ../mc/mc.make


mix-spec:	mix-spec.o	mix-util.o prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) mix-spec.o mix-util.o \
		  prior.o misc.o log.o rand.o -lm -o mix-spec

mix-spec.o:	mix-spec.c	misc.h prior.h model.h log.h data.h mix.h

mix-gen:	mix-gen.o	mix-util.o prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) mix-gen.o mix-util.o \
		  prior.o misc.o log.o rand.o -lm -o mix-gen

mix-gen.o:	mix-gen.c	misc.h prior.h model.h log.h data.h mix.h

mix-display:	mix-display.o	mix-util.o mix-print.o \
				prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) mix-display.o mix-util.o mix-print.o \
		  prior.o misc.o log.o rand.o -lm -o mix-display

mix-display.o:	mix-display.c	misc.h prior.h model.h log.h data.h mix.h

mix-cases:	mix-cases.o	mix-util.o mix-data.o mix-print.o prior.o \
				misc.o numin.o model.o data-trans.o log.o \
				rand.o 
		$(CC) $(LFLAGS) mix-cases.o mix-util.o mix-data.o mix-print.o \
		  prior.o misc.o numin.o model.o data-trans.o log.o rand.o \
		  -lm -o mix-cases

mix-cases.o:	mix-cases.c	misc.h prior.h model.h log.h data.h rand.h \
				mix.h mix-data.h

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


mix-data.o:	mix-data.c	misc.h log.h data.h prior.h model.h numin.h \
				mix.h mix-data.h
mix-util.o:	mix-util.c	misc.h log.h data.h prior.h model.h mix.h 
mix-print.o:	mix-print.c	misc.h log.h data.h prior.h model.h mix.h 
mix-quantities.o: mix-quantities.c misc.h log.h prior.h model.h mix.h data.h \
				mix-data.h quantities.h
