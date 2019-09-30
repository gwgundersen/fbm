# MAKEFILE FOR DIRICHLET DIFFUSION TREE PROGRAMS. 

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

programs:	dft-spec dft-gen dft-display dft-cases dft-mc dft-pred \
		dft-plt dft-tbl dft-hist dft-genp dft-dendrogram

clean:
	rm -f *.o dft-spec dft-gen dft-display dft-cases dft-mc dft-pred \
	          dft-plt dft-tbl dft-hist dft-genp dft-dendrogram


include ../util/util.make

Ponly=1		# Overrides setting of 0 in util.make

include ../mc/mc.make

dft-genp: 	# Just a stub, since standard mc methods aren't currently used
		echo >dft-genp; chmod 755 dft-genp

dft-spec:	dft-spec.o	dft-util.o prior.o misc.o log.o rand.o ars.o
		$(CC) $(LFLAGS) dft-spec.o dft-util.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o dft-spec

dft-spec.o:	dft-spec.c	misc.h prior.h model.h log.h data.h dft.h

dft-gen:	dft-gen.o	dft-util.o \
				prior.o misc.o log.o rand.o ars.o 
		$(CC) $(LFLAGS) dft-gen.o dft-util.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o dft-gen

dft-gen.o:	dft-gen.c	misc.h prior.h model.h log.h data.h dft.h

dft-display:	dft-display.o	dft-util.o dft-print.o \
				prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) dft-display.o dft-util.o dft-print.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o dft-display

dft-display.o:	dft-display.c	misc.h prior.h model.h log.h data.h dft.h

dft-dendrogram:	dft-dendrogram.o	dft-util.o dft-print.o \
				prior.o misc.o log.o rand.o
		$(CC) $(LFLAGS) dft-dendrogram.o dft-util.o dft-print.o \
		  prior.o misc.o log.o rand.o ars.o -lm -o dft-dendrogram

dft-dendrogram.o: dft-dendrogram.c misc.h prior.h model.h log.h data.h dft.h

dft-cases:	dft-cases.o	dft-util.o dft-data.o dft-print.o dft-div.o \
				dft-prob.o \
				prior.o misc.o numin.o model.o data-trans.o \
				log.o rand.o dft-sample.o ars.o uars.o phi.o
		$(CC) $(LFLAGS) dft-cases.o dft-util.o dft-data.o dft-print.o \
		  dft-div.o dft-prob.o prior.o misc.o numin.o model.o \
		  data-trans.o log.o rand.o dft-sample.o ars.o uars.o phi.o \
		  -lm -o dft-cases

dft-cases.o:	dft-cases.c	misc.h prior.h model.h log.h data.h rand.h \
				dft.h dft-data.h

dft-tbl:	dft-plt.o	dft-mc.o dft-util.o prior.o dft-data.o model.o \
				dft-quantities.o dft-sample.o ars.o \
				mc-quantities.o mc-util.o tbl.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				dft-prob.o dft-div.o uars.o phi.o
		$(CC) $(LFLAGS) dft-plt.o dft-util.o dft-mc.o prior.o \
		  dft-data.o model.o dft-quantities.o dft-sample.o ars.o \
		  mc-quantities.o mc-util.o tbl.o quantities.o uars.o phi.o \
		  misc.o log.o rand.o numin.o data-trans.o dft-prob.o dft-div.o\
		  -lm -o dft-tbl

dft-plt:	dft-plt.o	dft-mc.o dft-util.o prior.o dft-data.o model.o \
				dft-quantities.o dft-sample.o ars.o \
				mc-quantities.o mc-util.o plt.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				dft-prob.o dft-div.o uars.o phi.o
		$(CC) $(LFLAGS) dft-plt.o dft-util.o dft-mc.o prior.o \
		  dft-data.o model.o dft-quantities.o dft-sample.o ars.o \
		  mc-quantities.o mc-util.o plt.o quantities.o uars.o phi.o \
		  misc.o log.o rand.o numin.o data-trans.o dft-prob.o dft-div.o\
		  -lm -o dft-plt

dft-hist:	dft-plt.o	dft-mc.o dft-util.o prior.o dft-data.o model.o \
				dft-quantities.o dft-sample.o ars.o \
				mc-quantities.o mc-util.o hist.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o \
				dft-prob.o dft-div.o uars.o phi.o 
		$(CC) $(LFLAGS) dft-plt.o dft-mc.o dft-util.o prior.o \
		  dft-data.o model.o dft-quantities.o dft-sample.o ars.o \
		  mc-quantities.o mc-util.o hist.o quantities.o uars.o phi.o \
		  misc.o log.o rand.o numin.o data-trans.o dft-prob.o dft-div.o\
		  -lm -o dft-hist

dft-plt.o:	dft-plt.c	misc.h log.h quantities.h mc.h

dft-pred:	pred.o		dft-pred.o dft-data.o dft-util.o \
				dft-sample.o dft-prob.o dft-div.o \
				model.o misc.o log.o ars.o uars.o phi.o \
				rand.o numin.o data-trans.o prior.o 
		$(CC) $(LFLAGS) pred.o dft-pred.o dft-data.o dft-util.o \
		  dft-sample.o dft-prob.o dft-div.o \
		  model.o misc.o log.o rand.o numin.o data-trans.o prior.o \
		  ars.o uars.o phi.o -lm -o dft-pred

dft-pred.o:	dft-pred.c	misc.h prior.h model.h dft.h dft-data.h log.h \
				data.h numin.h rand.h mc.h pred.h phi.h

dft-mc:		mc.o dft-mc.o	misc.o log.o rand.o numin.o data-trans.o \
				dft-plt.o mc-iter.o mc-traj.o mc-util.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				prior.o dft-util.o dft-data.o dft-sample.o \
				model.o dft-quantities.o dft-prob.o dft-div.o \
				mc-quantities.o quantities.o mc-heatbath.o \
				ars.o uars.o phi.o
		$(CC) $(LFLAGS) mc.o dft-mc.o misc.o log.o rand.o numin.o\
		  data-trans.o mc-iter.o mc-traj.o mc-util.o dft-sample.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o ars.o \
		  prior.o dft-plt.o dft-quantities.o dft-prob.o dft-div.o \
		  mc-quantities.o quantities.o dft-util.o dft-data.o model.o \
		  uars.o phi.o -lm -o dft-mc

dft-mc.o:	dft-mc.c	misc.h rand.h log.h mc.h data.h prior.h \
				model.h dft.h dft-data.h


dft-data.o:	dft-data.c	misc.h log.h data.h prior.h model.h numin.h \
				dft.h dft-data.h
dft-util.o:	dft-util.c	misc.h log.h data.h prior.h model.h dft.h 
dft-prob.o:	dft-prob.c	misc.h log.h data.h prior.h model.h mc.h \
				rand.h dft.h dft-data.h
dft-sample.o:	dft-sample.c	misc.h log.h data.h prior.h model.h \
				rand.h ars.h uars.h phi.h dft.h dft-data.h
dft-div.o:	dft-div.c	misc.h log.h data.h prior.h model.h rand.h dft.h
dft-print.o:	dft-print.c	misc.h log.h data.h prior.h model.h dft.h 
dft-quantities.o: dft-quantities.c misc.h log.h prior.h model.h dft.h data.h \
				dft-data.h quantities.h
