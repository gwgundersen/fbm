# MAKEFILE FOR NEURAL NET PROGRAMS. 

# Copyright (c) 1995, 1996, 1997 by Radford M. Neal 
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

programs:	net-spec net-gen net-rej net-eval net-display net-pred net-mc \
		net-gd net-plt net-hist net-dvar net-grad-test net-stepsizes 

clean:
	rm -f *.o net-spec net-gen net-rej net-eval net-display net-pred net-mc\
	          net-gd net-plt net-hist net-dvar net-grad-test net-stepsizes 


include ../util/util.make
include ../mc/mc.make


net-func.o:	net-func.c	misc.h prior.h model.h net.h 
net-back.o:	net-back.c	misc.h prior.h model.h net.h 
net-grad.o:	net-grad.c	misc.h prior.h model.h net.h 
net-util.o:	net-util.c	misc.h prior.h model.h net.h data.h
net-setup.o:	net-setup.c	misc.h prior.h model.h net.h 
net-prior.o:	net-prior.c	misc.h prior.h model.h net.h rand.h
net-model.o:	net-model.c	misc.h prior.h model.h net.h rand.h
net-data.o:	net-data.c	misc.h prior.h model.h net.h data.h net-data.h \
				numin.h
net-print.o:	net-print.c	misc.h prior.h model.h net.h
net-quantities.o: net-quantities.c misc.h prior.h model.h net.h data.h \
				net-data.h quantities.h log.h


net-spec:	net-spec.o	net-util.o prior.o net-prior.o \
				misc.o log.o rand.o
		$(CC) $(LFLAGS) net-spec.o net-util.o prior.o net-prior.o \
		  misc.o log.o rand.o -lm -o net-spec

net-spec.o:	net-spec.c	misc.h prior.h model.h net.h log.h

net-gen:	net-gen.o	net-util.o prior.o net-prior.o net-setup.o \
				misc.o log.o rand.o
		$(CC) $(LFLAGS) net-gen.o net-util.o prior.o net-prior.o \
		  net-setup.o misc.o log.o rand.o -lm -o net-gen

net-gen.o:	net-gen.c	misc.h prior.h model.h net.h log.h rand.h

net-rej:	net-rej.o	net-util.o prior.o net-prior.o net-setup.o \
				net-data.o model.o net-func.o net-model.o \
				misc.o log.o numin.o data-trans.o rand.o
		$(CC) $(LFLAGS) net-rej.o net-util.o prior.o net-prior.o \
		  net-setup.o net-data.o model.o net-func.o net-model.o misc.o \
		  log.o numin.o data-trans.o rand.o -lm -o net-rej

net-rej.o:	net-rej.c	misc.h prior.h model.h net.h log.h rand.h

net-eval:	net-eval.o	net-util.o prior.o net-prior.o net-setup.o \
				net-func.o net-model.o model.o misc.o \
				log.o rand.o 
		$(CC) $(LFLAGS) net-eval.o net-util.o prior.o net-prior.o \
		  net-setup.o net-func.o net-model.o model.o misc.o \
		  log.o rand.o -lm -o net-eval

net-eval.o:	net-eval.c	misc.h prior.h model.h net.h log.h 

net-values:	net-values.o	net-util.o prior.o net-prior.o net-setup.o \
				net-func.o net-model.o misc.o log.o rand.o 
		$(CC) $(LFLAGS) net-values.o net-util.o prior.o net-prior.o \
		  net-setup.o net-func.o net-model.o \
		  misc.o log.o rand.o -lm -o net-values

net-values.o:	net-values.c	misc.h prior.h model.h net.h log.h 

net-display:	net-display.o	net-util.o prior.o net-prior.o net-setup.o \
				net-print.o misc.o log.o rand.o
		$(CC) $(LFLAGS) net-display.o net-util.o prior.o net-prior.o \
		  net-setup.o net-print.o misc.o log.o rand.o -lm -o net-display

net-display.o:	net-display.c	misc.h prior.h model.h net.h log.h

net-pred:	net-pred.o	net-setup.o net-func.o net-model.o net-data.o \
				model.o net-util.o misc.o log.o rand.o numin.o \
				data-trans.o 
		$(CC) $(LFLAGS) net-pred.o net-setup.o net-func.o net-model.o \
		  net-data.o model.o net-util.o misc.o log.o rand.o numin.o \
		  data-trans.o -lm -o net-pred

net-pred.o:	net-pred.c	misc.h prior.h model.h net.h net-data.h log.h \
				data.h numin.h rand.h mc.h

net-dvar:	net-dvar.o	net-util.o net-setup.o net-func.o net-model.o \
				misc.o log.o rand.o 
		$(CC) $(LFLAGS) net-dvar.o net-setup.o net-func.o net-model.o \
		  net-util.o misc.o log.o rand.o -lm -o net-dvar

net-dvar.o:	net-dvar.c	misc.h prior.h model.h net.h log.h 

net-grad-test:	mc-grad-test.o	net-mc.o ars.o net-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-metropolis.o mc-hybrid.o \
				mc-slice.o net-setup.o prior.o net-prior.o \
				net-model.o net-func.o net-back.o net-grad.o \
				net-data.o model.o net-util.o net-quantities.o \
				mc-quantities.o quantities.o mc-heatbath.o
		$(CC) $(LFLAGS) mc-grad-test.o net-mc.o ars.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  net-setup.o prior.o net-prior.o net-model.o net-func.o \
		  net-back.o net-grad.o net-data.o model.o net-util.o \
		  net-plt.o net-quantities.o mc-quantities.o quantities.o \
		  -lm -o net-grad-test

net-stepsizes:	mc-stepsizes.o	net-mc.o ars.o net-plt.o misc.o log.o rand.o \
				numin.o data-trans.o mc-iter.o mc-traj.o \
				mc-util.o mc-heatbath.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				net-setup.o prior.o net-prior.o net-model.o \
				net-func.o net-back.o net-grad.o net-data.o \
				model.o net-util.o net-quantities.o \
				mc-quantities.o quantities.o
		$(CC) $(LFLAGS) mc-stepsizes.o net-mc.o ars.o misc.o log.o \
		  rand.o numin.o data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  net-setup.o prior.o net-prior.o net-model.o net-func.o \
		  net-back.o net-grad.o net-data.o model.o net-util.o \
		  net-plt.o net-quantities.o mc-quantities.o quantities.o \
		  -lm -o net-stepsizes

net-mc:		mc.o net-mc.o	ars.o misc.o log.o rand.o numin.o data-trans.o \
				net-plt.o mc-iter.o mc-traj.o mc-util.o \
				mc-metropolis.o mc-hybrid.o mc-slice.o \
				net-setup.o prior.o net-prior.o net-model.o \
				net-func.o net-back.o net-grad.o net-data.o \
				model.o net-util.o net-quantities.o \
				mc-quantities.o quantities.o mc-heatbath.o
		$(CC) $(LFLAGS) mc.o net-mc.o ars.o misc.o log.o rand.o numin.o\
		  data-trans.o mc-iter.o mc-traj.o mc-util.o \
		  mc-metropolis.o mc-hybrid.o mc-slice.o mc-heatbath.o \
		  net-setup.o prior.o net-prior.o net-model.o net-func.o \
		  net-plt.o net-quantities.o mc-quantities.o quantities.o \
		  net-back.o net-grad.o net-data.o model.o net-util.o \
		  -lm -o net-mc

net-mc.o:	net-mc.c	misc.h rand.h ars.h log.h mc.h data.h prior.h \
				model.h net.h net-data.h

net-gd:		net-gd.o	net-setup.o net-func.o net-model.o net-data.o \
				net-back.o net-grad.o net-prior.o net-util.o \
				model.o misc.o log.o rand.o numin.o \
				data-trans.o prior.o 
		$(CC) $(LFLAGS) net-gd.o net-setup.o net-func.o net-model.o \
		  net-back.o net-grad.o net-data.o net-prior.o net-util.o \
		  model.o misc.o log.o rand.o numin.o data-trans.o prior.o \
		  -lm -o net-gd

net-gd.o:	net-gd.c	misc.h rand.h log.h mc.h data.h prior.h model.h\
				net.h net-data.h

net-plt:	net-plt.o	net-mc.o ars.o net-setup.o prior.o net-prior.o \
				net-func.o net-data.o model.o net-quantities.o \
				net-back.o net-grad.o net-util.o net-model.o \
				mc-quantities.o mc-util.o plt.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o
		$(CC) $(LFLAGS) net-plt.o net-mc.o ars.o net-setup.o prior.o \
		  net-model.o net-func.o net-data.o model.o net-quantities.o \
		  net-back.o net-grad.o net-util.o net-prior.o \
		  mc-quantities.o mc-util.o plt.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o -lm -o net-plt

net-hist:	net-plt.o	net-mc.o ars.o net-setup.o prior.o net-prior.o \
				net-func.o net-data.o model.o net-quantities.o \
				net-back.o net-grad.o net-util.o net-model.o \
				mc-quantities.o mc-util.o hist.o quantities.o \
				misc.o log.o rand.o numin.o data-trans.o
		$(CC) $(LFLAGS) net-plt.o net-mc.o ars.o net-setup.o prior.o \
		  net-prior.o net-model.o net-func.o net-data.o model.o \
		  net-quantities.o net-back.o net-grad.o net-util.o \
		  mc-quantities.o mc-util.o hist.o quantities.o \
		  misc.o log.o rand.o numin.o data-trans.o -lm -o net-hist

net-plt.o:	net-plt.c	misc.h log.h quantities.h mc.h
