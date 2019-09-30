# UTILITY MAKEFILE, FOR INCLUSION IN OTHER MAKEFILES.

# Copyright (c) 1995, 1996, 1998 by Radford M. Neal 
#
# Permission is granted for anyone to copy, use, or modify this program 
# for purposes of research or education, provided this copyright notice 
# is retained, and note is made of any changes that have been made. 
#
# This program is distributed without any warranty, express or implied.
# As this program was written for research purposes only, it has not been
# tested to the degree that would be advisable in any important application.
# All use of this program is entirely at the user's own risk.

log.h:
	ln -s ../util/log.h .
data.h:
	ln -s ../util/data.h .
misc.h: 
	ln -s ../util/misc.h .
matrix.h:
	ln -s ../util/matrix.h .
numin.h:
	ln -s ../util/numin.h .
quantities.h:
	ln -s ../util/quantities.h .
rand.h:
	ln -s ../util/rand.h .
ars.h:
	ln -s ../util/ars.h .
prior.h:
	ln -s ../util/prior.h .
model.h:
	ln -s ../util/model.h .
formula.h:
	ln -s ../util/formula.h .
pred.h:
	ln -s ../util/pred.h .

misc.c:
	ln -s ../util/misc.c .
matrix.c:
	ln -s ../util/matrix.c .
log.c:
	ln -s ../util/log.c .
quantities.c:
	ln -s ../util/quantities.c .
plt.c:
	ln -s ../util/plt.c .
tbl.c:
	ln -s ../util/tbl.c .
hist.c:
	ln -s ../util/hist.c .
rand.c:
	ln -s ../util/rand.c .
ars.c:
	ln -s ../util/ars.c .
numin.c:
	ln -s ../util/numin.c .
data-trans.c:
	ln -s ../util/data-trans.c .
prior.c:
	ln -s ../util/prior.c .
model.c:
	ln -s ../util/model.c .
formula.c:
	ln -s ../util/formula.c .
pred.c:
	ln -s ../util/pred.c .
digamma.c:
	ln -s ../util/digamma.c .


misc.o:		misc.c		misc.h
log.o:		log.c		log.h
quantities.o:	quantities.c	log.h quantities.h
plt.o:		plt.c		misc.h log.h quantities.h
tbl.o:		tbl.c		misc.h log.h quantities.h
hist.o:		hist.c		misc.h log.h quantities.h
ars.o:		ars.c		rand.h ars.h
numin.o:	numin.c		numin.h
data-trans.o:	data-trans.c	data.h
prior.o:	prior.c		rand.h prior.h 
matrix.o:	matrix.c	matrix.h
formula.o:	formula.c	formula.h rand.h
pred.o:		pred.c		misc.h rand.h log.h prior.h model.h data.h \
				mc.h pred.h
digamma.o:	digamma.c

rand.o:		rand.c		rand.h
	$(CC) $(CFLAGS) -DRAND_FILE=\"`pwd`/../util/randfile\" -c rand.c
