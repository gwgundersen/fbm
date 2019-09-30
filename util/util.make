# UTILITY MAKEFILE, FOR INCLUSION IN OTHER MAKEFILES.

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

log.h:
	ln -s ../util/log.h .
data.h:
	ln -s ../util/data.h .
misc.h: 
	ln -s ../util/misc.h .
numin.h:
	ln -s ../util/numin.h .
quantities.h:
	ln -s ../util/quantities.h .
rand.h:
	ln -s ../util/rand.h .

misc.c:
	ln -s ../util/misc.c .
log.c:
	ln -s ../util/log.c .
quantities.c:
	ln -s ../util/quantities.c .
plt.c:
	ln -s ../util/plt.c .
hist.c:
	ln -s ../util/hist.c .
rand.c:
	ln -s ../util/rand.c .
numin.c:
	ln -s ../util/numin.c .
data-trans.c:
	ln -s ../util/data-trans.c .


misc.o:		misc.c		misc.h
log.o:		log.c		log.h
quantities.o:	quantities.c	log.h quantities.h
plt.o:		plt.c		misc.h log.h quantities.h
hist.o:		hist.c		misc.h log.h quantities.h
rand.o:		rand.c		rand.h
numin.o:	numin.c		numin.h
data-trans.o:	data-trans.c	data.h
