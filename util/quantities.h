/* QUANTITIES.H - Interface to application modules that evaluate quantities. */

/* Copyright (c) 1995-2004 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided 
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */


/* This module establishes an interface though which an application-specific
   module can evaluate various quantities derived from the data stored in
   log files.  These values can then be used by various application-independent
   modules, such as one that plots traces of such quantities. 

   This header file must be included after log.h. 
*/


/* DESCRIPTION OF QUANTITIES WHOSE VALUES ARE WANTED. This array lists the
   quantities desired, identified by letter.  The 'available' flag is set
   if a quantity-evaluation module knows about the quantity so identified;
   the value +1 indicating that the module knows about the quantity and thinks
   it is meaningful, the value -1 indicating that the module knows about 
   the quantity but thinks it is not meaningful in the present context, and
   the value zero meaning that nobody has heard of this quantity.
   
   The 'modifier', 'low', and 'high' entries give the modifier for the
   quantity and the range of its values that are desired.  If the top-level 
   module sets the 'high' entry to -1, the actual value is filled in by the 
   application-specific module. */

#define Max_quantities 50 /* Max. number of quantities that can be specified */

typedef struct
{
  short letter;		/* Letter associated with quantity, zero for empty */
  short available;	/* Whether somebody knows how to get value of quantity*/
  int modifier;		/* Non-negative integer modifying value, -1 if none */
  int low;		/* Low end of range (non-negative), or -1 if no range */
  int high;		/* High end of range, -1 for indefinite */

} quantities_described[Max_quantities];


/* STORAGE FOR VALUES OF QUANTITIES WANTED.  Contains pointers to sufficient
   storage to hold the values of all quantities that are described in some
   quantities_described structure.  The 'updated' fields record whether the
   data stored has been brought up to date.  A 'next' field is provided to 
   facilitate linking together structures for successive iterations.  Freeing 
   this structure will free both it and the arrays of doubles it points to
   (since they are allocated as a single block). */

typedef struct quantities_held
{ 
  double *value[Max_quantities]; /* Storage for values of all quantities */
  short updated[Max_quantities]; /* Whether value stored is current */

  struct quantities_held *next;  /* Pointer to next set of values */

  double align;		/* Makes sure following storage is properly aligned */

} quantities_held;


/* APPLICATION-SPECIFIC MODULES.  This arrays are set up by the application
   to contain pointers to routines for looking at application specific
   arguments, for specifying record sizes, for initializing whatever needs it, 
   for determining whether quantities are known, for finding the values of
   quantities, and for cleaning up.  Each array is terminated by a zero entry. 
   See demo-quantities.c for an illustration of what these procedures do. */

extern void (*quant_app_arguments[])   (char ***);
extern void (*quant_app_record_sizes[])(log_gobbled *);
extern void (*quant_app_initialize[])  (log_gobbled *);
extern void (*quant_app_available[])   (quantities_described, log_gobbled *);
extern void (*quant_app_evaluate[])    (quantities_described, quantities_held *,
                                        log_gobbled *);
extern void (*quant_app_cleanup[])     (void);


/* PROCEDURES. */

int quantities_requested (quantities_described, char *, int);

void quantities_available (quantities_described, log_gobbled *);

quantities_held *quantities_storage (quantities_described);

void quantities_evaluate (quantities_described, quantities_held *, 
                          log_gobbled *);
