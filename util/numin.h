/* NUMIN.H - Interface to module for reading numeric input. */

/* Copyright (c) 1995 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


/* STRUCTURE HOLDING THE STATE OF A NUMERIC INPUT SOURCE.  This structure
   holds the specification of where to get the input, and the current status
   of reading the input.  The 'filename', 'first', 'last', 'complement', and
   'index' fields may be of significance as defaults even when a new source 
   is opened by the 'numin_spec' procedure; the other fields are initialized
   at this point. */

#define Max_items 130	/* Max number of items than can be required */

typedef struct
{
  int N_items;		/* Number of values needed */
  int last_index;	/* Last index used */

  char filename[100];	/* Name of data file */
  int first, last;	/* Range of lines to use */
  int complement;	/* Read all _but_ the range specified? */
  int index[Max_items];	/* Indexes for items needed */

  FILE *file;		/* File open for reading */
  int length;		/* Number of lines in file */
  int line;		/* Current line */

} numin_source;


/* PROCEDURES. */

void numin_spec  (numin_source *, char *, int);	/* Specify source of input */
int  numin_start (numin_source *);		/* Start reading input */
void numin_read  (numin_source *, double *);	/* Read next line of input */
void numin_close (numin_source *);		/* Close file being read */
