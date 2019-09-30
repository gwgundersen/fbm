/* NUMIN.H - Interface to module for reading numeric input. */

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


/* STRUCTURE HOLDING THE STATE OF A NUMERIC INPUT SOURCE.  This structure
   holds the specification of where to get the input, and the current status
   of reading the input.  The 'filename', 'first', 'last', 'complement', and
   'index' fields may be of significance as defaults even when a new source 
   is opened by the 'numin_spec' procedure; the other fields are initialized
   at this point. */

#define Max_items 10001	/* Max number of items than can be required */

typedef struct
{
  int N_items;		/* Number of values needed */
  int last_index;	/* Last index used */

  char filename[100];	/* Name of data file */
  int first, last;	/* Range of lines to use */
  int complement;	/* Read all _but_ the range specified? */
  int index[Max_items];	/* Indexes for items needed */

  /* The following has the information from "index" above in different form.  */

  int iused[Max_items+1]; /* Positions of items used, in order, with 0 at end */
  int ifor[Max_items];	  /* What they're used for */

  FILE *file;		/* File open for reading */
  int length;		/* Number of lines in file */
  int line;		/* Current line */

} numin_source;


/* PROCEDURES. */

void numin_spec  (numin_source *, char *, int);	/* Specify source of input */
int  numin_start (numin_source *);		/* Start reading input */
void numin_read  (numin_source *, double *);	/* Read next line of input */
void numin_close (numin_source *);		/* Close file being read */
