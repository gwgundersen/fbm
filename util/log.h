/* LOG.H - Interface to routines for reading and writing log files. */

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


/* LOG FILE HEADER/TRAILER.  Each record in a log file begins and ends 
   with the structure below.  In between is the data in the record. 
   The trailer at the end has the 'index' field complemented. */

#define Log_header_magic  0xf41a  /* Identifies header record */
#define Log_trailer_magic 0x9e06  /* Identifies trailer record */

typedef struct
{ unsigned short magic;	/* Magic number identifying header/trailer */
  short type;		/* Type of record, a character */
  int index;		/* Index of this record */
  int size;		/* Amount of data in record (in bytes) */
  int reserved;		/* Reserved for future extensions */
} log_header;


/* STATE OF A LOG FILE.  Used in connection with the low-level procedures 
   for reading and writing records, and maybe gobbling with the higher-level
   procedurs too.  The at_end flag is set when the file pointer is past
   the last record.  The at_beginning flag is set when the file pointer is
   at the start of the first record; it is not set when the file is empty. */

typedef struct
{ char *file_name;	/* Path name for log file */
  FILE *file_struct;	/* File structure for log file */
  log_header header;	/* Header for next record in log file */
  int at_end;		/* Whether we've read past the last record */
  int at_beginning;	/* Whether we're sitting at the first record */
  int last_index_known;	/* Is the index of the last record known? */
  int last_index;	/* Index of last record in log file, if known */
} log_file;


/* SET OF RECORDS GOBBLED UP FROM LOG FILE.  Used by the higher-level
   procedures for reading log files. */

typedef struct
{ int last_index;	/* Index of last record gobbled up */
  int req_size[128];	/* Required size for records of each type, -1 for any */
  int actual_size[128];	/* Actual size of current records */
  int index[128];	/* Index of current record stored for each type */
  void *data[128];	/* Data for last record of each type, 0 if none */
} log_gobbled;


/* Pointer to gobble structure containing records to compare against when
   records are written with log_append.  Set to zero when a comparison
   comes out unequal. */

extern log_gobbled *log_append_compare;


/* PROCEDURES. */

void log_file_create (log_file *);
void log_file_open   (log_file *, int);
void log_file_close  (log_file *);

void log_file_first    (log_file *);
void log_file_last     (log_file *);
void log_file_forward  (log_file *);
void log_file_backward (log_file *);

void log_file_read   (log_file *, void *, int);
void log_file_append (log_file *, void *);

void log_gobble_init (log_gobbled *, int);
void log_gobble      (log_file *, log_gobbled *);
int  log_gobble_last (log_file *, log_gobbled *);
