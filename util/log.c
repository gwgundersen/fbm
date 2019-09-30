/* LOG.C - Routines for handling log files. */

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "log.h"


/* This module contains routines for reading and writing log files.  All
   errors encountered are handled by displaying a message and terminating 
   the program. 

   Records from log files can be "gobbled up" using the somewhat higher-level
   log_gobble routines.
*/


static void read_header (log_file *, int);
static void check_header_trailer (log_header *, log_header *);


/* CREATE A NEW LOG FILE.  The name of the new file is taken from the
   log file state structure. */

void log_file_create
( log_file *logf	/* Log file state structure */
)
{ 
  logf->file_struct = fopen(logf->file_name,"w+b");

  if (logf->file_struct==NULL)
  { fprintf(stderr,"Can't create log file: %s\n",logf->file_name);
    exit(1);
  }

  logf->at_end = 1;
  logf->at_beginning = 0;
  logf->last_index_known = 0;
}


/* OPEN AN EXISTING LOG FILE.  The name of the file is taken from the
   log file state structure.  Also reads the header for the first 
   record.  If there is no first record, at_end is set. */

void log_file_open
( log_file *logf,	/* Log file state structure */
  int allow_append	/* Allow data to be appended? */
)
{
  logf->file_struct = fopen (logf->file_name, allow_append ? "r+b" : "rb");

  if (logf->file_struct==NULL)
  { fprintf(stderr,"Can't open log file: %s\n",logf->file_name);
    exit(1);
  }

  logf->at_end = 0;
  logf->last_index_known = 0;

  read_header(logf,0);
  logf->at_beginning = !logf->at_end;
}


/* CLOSE LOG FILE.  Should be called to check for any lingering errors. */

void log_file_close
( log_file *logf	/* Log file state structure */
)
{
  if (fclose(logf->file_struct)!=0)
  { fprintf(stderr,"Error closing log file\n");
    exit(1);
  }
}


/* MOVE TO FIRST RECORD OF LOG FILE.  Reads the header for the first 
   record of the log file.  If the log file is empty, at_end is set. */

void log_file_first
( log_file *logf	/* Log file state structure */
)
{ 
  if (fseek(logf->file_struct,0,0)!=0)
  { fprintf(stderr,"Error moving to first record of log file\n");
    exit(1);
  }
   
  logf->at_end = 0;

  read_header(logf,0);
  logf->at_beginning = !logf->at_end;
}


/* MOVE TO LAST RECORD OF LOG FILE.  Reads the header for the last record
   of the log file.  If the log file is empty, at_end is set. */

void log_file_last
( log_file *logf	/* Log file state structure */
)
{
  log_header trailer;

  if (fseek (logf->file_struct, 0, 2) != 0) 
  { fprintf(stderr,"Error moving to end of log file\n");
    exit(1);
  }

  if (ftell(logf->file_struct)==0)
  { logf->at_end = 1;
    logf->at_beginning = 0;
    return;
  }

  logf->at_end = 0;

  if (fseek (logf->file_struct, -sizeof(log_header), 1) != 0) 
  { fprintf(stderr,"Error moving to start of last trailer in log file\n");
    exit(1);
  }

  read_header(logf,1);

  if (logf->at_end)
  { fprintf(stderr,"Problem reading trailer at end of log file\n");
    exit(1);
  }

  trailer = logf->header;

  if (fseek (logf->file_struct, -2*sizeof(log_header) - trailer.size, 1) != 0) 
  { fprintf(stderr,"Error moving to start of last header in log file\n");
    exit(1);
  }

  logf->at_beginning = ftell(logf->file_struct)==0;

  read_header(logf,0);  

  if (logf->at_end)
  { fprintf(stderr,"Problem reading header at end of log file\n");
    exit(1);
  }

  check_header_trailer(&logf->header,&trailer);
}


/* MOVE FORWARD ONE RECORD.  Reads the header for the record following
   the current one, or sets at_end if there is no next record. */

void log_file_forward
( log_file *logf	/* Log file state structure */
)
{
  log_header trailer;

  if (logf->at_end)
  { fprintf(stderr,"Tried to move forward when at end of log file\n");
    exit(1);
  }

  trailer = logf->header;

  if (fseek (logf->file_struct, logf->header.size, 1)!=0)
  { fprintf(stderr,"Error skipping forward in log file\n");
    exit(1);
  }

  read_header(logf,1);

  if (logf->at_end)
  { fprintf(stderr,"Missing trailer in log file\n");
    exit(1);
  }

  check_header_trailer(&logf->header,&trailer);

  read_header(logf,0);

  logf->at_beginning = 0;
}


/* MOVE BACKWARD ONE RECORD.  Reads the header for the record preceding
   the current one.  It is an error to do this when past the last record 
   or at the first record. */

void log_file_backward
( log_file *logf	/* Log file state structure */
)
{
  log_header trailer;

  if (logf->at_end)
  { fprintf(stderr,"Tried to move backward when at end of log file\n");
    exit(1);
  }
  
  if (logf->at_beginning)
  { fprintf(stderr,"Tried to move backward when at beginning of log file\n");
    exit(1);
  }

  if (fseek (logf->file_struct, -2*sizeof(log_header), 1) != 0) 
  { fprintf(stderr,"Error moving to start of preceding trailer in log file\n");
    exit(1);
  }

  read_header(logf,1);

  if (logf->at_end)
  { fprintf(stderr,"Problem reading trailer going backwards in log file\n");
    exit(1);
  }

  trailer = logf->header;

  if (fseek (logf->file_struct, -2*sizeof(log_header) - trailer.size, 1) != 0) 
  { fprintf(stderr,"Error moving to start of preceding header in log file\n");
    exit(1);
  }

  logf->at_beginning = ftell(logf->file_struct)==0;

  read_header(logf,0);  

  if (logf->at_end)
  { fprintf(stderr,"Problem reading header going backwards in log file\n");
    exit(1);
  }

  check_header_trailer(&logf->header,&trailer);
}


/* READ DATA FROM CURRENT RECORD.  Also reads the header for the next record, 
   or sets at_end if there is no next record.  The caller must allocate 
   sufficient space in the location passed to hold the record, and pass the 
   expected size as the last parameter.  The actual size is available before 
   the call from logf->header.size. */

void log_file_read
( log_file *logf,	/* Log file state structure */
  void *data,		/* Place to store data */
  int size		/* Size of record expected. */
)
{
  log_header trailer;

  if (logf->at_end)
  { fprintf(stderr,"Tried to read data from log file when at end\n");
    exit(1);
  }

  if (logf->header.size!=size)
  { fprintf(stderr,
      "Record has wrong size: Type %c, Actual size %d, Required size %d\n",
      logf->header.type, logf->header.size, size);
    exit(1);
  }

  trailer = logf->header;

  if (fread(data, 1, logf->header.size, logf->file_struct) != logf->header.size)
  { fprintf(stderr,"Error reading data from log file\n");
    exit(1);
  }

  read_header(logf,1);

  if (logf->at_end)
  { fprintf(stderr,"Missing trailer in log file\n");
    exit(1);
  }

  check_header_trailer(&logf->header,&trailer);

  read_header(logf,0);

  logf->at_beginning = 0;
}


/* APPEND RECORD TO LOG FILE.  The type and index for the record are taken
   from the log file state structure, as is the size of the data block.  
   After the data has been appended, we will be positioned at the end of 
   the log file. 

   If the global variable log_append_compare is not zero, records written
   are compared with the record of the same type and index (if present) in the
   "gobble" structure it points to.  If the comparison comes out unequal,
   log_append_compare is set to zero to signal this (and to disable future
   comparisons).
*/

log_gobbled *log_append_compare;

void log_file_append 
( log_file *logf,	/* Log file state structure */
  void *data		/* Data to write */
)
{ 
  log_header trailer;

  if (logf->header.size<0)
  { fprintf(stderr,"Tried to write log record with negative size\n");
    exit(1);
  }
  if (logf->header.type<0 || logf->header.type>=128)
  { fprintf(stderr,"Tried to write log record with bad type: %d\n",
      logf->header.type);
    exit(1);
  }

  if (!logf->at_end)
  {
    if (fseek(logf->file_struct,0,2)!=0)
    { fprintf(stderr,"Error seeking to end of log file\n");
      exit(1);
    }

    logf->at_end = 1;
    logf->at_beginning = 0;
  }

  if (!logf->last_index_known)
  { 
    if (ftell(logf->file_struct)==0)
    { 
      logf->last_index = logf->header.index;  /* Pretend so check is disabled */
    }
    else
    { 
      if (fseek (logf->file_struct, -sizeof(log_header), 1) != 0) 
      { fprintf(stderr,"Error moving to start of last trailer in log file\n");
        exit(1);
      }

      if (fread (&trailer, 1, sizeof trailer, logf->file_struct) 
           != sizeof trailer
       || fseek(logf->file_struct,0,2)!=0)
      { fprintf(stderr,
          "Error reading last trailer before appending to log file\n");
        exit(1);
      }

      if (trailer.magic!=Log_trailer_magic)
      { fprintf(stderr,"End of log file does not look like a trailer\n");
        exit(1);
      }

      logf->last_index = trailer.index;
    }

    logf->last_index_known = 1;
  }

  if (logf->header.index<logf->last_index)
  { fprintf(stderr,
      "Tried to append record with index out of sequence (%d after %d)\n",
      logf->header.index, logf->last_index);
    exit(1);
  }

  logf->header.magic = Log_header_magic;
  logf->header.reserved = 0;	/* So future programs see this in old files */

  if (fwrite (&logf->header, sizeof logf->header, 1, logf->file_struct) != 1)
  { fprintf(stderr,"Error writing header to log file\n");
    exit(1);
  }

  if (fwrite(data, 1, logf->header.size, logf->file_struct)!=logf->header.size)
  { fprintf(stderr,"Error writing data to log file\n");
    exit(1);
  }

  trailer = logf->header;
  trailer.magic = Log_trailer_magic;

  if (fwrite (&trailer, sizeof trailer, 1, logf->file_struct) != 1)
  { fprintf(stderr,"Error writing trailer to log file\n");
    exit(1);
  }

  if (fflush(logf->file_struct)!=0)
  { fprintf(stderr,"Error writing record to log file\n");
    exit(1);
  }

  logf->last_index = logf->header.index;

  if (log_append_compare)
  { int i;
    i = logf->header.type;
    if (log_append_compare->index[i]==logf->header.index)
    { if (log_append_compare->actual_size[i]!=logf->header.size
       || memcmp(log_append_compare->data[i],data,logf->header.size)!=0)
      { log_append_compare = 0;
      }
    }
  }
}


/* READ HEADER/TRAILER FROM LOG FILE. */

static void read_header
( log_file *logf,	/* Log file state structure */
  int trailer		/* Is this supposed to be a trailer? */
)
{
  int n;

  n = fread (&logf->header, 1, sizeof logf->header, logf->file_struct);

  if (n==0)
  { logf->at_end = 1;
  }
  else if (n == sizeof logf->header)
  { if (logf->header.size<0)
    { fprintf(stderr,"Negative record size in log file header/trailer\n");
      exit(1);
    }
    if (logf->header.type<0 || logf->header.type>=128)
    { fprintf(stderr,"Bad type in log file header/trailer: %d\n",
        logf->header.type);
      exit(1);
    }
    if (logf->header.magic != (trailer ? Log_trailer_magic : Log_header_magic))
    { fprintf(stderr,"Bad magic number in header/trailer in log file (%04x)\n",
        logf->header.magic);
      exit(1);
    }
  }
  else
  { fprintf(stderr,"Error reading header/trailer from log file\n");
    exit(1);
  }
}


/* CHECK THAT HEADER AND TRAILER ARE CONSISTENT. */

static void check_header_trailer
( log_header *h,
  log_header *t
)
{
  if (h->type != t->type || h->index != t->index || h->size  != t->size)
  { fprintf(stderr,"Log record has inconsistent header/trailer\n");
    exit(1);
  }
}


/* INITIALIZE SET OF GOBBLED RECORDS TO BE EMPTY.  The 'fr' parameter controls
   whether the space currently in use should be freed.  If the structure being
   initialized currently contains garbage, 'fr' should be zero! */

void log_gobble_init
( log_gobbled *logg,	/* Records set of records gobbled up */
  int fr		/* Should old data be freed? */
)
{
  int c;

  logg->last_index = -1;

  for (c = 0; c<128; c++)
  { logg->req_size[c] = -1;
    logg->actual_size[c] = -1;
    logg->index[c] = -1;
    if (fr && logg->data[c]!=0) 
    { free(logg->data[c]);
    }
    logg->data[c] = 0;
  }
}


/* GOBBLE UP RECORDS WITH THE CURRENT INDEX.  Reads all consecutive records 
   that have the same index as the current record.  The data from these 
   records is stored in places pointed to by the 'data' pointers in the 
   log_gobbled structure passed.  The 'index' field of the log_gobbled
   structure is set to the index of the records read.  Note that data from
   records with previous indexes may still remain in the log_gobbled
   structure. 

   Data is stored by record type, with any previous record of the same type
   being discarded.  If the req_size field for the type in the log_gobbled
   structure is set to other than -1, the record must be of that size, otherwise
   any size is allowed.  The actual size of the currently residing record
   is stored in the actual_size field.  Space to hold the data is allocated 
   by this procedure if the current data pointer is null, or if the size of
   the current record is not the same as that of the new record. */

void log_gobble
( log_file *logf,	/* Log file state structure */
  log_gobbled *logg	/* Structure for holding gobbled data */
)
{ 
  int type;

  if (logf->at_end)
  { fprintf(stderr,"Trying to gobble records past end of log file\n");
    exit(1);
  }

  logg->last_index = logf->header.index;

  while (!logf->at_end && logf->header.index==logg->last_index)
  { 
    type = logf->header.type & 0177;

    if (logg->req_size[type]!=-1 && logg->req_size[type]!=logf->header.size)
    { fprintf(stderr,
        "Record has wrong size: Type %c, Actual size %d, Required size %d\n",
        type, logf->header.size, logg->req_size[type]);
      exit(1);
    }

    if (logg->data[type]==0 || logg->actual_size[type]!=logf->header.size)
    { if (logg->data[type]!=0) 
      { free(logg->data[type]);
      }
      logg->data[type] = malloc(logf->header.size);
      if (logg->data[type]==0)
      { fprintf(stderr,"Not enough memory!\n");
        exit(1);
      }
    }

    logg->actual_size[type] = logf->header.size;
    logg->index[type] = logf->header.index;

    log_file_read (logf, logg->data[type], logf->header.size);
  }
}


/* GOBBLE RECORDS AT END OF LOG FILE.  Gobbles all records with index 
   equal to the index of the last record in the log file, provided that
   index is non-negative.  Leaves the file pointer at the end.  Returns
   the index following that of the last record, or one if there are
   no records with non-negative indexes. */

int log_gobble_last
( log_file *logf,	/* Log file state structure */
  log_gobbled *logg	/* Structure for holding gobbled data */
)
{
  int index;

  log_file_last(logf);

  if (logf->at_end || logf->header.index<0)
  { 
    if (!logf->at_end)
    { log_file_forward(logf);
    }

    return 1;
  }

  index = logf->header.index;

  while (!logf->at_beginning && logf->header.index==index)
  { log_file_backward(logf);
  }

  if (logf->header.index!=index) 
  { log_file_forward(logf);
  } 

  log_gobble(logf,logg);

  return index+1;
}
