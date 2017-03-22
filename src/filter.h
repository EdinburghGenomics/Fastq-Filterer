//
//  filter.h
//  Fastq-Filterer
//
//  Created by Murray Wham on 22/03/2017.
//  Copyright (c) 2017 Edinburgh Genomics (see ../LICENCE)
//

#ifndef FastqFilterer_h
#define FastqFilterer_h

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "Fastq-Filterer"
#endif

#ifndef VERSION
#define VERSION "0.2"
#endif

#ifndef USAGE
#define USAGE "\
Fastq-Filterer\n\
Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> [--o1 <r1_filtered.fastq> --o2 <r2_filtered.fastq>] \
--threshold <filter_threshold>\n\
Fastq or fastq.gz files can be read in, but output will always be uncompressed.\n\
Options:\n\
--stats_file <fastq_filterer.stats> - write a file summarising the read pairs checked and removed\n\
--unsafe - use a simpler read function which is faster, but will chop lines over 4096 characters\n\
"
#endif

#endif
