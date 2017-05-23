#ifndef FastqFilterer_h
#define FastqFilterer_h

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "Fastq-Filterer"
#endif

#ifndef VERSION
#define VERSION "0.3.1"
#endif

#ifndef USAGE
#define USAGE "\
Fastq-Filterer\n\
Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> [--o1 <r1_filtered.fastq> --o2 <r2_filtered.fastq>] \
--threshold <filter_threshold>\n\
Fastq or fastq.gz files can be read in, but output will always be uncompressed.\n\
Options:\n\
--stats_file <stats_file> - write a file summarising the read pairs checked and removed\n\
--unsafe - use a simpler read function which is faster, but will chop lines over 4096 characters\n\
--remove_tiles <tile1,tile2,tile3...> - comma-separated list of tile ids to remove regardless of length\n\
--trim_r1 <max_len> - trim all reads for r1.fastq to a maximum length\n\
--trim_r2 <max_len> - as above for r2.fastq\n\
\n"
#endif

#endif
