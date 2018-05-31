#ifndef FastqFilterer_h
#define FastqFilterer_h

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "Fastq-Filterer"
#endif

#ifndef VERSION
#define VERSION "0.4"
#endif

#ifndef USAGE
#define USAGE "\
Fastq-Filterer\n\
Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> --threshold <filter_threshold>\n\
Fastq or fastq.gz files can be read in, but output will always be uncompressed.\n\
Options:\n\
--o1 <r1_filtered.fastq> - output file name for r1 (defaults to <input_path>_filtered.fastq)\n\
--o2 <r2_filtered.fastq> - as above for r2\n\
--f1 <r1_filtered_reads.fastq> - filtered reads file name for r1 (defaults to <input_path_filtered_reads.fastq)\n\
--f2 <r2_filtered_reads.fastq> - as above for r2\n\
--stats_file <stats_file> - write a file summarising the read pairs checked and removed\n\
--unsafe - use a simpler read function which is faster, but will chop lines over 4096 characters\n\
--remove_tiles <tile1,tile2,tile3...> - comma-separated list of tile ids to remove regardless of length\n\
--remove_reads <rm_reads.txt> - text file containing read names to filter out\n\
--trim_r1 <max_len> - trim all reads in the r1 output file to a maximum length\n\
--trim_r2 <max_len> - as above for r2\n\
\n"
#endif

#endif
