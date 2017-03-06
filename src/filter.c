//
//  filter.c
//  Fastq-Filterer
//
//  Created by Murray Wham on 05/09/2016.
//  Copyright (c) 2016 Edinburgh Genomics (see ../LICENCE)
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include <getopt.h>
#include <time.h>

#define block_size 2048
#define unsafe_block_size 4096

int threshold;

static void _log(char* msg) {
    
    time_t t = time(NULL);
    struct tm* now = localtime(&t);
    
    printf(
        "[%i-%i-%i %i:%i:%i][fastq_filterer] %s\n",
        now->tm_year + 1900, now->tm_mon, now->tm_mday, now->tm_hour, now->tm_min,  now->tm_sec,
        msg
    );
}


static void help_msg() {
    printf(
"Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> [--o1 <r1_filtered.fastq> --o2 <r2_filtered.fastq>] \
--threshold <filter_threshold> [--unsafe]\n\
Fastq or fastq.gz files can be read in, but output will always be uncompressed.\
--unsafe uses a faster read function, but will chop lines over 4096 characters.\
"
    );
}


static char* readln_unsafe(gzFile* f) {
    char* line = malloc(unsafe_block_size);
    line[0] = '\0';
    gzgets(f, line, unsafe_block_size);
    return line;
}


static char* readln(gzFile* f) {
    /*
     Read a line from a file. Each time this function is called on a file, the next
     line is read. Memory is dynamically allocated to allow reading of lines of any
     length.
     
     :input FILE* f: The file to read from
     :output: char* line (the line read), or null '\0' string if no input is available.
     */
    
    int _block_size = block_size;
    char* line;
    line = malloc(sizeof (char) * _block_size);
    line[0] = '\0';
    
    do {
        _block_size += block_size;
        line = realloc(line, _block_size + 1);
        char* _line_part = malloc(sizeof (char) * block_size);
        _line_part[0] = '\0';

        if (gzgets(f, _line_part, block_size)) {
            strcat(line, _line_part);
            free(_line_part);
        } else {
            free(_line_part);
            return line;
        }
    } while (line[strlen(line) - 1] != '\n');

    return line;
}

char* (*read_func)(gzFile*) = readln;


static int read_fastqs(char* r1i_path, char* r2i_path, char* r1o_path, char* r2o_path) {
    /*
     Read two fastqs (R1.fastq, R2.fastq) entry by entry, check whether the R1 and R2 for each read
     are both long enough, and output them to R1_filtered.fastq and R2_filtered.fastq if they are.
     
     :input char* r1_path: Path to R1.fastq input file
     :input char* r1_filtered: Path to R1_filtered.fastq output file
     :input char* r2_path: Path to R1.fastq input file
     :input char* r2_filtered: Path to R1_filtered.fastq output file
     */
    
    gzFile* r1i = gzopen(r1i_path, "r");
    gzFile* r2i = gzopen(r2i_path, "r");
    FILE* r1o = fopen(r1o_path, "w");
    FILE* r2o = fopen(r2o_path, "w");
    
    char* r1_header;
    char* r1_seq;
    char* r1_strand;
    char* r1_qual;

    char* r2_header;
    char* r2_seq;
    char* r2_strand;
    char* r2_qual;

    while (1) {
        r1_header = read_func(r1i);  // @read_1 1
        r1_seq = read_func(r1i);     // ATGCATGC
        r1_strand = read_func(r1i);  // +
        r1_qual = read_func(r1i);    // #--------

        r2_header = read_func(r2i);  // @read_1 2
        r2_seq = read_func(r2i);     // ATGCATGC
        r2_strand = read_func(r2i);  // -
        r2_qual = read_func(r2i);    // #--------
        
        if (*r1_header == '\0' || *r2_header == '\0') {
            if (*r1_header == *r2_header) {
                return 0;
            } else {
                return 1;
            }
        } else if ((strlen(r1_seq) > threshold) && (strlen(r2_seq) > threshold)) {
            fputs(r1_header, r1o);
            fputs(r1_seq, r1o);
            fputs(r1_strand, r1o);
            fputs(r1_qual, r1o);
            
            fputs(r2_header, r2o);
            fputs(r2_seq, r2o);
            fputs(r2_strand, r2o);
            fputs(r2_qual, r2o);
        }
        free(r1_header);
        free(r1_seq);
        free(r1_strand);
        free(r1_qual);
        
        free(r2_header);
        free(r2_seq);
        free(r2_strand);
        free(r2_qual);
    }
    return 0;
}


static char* build_output_path(char* input_path, char* file_ext) {
    /*
     Converts, e.g, basename.fastq to basename_filtered.fastq. Used when output fastq paths are not specified.
     */
    
    static char filtered[10] = "_filtered";
    size_t basename_len = strlen(input_path) - strlen(file_ext);
    char* basename = malloc(sizeof (char) * basename_len);
    strncat(basename, input_path, basename_len);
    char* new_output_path = malloc(sizeof (char) * (strlen(basename) + strlen(filtered) + strlen(file_ext)));
    strcpy(new_output_path, basename);
    strcat(new_output_path, filtered);
    strcat(new_output_path, file_ext);
    return new_output_path;
}


static int filter_fastqs(char* r1i_path, char* r2i_path, char* r1o_path, char* r2o_path) {
    if (r1o_path == NULL) {
        _log("No o1 argument given - deriving from i1");
        r1o_path = build_output_path(r1i_path, ".fastq");
    }
    if (r2o_path == NULL) {
        _log("No o2 argument given - deriving from i2");
        r2o_path = build_output_path(r2i_path, ".fastq");
    }
    
    char* msg1 = malloc(sizeof (char) * (strlen(r1i_path) + strlen(r1o_path) + 9));
    sprintf(msg1, "R1: %s -> %s", r1i_path, r1o_path);
    char* msg2 = malloc(sizeof (char) * (strlen(r2i_path) + strlen(r2o_path) + 9));
    sprintf(msg2, "R2: %s -> %s", r2i_path, r2o_path);
    
    _log(msg1);
    _log(msg2);
    
    return read_fastqs(r1i_path, r2i_path, r1o_path, r2o_path);
}


int main(int argc, char* argv[]) {
    int arg;
    char *r1i = NULL, *r2i = NULL, *r1o = NULL, *r2o = NULL;
    
    static struct option args[] = {
        {"help", no_argument, 0, 'h'},
        {"unsafe", no_argument, 0, 'f'},
        {"threshold", required_argument, 0, 't'},
        {"i1", required_argument, 0, 'r'},
        {"i2", required_argument, 0, 's'},
        {"o1", required_argument, 0, 'i'},
        {"o2", required_argument, 0, 'j'},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    
    while ((arg = getopt_long(argc, argv, "h:t:r:s:i:j:", args, &opt_idx)) != -1) {
        switch(arg) {
            case 'r':
                r1i = malloc(sizeof optarg);
                r1i = optarg;
                break;
            case 's':
                r2i = malloc(sizeof optarg);
                r2i = optarg;
                break;
            case 'i':
                r1o = malloc(sizeof optarg);
                r1o = optarg;
                break;
            case 'j':
                r2o = malloc(sizeof optarg);
                r2o = optarg;
                break;
            case 'h':
                help_msg();
                exit(0);
                break;
            case 'f':
                read_func = readln_unsafe;
                break;
            case 't':
                threshold = atoi(optarg);
                break;
            default:
                exit(1);
        }
    }
    if (r1i == NULL || r2i == NULL) exit(1);
    
    char* msg1 = malloc(sizeof (char) * (23));  // 19 + 4 chars, so can take a threshold up to 9999
    sprintf(msg1, "Filter threshold: %i", threshold);
    _log(msg1);
    
    int exit_status = filter_fastqs(r1i, r2i, r1o, r2o);
    char* msg2 = malloc(sizeof (char) * (31));  // 28 + 3 chars, so can take an exit status up to 999
    sprintf(msg2, "Completed with exit status %i", exit_status);
    _log(msg2);
    
    return exit_status;
}
