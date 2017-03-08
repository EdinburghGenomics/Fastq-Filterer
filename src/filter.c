//
//  filter.c
//  Fastq-Filterer
//
//  Created by Murray Wham on 05/09/2016.
//  Copyright (c) 2017 Edinburgh Genomics (see ../LICENCE)
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
char *r1i_path = NULL, *r2i_path = NULL, *r1o_path = NULL, *r2o_path = NULL, *stats_file = NULL;
int read_pairs_checked = 0, read_pairs_removed = 0, read_pairs_remaining = 0;


static void timestamp() {

    time_t t = time(NULL);
    struct tm* now = localtime(&t);
    
    printf(
        "[%i-%i-%i %i:%i:%i][fastq_filterer] ",
        now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min,  now->tm_sec
    );
}


static void help_msg() {
    printf(
"Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> [--o1 <r1_filtered.fastq> --o2 <r2_filtered.fastq>] \
--threshold <filter_threshold> [--stats_file <fastq_filterer.stats>] [--unsafe]\n\
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


static int filter_fastqs() {
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
            int ret_val = 0;
            if (*r1_header != *r2_header) {  // if either file is not finished
                timestamp();
                printf("Input fastqs have differing numbers of reads at line %i\n", read_pairs_checked * 4);
                ret_val = 1;
            }
            
            gzclose(r1i);
            gzclose(r2i);
            fclose(r1o);
            fclose(r2o);
            return ret_val;

        } else if ((strlen(r1_seq) > threshold) && (strlen(r2_seq) > threshold)) {
            read_pairs_checked++;
            read_pairs_remaining++;
            
            fputs(r1_header, r1o);
            fputs(r1_seq, r1o);
            fputs(r1_strand, r1o);
            fputs(r1_qual, r1o);
            
            fputs(r2_header, r2o);
            fputs(r2_seq, r2o);
            fputs(r2_strand, r2o);
            fputs(r2_qual, r2o);
            
        } else {
            read_pairs_checked++;
            read_pairs_removed++;
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
}


static char* build_output_path(char* input_path) {
    /*
     Convert, e.g, basename.fastq to basename_filtered.fastq. Used when output fastq paths are not specified.
     */
    
    int file_ext_len = 6;  // .fastq
    
    if (input_path[strlen(input_path) - 1] == 'z') {
        file_ext_len = 9;  // .fastq.gz
    }
    
    size_t basename_len = strlen(input_path) - file_ext_len;
    char* output_path = malloc(sizeof (char) * (basename_len + 15));  // _filtered.fastq
    strncpy(output_path, input_path, basename_len);
    output_path[basename_len] = '\0';
    strcat(output_path, "_filtered.fastq");
    return output_path;
}


static void check_file_paths() {
    if (r1i_path == NULL || r2i_path == NULL) exit(1);
    
    if (r1o_path == NULL) {
        timestamp();
        printf("No o1 argument given - deriving from i1\n");
        r1o_path = build_output_path(r1i_path);
    }
    if (r2o_path == NULL) {
        timestamp();
        printf("No o2 argument given - deriving from i2\n");
        r2o_path = build_output_path(r2i_path);
    }
    
    timestamp();
    printf("R1: %s -> %s\n", r1i_path, r1o_path);
    timestamp();
    printf("R2: %s -> %s\n", r2i_path, r2o_path);
    timestamp();
    printf("Filter threshold: %i\n", threshold);
}



static void output_stats() {
    FILE* f = fopen(stats_file, "w");
    
    char* stats = malloc(sizeof (char) * (83 + strlen(r1i_path) + strlen(r2i_path) + strlen(r1o_path) + strlen(r2o_path) + 24));
    sprintf(
        stats,
        "r1i %s\nr2i %s\nr1o %s\nr2o %s\nread_pairs_checked %i\nread_pairs_removed %i\nread_pairs_remaining %i\n",
        r1i_path, r2i_path, r1o_path, r2o_path, read_pairs_checked, read_pairs_removed, read_pairs_remaining
    );
    fputs(stats, f);
    free(stats);
    fclose(f);
}



int main(int argc, char* argv[]) {
    int arg;
    
    static struct option args[] = {
        {"help", no_argument, 0, 'h'},
        {"unsafe", no_argument, 0, 'f'},
        {"stats_file", required_argument, 0, 'w'},
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
                r1i_path = malloc(sizeof optarg);
                r1i_path = optarg;
                break;
            case 's':
                r2i_path = malloc(sizeof optarg);
                r2i_path = optarg;
                break;
            case 'i':
                r1o_path = malloc(sizeof optarg);
                r1o_path = optarg;
                break;
            case 'j':
                r2o_path = malloc(sizeof optarg);
                r2o_path = optarg;
                break;
            case 'h':
                help_msg();
                exit(0);
                break;
            case 'f':
                read_func = readln_unsafe;
                break;
            case 'w':
                stats_file = malloc(sizeof optarg);
                stats_file = optarg;
                break;
            case 't':
                threshold = atoi(optarg);
                break;
            default:
                exit(1);
        }
    }
    
    check_file_paths();
    int exit_status = filter_fastqs();
    
    timestamp();
    printf(
        "Checked %i read pairs, %i removed, %i remaining. Exit status %i\n",
        read_pairs_checked, read_pairs_removed, read_pairs_remaining, exit_status
    );

    if (stats_file != NULL) {
        timestamp();
        printf("Writing stats file %s\n", stats_file);
        output_stats();
    }
    
    return exit_status;
}
