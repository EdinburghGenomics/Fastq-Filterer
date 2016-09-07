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

#define initial_block_size 64
#define block_increment 64

int threshold;


static void _log(char* msg) {
    printf("[fastq_filterer] %s\n", msg);
}


static void help_msg() {
    printf(
"Usage: fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> [--o1 <r1_filtered.fastq> --o2 <r2_filtered.fastq>] \
--threshold <filter_threshold>\n\
.fastq.gz files can also be input, in which case zlib compression will be used (note: compression \
is decided based on the file extension of the input file, not the output)\n"
    );
}


static char* readln(FILE* f) {
    /*
     Read a line from a file. Each time this function is called on a file, the next
     line is read. Memory is dynamically allocated to allow reading of lines of any
     length.
     
     :input FILE* f: The file to read from
     :output: char* line (the line read), or null '\0' string if no input is available.
     */
    
    int _block_size = initial_block_size;
    char* line;
    line = malloc(sizeof (char) * _block_size);
    
    while (line[strlen(line) - 1] != '\n') {
        _block_size += block_increment;
        line = realloc(line, _block_size + 1);
        char _line_part[block_increment];
        
        char* p = fgets(_line_part, block_increment, f);
        if (p == NULL) {
            return line;
        }
        strcat(line, _line_part);
    }
    return line;
}


static char* readlngz(gzFile* f) {
    /*
     Equivalent of _readln for gz compressed files.
     
     :input gzFile* f: The gz file to read from
     :output: As _readln (decompressed)
     */
    
    int _block_size = initial_block_size;
    char* line;
    line = malloc(sizeof (char) * _block_size);
    
    while (line[strlen(line) - 1] != '\n') {
        _block_size += block_increment;
        line = realloc(line, _block_size + 1);
        char _line_part[block_increment];
        
        char* p = gzgets(f, _line_part, block_increment);
        if (p == NULL) {
            return line;
        }
        strcat(line, _line_part);
    }
    return line;
}


static char* find_file_ext(char* file_path) {
    /*
     Identify the extension of a file.
     
     :input char* file_path: The file path to identify
     :output: char* ext (the file extension) if it's known, else NULL
     */
    
    char* ext;
    char extensions[4][10] = {
        ".fastq",
        ".fastq.gz",
        ".fq",
        ".fq.gz"
    };
    
    for (int i=0; i<4; i++) {
        ext = strstr(file_path, extensions[i]);
        if (ext != NULL) {
            return ext;
        }
    }
    return NULL;
}


static int read_fastqs(char* r1i_path, char* r2i_path, char* r1o_path, char* r2o_path) {
    /*
     Read two fastqs (R1.fastq, R2.fastq) entry by entry, check whether the R1 and R2 for each read
     are both long enough, and output them to R1_filtered.fastq and R2_filtered.fastq if they are.
     
     :input char* r1_path: Path to R1.fastq input file
     :input char* r1_filtered: Path to R1_filtered.fastq output file
     :input char* r2_path: Path to R1.fastq input file
     :input char* r2_filtered: Path to R1_filtered.fastq output file
     */
    
    FILE* r1i = fopen(r1i_path, "r");
    FILE* r2i = fopen(r2i_path, "r");
    FILE* r1o = fopen(r1o_path, "w");
    FILE* r2o = fopen(r2o_path, "w");
    
    while (1) {
        char* r1_header = readln(r1i);  // @read_1 1
        char* r1_seq = readln(r1i);     // ATGCATGC
        char* r1_strand = readln(r1i);  // +
        char* r1_qual = readln(r1i);    // #--------
        
        char* r2_header = readln(r2i);  // @read_1 2
        char* r2_seq = readln(r2i);     // ATGCATGC
        char* r2_strand = readln(r2i);  // -
        char* r2_qual = readln(r2i);    // #--------
        
        if (*r1_header == '\0') {
            return 0;
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


static int read_gz_fastqs(char* r1i_path, char* r2i_path, char* r1o_path, char* r2o_path) {
    /*
     Equivalent of read_fastqs, but for gzipped fastqs. Outputs gz compressed fastqs.
     */
    gzFile* r1i = gzopen(r1i_path, "rb");
    gzFile* r2i = gzopen(r2i_path, "rb");
    gzFile* r1o = gzopen(r1o_path, "w");
    gzFile* r2o = gzopen(r2o_path, "w");
    
    while (1) {
        char* r1_header = readlngz(r1i);  // @read_1 1
        char* r1_seq = readlngz(r1i);     // ATGCATGC
        char* r1_strand = readlngz(r1i);  // +
        char* r1_qual = readlngz(r1i);    // #--------
        
        char* r2_header = readlngz(r2i);  // @read_1 2
        char* r2_seq = readlngz(r2i);     // ATGCATGC
        char* r2_strand = readlngz(r2i);  // -
        char* r2_qual = readlngz(r2i);    // #--------
        
        if (*r1_header == '\0') {
            gzflush(r1o, Z_FINISH);
            gzflush(r2o, Z_FINISH);
            return 0;
        } else if ((strlen(r1_seq) > threshold) && (strlen(r2_seq) > threshold)) {
            gzputs(r1o, r1_header);
            gzputs(r1o, r1_seq);
            gzputs(r1o, r1_strand);
            gzputs(r1o, r1_qual);
            
            gzputs(r2o, r2_header);
            gzputs(r2o, r2_seq);
            gzputs(r2o, r2_strand);
            gzputs(r2o, r2_qual);
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
    
    char* file_ext = find_file_ext(r1i_path);
    
    if (r1o_path == NULL) {
        _log("No o1 argument given - deriving from i1");
        r1o_path = build_output_path(r1i_path, file_ext);
    }
    if (r2o_path == NULL) {
        _log("No o2 argument given - deriving from i2");
        r2o_path = build_output_path(r2i_path, file_ext);
    }
    printf("[fastq_filterer] Input R1 file: %s\n", r1i_path);
    printf("[fastq_filterer] Input R2 file: %s\n", r2i_path);
    printf("[fastq_filterer] Output R1 file: %s\n", r1o_path);
    printf("[fastq_filterer] Output R2 file: %s\n", r2o_path);
    
    int gz;
    printf("[fastq_filterer] Using file extension '%s'\n", file_ext);
    
    char* gz_ext = strstr(file_ext, ".gz");
    if (gz_ext != NULL) {
        _log(".gz file extension detected - using gz file io");
        gz = 1;
    } else {
        _log("no .gz file extension detected - using basic file io");
        gz = 0;
    }
    
    if (gz == 1) {
        return read_gz_fastqs(r1i_path, r2i_path, r1o_path, r2o_path);
    } else {
        return read_fastqs(r1i_path, r2i_path, r1o_path, r2o_path);
    }
}


int main(int argc, char* argv[]) {
    
    int arg;
    char *r1i = NULL, *r2i = NULL, *r1o = NULL, *r2o = NULL;
    
    static struct option args[] = {
        {"help", no_argument, 0, 'h'},
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
            case 't':
                threshold = atoi(optarg);
                break;
            default:
                exit(1);
        }
    }
    if (r1i == NULL || r2i == NULL) exit(1);
    
    printf("[fastq_filterer] Filter threshold: %i\n", threshold);
    int exit_status = filter_fastqs(r1i, r2i, r1o, r2o);
    printf("[fastq_filterer] Completed with exit status %i\n", exit_status);
    return exit_status;
}
