#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <zlib.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
#include "uthash.h"
#include "filter.h"

#define block_size 2048
#define unsafe_block_size 4096

int threshold = -1;
bool quiet = false;
char *r1i_path = NULL, *r1o_path = NULL, *r1f_path = NULL;
char *r2i_path = NULL, *r2o_path = NULL, *r2f_path = NULL;
char *remove_reads_path = NULL;
int read_pairs_checked = 0, read_pairs_removed = 0, read_pairs_remaining = 0;
int trim_r1, trim_r2;
char* remove_tiles;
char** tiles_to_remove;


static void _log(char* fmt_str, ...) {
    if (quiet) {
        return;
    }
    
    time_t t = time(NULL);
    struct tm* now = localtime(&t);
    
    printf(
        "[%i-%i-%i %i:%i:%i][fastq_filterer] ",
        now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min,  now->tm_sec
    );
    va_list args;
    va_start(args, fmt_str);
    vprintf(fmt_str, args);
    va_end(args);
    
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
    char* line = malloc(sizeof (char) * _block_size);
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


typedef struct {
    char *header, *seq, *strand, *qual;
} FastqRead;


typedef struct {
    FastqRead r1, r2;
} FastqReadPair;


static char* get_tile_id(char* fastq_header) {
    char* field;
    field = strtok(fastq_header, ":");
    int i;
    for (i=0; i<4; i++) {
        field = strtok(NULL, ":");  // walk along the header 4 times to the tile ID
    }
    return field;
}


bool (**criteria)(FastqReadPair);
int ncriteria = 0;


static bool std_check_read(FastqReadPair read_pair) {
    if ((strlen(read_pair.r1.seq) > threshold) && (strlen(read_pair.r2.seq) > threshold)) {
        return true;
    } else {
        return false;
    }
}


static bool tile_check_read(FastqReadPair read_pair) {
    // get_tile_id will modify the string passed to it with strtok, so use a copy
    char* _read_id = malloc(sizeof (char) * (strlen(read_pair.r1.header) + 1));
    strcpy(_read_id, read_pair.r1.header);
    char* tile_id = get_tile_id(_read_id);
    int i = 0;
    
    char* comp = tiles_to_remove[i];
    
    bool ret_val = true;
    while (comp != NULL) {  // check for null terminator at end of remove_tiles
        if (strcmp(comp, tile_id) == 0) {
            ret_val = false;
        }
        i++;
        comp = tiles_to_remove[i];
    }
    
    free(_read_id);
    return ret_val;
}


typedef struct {
    const char* key;
    UT_hash_handle hh;
} HashTable;

HashTable* reads_to_remove;
HashTable* mask;


static bool id_check_read(FastqReadPair read_pair) {
    char* _read_id = malloc(sizeof (char) * (strlen(read_pair.r1.header) + 1));
    strcpy(_read_id, read_pair.r1.header);
    
    char* coord_information = strtok(_read_id, " ");
    
    HASH_FIND_STR(reads_to_remove, coord_information, mask);
    bool ret_val = true;
    if (mask) {
        ret_val = false;
    }
    free(_read_id);
    return ret_val;
}


static void std_include(FastqRead read, FILE* outfile) {
    fputs(read.header, outfile);
    fputs(read.seq, outfile);
    fputs(read.strand, outfile);
    fputs(read.qual, outfile);
}


static void _trim_include(FastqRead read, FILE* outfile, int trim_len) {
    if (strlen(read.seq) > trim_len + 1) {  // add 1 here to compensate for \n at end of line...
        read.seq[trim_len] = '\n';  // ...but 0-indexing means we don't need to add 1 here
        read.seq[trim_len + 1] = '\0';
        read.qual[trim_len] = '\n';
        read.qual[trim_len + 1] = '\0';
    }
    std_include(read, outfile);
}

static void trim_include_r1(FastqRead read, FILE* outfile) {
    _trim_include(read, outfile, trim_r1);
}


static void trim_include_r2(FastqRead read, FILE* outfile) {
    _trim_include(read, outfile, trim_r2);
}


void (*include_func_r1)(FastqRead, FILE*) = std_include;
void (*include_func_r2)(FastqRead, FILE*) = std_include;


static int filter_fastqs() {
    /*
     Read two fastqs, R1 and R2, entry by entry, checking whether the R1 and R2 for each read
     are both long enough, and output them to Rx_filtered.fastq if they are. If not, output them to
     Rx_filtered_reads.fastq.
     
     :input char* r1_path: Path to R1.fastq input file
     :input char* r1_filtered: Path to R1_filtered.fastq output file
     :input char* r2_path: Path to R1.fastq input file
     :input char* r2_filtered: Path to R1_filtered.fastq output file
     */
    
    gzFile* r1i = gzopen(r1i_path, "r");
    gzFile* r2i = gzopen(r2i_path, "r");
    FILE* r1o = fopen(r1o_path, "w");
    FILE* r2o = fopen(r2o_path, "w");
    FILE* r1f = fopen(r1f_path, "w");
    FILE* r2f = fopen(r2f_path, "w");
    
    FastqReadPair read_pair;
    
    while (true) {
        read_pair.r1.header = read_func(r1i);  // @read_1 1
        read_pair.r1.seq = read_func(r1i);     // ATGCATGC
        read_pair.r1.strand = read_func(r1i);  // +
        read_pair.r1.qual = read_func(r1i);    // #--------

        read_pair.r2.header = read_func(r2i);  // @read_1 2
        read_pair.r2.seq = read_func(r2i);     // ATGCATGC
        read_pair.r2.strand = read_func(r2i);  // -
        read_pair.r2.qual = read_func(r2i);    // #--------
        
        if (*read_pair.r1.header == '\0' || *read_pair.r2.header == '\0') {
            int ret_val = 0;
            if (*read_pair.r1.header != *read_pair.r2.header) {  // if either file is not finished
                _log("Input fastqs have differing numbers of reads, from line %i\n", read_pairs_checked * 4);
                ret_val = 1;
            }
            
            free(read_pair.r1.header);
            free(read_pair.r1.seq);
            free(read_pair.r1.strand);
            free(read_pair.r1.qual);
            
            free(read_pair.r2.header);
            free(read_pair.r2.seq);
            free(read_pair.r2.strand);
            free(read_pair.r2.qual);
            
            gzclose(r1i);
            gzclose(r2i);
            fclose(r1o);
            fclose(r2o);
            fclose(r1f);
            fclose(r2f);

            return ret_val;

        } else {
            bool read_included = true;
            int i;
            for (i=0; i<ncriteria + 1; i++) {
                bool (*func)(FastqReadPair) = criteria[i];
                //if (criteria[i](read_pair) == false) {
                if (func(read_pair) == false) {
                    read_included = false;
                    //break;
                }
            }
            
            read_pairs_checked++;
            if (read_included == true) {
                // include reads
                read_pairs_remaining++;
            
                include_func_r1(read_pair.r1, r1o);
                include_func_r2(read_pair.r2, r2o);
            } else {
                // exclude reads
                read_pairs_removed++;
                
                std_include(read_pair.r1, r1f);
                std_include(read_pair.r2, r2f);
            }
        }
        
        free(read_pair.r1.header);
        free(read_pair.r1.seq);
        free(read_pair.r1.strand);
        free(read_pair.r1.qual);
        
        free(read_pair.r2.header);
        free(read_pair.r2.seq);
        free(read_pair.r2.strand);
        free(read_pair.r2.qual);
    }
}


static char* build_output_path(char* input_path, char* new_extension) {
    /*
     Convert, e.g, basename.fastq to basename_filtered.fastq. Used when output fastq paths are not specified.
     */
    
    int file_ext_len = 6;  // .fastq
    
    if (input_path[strlen(input_path) - 1] == 'z') {
        file_ext_len = 9;  // .fastq.gz
    }
    
    size_t basename_len = strlen(input_path) - file_ext_len;
    char* output_path = malloc(sizeof (char) * (basename_len + strlen(new_extension) + 1));
    strncpy(output_path, input_path, basename_len);
    output_path[basename_len] = '\0';
    strcat(output_path, new_extension);
    return output_path;
}


static void build_remove_tiles() {
    if (remove_tiles == NULL) {  // no --remove_tiles argument
        return;
    }
    
    char* rm_tiles = malloc(sizeof (char) * (strlen(remove_tiles) + 1));
    strcpy(rm_tiles, remove_tiles);  // use a copy for strtok
    
    int ntiles = 1;
    char* comma = strchr(rm_tiles, ',');
    while (comma != NULL) {
        ntiles++;
        comma = strchr(comma+1, ',');
    }
    
    tiles_to_remove = malloc(sizeof (char*) * (ntiles + 1));
    int i = 0;
    char* field;
    field = strtok(rm_tiles, ",");
    while (field != NULL) {
        tiles_to_remove[i] = malloc(sizeof (char) * strlen(field) + 1);
        strcpy(tiles_to_remove[i], field);
        field = strtok(NULL, ",");
        i++;
    }
    tiles_to_remove[i] = NULL;  // set a null terminator
    free(rm_tiles);
}


static void build_remove_reads() {
    if (remove_reads_path == NULL) {
        return;
    }
    
    gzFile* rm_reads = gzopen(remove_reads_path, "r");
    reads_to_remove = NULL;
    mask = NULL;
    int i = 0;
    char* matchable_element;
    char* read_id;
    char* line;
    
    while (true) {
        line = readln(rm_reads);
        if (line == NULL || *line == '\0') {
            gzclose(rm_reads);
            return;
        }

        // match up to the first space
        char* space_match = strchr(line, ' ');
        if (space_match == NULL) {
            read_id = strtok(line, "\n");
        } else {
            read_id = strtok(line, " ");
        }
        
        if (read_id != NULL) {  // this can happen if there's a blank line in the file, i.e. '\n'
            matchable_element = malloc(sizeof (char) * strlen(read_id) + 2);  // include the @ and \0
            strcpy(matchable_element, "@");
            strcat(matchable_element, read_id);
            matchable_element[strlen(matchable_element)] = '\0';

            mask = malloc(sizeof (*mask));
            mask->key = matchable_element;
            HASH_ADD_KEYPTR(hh, reads_to_remove, mask->key, strlen(mask->key), mask);
            i++;
        }
        free(line);
    }
}


static void output_stats(char* stats_file) {
    FILE* f = fopen(stats_file, "w");
    
    fprintf(
        f,
        "r1i %s\nr1o %s\nr2i %s\nr2o %s\nr1f %s\nr2f %s\n",
        r1i_path, r1o_path, r2i_path, r2o_path, r1f_path, r2f_path
    );
    fprintf(
        f,
        "read_pairs_checked %i\nread_pairs_removed %i\nread_pairs_remaining %i\n",
        read_pairs_checked, read_pairs_removed, read_pairs_remaining
    );
    
    if (trim_r1) {
        fprintf(f, "trim_r1 %i\n", trim_r1);
    }
    if (trim_r2) {
        fprintf(f, "trim_r2 %i\n", trim_r2);
    }
    if (remove_tiles) {
        fprintf(f, "remove_tiles %s\n", remove_tiles);
    }
    if (remove_reads_path) {
        fprintf(f, "remove_reads %s\n", remove_reads_path);
    }
    
    fclose(f);
}


int main(int argc, char* argv[]) {
    
    /*char* s = malloc(sizeof (char) * 6);
    s = "this\n";
    char* t = malloc(sizeof (char) * 12);
    t = "that other\n";
    char* read_id = strtok(s, " ");
    printf("%s\n", read_id);
    char* read_id2 = strtok(t, " ");
    printf("%s\n", read_id2);
    
    return 0;*/
    
    int arg;
    
    static struct option args[] = {
        {"help", no_argument, 0, 1},
        {"version", no_argument, 0, 2},
        {"quiet", no_argument, 0, 3},
        {"unsafe", no_argument, 0, 4},
        {"stats_file", required_argument, 0, 5},
        {"threshold", required_argument, 0, 6},
        {"remove_tiles", required_argument, 0, 7},
        {"remove_reads", required_argument, 0, 8},
        {"trim_r1", required_argument, 0, 9},
        {"trim_r2", required_argument, 0, 10},
        {"i1", required_argument, 0, 11},
        {"i2", required_argument, 0, 12},
        {"o1", required_argument, 0, 13},
        {"o2", required_argument, 0, 14},
        {"f1", required_argument, 0, 15},
        {"f2", required_argument, 0, 16},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    char* stats_file = NULL;
    criteria = malloc(sizeof (bool(*)(FastqReadPair)) * (ncriteria + 1));
    criteria[ncriteria] = std_check_read;
    
    while ((arg = getopt_long(argc, argv, "", args, &opt_idx)) != -1) {
        switch(arg) {
            case 1:
                printf(USAGE);
                exit(0);
                break;
            case 2:
                printf("%s\n", VERSION);
                exit(0);
                break;
            case 3:
                quiet = true;
                break;
            case 4:
                read_func = readln_unsafe;
                break;
            case 5:
                stats_file = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(stats_file, optarg);
                break;
            case 6:
                threshold = atoi(optarg);
                break;
            case 7:
                remove_tiles = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(remove_tiles, optarg);
                build_remove_tiles();
                ncriteria++;
                criteria = realloc(criteria, sizeof (bool(*)(FastqReadPair)) * (ncriteria + 1));
                criteria[ncriteria] = tile_check_read;
                break;
            case 8:
                remove_reads_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(remove_reads_path, optarg);
                build_remove_reads();
                ncriteria++;
                criteria = realloc(criteria, sizeof (bool(*)(FastqReadPair)) * (ncriteria + 1));
                criteria[ncriteria] = id_check_read;
                break;
            case 9:
                trim_r1 = atoi(optarg);
                include_func_r1 = trim_include_r1;
                break;
            case 10:
                trim_r2 = atoi(optarg);
                include_func_r2 = trim_include_r2;
                break;
            case 11:
                r1i_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r1i_path, optarg);
                break;
            case 12:
                r2i_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r2i_path, optarg);
                break;
            case 13:
                r1o_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r1o_path, optarg);
                break;
            case 14:
                r2o_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r2o_path, optarg);
                break;
            case 15:
                r1f_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r1f_path, optarg);
                break;
            case 16:
                r2f_path = malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(r2f_path, optarg);
                break;
            default:
                exit(1);
        }
    }
    
    if (r1i_path == NULL || r2i_path == NULL || threshold == -1) {
        printf("Missing required arguments: r1i, r2i, threshold\n");
        exit(1);
    }
    
    if (r1o_path == NULL) {
        _log("No o1 argument given - deriving from i1\n");
        r1o_path = build_output_path(r1i_path, "_filtered.fastq");
    }
    if (r2o_path == NULL) {
        _log("No o2 argument given - deriving from i2\n");
        r2o_path = build_output_path(r2i_path, "_filtered.fastq");
    }
    
    if (r1f_path == NULL) {
        _log("No f1 argument given - deriving from i1\n");
        r1f_path = build_output_path(r1i_path, "_filtered_reads.fastq");
    }
    if (r2f_path == NULL) {
        _log("No f2 argument given - deriving from i2\n");
        r2f_path = build_output_path(r2i_path, "_filtered_reads.fastq");
    }

    _log("R1 input: %s\n", r1i_path);
    _log("R2 input: %s\n", r2i_path);
    _log("R1 output: %s\n", r1o_path);
    _log("R2 output: %s\n", r2o_path);
    _log("R1 filtered reads: %s\n", r1f_path);
    _log("R2 filtered reads: %s\n", r2f_path);
    _log("Filter threshold: %i\n", threshold);
    if (trim_r1) {_log("Trimming R1 to %i\n", trim_r1);}
    if (trim_r2) {_log("Trimming R2 to %i\n", trim_r2);}
    if (remove_tiles) {_log("Removing tiles: %s\n", remove_tiles);}
    if (remove_reads_path) {_log("Removing reads in: %s\n", remove_reads_path);}
    _log("Matching %i criteria\n", ncriteria + 1);
    
    int exit_status = filter_fastqs();
    
    _log("Checked %i read pairs, %i removed, %i remaining. Exit status %i\n",
         read_pairs_checked, read_pairs_removed, read_pairs_remaining, exit_status);

    if (stats_file != NULL) {
        _log("Writing stats file %s\n", stats_file);
        output_stats(stats_file);
    }
    
    return exit_status;
}
