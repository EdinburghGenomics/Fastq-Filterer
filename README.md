# Fastq-Filterer

This is a C script for filtering out short reads from paired-end fastq files.

The inputs are two (sorted!) paired-end fastqs (R1 and R2), which may be `fastq` or `fastq.gz` files. The
filterer then iterates pairwise through each read pair: if R1 and R2 are both longer than the minimum
specified, they are written to corresponding R1/R2 output files.

Output files will always be written uncompressed. Fastq-Filterer is intended to be used with
[pigz](https://github.com/madler/pigz) or similar multi-threaded compression tool, which is much faster than
compressing output on the fly in a single thread.

The filterer has two modes for reading in input files. By default, it uses dynamic memory rellocation and
`strcat`-ing to read in a file line of any length. In the second mode, a simpler reading function is used,
which is faster, but will also chop any lines longer than 4096 characters - you have been warned!

A file can also be output containing summary information on the input/output files and reads checked and
filtered.

Hash tables are implemented in this project via [uthash.h](https://github.com/troydhanson/uthash), an
unmodified copy of which is included in `src`.


## Installation
To set up the filterer, simply compile it in place via the Makefile:

- `make clean` to clean up previous builds
- `make` to compile
- `make check` to run the tests


## Usage
A minimum of three arguments are required:

    fastq_filterer --i1 <r1.fastq> --i2 <r2.fastq> --threshold <filter_threshold>

Running this will read in both files, apply the filter threshold, and output two files named after the input
files with the suffix '\_filtered.fastq'.

Other arguments can also be passed:
- `--o1 <r1_out.fastq>`: custom name for the R1 output file
- `--o2 <r2_out.fastq>`: custom name for the R2 output file
- `--stats_file <stats_file>`: write a file summarising the read pairs checked and removed
- `--unsafe`: use a simpler, faster but less safe read function
- `--remove_tiles <tile1,tile2,tile3...>`: comma-separated list of tile ids to remove regardless of length
- `--trim_r1 <max_len>`: trim all reads for r1.fastq to a maximum length
- `--trim_r2 <max_len>`: as above for r2.fastq


## Input files
A few assumptions are made about the input fastqs:
- It is assumed that both input fastqs have the same number of reads, and that they are both in the same
order. Therefore, a read starting on line `n` in r1.fastq should correspond to the read starting on line `n`
in r2.fastq
- If using `--remove_tiles`, Fastq-Filterer parses the flowcell tile ID from the fastq read headers, so it is
assumed that the read headers are in standard Illumina format:

    @instrument_id:run_id:flowcell_id:lane:tile_id:x:y read_number:filter_flag:0:idx_seq

For more information, see Illumina's bcl2fastq docs.
