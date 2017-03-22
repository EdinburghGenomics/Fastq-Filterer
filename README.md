# Fastq-Filterer
This is a small C script for filtering out short reads from paired-end fastq files.

This script accepts two (sorted!) paired-end fastq files (R1 and R2), which may be GZ compressed. It then
iterates pairwise through each pair of reads: if R1 and R2 are both longer than the minimum specified, they
are written to corresponding R1/R2 output files.

Output files are always written uncompressed. This script is intended to be used with an external compression
tool such as [pigz](https://github.com/madler/pigz), which is much faster than compressing output on the fly
in a single thread.

There are two functions that can be used to read in input files. The first one uses dynamic memory
reallocation and `strcat`-ing to read in a file line of any length, potentially sacrificing performance. Use
the argument `--unsafe` to use a function which is much simpler and faster, but will chop lines longer than
4096 characters - you have been warned!

The argument `--stats_file <fastq_filterer.stats>` can be used to output a file containing information on the
files input/output and reads checked/filtered.

## Installation
Installation is via the Makefile:

- `make clean` to clean up previous builds
- `make` to compile
- `make check` to run the tests
