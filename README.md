# Fastq-Filterer
This is a small C script for filtering out short reads from a pair of fastq files.

This script accepts two sorted paired-end fastq files (R1 and R2), which may be GZ compressed. It then
iterates pairwise through each pair of reads: if R1 and R2 are both longer than the minimum specified, they
are written to corresponding R1/R2 output files.

Output files are always written uncompressed. This script is intended to be used with an external compression
tool such as [pigz](https://github.com/madler/pigz), which is much faster than compressing output on the fly
in a single thread.

This script implements a function that can read in a file line of any length. This involves dynamic memory
reallocation and `strcat`-ing, which may sacrifice performance. There is an optional argument, `--unsafe`,
which uses a much simpler read function, but will chop lines longer than 4096 characters - you have been
warned!
