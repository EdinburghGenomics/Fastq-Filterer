#!/usr/bin/env bash

# Currently, this test script simply runs the filterer and compares the output files with expected output
# Ideally we should use something like check.h, but just implementing this for now.
# M.

scriptpath=$(dirname $0)
filterer="$scriptpath/../fastq_filterer --quiet"
exit_status=0

r1o=$scriptpath/R1_filtered.fastq
r2o=$scriptpath/R2_filtered.fastq

function _compare {
    diff $1 $2
    x=$?
    echo "Compared $1 and $2 with exit status $x"
    exit_status=$[$exit_status+$x]
    rm $1
}

function check_outputs {
    _compare $1 $scriptpath/R1_min_len_9.fastq
    _compare $2 $scriptpath/R2_min_len_9.fastq
    echo "______________________"
}


echo "Testing non-compressed"
$filterer --i1 $scriptpath/uncompressed_R1.fastq --i2 $scriptpath/uncompressed_R2.fastq --o1 $r1o --o2 $r2o --threshold 9
check_outputs $r1o $r2o


echo "Testing compressed"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $r1o --o2 $r2o --threshold 9
check_outputs $r1o $r2o


echo "Testing unsafe mode"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $r1o --o2 $r2o --threshold 9 --unsafe
check_outputs $r1o $r2o


echo "Testing implicit output paths"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --threshold 9
check_outputs $scriptpath/compressed_R1_filtered.fastq $scriptpath/compressed_R2_filtered.fastq


echo "Testing output stats"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $r1o --o2 $r2o --threshold 9 --stats_file $scriptpath/fastq_filterer.stats
echo
_compare $scriptpath/fastq_filterer.stats $scriptpath/expected_fastq_filterer.stats
check_outputs $r1o $r2o


echo "Testing tile removal"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $r1o --o2 $r2o --threshold 9 --remove_tiles 1102,2202
_compare $r1o $scriptpath/R1_remove_tiles.fastq
_compare $r2o $scriptpath/R2_remove_tiles.fastq
echo "______________________"


echo "Testing read trimming"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $r1o --o2 $r2o --threshold 9 --trim_r1 14 --trim_r2 16 --stats_file $scriptpath/fastq_filterer.stats
_compare $r1o $scriptpath/R1_trim_14.fastq
_compare $r2o $scriptpath/R2_trim_16.fastq
_compare $scriptpath/fastq_filterer.stats $scriptpath/trim_reads.stats
echo "______________________"

echo "Finished tests with exit status $exit_status"
