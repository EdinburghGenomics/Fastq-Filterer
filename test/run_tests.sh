#!/usr/bin/env bash

# Currently, this test script simply runs the filterer and compares the output files with expected output
# Ideally we should use something like check.h, but just implementing this for now.
# M.

scriptpath=$(dirname $0)
cd $scriptpath
r1o=R1_filtered.fastq
r2o=R2_filtered.fastq
r1f=R1_filtered_reads.fastq
r2f=R2_filtered_reads.fastq
filterer="../fastq_filterer --quiet --o1 $r1o --o2 $r2o --f1 $r1f --f2 $r2f --threshold 9"
exit_status=0

function compare {
    echo "Comparing $1 and $2"
    diff $1 $2
    x=$?
    echo "Exit status: $x"
    exit_status=$[$exit_status+$x]
    rm $1
}

function check_outputs {
    compare $r1o expected_outputs/${1}R1_filtered.fastq
    compare $r2o expected_outputs/${1}R2_filtered.fastq
    compare $r1f expected_outputs/${1}R1_filtered_reads.fastq
    compare $r2f expected_outputs/${1}R2_filtered_reads.fastq
    echo "______________________"
}


echo "Testing non-compressed"
$filterer --i1 inputs/R1.fastq --i2 inputs/R2.fastq
check_outputs


echo "Testing compressed"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz
check_outputs


echo "Testing unsafe mode"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz --unsafe
check_outputs


echo "Testing implicit output paths"
../fastq_filterer --quiet --threshold 9 --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz
mv inputs/R?_filtered*.fastq .
check_outputs


echo "Testing output stats"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz --stats_file inputs/fastq_filterer.stats
compare inputs/fastq_filterer.stats expected_outputs/fastq_filterer.stats
check_outputs


echo "Testing tile removal"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz --remove_tiles 1102,2202 --stats_file inputs/fastq_filterer.stats
compare inputs/fastq_filterer.stats expected_outputs/rm_tiles.stats
check_outputs rm_tiles_


echo "Testing read removal"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz --remove_reads inputs/rm_reads.txt --stats_file inputs/fastq_filterer.stats
compare inputs/fastq_filterer.stats expected_outputs/rm_reads.stats
check_outputs rm_reads_


echo "Testing read trimming"
$filterer --i1 inputs/R1.fastq.gz --i2 inputs/R2.fastq.gz --trim_r1 14 --trim_r2 16 --stats_file inputs/fastq_filterer.stats
compare inputs/fastq_filterer.stats expected_outputs/trim_reads.stats
check_outputs trim_reads_

echo "Finished tests with exit status $exit_status"
