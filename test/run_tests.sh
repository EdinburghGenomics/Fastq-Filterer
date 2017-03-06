#!/usr/bin/env bash

# Currently, this test script simply runs the filterer and compares the output files with expected output
# Ideally we should use something like check.h, but just implementing this for now.
# M.

scriptpath=$(dirname $0)
filterer=$scriptpath/../fastq_filterer
exit_status=0

function compare_files {
    echo "Comparing $1 and $2"
    diff $1 $2
    ex_st=$?
    echo "Finished comparison with exit status $ex_st"
    exit_status=$[$exit_status+$ex_st]
}


echo "Testing non-compressed"
$filterer --i1 $scriptpath/uncompressed_R1.fastq --i2 $scriptpath/uncompressed_R2.fastq --o1 $scriptpath/uncompressed_R1_filtered.fastq --o2 $scriptpath/uncompressed_R2_filtered.fastq --threshold 9
echo
compare_files $scriptpath/uncompressed_R1_filtered.fastq $scriptpath/R1_min_len_9.fastq
echo
compare_files $scriptpath/uncompressed_R2_filtered.fastq $scriptpath/R2_min_len_9.fastq
echo "______________________"
rm $scriptpath/*filtered.fastq

echo "Testing compressed"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $scriptpath/compressed_R1_filtered.fastq --o2 $scriptpath/compressed_R2_filtered.fastq --threshold 9
echo
compare_files $scriptpath/compressed_R1_filtered.fastq $scriptpath/R1_min_len_9.fastq
echo
compare_files $scriptpath/compressed_R2_filtered.fastq $scriptpath/R2_min_len_9.fastq
echo "______________________"
rm $scriptpath/*filtered.fastq

echo "Testing unsafe mode"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --o1 $scriptpath/compressed_R1_filtered.fastq --o2 $scriptpath/compressed_R2_filtered.fastq --threshold 9 --unsafe
echo
compare_files $scriptpath/compressed_R1_filtered.fastq $scriptpath/R1_min_len_9.fastq
echo
compare_files $scriptpath/compressed_R2_filtered.fastq $scriptpath/R2_min_len_9.fastq
echo "______________________"
rm $scriptpath/*filtered.fastq

echo "Testing implicit output paths and output stats"
$filterer --i1 $scriptpath/compressed_R1.fastq.gz --i2 $scriptpath/compressed_R2.fastq.gz --threshold 9 --stats_file $scriptpath/fastq_filterer.stats
echo
compare_files $scriptpath/compressed_R1_filtered.fastq $scriptpath/R1_min_len_9.fastq
echo
compare_files $scriptpath/compressed_R2_filtered.fastq $scriptpath/R2_min_len_9.fastq
echo
compare_files $scriptpath/fastq_filterer.stats $scriptpath/expected_fastq_filterer.stats
echo "______________________"
rm $scriptpath/*filtered.fastq

echo "Finished tests with exit status $exit_status"
