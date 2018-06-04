
0.4 (2018-06-04)
----------------
- Added `--remove_reads` option
- Redirecting filtered reads to `--f1` and `--f2` files
- Added uthash.h
- Test improvements

0.3.1 (2017-05-23)
------------------
- Fixed memory allocation bug in tile filtering

0.3 (2017-05-16)
----------------
- Added `--remove_tiles`, which filters out reads with specific tile ids (Illumina-formatted headers only)
- Added `--trim_r1` and `--trim_r2` for trimming reads down to a maximum length

0.2 (2017-03-22)
----------------
- Added `--stats_file`, which writes a record of read pairs checked and filtered
- Various refactoring and improvements, especially stdout

0.1 (2016-10-11)
----------------
- Initial release
