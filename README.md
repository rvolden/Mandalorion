# Mandalorion #
[![Github release](https://img.shields.io/github/tag/rvolden/Mandalorion.svg?label=Version)](https://github.com/rvolden/Mandalorion/tags)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

Episode III.5: Rogue Isoforms

Takes R2C2/C3POa data and defines high confidence isoforms.

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [medaka](https://github.com/nanoporetech/medaka)
- [abPOA](https://github.com/yangao07/abPOA)

The paths to these will need to be put into your config file [like this.](example_config) If you have the program installed or in your path already, replace the path with the name of the program.

## Usage ##
```bash
python3 Mandalorion.py [OPTIONS]
```

Running with default settings:
```bash
python3 Mandalorion.py -c config_file -p . -g gencodeV29.gtf -G hg38.fasta -a adapters.fasta -f R2C2_consensi.fasta -b R2C2_subreads.fastq
```

Required options:
```
-c  config file containing paths to required dependencies (above)
-p  output path
-g  annotation file (gtf)
-G  genome file (fasta)
-a  adapter file (fasta)
-f  R2C2 read file (fasta, can be comma separated list)
-b  R2C2 subread file (fastq, can be comma separated list)
```

Tweakable parameters:
```
-u  upstream buffer, defines leniency window for TSS and polyA definition (default 10)
-d  downstream buffer, defines leniency window for TSS and polyA definition (default 50)
-s  subsample consensus, defines how many random subreads are used to make isoforms (default 500)
-r  minimum ratio, proportion of reads that align to a locus required for an isoform (default 0.05)
-i  minimum internal ratio (default 0.125)
-R  minimum number of reads for an isoform (default 5)
-O  overhangs, defines bounds for unaligned bases on ends, format: min5',max5',min3',max3' (default 0,40,0,40)
-t  number of threads to use for minimap2 (default 4)
-e  ends, the ends of your sequences, format: 5prime,3prime (default 'ATGGG,AAAAA')
-I  minimum isoform length (default 500)
-n  minimum feature count (number of reads, default 2)
-w  splice site window (default 1)
-A  A proportion cutoff - isoforms with A content within a 30 nt window around the polyA higher than this will be discarded (default 0.5)
```

## Utils ##
These are the scripts used in the paper to do haplotype phasing and HLA analysis.
