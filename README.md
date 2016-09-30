# updownstream
Grab X bases upstream &amp; Y bases downstream based on PATRIC Genomic Features

## Usage
```
updownstream.py -i Genomicfeatures.tab -f Contigs.fna -u 100 -d 60 -o contigs
```

### Flags
#### `-i`
The Genomic Features file
#### `-f`
The Contigs fasta
#### `-u`
How many bases to grab upstream
#### `-d`
How many bases to grab downstream
#### `-o`
The **basename** for the output files


### Output
The script will write 8 files out. Hits are defined as
having a start codon found at the start position listed in the genomic features file.

* `basename_pos_up_hit.fa`
    * Positive-strand upstream hits
* `basename_pos_up_miss.fa`
    * Positive-strand upstream misses
* `basename_pos_down_hit.fa`
    * Positive-strand downstream hits
* `basename_pos_down_miss.fa`
    * Positive-strand downstream misses
* `basename_neg_up_hit.fa`
    * Negative-strand upstream hits
* `basename_neg_up_miss.fa`
    * Negative-strand upstream misses
* `basename_neg_down_hit.fa`
    * Negative-strand downstream hits
* `basename_neg_down_miss.fa`
    * Negative-strand downstream misses
