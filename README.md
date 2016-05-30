# fq-plots
Estimates and displays read depth and insert size plots

##Author

Jason Kwong (@kwongjc)

##Dependencies
* Python 2.7.x
* Pandas
* Samtools

##Usage

```
$ fq-plots.py -h
usage: 
  fq-plots.py BAMFILE

Estimates and plots insert sizes and read depth coverage of paired-end reads

positional arguments:
  BAMFILE              BAM file eg. snps.bam from Snippy

optional arguments:
  -h, --help           show this help message and exit
  --plot insert|depth  Plot insert sizes ("--plot insert") or read depth ("--plot depth")
  --centile %          Percentile filter for inserts eg. 95 = 95% most frequent insert sizes
  --locus LOCUS        Locus to display depth plots
  --coords START:END   Locus coordinates to display depth plots
  --interval LEN       Interval to draw depth plots (default=1000)
  --version            show program's version number and exit
```

##Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/fq-plots/issues).  

##Software Licence

[GPLv3](https://github.com/kwongj/fq-plots/blob/master/LICENCE)

