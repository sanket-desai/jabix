# Jabix - a tabix based Java classes for handling SAM/BAM, GTF and VCF files
## Version
0.1.0b
## Author
Sanket S Desai
## E-mail
desai\[.\]sanket12\[at\]gmail\[dot\]com
## Description
Jabix is an API initiative developed for handing and rapid extension of data from the most commonly used genomic NGS data formats like BAM, SAM, VCF, GTF, BED. Technically an extension of the TabixReader class created by Heng Li for rapid access of tabix [tabix](http://www.htslib.org/) indexed data files. The main intension of development was to provide a unified API to access all the data and then build applications on top of the same.

In the recent version of pysam [pysam](https://pysam.readthedocs.io/en/latest/api.html), it has integrated to the most awaited format parsers into their API, hence I would recommend use of PYSAM instead for the NGS downstream data handling. Pysam by all means perform most of the functions for which Jabix was planned, hence I do not plan to maintain and continue development of this code API further. However, I welcome use of this code in any application / tool if found suitable for.


Note: The TabixReader class depends on htslib (java) for gzip file access, hence have been provided along with the repository. Please refer HTSLIB [htslib](https://github.com/samtools/htsjdk) for more information and licence statement.

The code is released under MIT license.
