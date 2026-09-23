In order to generate the .cmp1 files QUILTS uses:

Contemporary Python 1 script instead of 2
-----------------------------------------
1. Download *chr*.fa.gz files from ftp://ftp.ensembl.org/pub/current_fasta/[your organism]/dna/
for human GRCh37, ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
for Gencode 42 ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/

2. Run “python gunzip_compress_dir.py --dir <path to folder with your .fa.gz files>”

So, if I had my files in references/genome/ensembl_hg19_GRCh37, I would run:

Z:\Projects\QUILTS>python pyQUILTS_KRC/prepare_genome/gunzip_compress_dir.py --dir references/genome/ensembl_hg19_GRCh37

Classic PERL
---------------

1. Download *chr*.fa.gz files from ftp://ftp.ensembl.org/pub/current_fasta/[your organism]/dna/ (for human GRCh37, ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/)
2. Run “perl gunzip_dir.pl <path to folder with your .fa.gz files>”, or simply “gunzip *.fa.gz”
3. Run “perl compress_chr_dir.pl <path to folder with your .fa files>”

So, if I had my files in /dir/hg38/, I would do:

perl gunzip_dir.pl /dir/hg38/
perl compress_chr_dir.pl /dir/hg38/


If you get an error saying “[…]/compress_chr: cannot execute binary file”, try recompiling compress_chr.c:

gcc -o compress_chr compress_chr.c

and then try running compress_chr_dir.pl again.