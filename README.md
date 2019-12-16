# ribotools
Analysis of Ribosome Footprints in Max-Planck Institute For Brain Research, Frankfurt am Main

## requirements
C library for high-throughput sequencing data formats [HTSLib](https://github.com/samtools/htslib)

C++17 compiler

## installation

```
git clone https://github.com/schumanlab/ribotools
cd ribotools
mkdir build
cd ./build
cmake ..
make
./ribotools --help
```

## tools

```
ribotools
a toolset to profile ana analyse ribosome footprint sequencing
version 0.0.1

usage: ribotools <subcommand> [-h|--help -v|--version -c|--contact]

[subcommands]
    basefreq       calculates base frequency per ORF
    cai            calculates codon adaptation index (CAI)
    codonfreq      calculates codon frequency per ORF
    codonrate      calculates codon decoding rate
    count          counts reads per ORF from BAM files
    depth          coverage per ORF from BAM files
    features       counts reads per gene features from BAM and BED files
    gcratio        calculates GC-ratio in reads from BAM file
    gcref          calculates GC-ratio in ORF from BED and FASTA files
    gtftobed       converts GTF file format to BED12 file format
    irate          calculates initiation rate based on BAM files from RNASeq and RiboSeq
    length         calculates reads length based on CIGAR string from BAM file
    metagene       project footprints to a metagene histogram
    mtdr           calculates mean transcript decoding rate (MTDR)
    pausing        calculates z-score pausing score per codon
    periodicity    creates P-site periodicity histogram
    poffset        calculates P-site offset per read length
    runoff         project runoff footprints normalised to transcript end
    structure      calculates A-site coverage over protein secondary structure provided by FASTA file
    translate      translates transcripts based on BED and FASTA files
    uorfs          screens for upstream open reading frame
    utrseq         exports UTR sequence from BED and FASTA files

[options]
    -h, --help     print this help message
    -v, --version  what major.minor.build version of ribotools is used
    -c, --contact  feature requests, bugs, mailing lists, etc.

```





