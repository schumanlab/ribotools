#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;
use Set::IntervalTree;
use Bed12;

sub loadBedTree($$);
sub processBams($$$);

MAIN:
{
    my $fileBed = shift;
    my @filesBam = @ARGV;
    my %bedTree = ();
    my %table = ();
    my $tic; 
    my $toc;

    # parse annotation
    $tic = time();
    #loadBedTree(\%bedTree, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(keys %bedTree), $toc - $tic);

    # process bam files
    $tic = time();
    processBams(\%table, \%bedTree, \@filesBam);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);


}


### PROCESSBAMS
sub processBams($$$)
{
    my $table = $_[0];
    my $bedTree = $_[1];
    my $filesBam = $_[2];

    foreach my $fileBam (@{$filesBam})
    {
        # parse file name
        my $fileName = fileparse($fileBam);
        my $readsUsed = 0;

        # start timer
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();

        # open Bam file
        my $hts = Bio::DB::HTSfile->open($fileBam);
        my $htsHeader = $hts->header_read;
        my $htsTargetNames = $htsHeader->target_name;
        while (my $read = $hts->read1($htsHeader)) {

            next if ($read->qual != 255);
            my $chrom = $htsTargetNames->[$read->tid];
            my $readStart = $read->pos;
            my $readEnd = $read->calend;
            my $readSpan = $readEnd - $readStart;
            my $readLength = $read->query->length;
            my $cigar = $read->cigar_str;

            print $readStart,"-",$readEnd,"\t",$readSpan,"/",$readLength,"\t",$cigar,"\n";


            #next if (!exists($bedTree->{$chrom}));
            #my $bed = $bedTree->{$chrom}->fetch($readStart, $readEnd);
            
            #print $chrom,"\t",$readStart,"\t",$readEnd,"\t",$readLength,"\n";


        }


        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

        #last;
    }
}


### LOADBEDTREE
sub loadBedTree($$)
{
    my $bedTree = $_[0];
    my $fileBed = $_[1];

    open(my $fh, "<", $fileBed) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);

        # allocate new tree
        if (!exists($bedTree->{$bed->chrom})) {
            $bedTree->{$bed->chrom} = Set::IntervalTree->new;
        }

        $bedTree->{$bed->chrom}->insert($bed, $bed->chromStart, $bed->chromEnd);
    }
    close($fh);
}
