#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub loadBedFile($$);
sub processBamFiles($$$);
sub printTable($);

MAIN:
{
    my $fileBed = shift;
    my @filesBam = @ARGV;
    my @dataBed = ();
    my %table = ();
    my $tic;
    my $toc;

    # parse annotation
    $tic = time();
    loadBedFile(\@dataBed, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(@dataBed), $toc - $tic);

    # process bam files
    $tic = time();
    processBamFiles(\%table, \@filesBam, \@dataBed);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);

    # print results
    printTable(\%table);
}

### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table})
    {
        foreach my $readLength (sort {$a <=> $b} keys %{$table->{$file}})
        {
            my $psiteOffsetBest = 0;
            my $psiteRepeatBest = 0;
            foreach my $psiteOffset (sort {$a <=> $b} keys %{$table->{$file}{$readLength}})
            {
                if ($psiteRepeatBest < $table->{$file}{$readLength}{$psiteOffset})
                {
                    $psiteRepeatBest = $table->{$file}{$readLength}{$psiteOffset};
                    $psiteOffsetBest = $psiteOffset;
                }
            }
            print $file,"\t",$readLength,"\t",$psiteOffsetBest,"\n";
        }
    }
}


### PROCESSBAMFILES
sub processBamFiles($$$)
{
    my $table = $_[0];
    my $filesBam = $_[1];
    my $dataBed = $_[2];

    foreach my $fileBam (@{$filesBam})
    {
        # parse file name
        my $fileName = fileparse($fileBam);
        my $readsUsed = 0;


        # start timer
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();

        # open Bam file
        my $hts = Bio::DB::HTS->new(-bam => $fileBam);

        # loop through bed records
        foreach my $bed (@{$dataBed})
        {
            my $queryStart = ($bed->strand eq "+") ? $bed->thickStart : ($bed->thickEnd - 1);
            my $queryEnd = ($bed->strand eq "+") ? ($bed->thickStart + 1) : $bed->thickEnd;

            my @reads = $hts->get_features_by_location(-seq_id => $bed->chrom,
                                                       -start  => $queryStart,
                                                       -end    => $queryEnd);
            
            next if(scalar(@reads) == 0);
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : $read->end;
                my $readSpan = $read->query->length;
                my $readLinear = $bed->toLinear($readPosition);
                next if($readLinear < 0);
                my $psiteOffset = $bed->txThickStart - $readLinear;

                if ((0 <= $psiteOffset) && ($psiteOffset <= $readSpan))
                {
                    $table->{$fileName}{$readSpan}{$psiteOffset}++;
                    $readsUsed++;
                }

                #print $bed->strand,"\t",$readPosition,"\t",$readSpan,"\t",$readLinear,"\t",$bed->txThickStart,"\t",$psiteOffset,"\n";
            }
            #last if($readsUsed == 100);
        }

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

        #last;
    }

}

### LOADBEDFILES
sub loadBedFile($$)
{
    my $dataBed = $_[0];
    my $fileBed = $_[1];

    open(my $fh, "<", $fileBed) or die $!;
    while(<$fh>)
    {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);
        push(@{$dataBed}, $bed);
    }
    close($fh);
}
