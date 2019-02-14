#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub loadBedFile($$);
sub processBamFiles($$$$$);
sub printTable($);

MAIN:
{
    my $readLength = shift;
    my $psiteOffset = shift;
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
    processBamFiles(\%table, \@filesBam, \@dataBed, $readLength, $psiteOffset);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);

    printTable(\%table);
}


### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table})
    {
        print $file,"\t",join(",", @{$table->{$file}[0]}),"\t",join(",", @{$table->{$file}[1]}),"\n";
    }
}


### PROCESSBAMFILES
sub processBamFiles($$$$$)
{
    my $table = $_[0];
    my $filesBam = $_[1];
    my $dataBed = $_[2];
    my $readLength = $_[3];
    my $psiteOffset = $_[4];

    foreach my $fileBam (@{$filesBam})
    {
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
            my @reads = $hts->get_features_by_location(-seq_id => $bed->chrom,
                                                       -start  => $bed->thickStart,
                                                       -end    => $bed->thickEnd);
            my $readsCount = scalar(@reads);
            
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : $read->end;
                my $readSpan = $read->query->length;
                next if($readSpan != $readLength);

                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                $readsUsed++;

                # relative offsets
                my $frame = ($readLinear + $psiteOffset - $bed->txThickStart) % 3;
                
            }
            last;
        }

        $table->{$fileName} = [\@coverage, \@counts];

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

        last;
    }

}


### LOADBEDFILE
sub loadBedFile($$)
{
    my $dataBed = $_[0];
    my $fileBed = $_[1];

    open(my $fh, "<", $fileBed) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);
        next if (($bed->lengthThick < 300) || (5000 < $bed->lengthThick));
        push(@{$dataBed}, $bed);
    }
    close($fh);
}

