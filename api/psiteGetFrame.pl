#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
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
    processBamFiles(\%table, \@dataBed, \@filesBam);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);

    printTable(\%table);
}


### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table}) {
        foreach my $span (sort {$a <=> $b} keys %{$table->{$file}}) {
            my $frame0 = exists($table->{$file}{$span}{"frame"}{0}) ? $table->{$file}{$span}{"frame"}{0} : 0;
            my $frame1 = exists($table->{$file}{$span}{"frame"}{1}) ? $table->{$file}{$span}{"frame"}{1} : 0;
            my $frame2 = exists($table->{$file}{$span}{"frame"}{2}) ? $table->{$file}{$span}{"frame"}{2} : 0;
            my $total = $frame0 + $frame1 + $frame2;
            my $count = $table->{$file}{$span}{"hist"};

            print $file,"\t",$span,"\t",$count,"\t",$frame0/$total,"\t",$frame1/$total,"\t",$frame2/$total,"\n";
        }
        
    }
}


### PROCESSBAMFILES
sub processBamFiles($$$)
{
    my $table = $_[0];
    my $dataBed = $_[1];
    my $filesBam = $_[2];
    
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
                                                       -start  => $bed->chromStart,
                                                       -end    => $bed->chromEnd);
            my $readsCount = scalar(@reads);
            next if ($readsCount < 1);
            
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : ($read->end - 1);
                my $readSpan = $read->query->length;
                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                my $offset = $bed->txThickStart - $readLinear;
                my $frame = abs($offset) % 3;
                $table->{$fileName}{$readSpan}{"frame"}{$frame}++;
                $table->{$fileName}{$readSpan}{"hist"}++;
                $readsUsed++;
            }
            
        }

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);
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
        push(@{$dataBed}, $bed);
    }
    close($fh);
}
