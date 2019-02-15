#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub loadPsiteTable($$);
sub loadBedFile($$);
sub processBamFiles($$$$);
sub printTable($);

MAIN:
{
    my $filePsite = shift;
    my $fileBed = shift;
    my @filesBam = @ARGV;
    my %offsets = ();
    my @dataBed = ();
    my %table = ();
    my $tic;
    my $toc;

    # parse psite table
    $tic = time();
    loadPsiteTable(\%offsets, $filePsite);
    $toc = time();
    printf(STDERR "Parsed %d rows in %.4f sec.\n", scalar(keys %offsets), $toc - $tic);

    # parse annotation
    $tic = time();
    loadBedFile(\@dataBed, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(@dataBed), $toc - $tic);

    # process bam files
    $tic = time();
    processBamFiles(\%table, \@filesBam, \%offsets, \@dataBed);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);

    printTable(\%table);
}


### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table}) {
        foreach my $label (sort keys %{$table->{$file}}) {
            foreach my $span (sort {$a <=> $b} keys %{$table->{$file}{$label}}) {

                next if (!exists($table->{$file}{$label}{$span}{"offset"}));
                next if (!exists($table->{$file}{$label}{$span}{"frame"}));
                my $frame0 = exists($table->{$file}{$label}{$span}{"frame"}{0}) ? $table->{$file}{$label}{$span}{"frame"}{0} : 0; 
                my $frame1 = exists($table->{$file}{$label}{$span}{"frame"}{1}) ? $table->{$file}{$label}{$span}{"frame"}{1} : 0;
                my $frame2 = exists($table->{$file}{$label}{$span}{"frame"}{2}) ? $table->{$file}{$label}{$span}{"frame"}{2} : 0;

                print $file,"\t",
                $label,"\t",
                $span,"\t",
                $table->{$file}{$label}{$span}{"offset"},"\t",
                $frame0,"\t",
                $frame1,"\t",
                $frame2,"\n";
            }
        }
    }
}


### PROCESSBAMFILES
sub processBamFiles($$$$)
{
    my $table = $_[0];
    my $filesBam = $_[1];
    my $offsets = $_[2];
    my $dataBed = $_[3];

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
            next if ($readsCount < 1);
            next if ($bed->strand eq "-");
            
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : ($read->end - 1);
                my $readSpan = $read->query->length;
                my $readLinear = $bed->toLinear($readPosition);

                next if ($readLinear < 0);
                next if (!exists($offsets->{$fileName}{$readSpan}));
                my $psiteOffset = $offsets->{$fileName}{$readSpan} + 1; # add one to get the position after the offset
                $readsUsed++;

                # frame
                my $frame = abs($readLinear + $psiteOffset - $bed->txThickStart + 1) % 3;
                my $label = "CDS";
                $label="5pUTR" if ($readLinear < $bed->txThickStart);
                $label="3pUTR" if ($readLinear > $bed->txThickEnd);
                $table->{$fileName}{$label}{$readSpan}{"frame"}{$frame}++;
                $table->{$fileName}{$label}{$readSpan}{"offset"} = $psiteOffset;

            }
            
        }

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
    my %genes = ();

    open(my $fh, "<", $fileBed) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);
        my ($transcript, $gene) = split(";", $bed->name, 2);
        next if (exists($genes{$gene}));
        $genes{$gene}++;
        push(@{$dataBed}, $bed);
    }
    close($fh);
}


### LOADPSITETABLE
sub loadPsiteTable($$)
{
    my $offsets = $_[0];
    my $filePsite = $_[1];

    open(my $fh, "<", $filePsite) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my ($fileName, $readLength, $offset, $score) = split("\t", $_, 4);
        $offsets->{$fileName}{$readLength} = $offset;
    }
    close($fh);
}