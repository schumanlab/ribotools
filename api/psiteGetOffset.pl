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
    my @dataBed = ();
    my %table = ();
    my %offsets = ();
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

    # print results
    printTable(\%table);
}

### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table}) {

        foreach my $readLength (sort {$a <=> $b} keys %{$table->{$file}}) {
            my @values = ();
            my $total = 0;
            my $count = $table->{$file}{$readLength}{"counts"};
            foreach my $frame (sort {$a <=> $b} keys %{$table->{$file}{$readLength}{"frame"}}) {
                my $val = $table->{$file}{$readLength}{"frame"}{$frame};
                push(@values, $val);
                $total += $val;
            }

            print $file,"\t",$readLength,"\t",$count,"\t",$values[0]/$total,"\t",$values[1]/$total,"\t",$values[2]/$total,"\n";

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
            my @reads = $hts->get_features_by_location(-seq_id => $bed->chrom,
                                                       -start  => $bed->thickStart,
                                                       -end    => $bed->thickEnd);
            
            my $readsCount = scalar(@reads);
            next if ($readsCount == 0);
            next if ($bed->lengthThick == 0);

            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : ($read->end - 1);
                my $readSpan = $read->query->length;

                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                next if (!exists($offsets->{$fileName}{$readSpan}));
                my $psiteOffset = $offsets->{$fileName}{$readSpan}; # add one to get the position after the offset
                $readsUsed++;

                # relative offsets
                my $relOffsetStart = $readLinear + $psiteOffset - $bed->txThickStart - 1;
                my $relOffsetCenter = $readLinear + $psiteOffset - $bed->txThickCenter - 1;
                my $relOffsetEnd = $readLinear + $psiteOffset - $bed->txThickEnd - 1;

                # accumulate frame
                my $frameStart = $relOffsetStart % 3;
                my $frameCenter = $relOffsetCenter % 3;
                my $frameEnd = $relOffsetEnd % 3;

                $table->{$fileName}{$readSpan}{"counts"}++;
                $table->{$fileName}{$readSpan}{"frame"}{$frameStart} += 1/3;
                $table->{$fileName}{$readSpan}{"frame"}{$frameCenter} += 1/3;
                $table->{$fileName}{$readSpan}{"frame"}{$frameEnd} += 1/3;
                
                $readsUsed++;
                
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
        #next if (($bed->lengthThick < 300) || (5000 < $bed->lengthThick));
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