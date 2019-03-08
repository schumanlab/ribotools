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

    foreach my $file (sort keys %{$table})
    {
        print $file,"\t",
              join(",", @{$table->{$file}[0]}),"\t",
              join(",", @{$table->{$file}[1]}),"\t",
              join(",", @{$table->{$file}[2]}),"\t",
              join(",", @{$table->{$file}[3]}),"\n";

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
        my @frame = (0) x 3;
        my @coverageStart = (0) x 101;
        my @coverageCenter = (0) x 101;
        my @coverageEnd = (0) x 101;
        
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
            
            my $depthScore = 1/($readsCount /$bed->lengthThick);
            
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : ($read->end - 1);
                my $readSpan = $read->query->length;
                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                my $frameTest = abs($bed->txThickStart - $readLinear) % 3;
                my $frame = 0 if($frameTest == 0);
                $frame = 1 if($frameTest == 2);
                $frame = 2 if($frameTest == 1);

                # debug
                #print $readSpan,"\t",$readLinear,"\t",$psiteOffset,"\t",$relOffsetStart,"\t",$relOffsetCenter,"\t",$relOffsetEnd,"\n";


                next if(!exists($offsets->{$fileName}{$readSpan}));
                my $psiteOffset =  exists($offsets->{$fileName}{$readSpan}{$frame}) ? $offsets->{$fileName}{$readSpan}{$frame} : 0; # add one to get the position after the offset
                $readsUsed++;

                # relative offsets
                my $relOffsetStart = $readLinear + $psiteOffset - $bed->txThickStart;
                my $relOffsetCenter = $readLinear + $psiteOffset - $bed->txThickCenter;
                my $relOffsetEnd = $readLinear + $psiteOffset - $bed->txThickEnd;

                
                # accumulate position
                if ((-25 <= $relOffsetStart) && ($relOffsetStart <= 75)) {
                    $coverageStart[$relOffsetStart + 25] += $depthScore;
                }

                if ((-50 <= $relOffsetCenter) && ($relOffsetCenter <= 50)) {
                    $coverageCenter[$relOffsetCenter + 50] += $depthScore;
                }

                if ((-75 <= $relOffsetEnd) && ($relOffsetEnd <= 25)) {
                    $coverageEnd[$relOffsetEnd + 75] += $depthScore;
                }

                # accumulate frame
                my $frameStart = $relOffsetStart % 3;
                #my $frameCenter = $relOffsetCenter % 3;
                #my $frameEnd = $relOffsetEnd % 3;
                $frame[$frameStart] += $depthScore;
                #$frame[$frameCenter] += $depthScore;
                #$frame[$frameEnd] += $depthScore;

            }
            #last;
        }

        $table->{$fileName} = [\@frame, \@coverageStart, \@coverageCenter, \@coverageEnd];

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

        #last;
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


### LOADPSITETABLE
sub loadPsiteTable($$)
{
    my $offsets = $_[0];
    my $filePsite = $_[1];

    open(my $fh, "<", $filePsite) or die $!;
    while (<$fh>)
    {
        chomp($_);
        
        next if($_ =~ m/^#/);
        my ($fileName, $readLength, $frame, $offset, $score) = split("\t", $_, 5);
        $offsets->{$fileName}{$readLength}{$frame} = $offset;
    }
    close($fh);
}