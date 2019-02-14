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
        print $file,"\t",join(",", @{$table->{$file}[0]}),"\t",join(",", @{$table->{$file}[1]}),"\n";
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
        my @coverage = (0) x 300;
        my @counts = (0) x 300;
        
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
            my $depth = 1;#($readsCount / $bed->lengthThick);
            next if ($readsCount < 128);
            
            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : $read->end;
                my $readSpan = $read->query->length;

                #next if($readSpan != 29);

                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                next if (!exists($offsets->{$fileName}{$readSpan}));
                my $psiteOffset = $offsets->{$fileName}{$readSpan};
                $readsUsed++;

                # relative offsets
                my $relOffsetStart = $readLinear + $psiteOffset - $bed->txThickStart;
                my $relOffsetCenter = $readLinear + $psiteOffset - $bed->txThickCenter;
                my $relOffsetEnd = $readLinear + $psiteOffset - $bed->txThickEnd;

                if ((-25 <= $relOffsetStart) && ($relOffsetStart <= 75))
                {
                    $relOffsetStart += 25;
                    $coverage[$relOffsetStart] += (1/$depth);
                    $counts[$relOffsetStart] += (1/$readsCount);
                }

                if ((-50 <= $relOffsetCenter) && ($relOffsetCenter <= 50))
                {
                    $relOffsetCenter += 151;
                    $coverage[$relOffsetCenter] += (1/$depth);
                    $counts[$relOffsetCenter] += (1/$readsCount);
                }

                if ((-75 <= $relOffsetEnd) && ($relOffsetEnd <= 25))
                {
                    $relOffsetEnd += 277;
                    $coverage[$relOffsetEnd] += (1/$depth);
                    $counts[$relOffsetEnd] += (1/$readsCount);
                }

            }
            #last;
        }

        $table->{$fileName} = [\@coverage, \@counts];

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
        next if (($bed->lengthThick < 300) || (5000 < $bed->lengthThick));
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