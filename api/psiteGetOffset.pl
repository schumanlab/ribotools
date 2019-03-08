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
    my %offsets = ();
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

    # print results
    printTable(\%table);
}

### PRINTTABLE
sub printTable($)
{
    my $table = $_[0];

    foreach my $file (sort keys %{$table}) {

        foreach my $readLength (sort {$a <=> $b} keys %{$table->{$file}}) {

            my %result = (
                0 => [0, 0],
                1 => [0, 0],
                2 => [0, 0]
            );
            foreach my $offset (sort {$a <=> $b} keys %{$table->{$file}{$readLength}}) {
                my $frame = ($offset % 3);
                my $counts = $table->{$file}{$readLength}{$offset};
                
                if ($result{$frame}[1] <= $counts) {
                    $result{$frame}[1] = $counts;
                    $result{$frame}[0] = $offset;
                }

                print $file,"\t",$readLength,"\t",$frame,"\t",$offset,"\t",$counts,"\n";
            }

            #print $file,"\t",$readLength,"\t","0","\t",$result{0}[0],"\t",$result{0}[1],"\n";
            #print $file,"\t",$readLength,"\t","1","\t",$result{1}[0],"\t",$result{1}[1],"\n";
            #print $file,"\t",$readLength,"\t","2","\t",$result{2}[0],"\t",$result{2}[1],"\n";
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
                                                       -start  => $bed->chromStart,
                                                       -end    => $bed->chromEnd);
            
            my $readsCount = scalar(@reads);
            next if ($readsCount == 0);
            next if ($bed->lengthThick == 0);

            foreach my $read (@reads)
            {
                my $readPosition = ($bed->strand eq "+") ? $read->start : ($read->end - 1);
                my $readSpan = $read->query->length;
                my $readLinear = $bed->toLinear($readPosition);
                next if ($readLinear < 0);
                my $shift = ($bed->txThickStart - $readLinear);
                my $steps = int($shift / $readSpan);
                my $delta = $readLinear + $steps * $readSpan;
                my $offset = abs($bed->txThickStart - $delta);
                
                
                #print $readLinear,"\t",$readSpan,"\t",$bed->txThickStart,"\t",$shift,"\t",$steps,"\t",$delta,"\t",$offset,"\n";

                if ((0 <= $offset) && ($offset <= $readSpan)) {
                    $table->{$fileName}{$readSpan}{$offset}++;
                    $readsUsed++;
                }

            }
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
        my ($fileName, $readLength, $offset, $score) = split("\t", $_, 4);
        $offsets->{$fileName}{$readLength} = $offset;
       
    }
    close($fh);
}