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

            foreach my $offset (sort {$a <=> $b} keys %{$table->{$file}{$readLength}}) {
                print $file,"\t",$readLength,"\t",$offset,"\t",$table->{$file}{$readLength}{$offset},"\n";
            }

=head
            my @values = ();
            my $total = 0;
            my $count = $table->{$file}{$readLength}{"counts"};
            foreach my $frame (sort {$a <=> $b} keys %{$table->{$file}{$readLength}{"frame"}}) {
                my $val = $table->{$file}{$readLength}{"frame"}{$frame};
                push(@values, $val);
                $total += $val;
            }

            print $file,"\t",$readLength,"\t",$count,"\t",$values[0]/$total,"\t",$values[1]/$total,"\t",$values[2]/$total,"\n";
=cut
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

                my $offsetStart = $bed->txThickStart - $readLinear;
                my $offsetEnd = $bed->txThickEnd - $readLinear - 6;
                
                if ((0 <= $offsetStart) && ($offsetStart <= $readSpan)) {
                    #print $readSpan,"\t",$offsetStart,"\t",($offsetStart%3),"\n";
                    $table->{$fileName}{$readSpan}{$offsetStart}++;
                    $readsUsed++;
                }

                if ((0 <= $offsetEnd) && ($offsetEnd <= $readSpan)) {
                    #print $readSpan,"\t",$offsetStart,"\t",($offsetStart%3),"\n";
                    #$table->{$fileName}{$readSpan}{$offsetEnd}++;
                    #$readsUsed++;
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