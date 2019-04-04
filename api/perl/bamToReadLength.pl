#!/usr/bin/perl

use warnings;
use strict;

use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;

MAIN:
{
    my @filesBam = @ARGV;
    my %table = ();
    my $tic;
    my $toc;

    foreach my $fileBam (@filesBam) {

        my $fileName = fileparse($fileBam);
        my $readsUsed = 0;

        # start timer
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();

        # open Bam file
        my $hts = Bio::DB::HTSfile->open($fileBam);
        my $header = $hts->header_read;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        while (my $read = $hts->read1($header)) {
            my $cigar = $read->cigar_str;
            my $readLength = $read->query->length;
            $table{$fileName}{$readLength}++;
        }

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

    }

    foreach my $fileBam (sort keys %table) {
        foreach my $readLength(sort {$a <=> $b} keys %{$table{$fileBam}}) {
            print $fileBam,"\t",$readLength,"\t",$table{$fileBam}{$readLength},"\n";
        }
    }

}
