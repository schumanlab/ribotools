#!/usr/bin/perl

use warnings;
use strict;

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS;
use Bio::DB::HTS::Tabix;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub processBamFiles($$);

MAIN:
{
    my $fileBed = shift;
    my @filesBam = @ARGV;
    my $tic;
    my $toc;

    # parse annotation
    $tic = time();
    processBamFiles(\@filesBam, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);

}

sub processBamFiles($$)
{
    my $filesBam = $_[0];
    my $fileBed = $_[1];

    # open Bed Index
    my $tabix = Bio::DB::HTS::Tabix->new( filename => $fileBed);

    foreach my $fileBam (@{$filesBam}) {

        # start timer
        my $fileName = fileparse($fileBam);
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();
        my $readsUsed = 0;
        my $test = 0;

        # open Bam file
        my $hts = Bio::DB::HTSfile->open($fileBam);
        my $header       = $hts->header_read;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        while (my $align = $hts->read1($header))
        {
            my $seqid     = $target_names->[$align->tid];
            my $start     = $align->pos+1;
            my $end       = $align->calend;

            my $iter = $tabix->query($seqid . ":" . $start . "-" . $end);
            next if (!defined($iter));
            my $i = 0;
            while (my $n = $iter->next) {
                $i++;
            }
            $test++ if ($i > 0);
            #my $cigar     = $align->cigar_str;
            $readsUsed++;

            last if ($readsUsed == 100000);
        }

        # stop timer
        my $toc = time();
        printf(STDERR "test %d - %d\n", $readsUsed, $test);
        printf(STDERR "done in %.4f sec. Reads used: %d\n", ($toc - $tic), $readsUsed);

    }

    $tabix->close;
}