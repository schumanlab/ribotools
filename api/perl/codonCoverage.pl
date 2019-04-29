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

MAIN:
{
    my $pathTrack = shift;
    my $fileBed = shift;
    my @filesBam = @ARGV;
    my @dataBed = ();
    my $tic;
    my $toc;

    # parse annotation
    $tic = time();
    loadBedFile(\@dataBed, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(@dataBed), $toc - $tic);

    # process bam files
    $tic = time();
    processBamFiles(\@filesBam, \@dataBed, $pathTrack);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesBam), $toc - $tic);





}

### PROCESSBAMFILES
sub processBamFiles($$$)
{
    my $filesBam = $_[0];
    my $dataBed = $_[1];
    my $pathTrack = $_[2];

    foreach my $fileBam (@{$filesBam}) {

        my $fileName = fileparse($fileBam);
        $fileName =~ s/_transcriptome_sorted_umi.bam//g;
        my $fileOut = $pathTrack . $fileName . ".txt.gz";
        open(my $fh, "| sort -k1,1 -k3,3n -k4,4n | bgzip > $fileOut") or die $!;
        my $readsUsed = 0;

        my @offsetNext = ([0, 2, 1], [1, 0, 2], [2, 1, 0]);
        my @offsetPrev = ([0, -1, -2], [-2, 0, -1], [-1, -2, 0]);

        # start timer
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();

        # open Bam file
        my $hts = Bio::DB::HTS->new(-bam => $fileBam);
        
        # loop through bed records
        foreach my $bed (@{$dataBed}) {
            my ($transcript, $gene) = split(";", $bed->name, 2);
            my @reads = $hts->get_features_by_location(-seq_id => $transcript,
                                                       -start  => 0,
                                                       -end    => $bed->lengthChrom);

            next if (scalar(@reads) == 0);

            my $indexCDS = ($bed->txThickStart % 3);
            my %track = ();
        
            foreach my $read (@reads) {

                my $readStart = $read->start + $offsetNext[$indexCDS][($read->start % 3)];
                my $readEnd = $read->end + $offsetPrev[$indexCDS][$read->end % 3];
                
                for (my $base = $readStart; $base <= $readEnd; $base++) {
                    $track{$base}++;
                }
            }

            print $fh $transcript,"\t",$gene,"\t",$bed->txThickStart,"\t",$bed->txThickEnd,"\t",$bed->lengthChrom,"\t";
            foreach my $base (sort {$a <=> $b} keys %track) {
                print $fh $base,",",$track{$base},",";
            }
            print $fh "\n";

            
        }

        close($fh);
        system("tabix --zero-based --sequence 1 --begin 3 --end 4 $fileOut");

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
    while (<$fh>) {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);
        push(@{$dataBed}, $bed);
    }
    close($fh);
}

sub closestDivisor($$)
{
    my $a = $_[0];
    my $b = $_[1];

    my $c_prev = $a - ($a % $b);
    my $c_next = ($a + $b) - ($a % $b);

    return (($a - $c_prev) > ($c_next - $a)) ? $c_next : $c_prev;
}