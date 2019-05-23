#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub loadBedFile($$);
sub processGbedFiles($$$);

MAIN:
{
    my $fileBed = shift;
    my $fileFasta = shift;
    my @filesGbed = @ARGV;
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
    processGbedFiles(\@filesGbed, \@dataBed, $fileFasta);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesGbed), $toc - $tic);
}

### PROCESSGBEDFILES
sub processGbedFiles($$$)
{
    my $filesBam = $_[0];
    my $dataBed = $_[1];
    my $fileFasta = $_[2];

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
