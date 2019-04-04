#!/usr/bin/perl

use warnings;
use strict;

use FindBin::libs;
use Bio::DB::HTS::Tabix;
use Time::HiRes qw(time);
use File::Basename;
use Bed12;

sub loadBedFile($$);
sub allocateTable($$$$);
sub processGbedFiles($$$$$);
sub printTable($);


MAIN:
{
    my $fileBed = shift;
    my $edgeMin = shift;
    my $edgeMax = shift;
    my @filesGbed = @ARGV;
    my @dataBed = ();
    my %table = ();
    my $tic;
    my $toc;

    # parse annotation
    $tic = time();
    loadBedFile(\@dataBed, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(@dataBed), $toc - $tic);

    # allocate table
    $tic = time();
    allocateTable(\%table, \@filesGbed, $edgeMin, $edgeMax);
    $toc = time();
    printf(STDERR "Allocated tables for %d GBED files in %.4f sec.\n", scalar(@filesGbed), $toc - $tic);

    # process gbed files
    $tic = time();
    processGbedFiles(\%table, \@filesGbed, \@dataBed, $edgeMin, $edgeMax);
    $toc = time();
    printf(STDERR "Parsed %d Gbed files in %.4f sec.\n", scalar(@filesGbed), $toc - $tic);

    # print table
    printTable(\%table);
    
}

### PRINT TABLE
sub printTable($)
{
    my $table = $_[0];
    foreach my $file (sort keys %{$table}) {
        my @listPosition = ();
        my @listDepth = ();
        my @listCount = ();
        foreach my $base (sort {$a <=> $b} %{$table->{$file}}) {

            next if(!exists($table->{$file}{$base}));
            push(@listPosition, $base);
            push(@listDepth, $table->{$file}{$base}[0]);
            push(@listCount, $table->{$file}{$base}[1]);
        }
        print $file,"\t",join(",", @listPosition),"\t",join(",", @listDepth),"\t",join(",", @listCount),"\n";
    }
}


### PROCESS GBED FILES
sub processGbedFiles($$$$$)
{
    my $table = $_[0];
    my $filesGbed = $_[1];
    my $dataBed = $_[2];
    my $edgeMin = $_[3];
    my $edgeMax = $_[4];

    foreach my $fileGbed (@{$filesGbed}) {
        
        # process current
        my $fileName = fileparse($fileGbed);
        printf(STDERR "Processing $fileName ... ");
        my $tic = time();

        # open GBED file
        my $hts = Bio::DB::HTS::Tabix->new(filename => $fileGbed);

        # loop through bed records
        foreach my $bed (@{$dataBed}) {

            my $query = $bed->chrom . ":" . $bed->chromStart . "-" . $bed->chromEnd;

            my $iter = $hts->query($query);
            next if(!defined($iter));

            while (my $gbed = $iter->next) {
                my ($chrom, $chromStart, $chromEnd, $depth) = split("\t", $gbed, 4);
                my $baseLinear = $bed->toLinear($chromStart);
                my $baseStart = $baseLinear - $bed->txThickStart + 1;
                my $baseEnd = $baseStart + ($chromEnd - $chromStart);

                my $condition = ($baseStart <= $edgeMax) & ($edgeMin <= $baseEnd);
                next if(!$condition);
                $baseStart = ($baseStart <= $edgeMin) ? $edgeMin : $baseStart;
                $baseEnd = ($baseEnd <= $edgeMax) ? $baseEnd : $edgeMax;

                for(my $k = $baseStart; $k < $baseEnd; $k++)
                {
                    $table->{$fileName}{$k}[0] += $depth;
                    $table->{$fileName}{$k}[1] += 1;
                }
            }

        }

        $hts->close();
        my $toc = time();
        printf(STDERR "done in %.4f sec.\n", ($toc - $tic));

        #last;

    }
}


### ALLOCATE TABLE
sub allocateTable($$$$)
{
    my $table = $_[0];
    my $filesGbed = $_[1];
    my $edgeMin = $_[2];
    my $edgeMax = $_[3];

    foreach my $file (@{$filesGbed}) {
        my $fileName = fileparse($file);
        for (my $base = $edgeMin; $base <= $edgeMax; $base++) {
            $table->{$fileName}{$base} = [0, 0];
        }
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
        #next if (($bed->lengthThick < 300) || (5000 < $bed->lengthThick));
        push(@{$dataBed}, $bed);
    }
    close($fh);
}