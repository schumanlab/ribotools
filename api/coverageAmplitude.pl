#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bio::DB::HTS::Tabix;
use Time::HiRes qw(time);
use List::Util qw(max sum);
use File::Basename;
use Bed12;

sub loadBedFile($$);
sub processGbedFiles($$$);
sub printTable($$);

MAIN:
{
    my $fileBed = shift;
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

    # process gbed files
    $tic = time();
    processGbedFiles(\%table, \@dataBed, \@filesGbed);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesGbed), $toc - $tic);

    # print table
    printTable(\%table, \@filesGbed);
}

### PRINTTABLE
sub printTable($$)
{
    my $table = $_[0];
    my $filesGbed = $_[1];

    my @filesName = ();
    foreach my $file (@{$filesGbed}) {
        my $name = fileparse($file);
        push(@filesName, $name);
    }

    print "# Name\tIndex\t",join("\t",@filesName),"\n";

    for (my $k = 0; $k < 3; $k++) {
        foreach my $row (sort keys %{$table}) {
            my @values = ();
            push(@values, $k);
            foreach my $key (@filesName) {
                my $v = exists($table->{$row}{$key}) ? $table->{$row}{$key}[$k] : 0;
                push(@values, $v);
            }
            print $row,"\t",join("\t",@values),"\n";
        }
    }

}

### PROCESSGBEDFILES
sub processGbedFiles($$$)
{
    my $table = $_[0];
    my $dataBed = $_[1];
    my $filesGbed = $_[2];

    foreach my $fileGbed (@{$filesGbed}) {
        # start timer
        my $fileName = fileparse($fileGbed);
        printf(STDERR "Processing $fileName ...");
        my $tic = time();

        # open tabix file
        my $tabix = Bio::DB::HTS::Tabix->new(filename => $fileGbed);

        # loop over bed records
        foreach my $bed (@{$dataBed}) {

            my $maxI = 0;
            my $maxE = 0;
            my $maxT = 0;
            my $sumE = 0;
            my $cntE = 1;
            my ($nameTx, $nameGn) = split(";", $bed->name, 2);
            my $query = $nameTx . ":0-" . $bed->lengthChrom;

            my $iter = $tabix->query($query);
            if (defined($iter)) {

                my $idxStart_I = $bed->txThickStart - 45;
                my $idxEnd_I = $bed->txThickStart + 45;

                my $idxStart_E = $bed->txThickStart + 45;
                my $idxEnd_E = $bed->txThickEnd - 15;

                my $idxStart_T = $bed->txThickEnd - 15;
                my $idxEnd_T = $bed->txThickEnd + 15;

                $cntE = $idxEnd_E - $idxStart_E;
                $cntE = 1 if($cntE < 0);

                while (my $line = $iter->next) {
                    my ($txName, $txStart, $txEnd, $depth) = split("\t", $line, 4);

                    # check initiation
                    if (($txStart <= $idxEnd_I) && ($idxStart_I <= $txEnd)) {
                        $maxI = $depth if ($maxI < $depth);
                    }

                    # check elongation
                    if (($txStart <= $idxEnd_E) && ($idxStart_E <= $txEnd)) {
                        $maxE = $depth if ($maxE < $depth);

                        my $isectStart = ($idxStart_E < $txStart) ? $txStart : $idxStart_E;
                        my $isectEnd = ($txEnd < $idxEnd_E) ? $txEnd : $idxEnd_E;
                        my $isectSize = ($isectEnd - $isectStart);
                        $sumE += ($isectSize * $depth);

                    }

                    # check termination
                    if (($txStart <= $idxEnd_T) && ($idxStart_T <= $txEnd)) {
                        $maxT = $depth if ($maxT < $depth);
                    }

                    
                } 
            }

            $maxI = ($maxI * $cntE) / ($sumE + 1);
            $maxE = ($maxE * $cntE) / ($sumE + 1);
            $maxT = ($maxT * $cntE) / ($sumE + 1);

            $table->{$bed->name}{$fileName} = [$maxI, $maxE, $maxT];

            #print $maxI,"\t",$maxE,"\t",$maxT,"\t",$sumE,"\t",$cntE,"\t",$bed->lengthThick,"\n";

            #last;
        }

        # close tabix file
        $tabix->close;

        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec.\n", ($toc - $tic));

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