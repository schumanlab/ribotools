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
sub printTable($$$);

MAIN:
{
    my $fileBed = shift;
    my @filesGbed = @ARGV;
    my %dataBed = ();
    my %table = ();
    my $tic;
    my $toc;

    # parse annotation
    $tic = time();
    loadBedFile(\%dataBed, $fileBed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", scalar(keys %dataBed), $toc - $tic);

    # process gbed files
    $tic = time();
    processGbedFiles(\%table, \%dataBed, \@filesGbed);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@filesGbed), $toc - $tic);

    # print table
    printTable(\%table, \%dataBed, \@filesGbed);
}

### PRINTTABLE
sub printTable($$$)
{
    my $table = $_[0];
    my $dataBed = $_[1];
    my $filesGbed = $_[2];

    my @filesName = ();
    foreach my $file (@{$filesGbed}) {
        my $name = fileparse($file);
        $name =~ s/\_transcriptome\_sorted\_umi\.gbed\.gz//g;
        push(@filesName, $name);
    }

    print "# Name\tcdsStart\tcdsEnd\t",join("\t",@filesName),"\t",join("\t",@filesName),"\n";

    
    foreach my $row (sort keys %{$table}) {

        my $bed = $dataBed->{$row};

        print $row,"\t",$bed->txThickStart,"\t",$bed->txThickEnd;
        foreach my $key (@filesName) {
            my $position = exists($table->{$row}{$key}) ? $table->{$row}{$key}[0] : 0;
            print "\t",$position;
        }

        foreach my $key (@filesName) {
            my $value = exists($table->{$row}{$key}) ? $table->{$row}{$key}[1] : 0;
            print "\t",$value;
        }

        print "\n";
        
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
        $fileName =~ s/\_transcriptome\_sorted\_umi\.gbed\.gz//g;
        printf(STDERR "Processing $fileName ...");
        my $tic = time();

        # open tabix file
        my $tabix = Bio::DB::HTS::Tabix->new(filename => $fileGbed);

        # loop over bed records
        foreach my $name (sort keys %{$dataBed}) {

            my $bed = $dataBed->{$name};
            # get current coverage
            my ($nameTx, $nameGn) = split(";", $bed->name, 2);
            my $query = $nameTx . ":0-" . $bed->lengthChrom;
            my $iter = $tabix->query($query);

            if (defined($iter)) {

                my @array = (0) x $bed->lengthChrom;
                my $cdsStart = $bed->txThickStart + 45;
                my $cdsEnd = $bed->txThickEnd;
                while (my $line = $iter->next) {

                    my ($txName, $txStart, $txEnd, $depth) = split("\t", $line, 4);
                    for (my $k = $txStart; $k < $txEnd; $k++) {
                        $array[$k] = $depth;
                    }
                }

                my $minPosition = 0;
                my $minValue = 0;
                my $cumSum = 0;
                my $total = 1;
                if ($cdsStart < $cdsEnd) {

                    for (my $k = $cdsStart; $k < $cdsEnd; $k++) {
                        $total += $array[$k];
                    }
                    $minValue = $total;

                    for (my $k = $cdsStart; $k < $cdsEnd; $k++) {
                        $cumSum += $array[$k];
                        my $tmpValue = abs(($cumSum/$total) - 0.5);
                        if ($tmpValue < $minValue) {
                            $minValue = $tmpValue;
                            $minPosition = $k;
                        }
                    }
                }

                #print $bed->name,"\t",$cdsStart,"\t",$cdsEnd,"\t",$minPosition,"\t",$minValue,"\t",$total,"\n";
                
                $table->{$bed->name}{$fileName} = [$minPosition, $total/($cdsEnd - $cdsStart)];

            }

            

            
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
        $dataBed->{$bed->name} = $bed;
    }
    close($fh);
}