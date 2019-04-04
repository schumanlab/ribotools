#!/usr/bin/perl

use warnings;
use strict;

sub parseOffsetsFile($$);
sub filterOffsets($$);

MAIN:
{
    my $fileOffsets = shift;
    my $filterBest = shift;
    my %table = ();

    parseOffsetsFile($fileOffsets, \%table);
    filterOffsets($filterBest, \%table);
}


sub filterOffsets($$)
{
    my $filterBest = $_[0];
    my $table = $_[1];

    foreach my $fileName (sort keys %{$table}) {
        my $k = 0;
        foreach my $score (sort {$b <=> $a} keys %{$table->{$fileName}}) {
            print $fileName,"\t",
            $table->{$fileName}{$score}[0],"\t",
            $table->{$fileName}{$score}[1],"\t",
            $score,"\n";
            last if ($k++ >= ($filterBest - 1));
        }
    }
}


sub parseOffsetsFile($$)
{
    my $fileOffsets = $_[0];
    my $table = $_[1];

    open(my $fh, "<", $fileOffsets) or die $!;
    while (<$fh>) {
        chomp($_);
        my ($fileName, $readLength, $frame, $offset, $score) = split("\t", $_, 5);
        $table->{$fileName}{$score} = [$readLength, $frame, $offset];
    }
    close($fh);
}