#!/usr/bin/perl

use warnings;
use strict;
use FindBin::libs;
use Bed12;

MAIN:
{
    my $fileBed = shift;
    my $offsetStart = shift;
    my $offsetEnd = shift;

    open(my $fh, "<", $fileBed) or die $!;
    while (<$fh>) {
        chomp($_);
        my $bed = Bed12->new();
        $bed->fromLine($_);

        my $elStart = $bed->txThickStart + $offsetStart;
        my $elEnd = $bed->txThickEnd - $offsetEnd;

        next if ($elStart >= $elEnd);

        my ($transcript, $gene) = split(";", $bed->name, 2);
        my $attribute = "gene_id \"$gene\"; transcript_id \"$transcript\";";

        # flip over strand
        if ($bed->strand eq "-") {
            $elStart = $bed->lengthChrom - $elStart;
            $elEnd = $bed->lengthChrom - $elEnd;
            ($elStart, $elEnd) = ($elEnd, $elStart);
        }

        my $offset = 0;
        my $thickStart = 0;
        my $thickEnd = 0;
        for (my $k = 0; $k < $bed->blocks; $k++) {
            my $exonStart = $bed->exonStart->[$k];
            my $exonEnd = $bed->exonEnd->[$k];
            my $offsetNext = $offset + ($exonEnd - $exonStart);

            if (($offset <= $elStart) && ($elStart <= $offsetNext)) {
                $thickStart = $exonStart + ($elStart - $offset);
            }

            if (($offset <= $elEnd) && ($elEnd <= $offsetNext)) {
                $thickEnd = $exonStart + ($elEnd - $offset);
            }

            $offset = $offsetNext;
        }


        for (my $k = 0; $k < $bed->blocks; $k++) {
            my $exonStart = $bed->exonStart->[$k];
            my $exonEnd = $bed->exonEnd->[$k];

            if (($exonStart <= $thickStart) && ($thickStart <= $exonEnd)) {
                $exonStart = $thickStart;
            }

            if (($exonStart <= $thickEnd) && ($thickEnd <= $exonEnd)) {
                $exonEnd = $thickEnd;
            }

            if (($exonStart <= $thickEnd) && ($thickStart <= $exonEnd)) {
                print $bed->chrom,"\t","ribotools","\t","thick","\t",$exonStart,"\t",$exonEnd,"\t",0,"\t",$bed->strand,"\t",".","\t",$attribute,"\n";
                #print $bed->chrom,"\t",$exonStart,"\t",$exonEnd,"\t",$bed->name,"\t",0,"\t",$bed->strand,"\n";
            }


        }


        #print $bed->name,"\t",$bed->txThickStart,"\t",$bed->txThickEnd,"\n" if($elStart >= $elEnd);

        #print $_,"\n";
        #last;

    }
    close($fh);
}