#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;
use RecordBed;


MAIN:
{
    my $bedFile = shift;
    my $gbedFile = shift;

    my $tabix = Bio::DB::HTS::Tabix->new(filename => $gbedFile);

    my @coverage = (0) x 100;
    my $lines = 0;
    open (my $fh, "<", $bedFile) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my $bed = RecordBed->new(line => $_);
        $bed->exons();
        next if(($bed->{spanCDS} < 100) || (5000 < $bed->{spanCDS}));

        next if($bed->{strand} eq "+");

        my @linear = (0) x $bed->{spanGene};
        my $query =  $bed->{chrom} . ":" . $bed->{chromStart} . "-" . $bed->{chromEnd};
        my $iter = $tabix->query($query);
        while (my $gbed = $iter->next)
        {
            my ($chrom, $chromStart, $chromEnd, $depth) = split("\t", $gbed, 4);
            my ($offset, $span) = $bed->find($chromStart, $chromEnd);

                if (($bed->{cdsStart}-25 <= $offset) && (($offset + $span) <= $bed->{cdsStart}+75))
                {
                    for (my $j = $offset; $j < ($offset + $span); $j++)
                    {
                        my $idx = $j - $bed->{cdsStart} - 25;
                        $coverage[$idx] += $depth;
                    }
                }
                
        }

        $lines++;
        last if($lines == 5000);
    }
    close($fh);

    my $idx = -25;
    for (my $k = 0; $k < 100; $k++)
    {
        print $idx+$k,"\t",$coverage[$k],"\n";
    }

    $tabix->close;
}