#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;

MAIN:
{
    my $file = shift;
    my $query = "chr18:56193977-56295869";

    my $idx = Bio::DB::HTS::Tabix->new(filename => $file);
    my $iter = $idx->query($query);

    my $lines = 0;
    while (my $gbed = $iter->next)
    {
        $lines++;
    }
    $idx->close;

    print $file,"\t",$lines,"\n";
}


