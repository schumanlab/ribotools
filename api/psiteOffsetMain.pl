#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;
use File::Basename;
use RecordBed;

MAIN:
{
    my $bedFile = shift;
    my $bamFile = shift;

    open(my $fh, "<", $bedFile) or die $!;
    while (<$fh>)
    {
        # parse bed line
        chomp($_);
        my $bed = RecordBed->new(line => $_);
        $bed->exons();
        next if ($bed->{spanCDS} <= 40);

        
    }
    close($fh);

}

