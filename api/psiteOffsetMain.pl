#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;
use File::Basename;
use RecordBed;

sub getOffsets($$);

MAIN:
{
    my $bedFile = shift;
    my $bamFile = shift;

    # bam handle
    my $hts = Bio::DB::HTS->new(-bam => $bamFile);

    open(my $fh, "<", $bedFile) or die $!;
    while (<$fh>)
    {
        # parse bed line
        chomp($_);
        my $bed = RecordBed->new(line => $_);
        $bed->exons();
        next if ($bed->{spanCDS} <= 40);

        my @left = $hts->get_features_by_location(-seq_id => $bed->{chrom},
                                                 -start  => $bed->{thickStart},
                                                 -end    => $bed->{thickStart}+1);
        
        my @ritght = $hts->get_features_by_location(-seq_id => $bed->{chrom},
                                                 -start  => $bed->{thickEnd}-1,
                                                 -end    => $bed->{thickEnd});

        getOffsets($bed, \@left);

    }
    close($fh);

}

sub getOffsets($$)
{
    my $bed = $_[0];
    my $aligns = $_[1];

    foreach my $bam (@{$aligns})
    {
        my $readStart = $bam->start;
        my $readEnd = $bam->end;
        my $readLength = $bam->query->length;
        my ($offset, $span) = $bed->find($readStart, $readEnd);
        #print $bed->{thickStart}," : ",
        #      $readStart, " - ",
        #      $readEnd, "\t",
        #      ($bed->{cdsStart} - $offset)," : ",
        #      ($offset + $span - $bed->{cdsStart}),"\n";
        print $readLength,"\t",($bed->{cdsStart}) - $offset,"\n";


    }

}

