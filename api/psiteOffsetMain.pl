#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;
use RecordBed;

sub getOffsets($$$);
sub printOffsets($);

MAIN:
{
    my $bedFile = shift;
    my $bamFile = shift;
    my %offsets = ();

    # bam handle
    my $hts = Bio::DB::HTS->new(-bam => $bamFile);

    my $lines = 0;
    open(my $fh, "<", $bedFile) or die $!;
    while (<$fh>)
    {
        # parse bed line
        chomp($_);
        my $bed = RecordBed->new(line => $_);
        $bed->exons();
        next if ($bed->{spanCDS} <= 40);

        my @alignments = ();
        if ($bed->{strand} eq "+")
        {
            @alignments = $hts->get_features_by_location(-seq_id => $bed->{chrom},
                                                            -start  => $bed->{thickStart},
                                                            -end    => $bed->{thickStart}+1);
        }
        else
        {
            @alignments = $hts->get_features_by_location(-seq_id => $bed->{chrom},
                                                            -start  => $bed->{thickEnd}-1,
                                                            -end    => $bed->{thickEnd});
            
        }
        getOffsets(\%offsets, \@alignments, $bed);

        $lines++;
        #last if($lines == 10000);
    }
    close($fh);

    printOffsets(\%offsets);
}

sub printOffsets($)
{
    my $offsets = $_[0];

    foreach my $readLength (sort keys %{$offsets})
    {
        my $dStepMax = 0;
        my $dStepValue = 0;
        foreach my $dStep (sort keys %{$offsets->{$readLength}})
        {
            if ($dStepValue < $offsets->{$readLength}{$dStep})
            {
                $dStepValue = $offsets->{$readLength}{$dStep};
                $dStepMax = $dStep;
            }
        }
        print $readLength,"\t",$dStepMax,"\n";
    }
}


sub getOffsets($$$)
{
    my $offsets = $_[0];
    my $alignments = $_[1];
    my $bed = $_[2];

    foreach my $bam (@{$alignments})
    {
        my $readStart = $bam->start + $bam->query->start - 1;
        my $readLength = $bam->query->length;
        my $step = ($bed->{strand} eq "+") ? $bed->linear($readStart) : $bed->linear($readStart + $readLength);
        next if($step == -1);
        my $dStep = ($bed->{cdsStart} - $step + 1);

        #print $readLength,"\t",$bed->{strand},"\t",$step,"\t",$bed->{cdsStart},"\t",$dStep,"\n";
        
        $offsets->{$readLength}{$dStep}++ if((0 <= $dStep) && ($dStep <= $readLength));
    }
}

