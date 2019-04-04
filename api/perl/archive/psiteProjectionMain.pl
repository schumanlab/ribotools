#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;

sub parsePSiteTable($$);

MAIN:
{
    my $psiteFile = shift;
    my $bamFile = shift;
    my %psites = ();

    parsePSiteTable(\%psites, $psiteFile);
    


    open(my $fh, "samtools view -h -q 255 $bamFile |") or die $!;
    while (<$fh>)
    {
        chomp($_);
        
        if ($_ =~ m/^@/)
        {
            print $_,"\n";
            next;
        }

        my $cigar = ($_ =~ /(([0-9]+[MIDNSHP=X])+)|\*/g) ? $1 : "<unknown>";
        #my $seq = ($_ =~ /([ACGTN]+)/g) ? $1 : "<unknown>";
        
        my @cigarNow = ();
        my $queryLength = 0;
        for my $op (grep defined, $cigar =~ /([0-9]+[MIDNSHP=X])/g)
        {
            my ($len, $type) = $op =~ /(\d*)(\D*)/a;
            $queryLength += $len if($type eq "M");
            #print $type,"\t",$len,"\n";
            push(@cigarNow, [$type, $len]);
        }

        if (!exists($psites{$queryLength}))
        {
            next;
        }

        my $dPsite = $psites{$queryLength} + 1;
        my $cigarNew = "";
        my $offset = 0;
        foreach my $temp (@cigarNow)
        {
            my $op = $temp->[0];
            my $span = $temp->[1];

            if ($op ne "M")
            {
                $cigarNew .= $span . $op;
                next;
            }

            if (($offset < $dPsite) && ($dPsite <= ($offset + $span)))
            {
                my $prevSpan = $dPsite - $offset - 1;
                my $nextSpan = ($offset + $span) - $dPsite;
                $cigarNew .= $prevSpan . "S" if ($prevSpan > 0);
                $cigarNew .= 1 . "M";
                $cigarNew .= $nextSpan . "S" if ($nextSpan > 0);
            }
            else
            {
                $cigarNew .= $span . "S";
            }

            $offset += $span;

            #print $op,"\t",$span,"\n";
        }

=head
        my $queryLengthNew = 0;
        for my $op (grep defined, $cigarNew =~ /([0-9]+[MIDNSHP=X])/g)
        {
            my ($len, $type) = $op =~ /(\d*)(\D*)/a;
            $queryLengthNew += $len if(($type eq "M") || ($type eq "S"));
            #print $type,"\t",$len,"\n";
        }
=cut


        #print $dPsite,"\t",$cigar,"\t",$cigarNew,"\n";
        #my $readLength = length($seq);
        #print $readLength,"\t",$queryLength,"\t",$queryLengthNew,"\t",$cigar,"\t",$cigarNew,"\t","\n" if ($readLength != $queryLengthNew);
        $_ =~ s/$cigar/$cigarNew/g;
        print $_,"\n";
    }
    close($fh);



}

sub parsePSiteTable($$)
{
    my $offsets = $_[0];
    my $psiteFile = $_[1];

    open(my $fh, "<", $psiteFile) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my ($readLength, $pSiteOffset) = split("\t", $_, 2);
        $offsets->{$readLength} = $pSiteOffset;
    }
    close($fh);
}
