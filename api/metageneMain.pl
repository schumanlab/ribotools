#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS::Tabix;
use File::Basename;
use RecordBed;
use MetaGene;

sub allocate_meta($$$$$);
sub tabix_open($$);
sub tabix_close($);
sub tabix_iter($$$);
sub bed_parser($$$);
sub averageCDSDepth($$$);
sub print_meta($$$$);

MAIN:
{
    my $bedFile = shift;
    my @gbedFiles = @ARGV;
    my @tabixHandles = ();
    my @metaHandles = ();

    my $edgeMin = -0.3;
    my $edgeMax = 1.3;
    my $binCounts = 100;
    my $binWidth = ($edgeMax - $edgeMin) / ($binCounts - 1);

    allocate_meta(\@metaHandles, \@gbedFiles, $edgeMin, $edgeMax, $binCounts);
    tabix_open(\@tabixHandles, \@gbedFiles);
    bed_parser($bedFile, \@tabixHandles, \@metaHandles);
    tabix_close(\@tabixHandles);
    print_meta(\@metaHandles, $binCounts, $binWidth, $edgeMin);
}


sub print_meta($$$$)
{
    my $metaHandles = $_[0];
    my $row_count = $_[1];
    my $bin_width = $_[2];
    my $edge_min = $_[3];
    my $col_count = scalar(@{$metaHandles});
    

    # print header
    print "xbin";
    foreach my $meta (@{$metaHandles})
    {
        print "\t",$meta->{name};
    }
    print "\n";
    
    for (my $r = 0; $r < $row_count; $r++)
    {
        my $xbin = $edge_min + ($r * $bin_width) + $bin_width/2;
        printf("%.6f",$xbin);
        foreach my $meta (@{$metaHandles})
        {

            my $value = ($meta->{bin_count}[$r] > 0) ? $meta->{bin_sum}[$r] / $meta->{bin_count}[$r] : 0.0;
            printf("\t%.6f", $value);
        }
        print "\n";

    }



}


sub allocate_meta($$$$$)
{
    my $metaHandles = $_[0];
    my $gbedFiles = $_[1];
    my $edgeMin = $_[2];
    my $edgeMax = $_[3];
    my $binCounts = $_[4];

    foreach my $file (@{$gbedFiles})
    {
        my $name = fileparse($file, ".gbed.gz");
        my $meta = MetaGene->new(name => $name,
                                 edge_min => $edgeMin,
                                 edge_max => $edgeMax,
                                 bin_count => $binCounts);
        push(@{$metaHandles}, $meta);
    }
}


sub bed_parser($$$)
{
    my $bedFile = $_[0];
    my $tabixHandles = $_[1];
    my $metaHandles = $_[2];

    my $loopSize = scalar(@{$tabixHandles});
    my $lines = 0;
    open(my $fh, "<", $bedFile) or die $!;
    while (<$fh>)
    {
        # parse bed line
        chomp($_);
        my $bed = RecordBed->new(line => $_);
        $bed->exons();
        next if(($bed->{spanCDS} < 100) || (5000 < $bed->{spanCDS}));

        # set iterator
        my $query =  $bed->{chrom} . ":" . $bed->{chromStart} . "-" . $bed->{chromEnd};
        my @listIterators = ();
        tabix_iter(\@listIterators, $tabixHandles, $query);

        # update histogram per file
        for (my $i = 0; $i < $loopSize; $i++)
        {
            my @linear = (0) x $bed->{spanGene};
            my $iter = $listIterators[$i];
            next if (!defined($iter));

            while (my $gbed = $iter->next)
            {
                my ($chrom, $chromStart, $chromEnd, $depth) = split("\t", $gbed, 4);
                my ($offset, $span) = $bed->find($chromStart, $chromEnd);

                for (my $j = $offset; $j < ($offset + $span); $j++)
                {
                    $linear[$j] = $depth;
                }
                #@linear[$offset..($offset+$span)] = ($depth) x $span;
            }

            my $depthCDS = averageCDSDepth(\@linear, $bed->{cdsStart}, $bed->{cdsEnd});
            next if($depthCDS < 2);

            for (my $j = 0; $j < $bed->{spanGene}; $j++)
            {
                my $value_x = ($j - $bed->{cdsStart}) / $bed->{spanCDS};
                my $value_y = $linear[$j]/$depthCDS;
                $metaHandles->[$i]->addValue($value_x, $value_y);
            }

        }

        $lines++;
        last if($lines == 10000);

    }
    close($fh);
    print STDERR "LINES: ",$lines,"\n";

}

sub tabix_iter($$$)
{
    my $iterators = $_[0];
    my $tabixHandles = $_[1];
    my $query = $_[2];

    foreach my $tabix (@{$tabixHandles})
    {
        my $iter = $tabix->query($query);
        push(@{$iterators}, $iter);
    }

}


sub tabix_open($$)
{
    my $tabixHandles =  $_[0];
    my $indexFiles = $_[1];

    foreach my $file (@{$indexFiles})
    {
        my $tabix = Bio::DB::HTS::Tabix->new(filename => $file);
        push(@{$tabixHandles}, $tabix);
    }
}

sub tabix_close($)
{
    my $tabixHandles =  $_[0];
    foreach my $tabix (@{$tabixHandles})
    {
        $tabix->close;
    }
}


sub averageCDSDepth($$$)
{
    my $array = $_[0];
    my $idxStart = $_[1];
    my $idxEnd = $_[2];

    my $value = 0;
    for (my $k = $idxStart; $k < $idxEnd; $k++)
    {
        $value += $array->[$k];
    }

    return ($value/($idxEnd - $idxStart));
}