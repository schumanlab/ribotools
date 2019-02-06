package RecordBed;

use strict;
use warnings;

sub new
{
    my ($class, %args) = @_;

    my ($chrom, 
        $chromStart, 
        $chromEnd,
        $name,
        $score,
        $strand,
        $thickStart,
        $thickEnd,
        $itemRGB,
        $blocks,
        $blockSizes,
        $blockStarts) = split("\t", $args{"line"}, 12);
    
    my @listBlockSizes = split(",", $blockSizes);
    my @listBlockStarts = split(",", $blockStarts);
    my @exonStart = (0) x $blocks;
    my @exonEnd = (0) x $blocks;
    my @offset = (0) x $blocks;

    return bless {
        chrom => $chrom,
        chromStart => $chromStart,
        chromEnd => $chromEnd,
        name => $name,
        score => $score,
        strand => $strand,
        thickStart => $thickStart,
        thickEnd => $thickEnd,
        itemRGB => $itemRGB,
        blocks => $blocks,
        blockSizes => \@listBlockSizes,
        blockStarts => \@listBlockStarts,
        exonStart => \@exonStart,
        exonEnd => \@exonEnd,
        offset => \@offset,
        cdsStart => -1,
        cdsEnd => -1,
        spanCDS => -1,
        spanGene => -1
    }, $class;
}

sub exons()
{
    my $self = $_[0];
    my $offset = 0;
    for (my $e = 0; $e < $self->{blocks}; $e++)
    {
        $self->{exonStart}[$e] = $self->{chromStart} + $self->{blockStarts}[$e];
        $self->{exonEnd}[$e] = $self->{exonStart}[$e] + $self->{blockSizes}[$e];
        $self->{offset}[$e] = $offset;

        if (($self->{exonStart}[$e] <= $self->{thickStart}) &&
            ($self->{thickStart} <= $self->{exonEnd}[$e]))
        {
            $self->{cdsStart} = $self->{thickStart} - $self->{exonStart}[$e] + $offset;
        }

        if (($self->{exonStart}[$e] <= $self->{thickEnd}) &&
            ($self->{thickEnd} <= $self->{exonEnd}[$e]))
        {
            $self->{cdsEnd} = $self->{thickEnd} - $self->{exonStart}[$e] + $offset;
        }
        
        $offset += $self->{blockSizes}[$e];
    }

    $self->{spanGene} = $offset;
    $self->{spanCDS} = $self->{cdsEnd} - $self->{cdsStart};

    if ($self->{strand} eq "-")
    {
        $self->{cdsStart} = $self->{spanGene} - $self->{cdsStart};
        $self->{cdsEnd} = $self->{spanGene} - $self->{cdsEnd};
        ($self->{cdsStart}, $self->{cdsEnd}) = ($self->{cdsEnd}, $self->{cdsStart});   
    }

}


sub linear()
{
    my $self = $_[0];
    my $qryPosition = $_[1];
    my $linearPosition = -1;

    for (my $e = 0; $e < $self->{blocks}; $e++)
    {
        if (($self->{exonStart}[$e] <= $qryPosition) &&
            ($qryPosition <= $self->{exonEnd}[$e]))
        {
            $linearPosition = $qryPosition - $self->{exonStart}[$e] + $self->{offset}[$e];
            last;
        }
    }

    $linearPosition = $self->{spanGene} - $linearPosition if($self->{strand} eq "-");
    $linearPosition = -1 if ($linearPosition > $self->{spanGene});
    return $linearPosition;
}

sub find()
{
    my $self = $_[0];
    my $qryStart = $_[1];
    my $qryEnd = $_[2];
    my $offsetPositive = $self->{spanGene};
    my $offsetNegative = 0;
    my $offset = -1;
    my $span = 0;

    for (my $e = 0; $e < $self->{blocks}; $e++)
    {
        if (($qryStart <= $self->{exonEnd}[$e]) && 
            ($self->{exonStart}[$e] <= $qryEnd))
        {
            my $isectStart = max($qryStart, $self->{exonStart}[$e]);
            my $isectEnd = min($qryEnd, $self->{exonEnd}[$e]);
            my $offsetStart = $isectStart - $self->{exonStart}[$e] + $self->{offset}[$e];
            my $offsetEnd = $isectEnd - $self->{exonStart}[$e] + $self->{offset}[$e];
            $offsetPositive = min($offsetPositive, $offsetStart);
            $offsetNegative = max($offsetNegative, $offsetEnd);
            $span += ($isectEnd - $isectStart);
        }
    }

    if ($self->{strand} eq "+")
    {
        $offset = $offsetPositive;
    }
    else
    {
        $offset = $self->{spanGene} - $offsetNegative;
    }

    return ($offset, $span);
}


sub min()
{
    my $x = $_[0];
    my $y = $_[1];
    return ($x <= $y) ? $x : $y;
}


sub max()
{
    my $x = $_[0];
    my $y = $_[1];
    return ($x <= $y) ? $y : $x;
}


sub swap()
{
    my $x = $_[0];
    my $y = $_[1];

    my $t = $x;
    $x = $y;
    $y = $t;

    return ($x, $y);
}
1;
