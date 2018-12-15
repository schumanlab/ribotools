package Bed12;

use Moose;
use warnings;
use strict;

has 'chrom' => (is => 'rw', isa => 'Str', default => '<unknown>');
has 'chromStart' => (is => 'rw', isa => 'Int', default => 0);
has 'chromEnd' => (is => 'rw', isa => 'Int', default => 0);
has 'name' => (is => 'rw', isa => 'Str', default => '<unknown>');
has 'score' => (is => 'rw', isa => 'Int', default => 0);
has 'strand' => (is => 'rw', isa => 'Str', default => '<unknown>');
has 'thickStart' => (is => 'rw', isa => 'Int', default => 0);
has 'thickEnd' => (is => 'rw', isa => 'Int', default => 0);
has 'itemRgb' => (is => 'rw', isa => 'Str', default => '<unknown>');
has 'blocks' => (is => 'rw', isa => 'Int', default => 0);
has 'blockSizes' => (is => 'rw', isa => 'ArrayRef');
has 'blockStarts' => (is => 'rw', isa => 'ArrayRef');
has 'txThickStart' => (is => 'rw', isa => 'Int', default => 0);
has 'txThickEnd' => (is => 'rw', isa => 'Int', default => 0);
has 'lengthThick' => (is => 'rw', isa => 'Int', default => 0);
has 'lengthChrom' => (is => 'rw', isa => 'Int', default => 0);
has 'exonStart' => (is => 'rw', isa => 'ArrayRef');
has 'exonEnd' => (is => 'rw', isa => 'ArrayRef');
has 'offsets' => (is => 'rw', isa => 'ArrayRef');

sub fromLine()
{
    my $self = $_[0];
    my $line = $_[1];

    my ($chrom, 
        $chromStart, 
        $chromEnd,
        $name,
        $score,
        $strand,
        $thickStart,
        $thickEnd,
        $itemRgb,
        $blocks,
        $blockSizes,
        $blockStarts) = split("\t", $line, 12);
    
    my @listBlockSizes = split(",", $blockSizes);
    my @listBlockStarts = split(",", $blockStarts);

    $self->chrom($chrom);
    $self->chromStart($chromStart);
    $self->chromEnd($chromEnd);
    $self->name($name);
    $self->score($score);
    $self->strand($strand);
    $self->thickStart($thickStart);
    $self->thickEnd($thickEnd);
    $self->itemRgb($itemRgb);
    $self->blocks($blocks);
    $self->blockSizes(\@listBlockSizes);
    $self->blockStarts(\@listBlockStarts);
    $self->exons();
}


sub exons()
{
    my $self = $_[0];
    my $txThickStart = 0;
    my $txThickEnd = 0;
    my @exonStart = (0) x $self->blocks;
    my @exonEnd = (0) x $self->blocks;
    my @offsets = (0) x $self->blocks;
    my $offsetNow = 0;

    for (my $e = 0; $e < $self->blocks; $e++)
    {
        $exonStart[$e] = $self->chromStart + $self->blockStarts->[$e];
        $exonEnd[$e] = $exonStart[$e] + $self->blockSizes->[$e];
        $offsets[$e] = $offsetNow;

        if (($exonStart[$e] <= $self->thickStart) && ($self->thickStart <= $exonEnd[$e]))
        {
            $txThickStart = $self->thickStart - $exonStart[$e] + $offsetNow;
        }

        if (($exonStart[$e] <= $self->thickEnd) && ($self->thickEnd <= $exonEnd[$e]))
        {
            $txThickEnd = $self->thickEnd - $exonStart[$e] + $offsetNow;
        }

        $offsetNow += $self->blockSizes->[$e];
    }

    if ($self->strand eq "-")
    {
        $txThickStart = $offsetNow - $txThickStart;
        $txThickEnd = $offsetNow - $txThickEnd;
        ($txThickStart, $txThickEnd) = ($txThickEnd, $txThickStart);   
    }

    $self->exonStart(\@exonStart);
    $self->exonEnd(\@exonEnd);
    $self->txThickStart($txThickStart);
    $self->txThickEnd($txThickEnd);
    $self->lengthChrom($offsetNow);
    $self->lengthThick($txThickEnd - $txThickStart);
    $self->offsets(\@offsets); 
}


sub toLinear()
{
    my $self = $_[0];
    my $chromQuery = $_[1];
    my $txQuery = -1;

    for (my $e = 0; $e < $self->blocks; $e++)
    {
        if (($self->exonStart->[$e] <= $chromQuery) &&
            ($chromQuery <= $self->exonEnd->[$e]))
        {
            $txQuery = $chromQuery - $self->exonStart->[$e] + $self->offsets->[$e];
            last;
        }
    }

    $txQuery = $self->lengthChrom - $txQuery if($self->strand eq "-");
    $txQuery = -1 if($txQuery > $self->lengthChrom);
    return $txQuery;
}

no Moose;
1;