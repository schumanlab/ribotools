#!/usr/bin/perl

use warnings;
use strict;
use Bio::DB::HTS;
use Set::IntervalTree;
use Time::HiRes qw(time);

use File::Basename;
use POSIX qw/ceil/;

sub parseConfigPSiteOffsets($$);
sub parseConfigScaleFactors($$);
sub bedToAnnotationTree($$);
sub bedToHash($$);
sub parseExons($);
sub processBams($$$$$);
sub printHistogram($);

MAIN:
{
    my $file_configPSiteOffsets = shift;
    my $file_configScaleFactors = shift;
    my $file_bed = shift;
    my @file_bam_list = @ARGV;
    
    my %configPSiteOffsets = ();
    my %configScaleFactors = ();
    my %bedHash = ();
    my %hist = ();
    
    my $tic;
    my $toc;
    
    # parse configuration file
    $tic = time();
    parseConfigPSiteOffsets(\%configPSiteOffsets, $file_configPSiteOffsets);
    $toc = time();
    printf(STDERR "Parsed config file in %.4f sec.\n", $toc - $tic);
    
    $tic = time();
    parseConfigScaleFactors(\%configScaleFactors, $file_configScaleFactors);
    $toc = time();
    printf(STDERR "Parsed config file in %.4f sec.\n", $toc - $tic);
    
    # parse annotation
    $tic = time();
    my $linesBed = bedToHash(\%bedHash, $file_bed);
    $toc = time();
    printf(STDERR "Parsed %d BED lines in %.4f sec.\n", $linesBed, $toc - $tic);
    
    # process list of bam files
    $tic = time();
    processBams(\%hist,
                \%configPSiteOffsets,
                \%configScaleFactors,
                \%bedHash,
                \@file_bam_list);
    $toc = time();
    printf(STDERR "Parsed %d BAM files in %.4f sec.\n", scalar(@file_bam_list), $toc - $tic);

    #printHistogram(\%hist);
}

sub printHistogram($)
{
    my $hist = $_[0];
    
    foreach my $rpos (sort {$a <=> $b} keys %{$hist})
    {
        foreach my $tag (sort keys %{$hist->{$rpos}})
        {
            my $fp_depth = exists($hist->{$rpos}{$tag}{"fpSide"}{"depth"}) ? $hist->{$rpos}{$tag}{"fpSide"}{"depth"} : 0;
            my $fp_cnt = exists($hist->{$rpos}{$tag}{"fpSide"}{"gene"}) ? scalar(keys %{$hist->{$rpos}{$tag}{"fpSide"}{"gene"}}) : 1;
            
            
            my $cnt_depth = exists($hist->{$rpos}{$tag}{"center"}{"depth"}) ? $hist->{$rpos}{$tag}{"center"}{"depth"} : 0;
            my $cnt_cnt = exists($hist->{$rpos}{$tag}{"center"}{"gene"}) ? scalar(keys %{$hist->{$rpos}{$tag}{"center"}{"gene"}}) : 1;
            
            
            my $tp_depth = exists($hist->{$rpos}{$tag}{"tpSide"}{"depth"}) ? $hist->{$rpos}{$tag}{"tpSide"}{"depth"} : 0;
            my $tp_cnt = exists($hist->{$rpos}{$tag}{"tpSide"}{"gene"}) ? scalar(keys %{$hist->{$rpos}{$tag}{"tpSide"}{"gene"}}) : 1;
            
            #print $tag,"\t",$rpos,"\t",$fp_cnt,"\t",$cnt_cnt,"\t",$tp_cnt,"\n";
            print $tag,"\t",$rpos,"\t",$fp_depth/$fp_cnt,"\t",$cnt_depth/$cnt_cnt,"\t",$tp_depth/$tp_cnt,"\n";
        }
    }
}


# PARSE CONFIG PSITE OFFSETS
sub parseConfigPSiteOffsets($$)
{
    my $config = $_[0];
    my $file_config = $_[1];
    
    open(my $fh, "<", $file_config) or die $!;
    while (<$fh>)
    {
        chomp($_);
        my ($tag, $span, $offset) = split("\t", $_, 3);
        $config->{$tag}{$span} = $offset;
    }
    close($fh);
}

# PARSE CONFIG SCALE FACTORS
sub parseConfigScaleFactors($$)
{
    my $config = $_[0];
    my $file_config = $_[1];
    
    open(my $fh, "<", $file_config) or die $!;
    while (<$fh>)
    {
        chomp($_);
        next if($_ =~ m/^#/);
        my ($tag, $group, $reads) = split("\t", $_, 3);
        $config->{$tag} = $reads;
    }
    close($fh);
}


# PROCESS BAMS
sub processBams($$$$$)
{
    my $hist = $_[0];
    my $config_offset = $_[1];
    my $config_factors = $_[2];
    my $data = $_[3];
    my $file_bam_list = $_[4];
    
    # loop through bam files
    foreach my $file_bam (@{$file_bam_list})
    {
        # start timer
        my $tic = time();
        
        # parse file name
        my $file_name = fileparse($file_bam);
        #my $file_tag = substr($file_name, 0, index($file_name,"_", index($file_name,"_") + 1));
        my $file_tag = substr($file_name, 0, index($file_name, "_genome"));
        printf(STDERR "Processing $file_tag ... ");
        
        # bam handle
        my $hts = Bio::DB::HTS->new(-bam => $file_bam);
        
        # loop through bed hash
        my $readsUsed = 0;
        my $bedLines = 0;
        foreach my $key (keys %{$data})
        {
            # current bed data
            my $bed = $data->{$key}[0];
            my $exonTree = $data->{$key}[1];
            my $cdsStart = $data->{$key}[2];
            my $cdsCenter = $data->{$key}[3];
            my $cdsEnd = $data->{$key}[4];
            my $span = $data->{$key}[5];
            
            # query bam file
            my @bam = $hts->get_features_by_location(-seq_id => $bed->[0],
                                                     -start  => $bed->[6],
                                                     -end    => $bed->[7]);
            
            # filter based on gene expression
            my $factor = 1000000 / $config_factors->{$file_tag};
            my $reads = scalar(@bam);
            my $rpm = ($factor * $reads);
            #next if(($rpm > 16) || ($reads == 0));
            next if($rpm < 16);
            
            # calculate incrementer
            #my $value = $factor / ($rpm/($cdsEnd - $cdsStart));
            my @values = (0) x 303;
            my $readsUsed = 0;
            my $readsRange = 0;
            
            foreach my $record (@bam)
            {
                my $readStart = $record->pos;
                my $readLength = $record->query->length;
                my $readPoffset = exists($config_offset->{$file_tag}{$readLength}) ? $config_offset->{$file_tag}{$readLength} : 0;
                my $info = $exonTree->fetch($readStart, $readStart+1);
                next if(!defined($info->[0][0]));
                #$hist->{"span"}{$readLength}++;
                $readsUsed++;
                
                # calculate read offset
                my $readOffset = ($readStart - $info->[0][0]) + $info->[0][1];
                
                # swap direction
                if ($bed->[5] eq "-")
                {
                    $readOffset = $span - $readOffset;
                }
                
                # add p-site offset
                $readOffset += $readPoffset;
                
                # relative position
                my $relOffsetStart = $readOffset - $cdsStart;
                my $relOffsetCenter = $readOffset - $cdsCenter;
                my $relOffsetEnd = $readOffset - $cdsEnd;
                
                # fill up histogram
                if ((-50 <= $relOffsetStart) && ($relOffsetStart <= 50))
                {
                    #$hist->{$relOffsetStart}{$file_tag}{"fpSide"}++;
                    #$hist->{$relOffsetStart}{$file_tag}{"fpSide"}{"depth"} += $value;
                    #$hist->{$relOffsetStart}{$file_tag}{"fpSide"}{"gene"}{$key}++;
                    $values[$relOffsetStart + 50]++;
                    $readsRange++;
                    
                }
                
                if ((-50 <= $relOffsetCenter) && ($relOffsetCenter <= 50))
                {
                    #$hist->{$relOffsetCenter}{$file_tag}{"center"}++;
                    #$hist->{$relOffsetCenter}{$file_tag}{"center"}{"depth"} += $value;
                    #$hist->{$relOffsetCenter}{$file_tag}{"center"}{"gene"}{$key}++;
                    $values[$relOffsetCenter + 151]++;
                    $readsRange++;
                }
                
                if ((-50 <= $relOffsetEnd) && ($relOffsetEnd <= 50))
                {
                    #$hist->{$relOffsetEnd}{$file_tag}{"tpSide"}++;
                    #$hist->{$relOffsetEnd}{$file_tag}{"tpSide"}{"depth"} += $value;
                    #$hist->{$relOffsetEnd}{$file_tag}{"tpSide"}{"gene"}{$key}++;
                    $values[$relOffsetEnd + 252]++;
                    $readsRange++;
                }
                
                
            }
            
            #print $file_tag,"\t",$bed->[3],"\t",$span,"\t",$cdsEnd-$cdsStart,"\t",$readsUsed,"\t",$readsRange,"\t",join(",",@values),"\n" if($readsRange > 0);
            

            
            $bedLines++;
            last if($bedLines == 100);
        }
        
        # stop timer
        my $toc = time();
        printf(STDERR "done in %.4f sec.\n", ($toc - $tic));
        #last;
        
    }
    
}


        
# BED TO HASH
sub bedToHash($$)
{
    my $data = $_[0];
    my $file_bed = $_[1];
    my $lines = 0;
    open(my $fh, "gunzip -c $file_bed|") or die $!;
    while(<$fh>)
    {
        # parse bed line
        chomp($_);
        my @bed = split(/\t/, $_, 12);
        
        # skip no CDS
        next if($bed[6] == $bed[7]);
        
        # parse exons
        my ($exonTree, $cdsStart, $cdsEnd, $span) = parseExons(\@bed);
        
        # swap direction
        if ($bed[5] eq "-")
        {
            my $temp = $cdsStart;
            $cdsStart = $span - $cdsEnd;
            $cdsEnd = $span - $temp;
        }
        
        # cds center
        my $cdsCenter = int(($cdsEnd - $cdsStart)/2);
        $cdsCenter = $cdsStart + ($cdsCenter - ($cdsCenter%3));
        
        # add to tree
        my $key = $bed[0] . ";" . $bed[3];
        $data->{$key} = [\@bed, $exonTree, $cdsStart, $cdsCenter, $cdsEnd, $span];
        
        # increment line counter
        $lines++;
    }
    close($fh);
    
    return $lines;
}


# BED TO ANNOTATION TREE
sub bedToAnnotationTree($$)
{
    my $data = $_[0];
    my $file_bed = $_[1];
    my $lines = 0;
    open(my $fh, "gunzip -c $file_bed|") or die $!;
    while(<$fh>)
    {
        # parse bed line
        chomp($_);
        my @bed = split(/\t/, $_, 12);
        
        # skip no CDS
        next if($bed[6] == $bed[7]);
        
        # parse exons
        my ($exonTree, $cdsStart, $cdsEnd, $span) = parseExons(\@bed);
        
        # swap direction
        if ($bed[5] eq "-")
        {
            my $temp = $cdsStart;
            $cdsStart = $span - $cdsEnd;
            $cdsEnd = $span - $temp;
        }
        
        # cds center
        my $cdsCenter = int(($cdsEnd - $cdsStart)/2);
        $cdsCenter = $cdsStart + ($cdsCenter - ($cdsCenter%3));
        
        # add to tree
        if (!exists($data->{$bed[0]}))
        {
            $data->{$bed[0]} = Set::IntervalTree->new;
        }
        $data->{$bed[0]}->insert([$exonTree, $cdsStart, $cdsCenter, $cdsEnd, $span, $bed[5]], $bed[1], $bed[2]);
        
        # increment line counter
        $lines++;
    }
    close($fh);
    
    return $lines;
}

# PARSE EXONS
sub parseExons($)
{
    my $bed = $_[0];
    my @blockSizes = split(/,/, $bed->[10]);
    my @blockStarts = split(/,/, $bed->[11]);
    
    my $exonTree = Set::IntervalTree->new;
    my $cdsStart = 0;
    my $cdsEnd = 0;
    my $span = 0;
    for (my $k = 0; $k < $bed->[9]; $k++)
    {
        my $exonStart = $bed->[1] + $blockStarts[$k];
        my $exonEnd = $exonStart + $blockSizes[$k];
        $exonTree->insert([$exonStart, $span], $exonStart, $exonEnd);
        
        if (($exonStart <= $bed->[6]) && ($bed->[6] <= $exonEnd))
        {
            $cdsStart = ($bed->[6] - $exonStart) + $span;
        }
        
        if (($exonStart <= $bed->[7]) && ($bed->[7] <= $exonEnd))
        {
            $cdsEnd = ($bed->[7] - $exonStart) + $span;
        }
        
        $span += $blockSizes[$k];
    }
    
    return ($exonTree, $cdsStart, $cdsEnd, $span);
}


=head
 # bam handle
 my $hbam = Bio::DB::HTSfile->open($file_bam);
 my $header = $hbam->header_read;
 my $target_names = $header->target_name;
 my $reads = 0;
 while (my $bam = $hbam->read1($header))
 {
 # current alignment
 my $chrom = $target_names->[$bam->tid];
 my $readStart = $bam->pos;
 my $readLength = $bam->query->length;
 my $readPoffset = exists($config->{$file_tag}{$readLength}) ? $config->{$file_tag}{$readLength} : 0;
 
 # find gene
 next if(!exists($bedTree->{$chrom}));
 my $hit = $bedTree->{$chrom}->fetch($readStart, $readStart+1);
 next if (!defined($hit->[0][0]));
 my $exonTree = $hit->[0][0];
 my $cdsStart = $hit->[0][1];
 my $cdsCenter = $hit->[0][2];
 my $cdsEnd = $hit->[0][3];
 my $span = $hit->[0][4];
 my $strand = $hit->[0][5];
 
 # find exon
 my $info = $exonTree->fetch($readStart, $readStart+1);
 next if(!defined($info->[0][0]));
 
 # calculate read offset
 my $readOffset = ($readStart - $info->[0][0]) + $info->[0][1];
 
 # swap direction
 if ($strand eq "-")
 {
 $readOffset = $span - $readOffset;
 }
 
 # add p-site offset
 #$readOffset += $readPoffset;
 
 # relative position
 my $relOffsetStart = $readOffset - $cdsStart;
 my $relOffsetCenter = $readOffset - $cdsCenter;
 my $relOffsetEnd = $readOffset - $cdsEnd;
 
 # fill up histogram
 if ((-50 <= $relOffsetStart) && ($relOffsetStart <= 50))
 {
 $hist->{$relOffsetStart}{$file_tag}{"fpSide"}++;
 }
 
 if ((-50 <= $relOffsetCenter) && ($relOffsetCenter <= 50))
 {
 $hist->{$relOffsetCenter}{$file_tag}{"center"}++;
 }
 
 if ((-50 <= $relOffsetEnd) && ($relOffsetEnd <= 50))
 {
 $hist->{$relOffsetEnd}{$file_tag}{"tpSide"}++;
 }
 
 
 $reads++;
 last if($reads == 1000000);
 }
=cut