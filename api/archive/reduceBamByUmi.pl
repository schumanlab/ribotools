#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long();

sub usage($);
sub parseBamHeader($$$);
sub parseBamData($$$$$);
sub printBamData($$);

MAIN:
{
    # define input selection
    my $version = sprintf("version 1.0, Aug 2018, Georgi Tushev\nbug report sciclist\@brain.mpg.de");
    my $file_bam;
    my $help;
    my @list_header = ();
    my $mapq = 0;
    my $read_length_min = 0;
    my $read_length_max = 1000;
    
    # set-up parameters
    Getopt::Long::GetOptions(
    "b|bam=s" => \$file_bam,
    "q|mapq=i" => \$mapq,
    "d|min=i" => \$read_length_min,
    "u|max=i" => \$read_length_max,
    "h|help" => \$help
    ) or usage("Error::invalid command line options");
    
    # parse inputs
    usage($version) if($help);
    usage("Error::bam file is required") unless defined($file_bam);
    
    # define output stream
    my $file_out = $file_bam;
    $file_out =~ s/\.bam$/\_umi\.bam/g;

    # open output handler
    open (my $fw, "|samtools view -bu -|samtools sort -@ 6 -o $file_out -") or die $!;
    
    # parse bam header
    parseBamHeader($file_bam, \@list_header, $fw);
    
    # parse bam data
    parseBamData($file_bam, $mapq, $read_length_min, $read_length_max, $fw);
    
    close($fw);
    
    # index result file
    system("samtools index $file_out");
    
    exit 0;
}

### print bam data
sub printBamData($$)
{
    my $data_ref = $_[0];
    my $fw = $_[1];
    foreach my $key (keys %{$data_ref})
    {
        my $tag = "XU:i:" . $data_ref->{$key}{"count"};
        my $line = "";
        my $freq_old = 0;
        foreach my $read (keys %{$data_ref->{$key}{"versions"}})
        {
            my $freq_new = $data_ref->{$key}{"versions"}{$read}{"count"};
            if ($freq_new > $freq_old)
            {
                $line = $data_ref->{$key}{"versions"}{$read}{"line"} if ($freq_new > $freq_old);
                $freq_old = $freq_new;
            }
        }
        
        print $fw $line,"\t",$tag,"\n";
    }
    
}

### parse bam data
sub parseBamData($$$$$)
{
    my $file_bam = $_[0];
    my $mapq = $_[1];
    my $read_length_min = $_[2];
    my $read_length_max = $_[3];
    my $fw = $_[4];
    
    my %data = ();
    my $chrom = "<unknown>";
    my $count_lines = 0;
    my $count_chroms = 0;
    open (my $fr, "samtools view -q $mapq $file_bam|") or die $!;
    while (<$fr>)
    {
        chomp($_);
        my @line = split('\t', $_, 12);
        
        # filter based on read length
        my $read_length = length($line[9]);
        next if($read_length_min > $read_length);
        next if($read_length_max < $read_length);
        $count_lines++;
        
        # extract umi
        my $umi_start = index($line[0],"_U:S:");
        my $umi_end = index($line[0],"_U:Q:");
        my $umi = substr($line[0],  $umi_start + 5, $umi_end - $umi_start - 5);
        
        # generate key
        my $key = $line[2] . ";" . $umi . ";" . $line[3] . ";" . $read_length;
        
        # check if hash is full
        if ($line[2] ne $chrom)
        {
            $count_chroms++;
            if (scalar(keys %data) >= 10000)
            {
                # print hash
                printBamData(\%data, $fw);
                
                # re-allocate hash
                %data = ();
            }
        }
        
        $data{$key}{"count"}++;
        $data{$key}{"versions"}{$line[9]}{"count"}++;
        $data{$key}{"versions"}{$line[9]}{"line"} = $_;
        
        # update chrom
        $chrom = $line[2];
    }
    close($fr);
    
    printBamData(\%data, $fw);
    
    print STDERR "LINES ",$count_lines,"\n";
    print STDERR "CHROMS ",$count_chroms,"\n";
}

### parse bam header
sub parseBamHeader($$$)
{
    my $file_bam = $_[0];
    my $list_header_ref = $_[1];
    my $fw = $_[2];
    
    open (my $fh, "samtools view -H $file_bam|") or die("samtools can't read bam header");
    while (<$fh>)
    {
        chomp($_);
        print $fw $_,"\n";
        if($_ =~ m/^\@SQ/)
        {
            my @line = split("\t", $_, 3);
            $line[1] =~ s/^SN\://g;
            push(@{$list_header_ref}, $line[1]);
        }
    }
    close($fh);
}


### usage
sub usage($)
{
    my $message = $_[0];
    if (defined $message && length $message)
    {
        $message .= "\n" unless ($message =~ /\n$/);
    }
    
    my $command = $0;
    $command =~ s#^.*/##;
    
    print STDERR (
    
    $message,
    "usage: $command \n" .
    "-b|--bam\n" .
    "\t/path/to/bams/file.bam (required)\n" .
    "-q|--mapq\n" .
    "\tmapping quality of a read (default 0)\n" .
    "-d|--min\n" .
    "\tminimum read length filter (default 0)\n" .
    "-u|--max\n" .
    "\tmaximum read length filter (default 1000)\n" .
    "-help\n" .
    "\tdefine usage\n"
    );
    
    die("\n");
}