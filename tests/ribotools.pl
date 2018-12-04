#!/usr/bin/perl

use warnings;
use strict;


use MetaGene;

my $VERSION = "1.0.0";
sub usage($);

MAIN:
{
    if (scalar(@ARGV) < 1)
    {
        usage("");
    }
    else
    {
        my $subcommand = shift @ARGV;
        if ($subcommand eq "metagene")
        {
            my $obj = MetaGene->new($subcommand);
        }
        else
        {
            usage("Error::ribotools: unknown subcommand.");
        }
    }
}


sub usage($)
{
    my $message = $_[0];
    if (defined $message && length $message)
    {
        $message .= "\n\n" unless($message =~ /\n\n$/);
        $message = "\n" . $message unless($message =~ /^\n/);
    }
    
    my $command = $0;
    $command =~ s#^.*/##;
    
    print STDERR (
        $message,
        "$command is a toolset for ribosome footprint analysis.\n" .
        "version $VERSION\n\n" .
        "About: developed in the Scientific Computing Facility at Max-Planck Institute For Brain Research\n" .
        "Docs: https://ribotools.github.molgen.mpg.de\n" .
        "Code: https://github.molgen.mpg.de/MPIBR-Bioinformatics/ribotools\n" .
        "Mail: sciclist\@brain.mpg.de\n\n" .
        "Usage: perl $command.pl <subcommand> [options]\n\n" .
        "The ribotools sub-commands include:\n" .
        "\tmetagene\tcreates a normalized histogram of footprint coverage\n" .
        "\n" .
        "General help:\n" .
        "\t-h | --help\tprint this help menu.\n" .
        "\t-v | --version\twhat version of ribotools are you using?\n" .
        "\t-c | --contact\tfeature requests, bugs, mailing lists, etc.\n" .
        "\n"
    );
    
    die("\n");
}


