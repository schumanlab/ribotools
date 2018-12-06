package MetaGene;

use strict;
use warnings;

sub new()
{
    my ($class, %args) = @_;

    my $name = $args{"name"};
    my $edge_min = $args{"edge_min"};
    my $edge_max = $args{"edge_max"};
    my $bin_count = $args{"bin_count"};
    my $bin_width = ($edge_max - $edge_min) / ($bin_count - 1);
    my @bin_sum = (0) x $bin_count;
    my @bin_freq = (0) x $bin_count; 

    return bless {
        name => $name,
        edge_min => $edge_min,
        edge_max => $edge_max,
        bin_count => $bin_count,
        bin_width => $bin_width,
        bin_sum => \@bin_sum,
        bin_count => \@bin_freq
    }, $class;

}


sub addValue()
{
    my $self = $_[0];
    my $value_x = $_[1];
    my $value_y = $_[2];

    # calculate bin index
    my $bin_index = int(($value_x - $self->{edge_min}) / $self->{bin_width});
    return if (($bin_index < 0) || ($self->{bin_count} <= $bin_index));
        
    $self->{bin_sum}[$bin_index] += $value_y;
    $self->{bin_count}[$bin_index]++;
}

1;