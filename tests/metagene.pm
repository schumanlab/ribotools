package MetaGene;

sub new
{
    my ($class, %args) = @_;
    print %args,"\n";
    return bless \%args, $class;
}

sub DESTROY
{
    my ($self) = @_;
}


sub open
{
    my $self = $_[0];
    my $argv = $_[1];

    print join(",",@{$argv}),"\n";
}

1;