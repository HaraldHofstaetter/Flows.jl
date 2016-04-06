
#$Rev: 35 $
#$Author: hofi $
#$Date: 2012-11-02 16:38:27 +0100 (Fri, 02 Nov 2012) $

use strict;

package TimeExpression;

our $t_zero = EmptyTimeLinearCombination->new();
our $t_dummy = TimeVariable->new("t_dummy");

use overload 
      "+" => "plus",
      "-" => "minus",
      "*" => "times",
      "fallback" => 1;

sub new {
    my $class = shift;
    my $self = { };
    return bless $self, $class;
}

sub plus {
    my ($self, $other) = @_;
    return TimeLinearCombination->new($self, 1, $other, 1);
}    

sub minus {
    my ($self, $other) = @_;
    return TimeLinearCombination->new($self, 1, $other, -1);
}    

sub times {
    my ($self, $other, $swap) = @_;
    return TimeLinearCombination->new($self, $other);
}

my %_expressions;   
my %_xxx;  

sub get_adr { $_[0] =~/HASH\(0x([0-9a-f]+)\)$/; return $1; }

sub _register {
    my $self = shift;
    my $key = $self->get_register_key();
    my $x;
    unless ($x = $_expressions{$key}[0]) {
       $_expressions{$key} = [$self, scalar(keys %_xxx)];
       $_xxx{$self} = scalar(keys %_xxx);
       $x = $self; 
    }
    return $x;
}

sub print_register {
    print "time_expressions:\n";
    for my $i (sort { $a->[1] <=> $b->[1] } values %_expressions){
        print '#' . $i->[1]. "\t". $i->[0]->_str_flat() ."\n";
    }   
}

sub _str_flat_arg_name {
    my $self = shift;
    return '#' . $_xxx{$self};
}

package TimeVariable;
our @ISA = qw( TimeExpression );

sub new {
    my $class = shift;
    my $self = TimeExpression->new();
    $self->{name} = shift;
    $self->{latex} = shift; #optional
    bless $self, $class;
    return $self;
}

sub str {
    my $self = shift;
    return $self->{name};
} 

sub latex {
    my $self = shift;
    return $self->{latex} || $self->{name};
}

sub coefficient {
    my ($self, $var) = @_;
    return ($self==$var ? 1 : 0);
}   

sub substitute_time_variable {
    my ($self, $this, $by) = @_;
    return ($self == $this ? $by : $self);
}
 

sub _str_flat_arg_name {
    my $self = shift;
    return $self->{name};
}


package TimeLinearCombination;
our @ISA = qw( TimeExpression );

use Scalar::Util qw( looks_like_number );

sub new {
    my $class = shift;
    my $self = TimeExpression->new();
    my @x=();
    while (@_) {
        my $expr = shift;
	die "TimeExpression expected, stopped" unless ref($expr) and $expr->isa("TimeExpression"); 
        my $coeff = shift;
        die "Number expected, stopped" unless looks_like_number($coeff);
    	push(@x, [$expr, $coeff] );
    }
    $self->{terms} =  \@x ;
    bless $self, $class;
    $self->simplify();
    return $self;
}

sub str {
    my $self = shift;
    my $s = join ('+', map { ($_->[1]==1 ? '' : ($_->[1]==-1 ? '-' :$_->[1] )) . $_->[0]->str() } 
            @{$self->{terms}});
    #$s = "(" . $s . ")"; # parentheses should not be necessary (no ambiguities expected)
    $s =~ s/\+\-/-/g;  # handle minus signs
    return $s
}

sub latex {
    my $self = shift;
    my $s = join ('+', map { ($_->[1]==1 ? '' : ($_->[1]==-1 ? '-' :$_->[1] )) . $_->[0]->latex() } 
            @{$self->{terms}});
    #$s = "(" . $s . ")"; # parentheses should not be necessary (no ambiguities expected)
    $s =~ s/\+\-/-/g;  # handle minus signs
    return $s
}



sub _str_flat {
    my $self = shift;
    my $s = join ('+', map { ($_->[1]==1 ? '' : ($_->[1]==-1 ? '-' :$_->[1] )) . 
            $_->[0]->_str_flat_arg_name() }  
            @{$self->{terms}});
    #$s = "(" . $s . ")"; # parentheses should not be necessary (no ambiguities expected)
    $s =~ s/\+\-/-/g;  # handle minus signs
    return $s
}

sub coefficient {
    my ($self, $var) = @_;
    my $c = 0;
    for my $i (@{$self->{terms}}) {
	$c += $i->[1] * $i->[0]->coefficient($var);
    }
    return $c;
}    

sub substitute_time_variable {
    my ($self, $this, $by) = @_;
    return TimeLinearCombination->new(map { ($_->[0]->substitute_time_variable($this, $by) , $_->[1]) } 
                                       @{$self->{terms}} );
}



sub _expand_lin_combs {
    my $self = $_[0];
    my @x=();
    for my $i (@{$self->{terms}}) {
        if (ref($i->[0]) and $i->[0]->isa("TimeLinearCombination")) {
	    $i->[0]->_expand_lin_combs();
            my $f = $i->[1];
            for my $j (@{$i->[0]->{terms}}) {
	        push(@x, [$j->[0], $f*$j->[1]]);
	    }
	}
	else {
	     push(@x, $i);
	}
    }
    $self->{terms} =  \@x ;
}

sub _collect {
    my $self = $_[0];
    my %x=();
    for my $i (@{$self->{terms}}) {
        my $r = $i->[0];
        if ($x{$r}) {
              $x{$r}[1] += $i->[1];
	}
	else {
              $x{$r} = [ $i->[0], $i->[1] ];
	}
    }
    $self->{terms} = [ grep { $_->[1]!=0 } values(%x) ];
    if (!@{$self->{terms}}) {
        # Each empty linear combination is replaced by $t_zero. 
        # Note: defined(...) necessary because otherwise $t_zero could never
	# get defined.
        $_[0] = $t_zero if defined($t_zero);
    }
    elsif ( @{$self->{terms}}==1 and $self->{terms}[0][1]==1 ) {
        # Each linear combination consisting of one term with coefficient 1 is replaced by this term.
        $_[0] = $self->{terms}[0][0];
    }
    else {
        $_[0] = $self->_register();
    }

}    

sub simplify {
    $_[0]->_expand_lin_combs();
    $_[0]->_collect();
}    

sub get_register_key {
    my $self = shift;
    return 'L'. join('|', map { TimeExpression::get_adr($_->[0]) . ':' . $_->[1] } @{$self->{terms}});
}

package EmptyTimeLinearCombination;
our @ISA = qw( TimeLinearCombination );

sub new {
    my $class = shift;
    my $self = TimeLinearCombination->new();
    return bless $self, $class;
}

sub str { return '0' }

sub _str_flat { return '0' }
sub _str_flat_arg_name { return '0' }




 
