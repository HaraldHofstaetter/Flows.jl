
#$Rev: 38 $
#$Author: hofi $
#$Date: 2012-12-05 11:25:47 +0100 (Wed, 05 Dec 2012) $

use strict;

use TimeExpression;

package Function;


sub new {
    my $class = shift;
    my $self = { };
    $self->{name} = shift;
    $self->{latex} = shift; #optional
    return bless $self, $class;
}

sub str {
    my $self = shift;
    return $self->{name};
}    

sub latex {
    my $self = shift;
    return $self->{latex} || $self->{name};
}


package AutonomousFunction;
our @ISA = qw( Function );

sub new {
    my $class = shift;
    my $self = Function->new( @_ );
    return bless $self, $class;
}


package SpaceExpression;

our $x_zero = EmptySpaceLinearCombination->new();

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
    return SpaceLinearCombination->new($self, 1, $other, 1);
}    

sub minus {
    my ($self, $other) = @_;
    return SpaceLinearCombination->new($self, 1, $other, -1);
}    

sub times {
    my ($self, $other, $swap) = @_;
    return SpaceLinearCombination->new($self, $other);
}    

our %_expressions;  
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
    print "space_expressions:\n";
    for my $i (sort { $a->[1] <=> $b->[1] } values %_expressions){
        print '#' . $i->[1]. "\t". $i->[0]->_str_flat() ."\n";
    }   
}

sub _str_flat_arg_name {
    my $self = shift;
    return '#' . $_xxx{$self};
}   

sub substitute {
    my $self = shift;
    my $this = shift;
    if ($this->isa("Function")) {
        my $by = shift;
        if ($by eq 0) {
           return $self->substitute_function_by_zero($this);
        }
        elsif (ref($by) && $by->isa("SpaceExpression")) {
           return $self->substitute_function_by_expression($this, $by, @_);
        }   
        else {
           return $self->substitute_function_by_function($this, $by, @_);
        };
    }
    elsif ($this->isa("SpaceVariable")) {
        return $self->substitute_space_variable($this, @_);
    }	
    elsif ($this->isa("TimeVariable")) {
        return $self->substitute_time_variable($this, @_);
    }	
    else {
        die "Function, SpaceVariable, or TimeVariable expected, stopped";
    };
}    


package SpaceVariable;
our @ISA = qw( SpaceExpression );

sub new {
    my $class = shift;
    my $self = SpaceExpression->new();
    $self->{name} = shift;
    $self->{latex} = shift; #optional
    return bless $self, $class;
}

sub str {
    my $self = shift;
    return $self->{name};
}    

sub latex {
    my $self = shift;
    return $self->{latex} || $self->{name};
}

sub differential {
    my ($self, $with_respect_to, $applied_to) = @_;
    return ($self==$with_respect_to ? $applied_to : $x_zero);
}

sub t_derivative {
    return $x_zero; #Time derivative of a space variable ...
}

sub substitute_space_variable {
    my ($self, $this, $by) = @_;
    return ($self==$this ? $by : $self);
}

sub substitute_time_variable {
    return $_[0];
}

sub substitute_function_by_expression {
    return $_[0];
}

sub substitute_function_by_function {
    return $_[0];
}

sub substitute_function_by_zero {
    return $_[0];
}

sub expand {
    return $_[0];
}

sub merge_flows {
    return $_[0];
}

sub commute_FE2DEF {
    return $_[0];
}

sub commute_DEF2FE {
    return $_[0];
}

sub reduce_order {
    return $_[0];
}

sub assume_linear {
    return $_[0];
}

sub _str_flat_arg_name {
    my $self = shift;
    return $self->{name};
}    



package SpaceLinearCombination;
our @ISA = qw( SpaceExpression );

use Scalar::Util qw( looks_like_number );


sub new {
    my $class = shift;
    my $self = SpaceExpression->new();
    my %x=();
    while (@_) {
        my $expr = shift;
	die ">>> ". $expr ."SpaceExpression expected, stopped" unless ref($expr) and $expr->isa("SpaceExpression"); 
        my $coeff = shift;
        die "Number expected, stopped" unless looks_like_number($coeff);
        if ($expr->isa("SpaceLinearCombination")) {
	    if ($coeff) {
                for my $j (@{$expr->{terms}}) {
                    my ($expr1, $coeff1) = @{$j};
		    my $r = SpaceExpression::get_adr($expr1);
                    if (exists($x{$r})) {
                          $x{$r}[1] +=  ($coeff * $coeff1);
	            }
	            else {
                          $x{$r} = [ $expr1 , $coeff * $coeff1 ];
	            }
	        }
	    }
	}
	else {
	     my $r = SpaceExpression::get_adr($expr);
             if (exists($x{$r})) {
                   $x{$r}[1] += $coeff;
	     }
	     else {
                   $x{$r} = [ $expr, $coeff ];
	     }
	}
    }

    $self->{terms} = [ grep { $_->[1]!=0 } values(%x) ];

    if (!@{$self->{terms}}) {
        # Each empty linear combination is replaced by $x_zero. 
        # Note: defined(...) necessary because otherwise $x_zero could never
	# get defined.
        return $x_zero if defined($x_zero);
        bless $self, $class;
        return $self;
    }
    elsif ( @{$self->{terms}}==1 and $self->{terms}[0][1]==1 ) {
        # Each linear combination consisting of one term with coefficient 1 is replaced by this term.
        return $self->{terms}[0][0];
    }
    else {
        bless $self, $class;
        return $self->_register();
    }   
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

sub differential {
    my ($self, $with_respect_to, $applied_to) = @_;
    return SpaceLinearCombination->new(map { ($_->[0]->differential($with_respect_to, $applied_to) , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub t_derivative {
    my ($self, $with_respect_to, $flag) = @_;
    my $r = SpaceLinearCombination->new(map { ($_->[0]->t_derivative($with_respect_to, $flag) , $_->[1]) } 
                                       @{$self->{terms}} );
    return $r;
}

sub substitute_time_variable {
    my ($self, $this, $by) = @_;
    return SpaceLinearCombination->new(map { ($_->[0]->substitute_time_variable($this, $by) , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub substitute_space_variable {
    my ($self, $this, $by) = @_;
    return SpaceLinearCombination->new(map { ($_->[0]->substitute_space_variable($this, $by) , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub substitute_function_by_expression {
    my ($self, $this, $by, $t, $x) = @_;
    if (ref($t) and $t->isa("SpaceVariable")) {
        $x = $t;
	$t = $TimeExpression::t_dummy;  
    }
    return SpaceLinearCombination->new(map { ($_->[0]->substitute_function_by_expression($this, $by, $t, $x) , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub substitute_function_by_function{
    my ($self, $this, $by) = @_;
    return SpaceLinearCombination->new(map { ($_->[0]->substitute_function_by_function($this, $by) , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub substitute_function_by_zero{
    my ($self, $this) = @_;
    return SpaceLinearCombination->new(map { ($_->[0]->substitute_function_by_zero($this) , $_->[1]) } 
                                       @{$self->{terms}} );
}


sub expand {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->expand() , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub merge_flows {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->merge_flows() , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub commute_FE2DEF {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->commute_FE2DEF() , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub commute_DEF2FE {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->commute_DEF2FE() , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub reduce_order {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->reduce_order() , $_->[1]) } 
                                       @{$self->{terms}} );
}

sub assume_linear {
    my $self = shift;
    return SpaceLinearCombination->new(map { ($_->[0]->assume_linear(@_) , $_->[1]) } 
                                       @{$self->{terms}} );
}


sub get_register_key {
    my $self = shift;
    return 'L'. join('|', map { SpaceExpression::get_adr($_->[0]) . ':' . $_->[1] } @{$self->{terms}});
}


package EmptySpaceLinearCombination;
our @ISA = qw( SpaceLinearCombination );

sub new {
    my $class = shift;
    my $self = SpaceLinearCombination->new();

    return bless $self, $class;
}

sub str { return '0' }
sub latex { return '0' }

sub _str_flat { return '0' }
sub _str_flat_arg_name { return '0' }



package FunctionExpression;

use Scalar::Util qw( looks_like_number );

our @ISA = qw( SpaceExpression );

sub new {
    my $class = shift;
    my $self = SpaceExpression->new();

    $self->{ident} = shift;
    $self->{t} = shift; 
    die "TimeExpression expected, stopped" unless ref($self->{t}) and $self->{t}->isa("TimeExpression");
    $self->{x} = shift;
    die "SpaceExpression expected, stopped" unless ref($self->{x}) and $self->{x}->isa("SpaceExpression");
    
    my $dt_order = $_[0]; #optional argument: order of time derivative
    if (looks_like_number($dt_order)) { 
        die "integer number >=0 expected, stopped" unless  $dt_order>=0 and $dt_order==int($dt_order);
        $self->{dt_order} = $dt_order;
	shift;
    }
    else {
        $self->{dt_order} = 0;
    }

    for my $i ( @_ ) {
        die " SpaceExpression expected, stopped" unless ref($i) and $i->isa("SpaceExpression");
	return $x_zero if $i eq  $x_zero;
        # Note: If the FunctionExpression $self is a derivative (i.e. a multilinear map), 
        # and one argument of this map is 0 then $self itself shall be zero.
    }
    $self->{d_args} = [ @_ ];

    bless $self, $class;

    $self=$self->_register();

    return $self;
}


sub differential {
    my ($self, $with_respect_to, $applied_to) = @_;
    my $dx =  $self->{x}->differential($with_respect_to, $applied_to);
    my @d_args = @{$self->{d_args}};
    my @terms;
    $terms[@d_args] = $x_zero; # preallocate array, size: length of @d_args + 1.
    $terms[0] =  ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args, $dx);
    # Note: ref($self)->new(...) generates Function of the same type as $self.
    for (my $i=0; $i<@d_args; $i++) {
       my @d_args1 = @d_args;
       $d_args1[$i] = $d_args[$i]->differential($with_respect_to, $applied_to);
       $terms[$i+1] = ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args1);
    }   
    return SpaceLinearCombination->new( map { ( $_, 1) } @terms );
}

sub t_derivative {
    my ($self, $with_respect_to, $flag) = @_;
    my $f = $self->{t}->coefficient($with_respect_to, $flag);
    my $dx =  $self->{x}->t_derivative($with_respect_to, $flag);
    my @d_args = @{$self->{d_args}};
    my @terms;
    $terms[@d_args+1] = $x_zero; # preallocate array, size: length of @d_args + 2.
    if ($f) {
    	$terms[0] =  ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}+1, @d_args);
    }
    else {	
        $terms[0] = $x_zero;
    }
    $terms[1] =  ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args, $dx);
    # Note: ref($self)->new(...) generates Function of the same type as $self.
    for (my $i=0; $i<@d_args; $i++) {
       my @d_args1 = @d_args;
       $d_args1[$i] = $d_args[$i]->t_derivative($with_respect_to, $flag);
       $terms[$i+2] = ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args1);
    }   
    my $r = SpaceLinearCombination->new( $terms[0] ,$f, map { ( $_, 1) } @terms[1..@d_args+1] );
    return $r;
}





sub substitute_space_variable {
    my ($self, $this, $by) = @_;
    return  ref($self)->new($self->{ident}, $self->{t}, $self->{x}->substitute_space_variable($this, $by), $self->{dt_order}, 
                              map { $_->substitute_space_variable($this, $by) } @{$self->{d_args}});
}	

sub substitute_time_variable {
    my ($self, $this, $by) = @_;
    return  ref($self)->new($self->{ident}, $self->{t}->substitute_time_variable($this, $by), 
                            $self->{x}->substitute_time_variable($this, $by), $self->{dt_order}, 
                            map { $_->substitute_time_variable($this, $by) } @{$self->{d_args}});
}	


sub substitute_function_by_expression {
    my ($self, $this, $by, $t, $x) = @_;
    my $x_var = SpaceVariable->new("x_var");
    my $t_var = TimeVariable->new("t_var");
    if (ref($t) and $t->isa("SpaceVariable")) {
        $x = $t;
	$t = $TimeExpression::t_dummy;  
    }
    my $xx = $self->{x}->substitute_function_by_expression($this, $by, $t, $x);
    my @d_args = map { $_->substitute_function_by_expression($this, $by, $t, $x) } @{$self->{d_args}};    
    if ($this eq $self->{ident}) {
        my $r = $by->substitute_space_variable($x, $x_var)
	            ->substitute_time_variable($t, $t_var);
        for (my $i=1; $i<=$self->{dt_order}; $i++) {
            $r = $r->t_derivative($t_var);
        }
        for my $i (@d_args) {
            $r = $r->differential($x_var, $i);
        }
        return $r->substitute_space_variable($x_var, $xx)
	         ->substitute_time_variable($t_var, $self->{t});
    }
    else {
        return  ref($self)->new($self->{ident}, $self->{t}, 
	                        $self->{x}->substitute_function_by_expression($this, $by, $t, $x), 
				$self->{dt_order}, 
				@d_args);
    }
}	

sub substitute_function_by_function {
    my ($self, $this, $by) = @_;
    if ($self->{ident}==$this) {
        return ref($self)->new($by, $self->{t}, $self->{x}->substitute_function_by_function($this, $by), $self->{dt_order}, 
                              map { $_->substitute_function_by_function($this, $by) } @{$self->{d_args}});
    }
    else {
        return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->substitute_function_by_function($this, $by), 
                               $self->{dt_order}, 
                               map { $_->substitute_function_by_function($this, $by) } @{$self->{d_args}});
    }
}	

sub substitute_function_by_zero {
    my ($self, $this) = @_;
    if ($self->{ident}==$this) {
        return $x_zero;
    }
    else {
        return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->substitute_function_by_zero($this), 
                               $self->{dt_order}, 
                               map { $_->substitute_function_by_zero($this) } @{$self->{d_args}});
    }
}	

sub _ddd {
   my ($self, $k) = @_;
   return $self if $k>=@{$self->{d_args}};
   return $self->_ddd($k+1) unless $self->{d_args}[$k]->isa("SpaceLinearCombination");
   my @d_args = @{$self->{d_args}};
   my $r = SpaceLinearCombination->new( map{ ( ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order},
      #@d_args[0..$k-1], $_->[0], @d_args[$k+1..(@d_args-1)])->_ddd($k), $_->[1]) } @{$self->{d_args}[$k]{terms}} );
      @d_args[0..$k-1], $_->[0], @d_args[$k+1..(@d_args-1)])->_ddd(0), $_->[1]) } @{$self->{d_args}[$k]{terms}} );
      #TODO: understand why here _ddd(0) is necessary ...
   return $r;   
}

sub expand {
   my $self = shift;
   my $x = ref($self)->new($self->{ident}, $self->{t}, $self->{x}->expand(), $self->{dt_order}, 
                              map { $_->expand() } @{$self->{d_args}});
   my $r = $x->_ddd(0);
   return $r;
}

sub merge_flows {
   my $self = shift;
   return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->merge_flows(), $self->{dt_order}, 
                              map { $_->merge_flows() } @{$self->{d_args}});
}

sub commute_FE2DEF {
   my $self = shift;
   return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->commute_FE2DEF(), $self->{dt_order}, 
                              map { $_->commute_FE2DEF() } @{$self->{d_args}});
}

sub commute_DEF2FE {
   my $self = shift;
   return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->commute_DEF2FE(), $self->{dt_order}, 
                              map { $_->commute_DEF2FE() } @{$self->{d_args}});
}

sub reduce_order {
   my $self = shift;
   return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->reduce_order(), $self->{dt_order}, 
                              map { $_->reduce_order() } @{$self->{d_args}});
}

sub assume_linear {
    my $self = shift;
    if (@{$self->{d_args}}>=2) {
        return $x_zero;
    }
    elsif ((@{$self->{d_args}}==1) && (grep { $self->{ident}==$_ } @_)) {
        return ref($self)->new($self->{ident}, $self->{t}, $self->{d_args}[0]->assume_linear(@_),  $self->{dt_order});
    }
    else {
        return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->assume_linear(@_), $self->{dt_order}, 
                               map { $_->assume_linear(@_) } @{$self->{d_args}});

    }
}



sub str {
    my $self = shift;
    my $s = $self->{ident}->str(); 
    my $n1 = $self->{dt_order};
    my $n2 = @{$self->{d_args}};
    if ($n1>0 or $n2>0) {
       $s .= '{' . $n1 . ',' . $n2 . '}';
    }
    $s .= '[' . $self->{t}->str().',' . $self->{x}->str() . ']';
    if ($n2>0) {
        $s .=  '('.  join(',', map { $_->str() } @{$self->{d_args}} ) .')';
    }
    return $s
}

sub latex {
    my $self = shift;
    my $n1 = $self->{dt_order};
    my $n2 = @{$self->{d_args}};
    my $s = '';
    if ($n1>=1) {
       $s .= '\partial_{1}';
    }   
    if ($n1>=2) {
       $s .= '^{' . $n1 . '}';
    }   
    if ($n2>=1) {
       $s .= '\partial_{2}';
    }   
    if ($n2>=2) {
       $s .= '^{' . $n2 . '}';
    }   
    $s .= $self->{ident}->latex(); 
    $s .= '(' . $self->{t}->latex().',' . $self->{x}->latex() . ')';
    if ($n2==1) { 
        if ($self->{d_args}[0]->isa("SpaceLinearCombination")) {
             $s .= '\cdot(' . $self->{d_args}[0]->latex() . ')';
        }
        else {
             $s .= '\cdot ' . $self->{d_args}[0]->latex();
        }
    }
    elsif ($n2>=2) {
        $s .=  '('.  join(',', map { $_->latex() } @{$self->{d_args}} ) .')';
    }
    return $s
}



sub _str_flat {
    my $self = shift;
    my $s = $self->{ident}->str(); 
    my $n1 = $self->{dt_order};
    my $n2 = @{$self->{d_args}};
    if ($n1>0 or $n2>0) {
       $s .= '{' . $n1 . ',' . $n2 . '}';
    }
    $s .= '[' . $self->{t}->_str_flat_arg_name().',' . $self->{x}->_str_flat_arg_name() . ']'; 
    if ($n2>0) {
        $s .=  '('.  join(',', map { $_->_str_flat_arg_name() } @{$self->{d_args}} ) .')';
    }
    return $s
}


sub get_register_key {
    my $self = shift;
    return 'F'. join(':', map { SpaceExpression::get_adr($_) } (  $self->{ident}, $self->{t}, $self->{x} )) . '|' .
               $self->{dt_order} . '|' .
               join(':', sort map { SpaceExpression::get_adr($_) } @{$self->{d_args}} );
}


package AutonomousFunctionExpression;

our @ISA = qw( FunctionExpression );

sub new {
    my $class = shift;
    my $ident = shift;
    shift if ref($_[0]) and $_[0]->isa("TimeExpression"); #ignore time if present
    my $self = $class->SUPER::new($ident, $TimeExpression::t_zero, @_);
    return $x_zero if $self eq $x_zero;
    die "AutonomousFunction expected, stopped" unless $self->{ident}->isa("AutonomousFunction");
    return  $x_zero if $self->{dt_order}>0; #Time derivative of a time-independent function ...
    return bless $self, $class;
}

sub str {
    my $self = shift;
    my $s = $self->SUPER::str();
    $s =~ s/\[0,/[/;
    $s =~ s/\{0,/{/;
    return $s
}

sub latex {
    my $self = shift;
    my $s = $self->{ident}->latex(); 
    my $n = @{$self->{d_args}};
    if ($n==1) {
       $s .= "'";
    }   
    elsif ($n==2) {
       $s .= "''";
    }   
    elsif ($n==3) {
       $s .= "'''";
    }   
    elsif ($n>=4) {
       $s .= '^{(' . $n . ')}';
    }
    $s .= '(' . $self->{x}->latex() . ')';
    if ($n==1) { 
        if ($self->{d_args}[0]->isa("SpaceLinearCombination")) {
             $s .= '\cdot(' . $self->{dargs}[0]->latex() . ')';
        }
        else {
             $s .= '\cdot ' . $self->{d_args}[0]->latex();
        }
    }
    elsif ($n>=2) {
        $s .=  '('.  join(',', map { $_->latex() } @{$self->{d_args}} ) .')';
    }
    return $s
}

sub _str_flat {
    my $self = shift;
    my $s = $self->SUPER::_str_flat();
    $s =~ s/\[0,/[/;
    $s =~ s/\{0,/{/;
    return $s
}


sub get_register_key {
    my $self = shift;
    return 'A' . join(':', map { SpaceExpression::get_adr($_) } (  $self->{ident}, $self->{x} )) . '|' .
                 join(':', sort map { SpaceExpression::get_adr($_) } @{$self->{d_args}} );
}

sub commute_FE2DEF {
   my $self = shift;
   if ( $self->{x}->isa("FlowExpression") && ($self->{ident}==$self->{x}{ident}) &&
       (@{$self->{d_args}}==0) && (@{$self->{x}{d_args}}==0)) {
       my $xx = $self->{x}{x}->commute_FE2DEF();
       return ref($self->{x})->new($self->{ident}, $self->{x}{t}, $xx, $self->{dt_order},
                                  ref($self)->new($self->{ident}, $self->{t}, $xx));
   }
   elsif ( $self->{x}->isa("FlowExpression") && ($self->{ident}==$self->{x}{ident}) &&            
           (@{$self->{x}{d_args}}==0) && (@{$self->{d_args}}==1) && 
           $self->{d_args}[0]->isa("FlowExpression") && ($self->{d_args}[0]->{ident}==$self->{ident}) &&
           (@{$self->{d_args}[0]{d_args}}==1) &&  ($self->{x}{x}==$self->{d_args}[0]{x}) ) {
           my $xx = $self->{x}{x}->commute_FE2DEF();
           my $v = $self->{d_args}[0]{d_args}[0]->commute_FE2DEF();
           my $r = (ref($self->{x})->new($self->{ident}, $self->{x}{t}, $xx, $self->{dt_order},
                                  ref($self)->new($self->{ident}, $self->{t}, $xx, $v))
                  +ref($self->{x})->new($self->{ident}, $self->{x}{t}, $xx, $self->{dt_order}, 
                                  ref($self)->new($self->{ident}, $self->{t}, $xx), $v));

           return $r;
   }
   else {
       return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->commute_FE2DEF(), $self->{dt_order}, 
                              map { $_->commute_FE2DEF() } @{$self->{d_args}});
   }
}



package FlowExpression;

our @ISA = qw( FunctionExpression );

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $x_zero if $self eq $x_zero;
    die "AutonomousFunction expected, stopped" unless $self->{ident}->isa("AutonomousFunction");
    die "time derivative for FlowExpression constructor not supported, stopped" unless $self->{dt_order}==0;
    return bless $self, $class;
}

sub str {
    my $self = shift;
    return "E_". $self->SUPER::str();
}


sub latex {
    my $self = shift;
    my $n1 = $self->{dt_order};
    my $n2 = @{$self->{d_args}};
    my $s = '';
    if ($n1>=1) {
       $s .= '\partial_{1}';
    }   
    if ($n1>=2) {
       $s .= '^{' . $n1 . '}';
    }   
    if ($n2>=1) {
       $s .= '\partial_{2}';
    }   
    if ($n2>=2) {
       $s .= '^{' . $n2 . '}';
    }   
    #$s .= '\mathcal{E}_{' . $self->{ident}->latex() . '}' ; 
    $s .= '\niceE_{' . $self->{ident}->latex() . '}' ; 
    $s .= '(' . $self->{t}->latex().',' . $self->{x}->latex() . ')';
    if ($n2==1) { 
        if ($self->{d_args}[0]->isa("SpaceLinearCombination")) {
             $s .= '\cdot(' . $self->{d_args}[0]->latex() . ')';
        }
        else {
             $s .= '\cdot ' . $self->{d_args}[0]->latex();
        }
    }
    elsif ($n2>=2) {
        $s .=  '('.  join(',', map { $_->latex() } @{$self->{d_args}} ) .')';
    }
    return $s
}





sub _str_flat {
    my $self = shift;
    return "E_". $self->SUPER::_str_flat();
}

sub get_register_key {
    my $self = shift;
    return 'E' . join(':', map { SpaceExpression::get_adr($_) } (  $self->{ident}, $self->{t}, $self->{x} )) . '|' .
               $self->{dt_order} . '|' .
              join(':', sort map { SpaceExpression::get_adr($_) } @{$self->{d_args}} );

}

sub merge_flows {
   my $self = shift;
   my $x =  $self->{x}->merge_flows();
   if (($self->{t} eq $TimeExpression::t_zero) and (@{$self->{d_args}}==0)) {
      return $x;
   }
   elsif ((@{$self->{d_args}}==0) and ref($x) and $x->isa("FlowExpression") 
      and ($self->{ident} eq $x->{ident}) and (@{$x->{d_args}}==0)) {
      my $t = $self->{t} + $x->{t};
      return $x->{x} if $t eq  $TimeExpression::t_zero;
      return ref($self)->new($self->{ident}, $t, $x->{x});
   }
   else {
      return ref($self)->new($self->{ident}, $self->{t}, $x, 
                             map { $_->merge_flows() } @{$self->{d_args}});
   }
}

sub commute_DEF2FE {
   my $self = shift;
   if ((@{$self->{d_args}}==1) && ($self->{ident}==$self->{d_args}[0]{ident}) &&
       $self->{d_args}[0]->isa("AutonomousFunctionExpression") &&
       ($self->{x}==$self->{d_args}[0]{x}) && (@{$self->{d_args}[0]{d_args}}==0)) {
       my $r =  $self->{d_args}[0];
       return ref($r)->new($r->{ident}, $r->{t}, 
                    ref($self)->new($self->{ident}, $self->{t}, $self->{x}->commute_DEF2FE(), $self->{dt_order}),
                    $r->{dt_order});
   }                 
   elsif ((@{$self->{d_args}}==1) && ($self->{ident}==$self->{d_args}[0]{ident}) &&
       $self->{d_args}[0]->isa("AutonomousFunctionExpression") &&
       ($self->{x}==$self->{d_args}[0]{x}) && (@{$self->{d_args}[0]{d_args}}==1)) {
       my $r = $self->{d_args}[0];
       my $v = $self->{d_args}[0]{d_args}[0]->commute_DEF2FE();
       my $x = $self->{x}->commute_DEF2FE();
       return ref($r)->new($r->{ident}, $r->{t}, 
                    ref($self)->new($self->{ident}, $self->{t}, $x, $self->{dt_order}),
                    ref($self)->new($self->{ident}, $self->{t}, $x, $self->{dt_order}, $v))
                    -ref($self)->new($self->{ident}, $self->{t}, $x, $self->{dt_order},
                             ref($r)->new($r->{ident}, $r->{t}, $x, $r->{dt_order}), $v);
   }
   else {
       return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->commute_DEF2FE(), $self->{dt_order}, 
                              map { $_->commute_DEF2FE() } @{$self->{d_args}});
   }
}

my $x_var = SpaceVariable->new("x_var");
my $ro_F = AutonomousFunction->new("F");
my $ro_u = SpaceVariable->new("u");
my $ro_t = TimeVariable->new("t");
my $ro_m = 5;
my @ro_vars = map {SpaceVariable->new("v" . $_); } (0 .. $ro_m-1);
my @ro_expr;

for (my $i=0; $i<$ro_m; $i++) {
    $ro_expr[$i] = FlowExpression->new($ro_F, $ro_t, $ro_u, @ro_vars[0..($i-1)] )->t_derivative($ro_t, 0) 
           -FlowExpression->new($ro_F, $ro_t, $ro_u, @ro_vars[0..($i-1)] )->t_derivative($ro_t, 1)
           +FlowExpression->new($ro_F, $ro_t, $ro_u, AutonomousFunctionExpression->new($ro_F, $ro_t, $ro_u), @ro_vars[0..($i-1)] );
}


sub reduce_order {
    my $self = shift;
    my $Fu = AutonomousFunctionExpression->new($self->{ident}, $self->{t}, $self->{x});
    my $j = -1;
    my $m = @{$self->{d_args}};
    for (my $i=0; $i<$m; $i++) {
        if ($self->{d_args}[$i]==$Fu) {
            $j=$i;
            last;
        }
    }    
    if (($j>=0) && ($m<=$ro_m)) {
        my @other_args = (@{$self->{d_args}}[0..$j-1], @{$self->{d_args}}[$j+1..$m-1]);
        my $r = $ro_expr[$m-1]
                ->substitute($ro_F, $self->{ident})
                ->substitute($ro_t, $self->{t})
                ->substitute($ro_u, $self->{x}->reduce_order());
        for (my $i=0; $i<$m-1; $i++) {
           $r = $r->substitute($ro_vars[$i], $other_args[$i]->reduce_order());
        }
        return $r->reduce_order();
    }
    else {
         return $self->SUPER::reduce_order();
    }
}

#my $x_var = SpaceVariable->new("x_var");
sub t_derivative {
    my ($self, $with_respect_to, $flag) = @_;
    my $f = $self->{t}->coefficient($with_respect_to);
    my $dx =  $self->{x}->t_derivative($with_respect_to, $flag);
    my @d_args = @{$self->{d_args}};
    my @terms;
    $terms[@d_args+1] = $x_zero; # preallocate array, size: length of @d_args + 2.
    if ($f) {
        my $T;
        if (!$flag) {
            $T = AutonomousFunctionExpression->new($self->{ident},
                  FlowExpression->new($self->{ident}, $self->{t}, $x_var));
        }
        else {
            $T = FlowExpression->new($self->{ident}, $self->{t}, $x_var, 
                  AutonomousFunctionExpression->new($self->{ident}, $x_var));
        }
        for my $d_arg (@d_args) {
            $T = $T->differential($x_var, $d_arg);        
        }
        $terms[0] = $T->substitute_space_variable($x_var, $self->{x});
    }
    else {	
        $terms[0] = $x_zero;
    }
    $terms[1] =  ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args, $dx);
    # Note: ref($self)->new(...) generates Function of the same type as $self.
    for (my $i=0; $i<@d_args; $i++) {
       my @d_args1 = @d_args;
       $d_args1[$i] = $d_args[$i]->t_derivative($with_respect_to, $flag);
       $terms[$i+2] = ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args1);
    }   
    $terms[1] =  ref($self)->new($self->{ident}, $self->{t}, $self->{x}, $self->{dt_order}, @d_args, $dx);
    my $r = SpaceLinearCombination->new( $terms[0] ,$f, map { ( $_, 1) } @terms[1..@d_args+1] );
    return $r;
}


sub substitute_function_by_expression {
    my ($self, $this, $by, $t, $x) = @_;
    if (ref($t) and $t->isa("SpaceVariable")) {
        $x = $t;
	$t = $TimeExpression::t_dummy;  
    }
    my @d_args = map { $_->substitute_function_by_expression($this, $by, $t, $x) } @{$self->{d_args}};    
    return  ref($self)->new($self->{ident}, $self->{t}, 
                            $self->{x}->substitute_function_by_expression($this, $by, $t, $x), 
          		    $self->{dt_order}, 
			    @d_args);
}

sub substitute_function_by_zero {
    my ($self, $this) = @_;
    if ($self->{ident}==$this) {
        if (@{$self->{d_args}}==0) {
            return $self->{x};
        }
        elsif (@{$self->{d_args}}==1) {     
            return $self->{d_args}[0];
        }
        else {
            return $x_zero;
        }
    }
    else {
        return ref($self)->new($self->{ident}, $self->{t}, $self->{x}->substitute_function_by_zero($this), 
                               $self->{dt_order}, 
                               map { $_->substitute_function_by_zero($this) } @{$self->{d_args}});
    }
}	

1;

