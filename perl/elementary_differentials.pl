#$Rev: 27 $
#$Author: hofi $
#$Date: 2012-10-24 12:58:50 +0200 (Wed, 24 Oct 2012) $

use strict;

use TimeExpression;
use SpaceExpression;


my $t = TimeVariable->new("t");
my $x = SpaceVariable->new("x");

my $F = new AutonomousFunction('F');
my $u = FlowExpression->new($F, $t, $x);

sub uuu { $_[0] =~ s/E_F\[t,x\]//g; $_[0] =~ s/\[\]//g; return $_[0] }
for (my $i=1; $i<=11; $i++) {
   $u = $u->t_derivative($t)->expand();
   print $i . "\t" . uuu($u->str())  . "\t(" . 
     ($u->isa("SpaceLinearCombination")  ? scalar(@{$u->{terms}}) : 1 )
          ." terms)\n";
}

#TimeExpression::print_register();
#SpaceExpression::print_register();
