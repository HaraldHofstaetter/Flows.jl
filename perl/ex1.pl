use strict;
use TimeExpression; 
my $t = TimeVariable->new('t'); 
my $tau = TimeVariable->new('tau', '\tau'); 
my $tt1 = TimeLinearCombination->new($t, 2, $tau, 3); 
my $tt2 = 2*$t+3*$tau;     #overloaded operators
my $tt3 = TimeLinearCombination->new($tt1, 2, $t, -1); 
my $tt4 = $tt1 - $tt2;
print $tt1->str() . "\n";
print $tt1->latex() . "\n";
print $tt2->str . "\n";
print $tt3->str . "\n";
print $tt4->str . "\n";
