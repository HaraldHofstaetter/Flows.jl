use strict;
use TimeExpression; 
use SpaceExpression;
my $t = TimeVariable->new('t');
my $u = SpaceVariable->new('u');
my $A = AutonomousFunction->new('A');
my $B = AutonomousFunction->new('B');
my $Stu = FlowExpression->new($B, $t, 
                 FlowExpression->new($A, $t, $u));
my $Cu = AutonomousFunctionExpression->new($A, $u, 
                 AutonomousFunctionExpression->new($B, $u))
        -AutonomousFunctionExpression->new($B, $u, 
                 AutonomousFunctionExpression->new($B, $u));    
print $Stu->str() . "\n";
print $Cu->str() . "\n";

my $v = SpaceVariable->new('v');
my $E_Atu = FlowExpression->new($A, $t, $u);
my $E_Btv = FlowExpression->new($B, $t, $v);
my $Stu = $E_Btv->substitute($v, $E_Atu);
print $Stu->str() . "\n";

# $Cu already defined
my $Du = $Cu->substitute($B, $Cu, $u);
print $Du->str() . "\n";

my $S1tu = $Stu->t_derivative($t)
               -AutonomousFunctionExpression->new($A, $Stu)
               -AutonomousFunctionExpression->new($B, $Stu);
print $S1tu->str() . "\n"; 

my $F = AutonomousFunction->new('F');
my $Zu = AutonomousFunctionExpression->new($F,
                FlowExpression->new($F, $t, $u))
        -FlowExpression->new($F, $t, $u,
		AutonomousFunctionExpression->new($F, $u));
my $Z1uv = $Zu->differential($u, $v);
my $w = SpaceVariable->new('w');
my $Z2uvw = $Z1uv->differential($u, $w);
print $Zu->str() . "\n";
print $Z1uv->str() . "\n";
print $Z2uvw->str() . "\n";
print $Zu->latex() . "\n";
print $Z1uv->latex() . "\n";
print $Z2uvw->latex() . "\n";

