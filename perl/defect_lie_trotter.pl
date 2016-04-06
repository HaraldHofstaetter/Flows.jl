use strict;
use TimeExpression; 
use SpaceExpression;
my $t = TimeVariable->new('t');
my $u = SpaceVariable->new('u');
my $v = SpaceVariable->new('v');
my $A = AutonomousFunction->new('A');
my $B = AutonomousFunction->new('B');
my $Av = AutonomousFunctionExpression->new($A, $v);
my $Bv = AutonomousFunctionExpression->new($B, $v);
my $E_Atu = FlowExpression->new($A, $t, $u);
my $E_Btv = FlowExpression->new($B, $t, $v);

# --- Check (1) ------------------------------------
my $Stu = $E_Btv->substitute($v, $E_Atu);
my $S1tv = $E_Btv->differential($v, $Av)
          -$Av->substitute($v, $E_Btv);
my $S1tu = $Stu->t_derivative($t)
          -$Av->substitute($v, $Stu)
          -$Bv->substitute($v, $Stu);
my $h = $S1tu - $S1tv->substitute($v, $E_Atu);
print "(1): " . $h->expand()->reduce_order()->expand()->str() . "\n";

# --- Check (2) ------------------------------------
# commutator [B,A](v):
my $C_BAv = $Bv->differential($v, $Av) - $Av->differential($v, $Bv);
my $dS1tv = AutonomousFunctionExpression->new($B, $E_Btv, $S1tv)
           +$C_BAv->substitute($v, $E_Btv);
my $h = $S1tv->t_derivative($t) - $dS1tv;
print "(2): " . $h->expand()->reduce_order()->expand()->str() . "\n";

# --- Check (3) ------------------------------------
# should hold for generic H and S
{
my $H = AutonomousFunction->new('H');
my $S = Function->new('S');
my $Stu = FunctionExpression->new($S, $t, $u);
my $S1tu = $Stu->t_derivative($t)
          -AutonomousFunctionExpression->new($H, $Stu);
my $S2tu = $S1tu->t_derivative($t)
          -AutonomousFunctionExpression->new($H, $Stu, $S1tu);
my $T = TimeVariable->new('T');
my $FtTu = FlowExpression->new($H, $T-$t, $Stu, $S1tu);
my $dFtTu = FlowExpression->new($H, $T-$t, $Stu, $S2tu)
           +FlowExpression->new($H, $T-$t, $Stu, $S1tu, $S1tu);
print "FtTu->t_derivative(t): " . $FtTu->t_derivative($t)->str() ."\n"; 
print "FtTu->t_derivative(t): " . $FtTu->t_derivative($t)->expand()->str() ."\n"; 
print "dFtTu: " . $dFtTu->str() ."\n"; 
print "dFtTu: " . $dFtTu->expand()->str() ."\n"; 
my $h = $FtTu->t_derivative($t) - $dFtTu;
print "h: " . $h->str() ."\n"; 
print "h: " . $h->expand()->str() ."\n"; 
print "(3): " . $h->expand()->reduce_order()->expand()->str() . "\n";
}

# --- Check (4) ------------------------------------
my $S2tv = $S1tv->differential($v, $Av)
          -AutonomousFunctionExpression->new($A, $E_Btv, $S1tv)
          +$C_BAv->substitute($v, $E_Btv);
my $S2tu = $S1tu->t_derivative($t)
          -AutonomousFunctionExpression->new($A, $Stu, $S1tu)
          -AutonomousFunctionExpression->new($B, $Stu, $S1tu);
my $h = $S2tu - $S2tv->substitute($v, $E_Atu);
print "(4): " . $h->expand()->reduce_order()->expand()->str() . "\n";

# --- Check (5) ------------------------------------
# commutator [B,[B,A]](v):
my $C_BBAv = $Bv->differential($v, $C_BAv) - $C_BAv->differential($v, $Bv);
# commutator [A,[B,A]](v):
my $C_ABAv = $Av->differential($v, $C_BAv) - $C_BAv->differential($v, $Av);
# derivative of commutator [B,A]'(v)(w):
my $w = SpaceVariable->new('w');
my $C_BAvw = $C_BAv->differential($v, $w);
my $dS2tv = AutonomousFunctionExpression->new($B, $E_Btv, $S2tv)
           +AutonomousFunctionExpression->new($B, $E_Btv, $S1tv, $S1tv)
           -$C_BBAv->substitute($v, $E_Btv)
           -$C_ABAv->substitute($v, $E_Btv)
           +2*$C_BAvw->substitute($v, $E_Btv)->substitute($w, $S1tv);
#print "dS2tv: " . $dS2tv->str() ."\n"; 
my $h = $S2tv->t_derivative($t) - $dS2tv;
print "(5): " . $h->expand()->reduce_order()->expand()->str() . "\n";

# --- Check (6) ------------------------------------
my $H = AutonomousFunction->new('H');
my $T = TimeVariable->new('T');
my $Xu = FlowExpression->new($H, $T-$t, $Stu,
           FlowExpression->new($B, $t, $E_Atu, 
           $C_BAv->substitute($v, $E_Atu)));
my $dXv = FlowExpression->new($H, $T-$t, $E_Btv,
           $S1tv->differential($v, $C_BAv))
         +FlowExpression->new($H, $T-$t, $E_Btv, $S1tv,
           $E_Btv->differential($v, $C_BAv))
         -FlowExpression->new($H, $T-$t, $E_Btv,
           $E_Btv->differential($v, $C_ABAv));
my $h = $Xu->t_derivative($t) - $dXv->substitute($v, $E_Atu);
print "(6): " . $h->commute_FE2DEF()->substitute($H, $Av+$Bv, $v)
                  ->expand()->reduce_order()->expand()->str() . "\n";
         
