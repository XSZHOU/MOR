Demo_heat: visualization of solving the PDE
Demo_abcd: Homework parts a,b,c,d
Demo_efgh: Homework parts e,f,g,h


->for the recduced model. we cannot apply the usual
boundary condition procedure to P,Q,S.in 'PlotTransient.m'
%P(1,:)=0;P(1,1)=1;Q(1,:)=0;S(1)=B(1); should be commented

although in real space, without this line, the transient
temperature distribution would exhibit spurious solution,
in the case"heat @left, measure @right configuration",
but notice that:
-1 the reduced system do not have physical meaning in real
 space any more, we cannot apply B.C to a q x q matrix

-2 the set of equation use oberservable y=Cx as an estimator,
 in such a way, the spurious solution in real space(x) may
 get circumvented somehow. 

