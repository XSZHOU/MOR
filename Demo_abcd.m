clear all;clc
n = 501;  %use nomalized length in assem
dt =0.05; nt=100;
heat_in='gaussian'; heat_out='right';
%--------------------------------------------------------------------------
% part(a)(b)
[b,C]=heat_source(n,heat_in,heat_out);
%
[A,B,E,P,Q,S,x,u]=assem_heat_1D(n,dt,nt,b,heat_in);
u_steady = -A\B;
%
figure(1),hold on
title('heat equation space-time discretization: Galerkin based \theta MOL')
plot(x,u_steady,'*b'),xlim([0 x(end)])
for i=1:nt
y=P\(Q*u(:,i)+S);
plot(x,u(:,i),'-.r');
drawnow
u(:,i+1)=y;
end
p1=plot([2 2],[1 2],'*b');
p2=plot([2 2],[1 2],'-.r');
legend([p1 p2],'steady state','u(t)')
%--------------------------------------------------------------------------
% part(c)
nt=100;
[A,B,E,P,Q,S,x,u]=assem_heat_1D(n,dt,nt,b,heat_in);
figure(2),hold on,xlim([0 dt*nt])
title('output of the dynamic system: time domain')
yy=zeros(1,nt);t_now=0;
for i=1:nt-1   
plot(t_now,yy(i),'*b'),drawnow
y=P\(Q*u(:,i)+S);
yy(i+1)=C*y;
u(:,i+1)=y;
t_now=t_now+dt;
end
%--------------------------------------------------------------------------
% part (d)
w=logspace(-8,4,n);
H=PlotFreqResp(w,E,A,B,C);
figure(3), grid on
title('frequency response of the dynamic system')
p1=semilogx(w,real(H),'-.r','linewidth',2);
hold on;
p2=semilogx(w,imag(H),'-.b','linewidth',2);
legend([p1 p2],'real','imaginary')