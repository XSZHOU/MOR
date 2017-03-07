clear all; close all;clc
n = 501;  %use nomalized length in assem
dt =0.025; nt=200;
h = 1/(n-1); t=0:dt:dt*(nt-1);
heat_in='left'; heat_out='right';
[b,C]=heat_source(n,heat_in,heat_out);
[A,B,E,P,Q,S,x,u]=assem_heat_1D(n,dt,nt,b,heat_in);
w=logspace(-8,4,n);H=PlotFreqResp(w,E,A,B,C);%full system
%--------------------------------------------------------------------------
%part(e), eigen-mode truncation: slide_07
yy=zeros(1,nt);
for i=1:nt-1
y=P\(Q*u(:,i)+S);
yy(i+1)=C*y;
u(:,i+1)=y;
end
figure(1),hold on,xlim([0 dt*nt])
title('Model truncation, reasonably accurate in time domain with small q')
p1=plot(t,yy,'*r');
%
[V,D]=eig(full(E\A));
[dd,ind]=sort(diag(D),'ascend');V=V(:,ind); D(logical(eye(size(A))))=dd;
Anew=D;
Bnew=V\(E\B);
Cnew=C*V;
y_reduc=zeros(1,nt);
q_range=(n-11):n;
for q= q_range
    y_reduc=y_reduc+Cnew(q)*Bnew(q)/D(q,q).*(exp(D(q,q)*t)-1);
end
p2=plot(t,y_reduc,'og');
legend([p1 p2],'full','reduced')
%
q_range=51:501; %use 450 q, model truncation: ineffective
Hr=PlotFreqResp(w,E(q_range,q_range),A(q_range,q_range),B(q_range),C(q_range));
%
figure(2), grid on
p1=semilogx(w,real(H),'*r');hold on;
p2=semilogx(w,imag(H),'*g');
p3=semilogx(w,real(Hr),'ob');
p4=semilogx(w,imag(Hr),'om');hold off;
title('Model tuncation, ineffective in frequency domain even for large q')
legend([p1 p2 p3 p4],'real,full','imaginary,full','real,reduced','imaginary,reduced')
%--------------------------------------------------------------------------
%part(f), transient signal
n = 501;
dt =10; nt=1000;t=0:dt:dt*(nt-1);
[B,C]=heat_source(n,heat_in,heat_out);
[A,B,E,P,Q,S,x,~]=assem_heat_1D(n,dt,nt,B,heat_in);
u=sin(0.01*t);[y,U]= PlotTransient(n,dt,nt,u,A,B,C,E);
figure(3),hold on
p1=plot(t,y,'*r');
%
[V,D]=eig(full(E\A));
[dd,ind]=sort(diag(D),'ascend');V=V(:,ind); D(logical(eye(size(A))))=dd;
y_reduc=zeros(1,nt);
Anew=D;
Bnew=V\(E\B);
Cnew=C*V;
q_range=n-11:n;
for q= q_range
    y_reduc=y_reduc+Cnew(q)*Bnew(q)/D(q,q).*(exp(D(q,q)*t)-1).*sin(t*0.01);
end
p2=plot(t,y_reduc,'o-b');legend([p1 p2],'full','reduced');
%--------------------------------------------------------------------------
%part(g) balanced truncation, not available in Matlab provided by Polito
%--------------------------------------------------------------------------
%part(h) moment matching
v_expan=[10^-2,1; 10^3,1]; % match s^th order moments at two freq-points
[Hr,Er,Ar,Br,Cr,V]=PlotMORmm(w,v_expan,E,A,B,C);
disp(['N x q: ',num2str(size(V))]);%tall and thin
hmin=min([real(H(:));imag(H(:))]);
figure(4);grid on;
semilogx(w,real(H),'*r',w,imag(H),'*g',w,real(Hr),'bo',w,imag(Hr),'mo');
hold on
plot(v_expan(:,1),[hmin hmin],'rp','MarkerSize',12);
legend('real,full','imag,full','real,reduced','imag,reduced','mm points');
title('Moment matching: two frequency points with s^{th} moment')
%
%repeat (e)
dt =0.025; nt=200;t=0:dt:dt*(nt-1);
u=ones(length(t),1);
[y,U]= PlotTransient(n,dt,nt,u,A,B,C,E);
[yr,Ur]= PlotTransient(size(V,2),dt,nt,u,Ar,Br,Cr,Er);
figure(5),hold on
p1=plot(t,y,'*r');
p2=plot(t,yr,'o-b');legend([p1 p2],'full','reduced');
%
%repeat (f)
dt =10; nt=1000;t=0:dt:dt*(nt-1);
u=sin(0.01*t);
[y,U]= PlotTransient(n,dt,nt,u,A,B,C,E);
[yr,Ur]= PlotTransient(size(V,2),dt,nt,u,Ar,Br,Cr,Er);
figure(6),hold on
p1=plot(t,y,'*r');
p2=plot(t,yr,'o-b');legend([p1 p2],'full','reduced');
