function[y,U]= PlotTransient(n,dt,nt,u,A,b,C,E)
yy=zeros(1,nt);t_now=0;
U=NaN*zeros(n,nt);U(:,1)=0;
%
%theta =0;  % explicit Euler 
theta =0.5; % trapezoidal
%theta=1;   % implicit(backward) Euler
P = E   -theta *dt*A;                         %Dirichlet BC of E irrelevant                    
Q = E+(1-theta)*dt*A;                                    %overwritten by PQ
for ii=1:nt-1
B=b*u(ii);
S = dt*B;
%P(1,:)=0;P(1,1)=1;Q(1,:)=0;S(1)=B(1);                      %Dirichlet @left
%
y=P\(Q*U(:,ii)+S);
yy(ii+1)=C*y;
U(:,ii+1)=y;
t_now=t_now+dt;
end

y=C*U;