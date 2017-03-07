%assemble the time-dependent pde, space->Galerkin, time->theta
function [A,B,E,P,Q,S,x,u]=assem_heat_1D(n,dt,nt,B,heat_in)
h = 1/(n-1);
x = 0:h:1;
%
% A refers to -d^2/dx^2, the standard diffusion
if(strcmpi(heat_in,'left'))
bnd=B(1); %Dirichlet @left
else
bnd=0.0;  %Dirichlet @left
end       
A = 1/h*spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);          %Stiffness
A(1,:)=0;A(1,1)=1; A(end,end)=1/h;                      %Dirichlet&Neumann
B = h*B; B(1)=bnd ; B(end)=1/2*B(end);                  %Dirichlet&Neumann
%
E = h/6*spdiags(ones(n,1)*[1 4 1],-1:1,n,n);            %Mass
E(1,1)=h/3;E(end,end)=h/3;                              %Neumann,see below 
%
% A refers to +d^2/dx^2, heat equation is not effective mass equation!
A=-A;
%theta =0;  % explicit Euler 
theta =0.5; % trapezoidal
%theta=1;   % implicit(backward) Euler
%
P = E   -theta *dt*A;                         %Dirichlet BC of E irrelevant                    
Q = E+(1-theta)*dt*A;                                    %overwritten by PQ
S = dt*B;
P(1,:)=0;P(1,1)=1;Q(1,:)=0;S(1)=bnd;                       %Dirichlet @left
u=NaN*zeros(n,nt);u(:,1)=0;                                  %initial state
end