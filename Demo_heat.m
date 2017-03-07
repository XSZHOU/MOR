%heat_1D, original implementation S.Z@Polito
clear all;clc
dt =0.05; nt=200;
n = 501;
h = 1/(n-1);
x = 0:h:1;
%
d=0.0;                                                  % Dirichlet val
A = 1/h*spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);          % Stiffness
A(1,:)=0;A(1,1)=1; A(end,end)=1/h;                      % Dirichlet/Neumann
B = 1/(0.1*sqrt(2*pi))*exp( -(x-0.5).^2/2/0.1^2)';      % Gaussian
B = h*B; B(1)=d ; B(end)=1/2*B(end);                    % Dirichlet/Neumann
u_steady = A\B;
%
figure(1),hold on
plot(x,u_steady,'*b'),xlim([0 x(end)])
p1=plot([2 2],[1 2],'*b');
p2=plot([2 2],[1 2],'-.r');
legend([p1 p2],'steady state','u(t)')
title('heat equation space-time discretization: Galerkin based \theta MOL')
%
M = h/6*spdiags(ones(n,1)*[1 4 1],-1:1,n,n);            % mass
M(1,1)=h/3;M(end,end)=h/3;                              % Neumann
%                                       
A =-A;   %A is +d^2x/dx^2 now
theta =0.5;   
P = M   -theta *dt*A;                               
Q = M+(1-theta)*dt*A;
S = dt*B;       
P(1,:)=0;P(1,1)=1;Q(1,:)=0;S(1)=d;                      % Dirichlet
u=NaN*zeros(n,nt);u(:,1)=0;                             % initial state
%
for i=1:nt
y=P\(Q*u(:,i)+S);
plot(x,u(:,i),'-.r');
if i<nt/5; drawnow, end
u(:,i+1)=y;
end
display(['residul err: ' num2str( norm(u(:,end)-u_steady))])
%
figure(2)
surf(x,dt:dt:(dt*nt/5),u(:,1:nt/5)','LineStyle','none'); pbaspect([1 3 1])
xlabel x; ylabel time; zlabel Temperature;
title('heat equation : temperature field marching in time ')