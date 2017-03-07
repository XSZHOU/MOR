function H = PlotFreqResp(w,E,A,B,C)
%plot the frequency response of the dynamic system, c.f. slide_06<-2014_MOR
%  E*sX(s)=A*X(s)+B*U(s)
%  Y(s)=C'*X(s)
% we can derive the transfer function:
%  X(s)= (E*s-A)^-1 *BU(s)
%  Y(s)= H(s)U(s)= C'*X(s)
% so we have: H=Y/U
%  H(s)= C'*(E*s-A)^-1 *B

s=1i*w;
H=zeros(length(s),1);
for j=1:length(s)
H(j)= C*((E*s(j)-A)\B);  %notice here A refer to d^2/dx2
end

end