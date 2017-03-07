function [B,C] = heat_source(n,heat_in,heat_out)
%generate heat input confuguration and output coupling
% B
if(strcmpi(heat_in,'left'))
B = zeros(n,1);B(1)=1;%B(1) can be overwritten in assem_heat_1D
elseif(strcmpi(heat_in,'uniform'))
B = ones(n,1);  
elseif(strcmpi(heat_in,'Gaussian'))
h = 1/(n-1);
B = 1/(0.1*sqrt(2*pi))*exp( -((0:h:1) - 0.5).^2/2/0.1^2)'; 
else
error('error heat configuration')
end
% C
if(strcmpi(heat_out,'right'))
C=zeros(1,n);C(end)=1;  
elseif(strcmpi(heat_out,'average'))
C=ones(1,n)/n;            
else
error('error heat configuration')
end

end