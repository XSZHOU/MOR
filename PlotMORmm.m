function [Hr,Er,Ar,Br,Cr,V]=PlotMORmm(w,vexpan,E,A,B,C)
% Generate moment matching Krylov vectors, use W=V
V = mmKrylov(vexpan,E,A,B);W=V; 

% Do projection
[Er,Ar,Br,Cr]=DoProj(W,V,E,A,B,C);

% MOR freq-response
Hr=PlotFreqResp(w,Er,Ar,Br,Cr);


function V = mmKrylov(vexpan,E,A,B)
% Cousebook: pp.174~180 | slide_08: s5@pp.2 && s3@pp.6
V=[];
for iw=1:size(vexpan,1) 
    M= (A-1i*vexpan(iw,1)*E); % we have E instead of I
    Akrylov=M\E; Vk =M\B;
    V=[V real(Vk)];
    V=orth(V);   
    V=[V imag(Vk)];
    V=orth(V);         
    for ik=1:vexpan(iw,2)     % match higher moments
        Vk=Akrylov*Vk;
        V=[V real(Vk)]; 
        V=orth(V);   
        V=[V imag(Vk)];  
        V=orth(V);
    end
end

function [Er,Ar,Br,Cr]= DoProj(W,V,E,A,B,C)
% name convention -> Course book: pp.173
%W - N x q
%V - N x q 
%E - N x N  Er- q x q
%A - N x N  Ar- q x q 
%B - N x 1  Br- q x 1
%C - 1 x N  Cr- 1 x q
%
Er=W'*E*V;
Ar=W'*A*V;
Br=W'*B;
Cr=C* V;