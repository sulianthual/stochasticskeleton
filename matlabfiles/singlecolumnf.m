function f=singlecolumnf(fm,psim)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% computes f(x,yk) in physical space from fm(x,m) the projection on psim functions
% 
% Inputs:
% - fm(nx,m)
% - psim(nyk,nym)
% 
% Outputs
% - f(x,nyk) in physical space
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(fm);
nx=term(1);
nym=term(2);
term=size(psim);% beware they can be more psim than nyk, that are not to be used
nyk=term(1);
f=zeros(nx,nyk);
for ik=1:nyk;
for im=1:nym;
f(:,ik)=f(:,ik)+fm(:,im)*psim(ik,im);
end
end


