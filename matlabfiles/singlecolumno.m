function o=singlecolumno(K,Rm,psim)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% o= - K*psi0/sqrt2 - sum Rm/4*[psim+1/sqrt(m+1)+ psim-1/sqrt(m)]
%
% Inputs:
% - K(x) and Rm(x,m) in spectral space (cylinder functions)
% - psim(y,m) values at locations of physical grid
%
% Outputs
% - o on physical grid (x,y)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(K);
nx=term(1);
term=size(psim);
nyk=term(1);
nym=term(2);
o=zeros(nx,nyk);
%
for im=1:nym
o=o+singlecolumnom(K,Rm,im-1)*psim(:,im)';
end
