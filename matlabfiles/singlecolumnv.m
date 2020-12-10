function v=singlecolumnv(Rm,Sm,psim,dx)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% v= sum [ dxRm +sqrt(m+1)*Sm+1 - sqrt(m)*Sm-1 ]* psim/sqrt(2)/(2m+1)
%
% Inputs:
% - Rm(x,m) in spectral space (cylinder functions)
% - Sm(x,m) in spectral space
% - psim(y,m) values at locations of physical grid
% - dx zonal timestep to compute zonal derivative
%
% Outputs
% - v on physical grid (x,y)
%  
% Rq:
% with psim of size nm, there are nm-2 Rossby waves, yet nm values Sm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(Rm);
nx=term(1);
term=size(psim);
nyk=term(1);
nym=term(2);
v=zeros(nx,nyk);
%
for im=1:nym
v=v+singlecolumnvm(Rm,Sm,dx,im-1)*psim(:,im)';
end

































