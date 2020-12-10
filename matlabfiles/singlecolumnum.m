function um=singlecolumnum(K,Rm,im)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% u=  K*psi0/sqrt2 + sum Rm/4*[psim+1/sqrt(m+1) - psim-1/sqrt(m)]
% BEWARE OF STARTS AND ENDS OF SUM !!!!
%
% Inputs:
% - K(x) and Rm(x,m) in spectral space (cylinder functions)
% - psim(y) values at locations of physical grid
%
% Outputs
% - um(x) contribution on hermite im : BEWARE im is not the indice of psim,
%   (that starts from 1) but of the Hermites (thats starts from zero)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(K);
nx=term(1);
term=size(Rm);
nrm=term(2);
um=zeros(nx,1);
if im==0;
um(:,1)= K/sqrt(2) - Rm(:,1)/4;
else
if im>=2;
um(:,1)=um(:,1) + Rm(:,im-1)/4/sqrt(im);
end
if im<=nrm-1;
um(:,1)=um(:,1)- Rm(:,im+1)/4/sqrt(im+1);
end
end
