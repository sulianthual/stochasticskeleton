function vm=singlecolumnvm(Rm,Sm,dx,im)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% v= sum [ dxRm +sqrt(m+1)*Sm+1 - sqrt(m)*Sm-1 ]* psim/sqrt(2)/(2m+1)
%
% Inputs:
% - Rm(x,m) in spectral space (cylinder functions)
% - Sm(x,m) in spectral space
% - dx zonal timestep to compute zonal derivative
%  - im desired contribution
%
% Outputs
% - vm(x) contribution on hermite im : BEWARE im is not the indice of psim,
%   (that starts from 1) but of the Hermites (thats starts from zero)
%  Beware also Sm starts at m=0 like psim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(Rm);
nx=term(1);
nrm=term(2);
term=size(Sm);
nym=term(2);
vm=zeros(nx,1);
kf=fftkspe(nx,dx); % wavenumber=2*pi*f (nk=nx everywhere)
ii=complex(0,1);
%
if im==0;
if nym>=2; vm(:,1)= Sm(:,2)/sqrt(2); end
else
if im<=nrm;
term=fftspe(Rm(:,im)); term=term.*(-ii*kf)'; dxRm=fftispe(term);
vm(:,1)=vm(:,1)+dxRm/sqrt(2)/(2*im+1);
end
if im>=1;
vm(:,1)=vm(:,1) - sqrt(im)*Sm(:,im)/sqrt(2)/(2*im+1);
end
if im<=nym-2;
vm(:,1)=vm(:,1) + sqrt(im+1)*Sm(:,im+2)/sqrt(2)/(2*im+1);
end
end
