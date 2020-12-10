function k=fftkspe(nx,dx)
% Computes the pulsation of Fourier transform (arranged)
%  to work with fftspe.m and fftispe.m (that use fftshift)
% by Sulian Thual
%
% Input:
% - nx: size of initial field
% - dx : step
% Output:
% - pulsations k corresponding to fftspe.m

% NOTE: to get the frequency (i.e. frequency wavenumber), do kf=k/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dk=1/(nx*dx);
kf=dk*(-nx/2:nx/2-1);% frequency wavenumber
k=2*pi*kf;% pulsation wavenumber

