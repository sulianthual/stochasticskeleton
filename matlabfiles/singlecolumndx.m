function b=singlecolumndx(a,x)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% compute zonal derivative using Fourier series
%  input=u(x,y) periodic and grid x
% output=dxu(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ii=complex(0,1); % the complex number
%k=2*pi*k; % not necessary if uses fftkspe.m
pass=size(a); nx=pass(1); ny=pass(2);
dx=x(2)-x(1); kg=fftkspe(nx,dx)/(2*pi);
b=zeros(nx,ny);
for ik=1:ny
af=fftshift(fft(a(:,ik)));
af=2*pi*ii*af.*kg';
b(:,ik)=ifft(ifftshift(af));
end
b=real(b);
