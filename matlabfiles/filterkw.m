function x=filterkw(xin,dx,dt,krange,wrange)
% Filter x(x,t) using box with k-w space
% by Sulian Thual
%
% Input:,
% - x(x,t)
% - dx spacing and dt spacing of data
% - krange=[kmin,kmax] and wrange=[wmin,wmax]
% Output:
% - xf(x,t) filtered
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
x=xin;
[nx nt]=size(x);
kg=fftkspe(nx,dx)/(2*pi); wg=fftkspe(nt,dt)/(2*pi);

kmin=krange(1); kmax=krange(2);
wmin=wrange(1); wmax=wrange(2);

% box 0:
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
%
% box 1: double symetric of box0
iwmin1=searchclosest(wg,-wmin); iwmax1=searchclosest(wg,-wmax); 
ikmin1=searchclosest(kg,-kmin); ikmax1=searchclosest(kg,-kmax); 
if iwmin1>=iwmax1; term=iwmax1; iwmax1=iwmin1; iwmin1=term; end
if ikmin1>=ikmax1; term=ikmax1; ikmax1=ikmin1; ikmin1=term; end
%
% keep cc component (do in two parts)
x0=x;
for kts=1:nt; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
xf=x*0; xf(ikmin:ikmax,:)=x(ikmin:ikmax,:); 
x=xf; xf=0;
for kts=1:nt; x(:,kts)=ifft(ifftshift(x(:,kts))); end; % zonal
for i=1:nx; x(i,:)=fftshift(fft(x(i,:)')); end % temporal
xf=x*0; xf(:,iwmin:iwmax)=x(:,iwmin:iwmax); 
x=xf; xf=0;
for i=1:nx; x(i,:)=ifft(ifftshift(x(i,:).')); end % temporal
%
x1=x0;
for kts=1:nt; x1(:,kts)=fftshift(fft(x1(:,kts))); end; % zonal
xf=x1*0; xf(ikmin1:ikmax1,:)=x1(ikmin1:ikmax1,:); 
x1=xf; xf=0;
for kts=1:nt; x1(:,kts)=ifft(ifftshift(x1(:,kts))); end; % zonal
for i=1:nx; x1(i,:)=fftshift(fft(x1(i,:)')); end % temporal
xf=x1*0; xf(:,iwmin1:iwmax1)=x1(:,iwmin1:iwmax1); 
x1=xf; xf=0;
for i=1:nx; x1(i,:)=ifft(ifftshift(x1(i,:).')); end % temporal
% add both components
x=x+x1;


