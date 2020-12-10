% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% Sulian Thual
% 
% runs the skeleton (graphs are separated)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving
%
%
if dorun==1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of Run
%
% Prepare Netcdf (Saved Variables): (note o and u are non-necessary for restarts)
term=zeros(nx,nts)+NaN; term(:,1)=Kini;
ncdfmakevar(fileout,'Ks',{'X','T'},term,NaN,2);
term=zeros(nx,nrm,nts)+NaN; term(:,:,1)=Rmini;
ncdfmakevar(fileout,'Rms',{'X','M','T'},term,NaN,1);
term=zeros(nx,nyk,nts)+NaN; term(:,:,1)=etaini;
ncdfmakevar(fileout,'etas',{'X','Y','T'},term,NaN,1);
term=zeros(nx,nyk,nts)+NaN; term(:,:,1)=Zini;
ncdfmakevar(fileout,'Zs',{'X','Y','T'},term,NaN,1);
term=zeros(nx,nyk,nts)+NaN; term(:,:,1)=Zini*0;
ncdfmakevar(fileout,'ntaus',{'X','Y','T'},term,NaN,1);
ncdfmakevar(fileout,'ts',{'one','T'},ts,NaN,1); term=0;
%Initialise run and increment variables
eta=etaini; Z=Zini; K=Kini; Rm=Rmini; o=singlecolumno(Kini,Rmini,psim);
etae=eta*0; Ze=Z*0; Ke=K*0; Rme=Rm*0; Kep=K*0; Rmep=Rm*0; oe=o*0; ntaue=Z*0;
% Zonal Fourier for run
kf=fftkspe(nx,dx); % wavenumber=2*pi*f (nk=nx everywhere)
K=fftspe(K); for im=1:nrm; Rm(:,im)=fftspe(Rm(:,im)); end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop of Run
%
for kt=1:nt-1
{'tend=',num2str(tlength+tini),'t=',num2str(tg(kt)),'Z=',num2str(Ze(1,1)),'eta=',num2str(etae(1,1))}
%
%%%%%Seasonal cycle
if doseaso==1; 
so=seasowarmpool(tg(kt),yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,1); %tg(kt) is current time in days
sq=so;%(or could recall using sqref if unbalanced)
end% so for a seasonal cycle

%%%%% Solve Equatorial Waves
% (K and R can be either in Fourier/Zonal, from dofourier)
% (Z,a,sq,so are in physical space=strips)
% Forcing
term=(HH*DDA*eta-so); % nonlinear
Sm=zeros(nx,nrm+2);% size nrm+2 is to solve all Rossby (but m>nym is null forcing)
Sm(:,1:nyk)=singlecolumnfm(term,psim,Hyk); % to Hermite coefficients
%
% Kelvin
term=-Sm(:,1)/sqrt(2);
term=fftspe(term);
Ke=longwavefana(K,ekr,1,term,kf,dt);
% Rossbys 
for im=1:nrm
term=-2*sqrt(im*(im+1))/(2*im+1)*( sqrt(im)*Sm(:,im+2)+sqrt(im+1)*Sm(:,im) );% GOOD !
term=fftspe(term);
Rme(:,im)=longwavefana(Rm(:,im),ekr,-1/(2*im+1),term,kf,dt);
end% loop im
% K-Rm (without zonal zonalFourier)
Kep=fftispe(Ke); for im=1:nrm; Rmep(:,im)=fftispe(Rme(:,im)); end;
%%%%% Solve Single-Column
% o, u (u is only for ncdf output)
oe=singlecolumno(Kep,Rmep,psim);
%%%%%%
% Z and eta (physical space)
% (beware, in spectral space dtam=GG*<q*a.psim> would be non trivial)
for i=1:nx; for j=1:nyk; % stochastic (with inner stochastic timestep)
[term1, term2, term3]=...% split method requires oe
singlecolumnstoch(eta(i,j),Z(i,j),oe(i,j),dt,rr,rr0,DDA,GG,HH,QQ,sq(i,j),so(i,j));
etae(i,j)=term1; Ze(i,j)=term2; ntaue(i,j)=term3;
end; end; % loop i,j
%
%%%%% Save external variables (physical space)
term=(kt+1)/mts;%modulo mts for saving (kt=tis)
if term==fix(term) & term<=nts;
ncdfmakevar(fileout,'Ks',{'X','T'},Kep,[1,term],0);
ncdfmakevar(fileout,'Rms',{'X','M','T'},Rmep,[1,1,term],0);
ncdfmakevar(fileout,'etas',{'X','Y','T'},etae,[1,1,term],0);
ncdfmakevar(fileout,'Zs',{'X','Y','T'},Ze,[1,1,term],0);
ncdfmakevar(fileout,'ntaus',{'X','Y','T'},ntaue,[1,1,term],0);
ncdfmakevar(fileout,'ts',{'one','T'},ts(1,term),[1,term],0);
end
%
%%%%% Increment variables
K=Ke; Rm=Rme; Z=Ze; eta=etae; o=oe;
end%loop kt
%
end; % dorun
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
