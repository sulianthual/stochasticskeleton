% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% : hovmuller of projection on one eigenmode (its amplitude), to evidence low-frequency modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
dorecalc=       1          ; % recompute
%
% Initial Filtering 
dofilterkw=     1          ; % filter within k-w the signal before projection (ref=1)
kmin=1; kmax=3;% in (2pi/40000km)
wmin=1/90; wmax=1/30; %
%
% Projection
imodshow=        2        ; % eigenmode to project on
inix=           10      ; % where to compute linear stability (in 1000 km)
donullmean=     0          ; % put k=0 to zero when projecting (ref=0) 
donullnyquist=  0          ; % put nyquist k to zero (ik=1) when projecting (ref=0)
dorefphase=     1          ; % impose reference phase (ref=1)
%
% smooth for graph
dosmooth=   1; % simple smooth for graphs (faster contour)
nxsmoo=     5; % x size of smooth (odd)
ntsmoo=     5; % t size of smooth (odd)
ismoo=6; % number of time smoothed
%
% graph
showampl=        1        ; % show abs (amplitude) of projection rather than real part (ref=1)
dolevsproj=       1       ; % do levels
levels1proj=(0:0.1:1)*0.2  ; 
timeunits=10000; % modify time units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ AND PREPARE OUPUTS FIRST

% Recompute

if dorecalc==1
'read'
% Read variable to concatenate
xconK=((1:nx)*0)'; 
xconR=((1:nx)*0)'; 
if nrm>1;xconR2=((1:nx)*0)'; end% for the Rossby, need to do each by hand
if nrm>2; xconR3=((1:nx)*0)'; end 
xconZ=((1:nx)*0)'; 
xcona=((1:nx)*0)'; 
tcon=[0,0];% first timestep is dummy
for ilpps=ilpmin: ilpmax
[ilpmax ilpps ]
%
indexrestart=ilpps; run(fileini); % get the infos
%
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas'); 
Zs=ncdfgetvar(fileout,'Zs');
% Zs and ts have to be converted to spectral-y coefficients 
for kts=1:nts; 
etas(:,:,kts)=singlecolumnfm(etas(:,:,kts),psim,Hyk);
Zs(:,:,kts)=singlecolumnfm(Zs(:,:,kts),psim,Hyk);
end
term=Ks; xconK=[xconK,term];  
term=squeeze(Rms(:,1,:)); xconR=[xconR,term]; 
if nrm>1;term=squeeze(Rms(:,2,:)); xconR2=[xconR2,term]; end
if nrm>2; term=squeeze(Rms(:,3,:)); xconR3=[xconR3,term]; end
term=squeeze(Zs(:,1,:)); xconZ=[xconZ,term]; 
term=squeeze(etas(:,1,:)); xcona=[xcona,term]; 
tcon=[tcon,ts]; % concatenate
end
%xcon=xcon(:,2:end); 
mrr=2;
tcon=tcon(mrr+1:end);
xconK=xconK(:,mrr:end);
xconR=xconR(:,mrr:end); 
if nrm>1; xconR2=xconR2(:,mrr:end); end
if nrm>2; xconR3=xconR3(:,mrr:end); end
xconZ=xconZ(:,mrr:end);
xcona=xcona(:,mrr:end);
%
nts=length(tcon); ts=tcon; tcon=0;
Ks=xconK; xconK=0;
Rms=zeros(nx,nrm,nts); Rms(:,1,:)=xconR; xconR=0;
if nrm>1; Rms(:,2,:)=xconR2; xconR2=0;  end
if nrm>2; Rms(:,3,:)=xconR3; xconR3=0; end
Zs=zeros(nx,1,nts); Zs(:,1,:)=xconZ; xconZ=0;
etas=zeros(nx,1,nts); etas(:,1,:)=xcona; xcona=0;
%
% Filter kw (only modify at jshow) 
'filter'
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
Ks=filterkw(Ks,dxg,dtg,[kmin,kmax],[wmin,wmax]);
for irm=1:nrm
Rms(:,irm,:)=filterkw(squeeze(Rms(:,irm,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
end
for iyk=1:nyk
Zs(:,iyk,:)=filterkw(squeeze(Zs(:,iyk,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
etas(:,iyk,:)=filterkw(squeeze(etas(:,iyk,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
end
end% end filter

% Compute FFT zonal
'fftzonal'
x=Ks;
for kts=1:nts; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
if donullmean==1; x(nx/2+1,:)=0; end
if donullnyquist==1; x(1,:)=0; end
Ks=x;
for j=1:nrm
x=squeeze(Rms(:,j,:));
for kts=1:nts; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
if donullmean==1; x(nx/2+1,:)=0; end
if donullnyquist==1; x(1,:)=0; end
Rms(:,j,:)=x;
end
for j=1:nyk
x=squeeze(Zs(:,j,:));
for kts=1:nts; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
if donullmean==1; x(nx/2+1,:)=0; end
if donullnyquist==1; x(1,:)=0; end
Zs(:,j,:)=x;
end
for j=1:nyk
x=squeeze(etas(:,j,:));
for kts=1:nts; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
if donullmean==1; x(nx/2+1,:)=0; end
if donullnyquist==1; x(1,:)=0; end
etas(:,j,:)=x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute projection
%
xxstab=searchclosest(xg,inix); % where to compute stability
'project1'
% Projection as in MS2011
% do zonal fft, then project, then come back (coeffs should be reals)
% Compute stability and project
kk=fftkspe(nx,dx); 
xproj=zeros(nx,nts);
% loop 
for ik=1:length(kk)
% compute stability
%[wr,wi,vects]=stabskel(kk(ik),ekr,QQ,GG,HH,soref,nyk,nrm);
[wr,wi,vects]=stabskelnew(kk(ik),ekr,QQ,GG,HH,squeeze(so(xxstab,:)),nyk,nrm,psim,Hyk);% new is with better parameters
wr=wr(imodshow); 
wrtest(ik)=wr;
vects=squeeze(vects(:,imodshow)); % X=K,Rm,Am,Zm
% CONVENTION: modify vect phase such that phase(K)=0
if dorefphase==1
term=vects(1); phaseK=phase(term);
vects=vects.*exp(-ii*phaseK);
end
% project
for kts=1:nts
% vector from simu (with zonal FFT)
Kp=Ks(ik,kts); Rmp=squeeze(Rms(ik,:,kts)); Ap=DDA*squeeze(etas(ik,:,kts)); Zp=squeeze(Zs(ik,:,kts)); 
nxs=1+nrm+nyk+nyk;
Xsim=[Kp,Rmp,Ap,Zp].'; 
xproj(ik,kts)=Xsim.'*vects; 
end% kts
end% ik
for kts=1:nts; xproj(:,kts)=ifft(ifftshift(xproj(:,kts))); end; % zonal
tss=ts/timeunits;
end% end of dorecalc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
term=abs(xproj); % take amplitude of projection rather !
if dosmooth==1; for ii=1:ismoo; term=calc_smoothn(term,[nxsmoo,ntsmoo]);  end; end
%
% Graph
figure(iwind); iwind=iwind+1;, clf;
subplot(1,2,1)
if dolevsproj; contourf(xg,tss,term',levels1proj,'LineStyle','none');
else contourf(xg,tss,term','LineStyle','none'); end
caxis([levels1proj(1) levels1proj(end)])
colorbar('location','southoutside');
xlabel('x (1000km)')
ylabel('time (days)')
title('Re')
%ylim(tlim);
%set(gca,'XTick',[]); set(gca,'YTick',[]);

