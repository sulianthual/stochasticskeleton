% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Graphi 2: hovmuller quick
% does hovmuller U,o,q,a non filtered, and MJO projection next to it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hovmullers of variables
showhermite=     0        ; % Show strips values/spectral coefficients on all graphs
yshow=           0        ; %1.4 for north% strip location (yg-1000km) to show (for showhermite=1)
mshow=           0        ; % Hermite to show (0<m<nyk-1) (for showhermite=0)
rshow=           1        ; % Rossby index (1<=m<=nrm) to show 
tlim=        [0,tlength]  ; % range on t for graphs (adding start time from netcdf)
%tlim=        [0, 1000]  ; 
%
dofilterkw=     1          ; % filter within k,w, for the hovmullers (ref=0) 
kmin=1; kmax=3; wmin=1/90; wmax=1/30; % in cpd, for hovmullers filtering
%
dosmooth=   1; % simple smooth for graphs (faster contour)
nxsmoo=     1; % x size of smooth (odd), ref=1
ntsmoo=     5; % t size of smooth (odd), ref=5
%
levelsu=     (-1:0.1:1)*4  ; % levels a for contours (m.s-1), dograph2
levelso=     (-1:0.1:1)*1  ; % levels a for contours (K), dograph2
levelsq=     (-1:0.1:1)*1  ; % levels a for contours (K), dograph2
levelsHa=    (0:0.05:1)*2  ; % levels a for contours (K.d-1), dograph2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data projection on mode
%
iimodshow= 2; % eignemode from linear stability on which to project (e.g., 2, or [2,3,4] for multiples)
inix=    10; % in 1000km, where to take so(y) profile for computing initial condition wave
%
dofilterkw2=     1          ; % filter within k-w the signal before projection (ref=1)
kmin2=1; kmax2=3;% in (2pi/40000km)
wmin2=1/90; wmax2=1/30; %
%
donullmean=     0          ; % put k=0 to zero when projecting (ref=0) 
donullnyquist=  0          ; % put nyquist k to zero (ik=1) when projecting (ref=0)
dorefphase=     1          ; % impose reference phase (ref=1)
%
dolevsproj=       1       ; % levels for data projection
levels1proj=(-1.2:0.1:1.2)*0.1  ; % real
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change time axis
%
modtime=0  ; % Modify time axis (ref=0) 
timeunits=1; % (scaling factor of time units, ref=1)
%
ngraphs=4+length(iimodshow)+1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs of run: hovmullers
%
% Read Outputs 
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas');
Zs=ncdfgetvar(fileout,'Zs');
ts=ncdfgetvar(fileout,'ts');
tlim=tlim+ts(1);
%
% Compute os and us
os=zeros(nx,nyk,nts); us=zeros(nx,nyk,nts);
for kts=1:nts
os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim);
us(:,:,kts)=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim);
end

% Rather show spectral coefficients
if showhermite==1
{'showhermite=1: graphs are of spectral coefficients'}
for kts=1:nts; 
os(:,:,kts)=singlecolumnfm(os(:,:,kts),psim,Hyk);
us(:,:,kts)=singlecolumnfm(us(:,:,kts),psim,Hyk);
etas(:,:,kts)=singlecolumnfm(etas(:,:,kts),psim,Hyk);
Zs(:,:,kts)=singlecolumnfm(Zs(:,:,kts),psim,Hyk);
end
end
%
% meridional location to show
if showhermite==1
jshow=mshow+1; % because 0<=m<=nyk while matlab array starts at =1
else
[jshow,term]=searchclosest(yg,yshow);% selected location for graphs-print
end
us=squeeze(us(:,jshow,:));
os=squeeze(os(:,jshow,:));
Zs=squeeze(Zs(:,jshow,:));
etas=squeeze(etas(:,jshow,:));
qs=Zs-QQ*os;
%
% Filter kw 
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
us=filterkw(us,dxg,dtg,[kmin,kmax],[wmin,wmax]);
os=filterkw(os,dxg,dtg,[kmin,kmax],[wmin,wmax]);
qs=filterkw(qs,dxg,dtg,[kmin,kmax],[wmin,wmax]);
etas=filterkw(etas,dxg,dtg,[kmin,kmax],[wmin,wmax]);
end% end filter
%
% Smooth 
if dosmooth==1
us=calc_smoothn(us,[nxsmoo,ntsmoo]); 
os=calc_smoothn(os,[nxsmoo,ntsmoo]); 
Zs=calc_smoothn(Zs,[nxsmoo,ntsmoo]); 
etas=calc_smoothn(etas,[nxsmoo,ntsmoo]); 
qs=calc_smoothn(qs,[nxsmoo,ntsmoo]); 
end
%
% Modify time axis
if modtime==1
ts=ts-tref;
tlim=tlim2;
end
%
% Hovmuller Contour Physical (on one hermite): u,o,q,Ha (with dims)
figure(iwind); iwind=iwind+1;, clf;
subplot(1,ngraphs,1)
term=ua*us;
if dolevs; contourf(xg,ts/timeunits,term',levelsu,'LineStyle','none');
else contourf(xg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsu(1) levelsu(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits); 
colorbar('location','southoutside');
xlabel('x(1000km)')
ylabel('time (days)')
title('u')
%
subplot(1,ngraphs,2)
term=oa*os;
if dolevs; contourf(xg,ts/timeunits,term',levelso,'LineStyle','none');
else contourf(xg,ts/timeunits,term','LineStyle','none'); end
caxis([levelso(1) levelso(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('t (days)')
title('theta')
%
subplot(1,ngraphs,3)
term=oa*qs; 
if dolevs; contourf(xg,ts/timeunits,term',levelsq,'LineStyle','none');
else contourf(xg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsq(1) levelsq(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('t (days)')
title('q')
%axis off;
%
subplot(1,ngraphs,4)
term=HH*DDA*oa/(ta/oneday)*etas; 
if dolevs; contourf(xg,ts/timeunits,term',levelsHa,'LineStyle','none');
else contourf(xg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsHa(1) levelsHa(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('t (days)')
title('Ha')
%axis off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data projection

% Read Outputs 
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas'); 
Zs=ncdfgetvar(fileout,'Zs');
ts=ncdfgetvar(fileout,'ts');
ii=complex(0,1);
% Zs and ts have to be converted to spectral-y coefficients 
if 1==1
for kts=1:nts; 
etas(:,:,kts)=singlecolumnfm(etas(:,:,kts),psim,Hyk);
Zs(:,:,kts)=singlecolumnfm(Zs(:,:,kts),psim,Hyk);
end
end

% Filter kw (only modify at jshow) 
if dofilterkw2==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin2); iwmax=searchclosest(wg,wmax2); 
ikmin=searchclosest(kg,kmin2); ikmax=searchclosest(kg,kmax2); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
Ks=filterkw(Ks,dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
for irm=1:nrm
Rms(:,irm,:)=filterkw(squeeze(Rms(:,irm,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
end
for iyk=1:nyk
Zs(:,iyk,:)=filterkw(squeeze(Zs(:,iyk,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
etas(:,iyk,:)=filterkw(squeeze(etas(:,iyk,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
end
end% end filter

% Compute FFT zonal
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
%
% Compute stability
kk=fftkspe(nx,dx); 
xproj=zeros(nx,nts);
% loop 
icount=1;
for imodshow=min(iimodshow):max(iimodshow)
for ik=1:length(kk)
% compute stability
xxstab=searchclosest(xg,inix); % where to compute stability
[wr,wi,vects]=stabskelnew(kk(ik),ekr,QQ,GG,HH,squeeze(so(xxstab,:)),nyk,nrm,psim,Hyk);% new is with better parameters
wr=wr(imodshow); 
wrtest(ik)=wr;
vects=squeeze(vects(:,imodshow)); % X=K,Rm,Am,Zm
%
% Modifiy phase (convention here is that phase(K)=0)
if dorefphase==1
term=vects(1); phaseK=phase(term);
vects=vects.*exp(-ii*phaseK);
end
%
% compute data projection
for kts=1:nts
% vector from simu (with zonal FFT)
Kp=Ks(ik,kts); Rmp=squeeze(Rms(ik,:,kts)); Ap=DDA*squeeze(etas(ik,:,kts)); Zp=squeeze(Zs(ik,:,kts)); 
nxs=1+nrm+nyk+nyk;
Xsim=[Kp,Rmp,Ap,Zp].'; 
xproj(ik,kts)=Xsim.'*vects; 
end% kts
end% ik
%
% Inverse Fourier
for kts=1:nts; xproj(:,kts)=ifft(ifftshift(xproj(:,kts))); end; 

% Modify time axis
if modtime==1
ts=ts-tref;
tlim=tlim2;
end

% Graph data projection 
if 1==1
%figure(iwind); iwind=iwind+1;, clf;
subplot(1,ngraphs,icount+4); icount=icount+1
term=real(xproj);
if dolevsproj; contourf(xg,ts/timeunits,term',levels1proj,'LineStyle','none');
else contourf(xg,ts/timeunits,term','LineStyle','none'); end
caxis([levels1proj(1) levels1proj(end)]);
colorbar('location','southoutside');
xlabel('x(1000km)')
%ylabel('t (days)')
pass=['m=',num2str(imodshow)];
title(pass)
ylim(tlim/timeunits);
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
end
%
end% loop on iimodshow

