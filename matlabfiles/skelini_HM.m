% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% skelini: initialises the skeleton 
%
% Input (to define prior): 
% - indexrestart. indexrestart=1 does initial conditions
%                 else a restart, from the file dfolder/filename_indexrestart.nc
% - setup: name of setup
% 
% Output: 
% - all parameters for one simulation (to initialise, or do to post-processing)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Parameters
%
% Input parameter (elsewhere)
%indexrestart=    1       ; % index of file of restart, defined in skelmain_all
%
% Meridional Truncation
nyk=             1        ; % Meridional Truncation (strips/Hermites)
nrm=             1        ; % Rossby solved (can have nrm>=nyk-2,for which Sm=psim=0)
%
% Time Grid
tlength=         2000     ; % length of each restart (in days)
dt=              NaN      ; % timestep for solving (adim)(dt=NaN puts dt=dx/2)
mts=             10       ; % modulo of timestep for saving
%
% Parameters
DDA=             0.001    ; % ref=0.001
rr=              0        ; % ref=0
rr0=             1        ; % ref=1 (is the kronecker delta for boundary condition)
QQ=              0.9      ; % mean low-level moisture adim (=0.9);
GG=              1.66     ; % Gamma adim (=1.66)
HH=              0.22     ; % overbar H adim (=0.22)
ekr=             0        ; % friction on longwaves adim (=0)
%
% Background state of external forcings
backgroundtype=  0        ; % type of background to use, 0=homogeneous, 1=gaussian warmpool
sqref=           0.022    ; % reference sq external forcing adim (in physical space)
soref=           sqref    ; % reference so external forcing adim (in physical space)
wpint=           2        ; % intensity of warm pool (when used)
yoo=             0        ; % northward displacement of warm pool (when used)
doseaso=         0        ; % do a seasonal cycle (where the WP displaces from +-yo each year)
pseaso=          365       ; % seasonal cycle period (in days)
yseaso=          0;%0.6       ; % maximal yoo of seasonal cycle (otherwise yoo is used, e.g. for initiation/diagnostics)
%
% Files
dfolder=strcat('data',setup); % setup is defined in main
filein=strcat(dfolder,'/skeleton_',num2str(indexrestart-1), '.nc');
fileout=strcat(dfolder,'/skeleton_',num2str(indexrestart), '.nc');
% 
% Initial conditions
if indexrestart==1; initype=2; else initype= 1; end; % initial conditions(=1 file, =2 RCE+wave) 
inikwave=        2        ; % for initype=2, integer wavenumber(s) (inikwave=[1,2]... sums components)
iniampl=         0.05     ; % for initype=2,amplitude factor of eigenmodes (put 0.05 for amplitude 2 of WP !)
inimod=          2        ; % for initype=2, eigenmode (1-4, see corresponding row in dograph5)
inix=            10       ; % in 1000km, where to take so(y) profile for computing initial condition wave
%
% Specific Keys for solving method (do not touch)
dodefault=       1        ; % do default (faster) stochastic run with no keys
dorun =          1        ; % recompute skeleton (security)
dostocha=        1        ; % do stochastic/deterministic skeleton
dofourier=       1        ; % define and do waves with zonalfourier/finitediff (ref=1)
dofreewaves=     0        ; % remove forcing of long-waves (& add ideal Cini) (ref=0)
dofreecolumn=    0        ; % remove forcing of single-column (QQ=0) (ref=0)
looprun=         0        ; % restart from same file to loop runs (also sets initype=1) (ref=0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%
% Keys dependencies
if dostocha~=1; DDA=1; end; % for deterministic A=eta
if looprun==1; initype=1; end
if looprun==1; filein=fileout; end
if dofreecolumn; QQ=0; end; % single-columns
if isfinite(rr0)==0; rr0=rr; end 
%
%%%%%%%%
% Constants and Dimensionalisation
onehour=3600.;% one hour in seconds
oneday=onehour*24.;% one day in seconds
onekm=1000.;% one km in meters
iwind=2; % first window number
ii=complex(0,1);
xa=1500.*onekm; %x adim (meters),
ya=1500.*onekm; %y adim (meters)
ta=8*onehour; %t adim (seconds),
oa=15; % temperature adim (deg Kelvin)
ua=50; % velocity adim (m.s-1)
%
%%%%%%%%
% Grids
nx=64;
dx=625.*onekm/xa;% dx grid step (adim)
xx=(0:nx-1)*dx; % adim x axis
LL=nx*dx;
xg=(0:nx-1)*dx*xa/1000/1000; % axis x for graphs, in 1000km
yk=hermiteroots(nyk); % adim y axis (nyk in total)
yg=yk*ya/1000/1000; % axis y for graphs, in 1000km 
Hyk=hermitegaussw(yk);% Hermite Gauss coefficients
nym=nyk;% number of Hermites 
psim=zeros(nyk,nym); % Hermite functions values at strips yk
for im=1:nyk; psim(:,im)=hermitefunc(im-1,yk); end;
%
if isfinite(dt)==0; dt=dx/2; end 
nt=fix(tlength*oneday/ta/dt); % for solving
tt=(0:nt-1)*dt; % adim t axis 
tg=tt*ta/oneday; % axis t for graphs, in days (all computed steps)
nts=fix(tlength*oneday/ta/(mts*dt)); % for saving (is a guess)
tis=1:mts:1+nts*mts; tis=tis(tis<=nt); nts=length(tis)-1;
tts=(0:nts-1)*dt*mts; % adim t axis (of saves)
ts=tts*ta/oneday;% axis for saving and graphs, in days (restricted to saving steps)
%
%%%%%%%%
% RCE From External forcings (in physical space=strips)
% homogeneous
if backgroundtype==0
sq=zeros(nx,nyk)+sqref; so=zeros(nx,nyk)+soref; 
end
% warm pool
if backgroundtype==1
sq=zeros(nx,nyk)+sqref; so=zeros(nx,nyk)+soref; 
psimnorth=zeros(nyk,1); psimnorth(:,1)=hermitefunc(0,yk-yoo);
for ik=1:nyk; sq(:,ik)=sqref*(1-wpint*0.3*cos(2*pi/LL*xx))/0.7511*psimnorth(ik,1); end; so=sq; % zonal warm pool
end
% RCE state (in physical space=strips)
MMeta=so/(HH*DDA); % MMeta(nx,nyk). also RCE is with K,Rm,Zm=0
%
%%%%%%%%
% Initial conditions 
Kini=zeros(nx,1); Rmini=zeros(nx,nrm);
etaini=zeros(nx,nyk); Zini=zeros(nx,nyk);
%
% Initial conditions from (restart) file
if initype==1
term=ncdfgetvar(filein,'Ks'); Kini(:,1)=term(:,end); 
term=ncdfgetvar(filein,'Rms'); Rmini(:,:)=term(:,:,end); 
term=ncdfgetvar(filein,'etas'); etaini(:,:)=term(:,:,end); 
term=ncdfgetvar(filein,'Zs'); Zini(:,:)=term(:,:,end); 
term=ncdfgetvar(filein,'ts'); tini=term(1,end);
tg=tg+tini; ts=ts+tini; % HERE MANAGES THE UPDATE OF TIME IN RESTARTS !
end%initype=1
%
% Initial conditions from stability +RCE state
if initype==2
passx=searchclosest(xg,inix);
for ikk=1:length(inikwave)
kk=(2*pi/LL)*inikwave(ikk); 
[wr,wi,vects]=stabskelnew(kk,ekr,QQ,GG,HH,squeeze(so(passx,:)),nyk,nrm,psim,Hyk);
Kini(:,1)=Kini(:,1)+iniampl*real(vects(1,inimod)*exp(-ii*kk*xx))';% Beware of sign -kk !
for im=1:nrm
Rmini(:,im)=Rmini(:,im)+iniampl*real(vects(1+im,inimod)*exp(-ii*kk*xx))';
end
for ik=1:nyk; % eta and Z are here in spectral Hermite space
etaini(:,ik)=etaini(:,ik)+(iniampl/DDA)*real(vects(1+nrm+ik,inimod)*exp(-ii*kk*xx))';
Zini(:,ik)=Zini(:,ik)+iniampl*real(vects(1+nrm+nyk+ik,inimod)*exp(-ii*kk*xx))';
end;
end%ikk
etaini=singlecolumnf(etaini,psim); 
Zini=singlecolumnf(Zini,psim); 
etaini=etaini+MMeta;% add RCE state to the perturbation
tini=0;
end%initype=2
%
% Check on initial conditions
if dostocha==1;etaini=round(etaini); end; % eta stochastic must be integer
if min(etaini(:))<0;% eta must be positive
{'BEWARE: eta is negative=', min(etaini(:))}
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
