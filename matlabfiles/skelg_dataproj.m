% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% script to compute data projections (write into files): 
% on files from ilpmin to ilpmax (in skelmain_all.m)
% on data projection from imodmin to imodmax
% NOTE: the output file where proj is put depends also on the ini file
%filecalc=strcat(dfolder,'/lagcodatabs1_',num2str(imodshow),'_',num2str(imodshow1),'.nc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% min an max range of modes to compute
imodmin=2;
imodmax=3; 
%
% filtering and other parameters of the projection
dofilterkw2=     1          ; % filter within k-w the signal before projection (ref=1)
kmin2=1; kmax2=3; wmin2=1/90; wmax2=1/30;% in (2pi/40000km) and in cpd
%
donullmean=     0          ; % put k=0 to zero when projecting (ref=0) 
donullnyquist=  0          ; % put nyquist k to zero (ik=1) when projecting (ref=0)
dorefphase=     1          ; % impose reference phase (ref=1)
dopropaphase=    0          ; % further add propagation phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection
%
for imodshow=imodmin:imodmax;% loop on modes to project on
for ilpps=ilpmin:ilpmax% loop on restart files
%
{imodshow,imodmax,ilpps,ilpmax}% print current
indexrestart=ilpps; run(fileini);% get all infos, notably the restart file
fileproj=strcat(dfolder,'/dataproj_',num2str(indexrestart),'_',num2str(imodshow), '.nc');% file to write 
%
% Read Outputs 
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas'); 
Zs=ncdfgetvar(fileout,'Zs');
ts=ncdfgetvar(fileout,'ts');
ii=complex(0,1);
% Zs and ts have to be converted to spectral-y coefficients 
for kts=1:nts; 
etas(:,:,kts)=singlecolumnfm(etas(:,:,kts),psim,Hyk);
Zs(:,:,kts)=singlecolumnfm(Zs(:,:,kts),psim,Hyk);
end

% Filter kw (only modify at jshow): filtering is done independently for each restart 
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
% Compute stability and project
kk=fftkspe(nx,dx); 
xproj=zeros(nx,nts);
for ik=1:length(kk)
% compute stability
xxstab=searchclosest(xg,inix); % where to compute stability
[wr,wi,vects]=stabskelnew(kk(ik),ekr,QQ,GG,HH,squeeze(so(xxstab,:)),nyk,nrm,psim,Hyk);% new is with better parameters
vects=squeeze(vects(:,imodshow)); % X=K,Rm,Am,Zm
%
% modify vect phase such that phase(K)=0
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
%
% write down
ncdfmakevar(fileproj,'xproj',{'X','T'},xproj,NaN,2);
%
end% loop on modes to project on
end% loop on restart files
