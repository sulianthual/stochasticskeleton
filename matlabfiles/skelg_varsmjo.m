% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Compute all variables (u,v,o,q,a) in non-dimensional units, spectral coefficients,
% HERE FILTERED IN MJO BAND
% and write down in files
%
% Input:
% - indexrestart: the restart file
% - fileini: the ini file with all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get all infos
run(fileini); % not indexrestart must also given for the fileini
%
% Get filtering box (MJO band)
kmin2=1; kmax2=3;% in (2pi/40000km)
wmin2=1/90; wmax2=1/30; %
%
% determine input/output file for here
filevars=strcat(dfolder,'/vars_',num2str(indexrestart), '.nc'); % input file
filevarsmjo=strcat(dfolder,'/varsmjo_',num2str(indexrestart), '.nc'); % input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up filtering
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin2); iwmax=searchclosest(wg,wmax2); 
ikmin=searchclosest(kg,kmin2); ikmax=searchclosest(kg,kmax2); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
%
% Compute and write down each variable (those are spectral coefficients)
%
% u 
passf=ncdfgetvar(filevars,'ums');
for iy=1:nyk;
pass=filterkw(squeeze(passf(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
passf(:,iy,:)=real(pass);
end; pass=0;
ncdfmakevar(filevarsmjo,'ums',{'x','m','t'},passf,NaN,2);
%
% v 
passf=ncdfgetvar(filevars,'vms');
for iy=1:nyk;
pass=filterkw(squeeze(passf(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
passf(:,iy,:)=real(pass);
end; pass=0;
ncdfmakevar(filevarsmjo,'vms',{'x','m','t'},passf,NaN,1);
%
% o
passf=ncdfgetvar(filevars,'oms');
for iy=1:nyk;
pass=filterkw(squeeze(passf(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
passf(:,iy,:)=real(pass);
end; pass=0;
ncdfmakevar(filevarsmjo,'oms',{'x','m','t'},passf,NaN,1);
%
% q 
passf=ncdfgetvar(filevars,'qms');
for iy=1:nyk;
pass=filterkw(squeeze(passf(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
passf(:,iy,:)=real(pass);
end; pass=0;
ncdfmakevar(filevarsmjo,'qms',{'x','m','t'},passf,NaN,1);
%
% eta 
passf=ncdfgetvar(filevars,'etams');
for iy=1:nyk;
pass=filterkw(squeeze(passf(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
passf(:,iy,:)=real(pass);
end; pass=0;
ncdfmakevar(filevarsmjo,'etams',{'x','m','t'},passf,NaN,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
