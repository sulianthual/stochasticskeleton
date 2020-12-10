% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Compute all variables (u,v,o,q,a) in non-dimensional units, spectral coefficients,
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
% determine output file for here
filevars=strcat(dfolder,'/vars_',num2str(indexrestart), '.nc');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and write down each variable (those are spectral coefficients)
%
% Read (some) ncdf input
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
%
% u
ums=zeros(nx,nyk,nts);
for kts=1:nts
pass=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim);
ums(:,:,kts)=singlecolumnfm(pass,psim,Hyk);
end
ncdfmakevar(filevars,'ums',{'x','m','t'},ums,NaN,2); 
ums=0; % cleaning
%
% o
oms=zeros(nx,nyk,nts);
for kts=1:nts
pass=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim);
oms(:,:,kts)=singlecolumnfm(pass,psim,Hyk);
end
ncdfmakevar(filevars,'oms',{'x','m','t'},oms,NaN,1); 
Ks=0; % cleaning
%
% Read (some) ncdf input
Zs=ncdfgetvar(fileout,'Zs');
%
% q (after o)
qms=zeros(nx,nyk,nts);
for kts=1:nts
qms(:,:,kts)=singlecolumnfm(Zs(:,:,kts),psim,Hyk)-QQ*oms(:,:,kts);
end
ncdfmakevar(filevars,'qms',{'x','m','t'},qms,NaN,1); 
oms=0; qms=0; Zs=0; % cleaning
%
% Read (some) ncdf input
etas=ncdfgetvar(fileout,'etas');
%
% etas
etams=zeros(nx,nyk,nts);
for kts=1:nts
etams(:,:,kts)=singlecolumnfm(etas(:,:,kts),psim,Hyk);
end
ncdfmakevar(filevars,'etams',{'x','m','t'},etams,NaN,1); 
%
% Read (some) ncdf input
ts=ncdfgetvar(fileout,'ts');
%
% v (after eta, needs so=sq also)
vms=zeros(nx,nyk,nts);
for kts=1:nts;
if doseaso==1
spass=seasowarmpool(ts(kts),yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,1);
Sm=HH*DDA*etams(:,:,kts)-singlecolumnfm(spass,psim,Hyk);
else
Sm=HH*DDA*etams(:,:,kts)-singlecolumnfm(so,psim,Hyk);
end
pass=singlecolumnv(Rms(:,:,kts),Sm,psim,dx); 
vms(:,:,kts)=singlecolumnfm(pass,psim,Hyk);
end
ncdfmakevar(filevars,'vms',{'x','m','t'},vms,NaN,1); 
etams=0; vms=0; % cleaning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

