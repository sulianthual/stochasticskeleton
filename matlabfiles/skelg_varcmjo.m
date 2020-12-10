% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Compute circulation variables (dxu,dyu, dxv, dyv) in non-dimensional, spectral coefficients,
% and write down
% This is done from variables filtered in the MJO band
% from then can compute div=dxu+dyv and curl=dxy-dyu
%
% Input:
% - indexrestart: the restart file
% - fileini: the ini file with all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get all infos
run(fileini); % not indexrestart must also given for the fileini
%
% determine input/output file for here
filevarsmjo=strcat(dfolder,'/varsmjo_',num2str(indexrestart), '.nc'); % input file
filevarcmjo=strcat(dfolder,'/varcmjo_',num2str(indexrestart), '.nc'); % input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute and write down each variable (those are spectral coefficients)
%
% Get u (filtered mjo band)
ums=ncdfgetvar(filevarsmjo,'ums');
% Compute dxu, dyu
dxums=zeros(nx,nyk,nts); dyums=zeros(nx,nyk,nts);
for kts=1:nts
pass=squeeze(ums(:,:,kts));
dxums(:,:,kts)=singlecolumndx(pass,xx); 
dyums(:,:,kts)=singlecolumndy(pass); 
% note for dy the spectral basis remains M-truncated, such that appearing higher-terms are omitted !  
end
ncdfmakevar(filevarcmjo,'dxums',{'x','m','t'},dxums,NaN,2);
ncdfmakevar(filevarcmjo,'dyums',{'x','m','t'},dyums,NaN,1);
%
% Get v (filtered mjo band)
vms=ncdfgetvar(filevarsmjo,'vms');
% Compute dxu, dyu
dxvms=zeros(nx,nyk,nts); dyvms=zeros(nx,nyk,nts);
for kts=1:nts
pass=squeeze(vms(:,:,kts));
dxvms(:,:,kts)=singlecolumndx(pass,xx); 
dyvms(:,:,kts)=singlecolumndy(pass); 
% note for dy the spectral basis remains M-truncated, such that appearing higher-terms are omitted !  
end
ncdfmakevar(filevarcmjo,'dxvms',{'x','m','t'},dxvms,NaN,1);
ncdfmakevar(filevarcmjo,'dyvms',{'x','m','t'},dyvms,NaN,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
