function [fout]=ncdfmakevar(filenc,varnc,dimnc,inputnc,startnc,erasenc)
% by Sulian Thual

% entries
%inputnc=ones(1,2,3);
%filenc='test.nc';
%varnc='a';
%dimnc={'X','Y','Z'};
% startnc=array of starting locations (NaN for regular)
% erasenc= (change is for create or overwrite)
%        0 : not change file, not change var, change data
%        1 : not change file, change var, change data
%        2 : change file, change var, change data
%
% examples :
% ncdfmakevar('test.nc','var',{'X','Y','Z'},ones(1,2,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Change file
if erasenc==2;
delete(filenc);
end
%
% Change Var
if erasenc>=1;
sdimnc=size(inputnc);
%dimnc=fliplr(dimnc); % indeed is writen in wrong sense for matlab
ndimnc=length(sdimnc);
cellnc={};
for inc=1:ndimnc 
cellnc=[cellnc,dimnc(inc),sdimnc(inc)];
end
nccreate(filenc,varnc,'Dimensions', cellnc, 'Format', 'classic');
end
% (format=classic allows to use ncdump)
%
% Write
if min(isfinite(startnc))==1;% write with start location
%stridenc=startnc*0+1;
ncwrite(filenc,varnc,inputnc,startnc);
else
ncwrite(filenc,varnc,inputnc); % write classic
end

% Display
displaync=0;
if displaync, ncdisp(filenc);, end

% Return
fout=1;

