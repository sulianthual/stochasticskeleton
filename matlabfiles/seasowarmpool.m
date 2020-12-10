function so = seasowarmpool(time,yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,type)
% compute warm pool position
% by Sulian Thual

% if type=1, return so
% if type=2, return ynow the warm pool center (and can be many times all at once)
%
ynow=yseaso*sin(2*pi*time/pseaso);
%
if type==1
% compute new warm pool
so=zeros(nx,nyk)+soref;
psimnorth=zeros(nyk,1); psimnorth(:,1)=hermitefunc(0,yk-ynow);
for ik=1:nyk; 
so(:,ik)=soref*(1-wpint*0.3*cos(2*pi/LL*xx))/0.7511*psimnorth(ik,1); 
end;
end
%
if type==2
so=ynow;
end
% so=sq; % zonal warm pool
