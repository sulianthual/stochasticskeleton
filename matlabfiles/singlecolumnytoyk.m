function g=ytoyk(f,yk,psim,y)

% Single column skeleton model (deterministic)
% by Sulian Thual
%
% given f(x,yk), gives same field g(x,y) with y as wanted
% with psim of size nm, there are nm-2 Rossby waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nyk=length(yk);
term=size(f); nx=term(1);
ny=length(y);
term=size(psim);
nym=term(2);
g=zeros(nx,ny);
%
Hyk=hermitegaussw(yk);
%
for im=1:nym
for i=1:nx
fm=hermitegauss(f(i,:)',psim(:,im),Hyk');
g(i,:)=g(i,:)+fm*hermitefunc(im-1,y);
end
end
