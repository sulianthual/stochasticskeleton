function b=singlecolumndy(a)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% compute meridional derivative using Hermite series
%  input=u(x,y) in Hermite spectral coefficients
% output=dyu(x,y) in Hermite spectral coefficients

% use BM2006 convention for Hermites:
% sqrt2 dypsim = -sqrt(m+1)*psim+1 +sqrt(m)*psim-1
% so dxum= - um-1 * sqrtm/sqrt2 + um+1*sqrt(m+1)/sqrt2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
pass=size(a); nx=pass(1); nyk=pass(2);
b=zeros(nx,nyk);
for m=0:nyk-1
mi=m+1;
if mi>1; b(:,mi)=b(:,mi) - a(:,mi-1)*sqrt(m)/sqrt(2); end
if mi<nyk; b(:,mi)=b(:,mi) + a(:,mi+1)*sqrt(m+1)/sqrt(2); end
end

