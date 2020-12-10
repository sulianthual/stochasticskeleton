function hg = hermitegauss(f,g,hgw)
% HERMITE: compute the inner product of f and g using Gauss-Hermite quadrature
% by Sulian Thual
% 
%   h = hermitegauss(f,g)
% 
% Inputs:
%  - f and g must be expressed at the roots hroots of the Hermite polynomial considered (vectors)
%  - hweigths are the Gauss Hermite quadrature coefficients (vector)
% 
% Outputs:
% will return the inner product <f.g> using gauss hermite quadrature
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%hg=sum((f.*g).*hgw);
%hg=sum((f'.*g').*hgw');
hg=f'*(g'.*hgw')';
