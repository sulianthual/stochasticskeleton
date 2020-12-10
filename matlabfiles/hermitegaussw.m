function Hyk = hermitegaussw(hroots)
% HERMITE: computes the weigth coefficients of Gauss-Hermite quadrature
% by Sulian Thual
% 
%   h = hermitegaussweigths(hroots)
% 
% Inputs:
%  - hroots the roots of the Hermite polynomial considered 
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% deduce order N of Hermite polynomial
n=length(hroots);
% compute Gauss-Hermite weighting coefficients
Hyk=1./(n*hermitefunc(n-1,hroots).^2);
