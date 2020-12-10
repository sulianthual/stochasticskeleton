function hf = hermitefunc(n,x)
% HERMITE: compute the Hermite functions (i.e. parabolic cylinder functions)
% by Sulian Thual
% 
%   h = hermitefunc(n)
%   h = hermitefunc(n,x)
% 
% Inputs:
%   - n is the order of the Hermite function i.e. polynomial (n>=0).
%   - x is values to be evaluated on the resulting Hermite
%     polynomial function.
% 
% Additional infos:
% uses hermitepoly.m (initial name is hermite.m) to compute hermite polynomials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute hermite polynomial
hpoly=hermitepoly(n,x);
% compute scaling factor
sfac=1./sqrt(2^n*factorial(n)*sqrt(pi));
% compute hermite function
hf=hpoly.*exp(-x.^2/2)*sfac;
