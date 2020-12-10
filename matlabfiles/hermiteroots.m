function hr = hermiteroots(n)
% HERMITE: compute the roots of an Hermite function (i.e. parabolic cylinder function)
%          this is equivalent to the roots of its Hermite polynomial
% by Sulian Thual
% 
%   h = hermiteroots(n)
% 
% Inputs:
%   - n is the order of the Hermite function i.e. polynomial (n>=0).
% Outputs:
% returns array with all roots
% 
% Additional infos:
% uses hermitepoly.m (initial name is hermite.m) to compute hermite polynomials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute hermite polynomial coefficients
hpolycoeffs=hermitepoly(n);
% compute roots of polynomial
hr=roots(hpolycoeffs)';
[hr,hrindex]=sort(hr);
