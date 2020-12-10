function F=fftspe(f)
% Computes Fast Fourier Transform of f (arranged version)
% by Sulian Thual
%
% Input:
% - field f(x) (supposed periodic, and size nx MUST BE EVEN (pair)
% Output:
% - field F with Fast Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=fftshift(fft(f));
