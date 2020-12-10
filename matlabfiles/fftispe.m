function f=fftispe(F)
% Computes Inverse Fast Fourier Transform of f 
% by Sulian Thual
%
% Input:
% - field F(k) (in Fourier space)
% - dx : stepping of field
% Output:
% - field f 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=ifft(ifftshift(F)); 
f=real(f); % remove potential round-off imaginary part


