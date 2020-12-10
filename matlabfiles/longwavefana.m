function Ke = longwavefana(K,Ek,Ck,Pk,k,dt)
% Solves evolution of a long-wave using analytic solution in zonal Fourier space
% by Sulian Thual
%
% solves: 
% dtK+Ek*K+Ck*dxK=Pk with periodic boundary condition
% which gives in zonal Fourier space (dx=i*k): 
% dtK+ik*Ck*K=Pk
% K(k,t+dt)=K(k,t)exp(-i*k*Ck*dt)+Pk/(i*k) FOR k ne 0
% K(0,t+dt)=K(0,t) + Pk*dt
% 
% Inputs :
% - K(k,t), in zonal Fourier space  at time t 
% - speed Ck 
% - forcing Pk(k,t) in zonal Fourier space at time t
% - k: zonal wavenumber (give row of values) (in Fourier it can also be 2*pi*k)
% - dt: timestep
%
% Outputs:
% K(k,t+dt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
ii=complex(0,1); % the complex number
%k=2*pi*k; % not necessary if uses fftkspe.m
%
Ke=K.*exp(-Ek*dt-ii*k'*Ck*dt) +(Pk./(Ek+ii*k'*Ck)).*(1-exp(-Ek*dt-ii*k'*Ck*dt));

if Ek==0;% Ek is the damping
pass=k==0; pass=max(pass);
if pass==1
[i0,term]=searchclosest(k,0); %k(i0)
Ke(i0)=K(i0) + Pk(i0)*dt;
end
end



