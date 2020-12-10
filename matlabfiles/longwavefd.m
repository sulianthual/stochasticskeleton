function Ke = longwavefd(K,Ek,Ck,Pk,dx,dt)
% by Sulian Thual
% Solves evolution of a long-wave using finite difference
% at 1st order: upward sheme par derivatives and Euler in time
%
% solves: 
% dtK+Ek*K+Ck*dxK=Pk with periodic boundary condition
% 
% Inputs :
% - K(x,t), on grid x and  at time t 
% - speed Ck (note if Ck<0 the upward scheme is inverted)
% - forcing Pk(x,t) on grid x and at time t
% - spatial step dx of grid x
% - temporal step dt
%
% Outputs:
% K(x,t+dt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% compute spatial derivative with 1st order finite difference
nx=length(K);
if Ck>=0
for i=2:nx; dxK(i)=(K(i)-K(i-1))/dx; end
dxK(1)=(K(1)-K(nx))/dx;% periodic boundary
else
for i=1:nx-1; dxK(i)=(K(i+1)-K(i))/dx; end
dxK(nx)=(K(1)-K(nx))/dx;% periodic boundary
end
%
% Compute time evolution with 1st order finite difference (Euler)
Ke=K*0;
for i=1:nx
Ke(i)=K(i) + dt*(-Ek*K(i) -Ck*dxK(i) + Pk(i) );
end

