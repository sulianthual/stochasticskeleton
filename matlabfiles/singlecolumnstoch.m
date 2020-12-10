function [etae,Ze,ntaue]=singlecolumnstoch(eta,Z,o,dt,rr,rr0,DDA,GG,HH,QQ,sq,so)
% Single column skeleton model (stochastic)
% by Sulian Thual
%
% Solves from t to td+dt the coupled system:
%dteta=GG*(Z-QQ*o)*eta but with Master equation (with random walk r)
%dtZ=(QQ-1)H*DDA*eta+(sq-QQ*so)
%
% ntaue is the number of jumps that have occured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gillespie algorithm
t=0; etae=eta; Ze=Z;
ntaue=0; % count number of jumps
while t<dt;
%
% Compute transition rates
qe=Ze-QQ*o;
if etae>0; 
krod=0; rrput=rr; 
else 
krod=1; rrput=rr0; 
end
if qe>=0; 
rup=GG*qe*etae+rrput; rdown=rrput*(1-krod);
else 
rup=rrput; rdown=-GG*qe*etae+rrput*(1-krod); 
end
% Compute time until next jump
tau=1/(rup+rdown)*log(1/rand);
% Jump: 
% compute if is not last/last(=omitted)
if t+tau<=dt;% not last jump
% Z update (not last jump)
Ze=Ze+( (QQ-1)*HH*DDA*etae+(sq-QQ*so) )*tau;
% eta update (not last jump)
if rand<=rup/(rup+rdown);
etae=etae+1;
else; 
etae=etae-1; 
end
if etae<0; etae=0; end% force if error eta<0
ntaue=ntaue+1; 
else% last jump
% Z update (last jump=omitted)
Ze=Ze+( (QQ-1)*HH*DDA*etae+(sq-QQ*so) )*(dt-t);
end
% Time update
t=t+tau;
end%while
