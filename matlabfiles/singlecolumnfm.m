function fm=singlecolumnfm(f,psim,Hyk)
% Single column skeleton model (deterministic)
% by Sulian Thual
%
% computes projection of f(x,yk) on psim functions
% 
% Inputs:
% - f(nx,nyk)
% - psim(nyk,nym)
% - Hyk(nyk)
% 
% Outputs
% - fm(x,nym) contribution on psim
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
term=size(f);
nx=term(1);
term=size(psim);
nym=term(2);
fm=zeros(nx,nym);
for i=1:nx;
for im=1:nym;
fm(i,im)=hermitegauss(f(i,:)',psim(:,im),Hyk');
end
end


