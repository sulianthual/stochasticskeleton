function [wr,wi,vects]=stabskelnew(kk,ekr,QQ,GG,HH,soref,nyk,nrm,psim,Hyk)
%
% Skeleton model (deterministic or stochastic)
% by Sulian Thual
% 
% computes stability on K,Rm,Zm,Am (defined on Hermite functions)
% Non dimensional units from Majda Stechmann2011

% Input:
% - kk=desired wavenumber, kk=j*(2pi/L) with j integer (adim)
% - parameters...
% - nyk=numbers of Hermites (for Am,Zm), from m=0
% - nrm=number of Rossbys, from m=1
% - Here it is assumed so(x,y)=so(y)  i.e. homogeneous zonally
% - 
% Output:
% [wr,wi,Vects]:
% - eigenvalues(nmods)
% - vects(nxs,nmods)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RCE state
MMa=zeros(nyk,1); 
for ik=1:nyk; MMa(ik)=soref(ik)/HH; end
%
% Compute Matrix for : dtX + AA X + CC dxX=0
% X=(K, Rm, Am, Zm ) with all hermite-strips concatenation   
% Rm starts from m=1
% Am, Zm starts from m=0
nxs=1 + nrm + nyk + nyk; % size of X
ii=complex(0,1);
MAm=singlecolumnfm(MMa',psim,Hyk); % to Hermite coefficients
%
% Compute CC:
CC=zeros(nxs,nxs);
CC(1,1)=1;% K
for im=1:nrm
CC(1+im,1+im)=-1/(2*im+1);% Rm
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute AA:
AA=zeros(nxs,nxs);
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Kelvin: (dt+e)K+dxK+So/sqrt(2)=0, So=<HH*a.psim>
AA(1,1)=ekr;%K:K
AA(1,1+nrm+1)=HH/sqrt(2);%K:A0 (Ok)
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Rossby-m: Rm, for m=1:nrm
% dtRm-dxRm/(2m+1)+ 2*sqrt(m(m+1))/(2m+1)(  delta(m<=M-2)*sqrt(m)*HH*Am+1 + sqrt(m+1)*HH*Am-1) =0, for 1<=m<=nrm
for im=1:nrm
AA(1+im,1+im)=ekr;%Rm:Rm
end
for im=1:nrm
AA(1+im,1+nrm+1+(im-1))=2*sqrt(im*(im+1))/(2*im+1)*sqrt(im+1)*HH;%Rm:Am-1 (Ok)
if im<=nyk-2;
AA(1+im,1+nrm+1+(im+1))=2*sqrt(im*(im+1))/(2*im+1)*sqrt(im)*HH; %Rm:Am+1 (Ok)
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Zm: dtZm -(QQ-1)*HH*Am=0
for ik=0:nyk-1
AA(1+nrm+nyk+(1+ik),1+nrm+(1+ik))=-(QQ-1)*HH;%Zm:Am (Ok)
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Am: Am, for m=0:nyk-1
% dtAm - GG* [sum(i=0,nyk-1) Qi *[ sum(j=0,nyk-1) MAj * gamma_ijm]  ] =0
% where:
% gamma_ijm= < phi_i*phi_j*phi_m> (scalar product <>)
% Qi= Zi - QQ * Oi
% Oi=-delta(i=0)*K/sqrt(2) - delta(i>=2)*R_i-1/4/sqrt(i) - delta(i<=nrm-1)*R_i+1/4/sqrt(i+1)
%
% compute gammaij
gammaijm=zeros(nyk,nyk,nyk);
for ik=0:nyk-1
% compute the gamma_ik_jm
for jk=0:nyk-1
for mk=0:nyk-1
f1=psim(:,ik+1); f2=psim(:,jk+1); f3=psim(:,mk+1); 
hg=((f1'.*f2').*f3')*Hyk';
gammaijm(ik+1,jk+1,mk+1)=hg;
end
end
end
%Am:Zi
% % dtAm -  GG* sum (i=0,nyk-1) * Zi * [ sum (j=0,nyk-1) MAj * gamma_ijm ] + .... =0
for mk=0:nyk-1% Am
for ik=0:nyk-1% Zi
% here has to sum over all jk
for jk=0:nyk-1% j 
AA(1+nrm+(1+mk),1+nrm+nyk+(1+ik))=AA(1+nrm+(1+mk),1+nrm+nyk+(1+ik)) -GG* MAm(jk+1)*gammaijm(ik+1,jk+1,mk+1); % (Ok)
end
end
end
% 
% Am: K and Rm (in the Oi term)
% % dtAm + GG *QQ * sum (i=0,nyk-1)  *  Oi *  [ sum (j=0,nyk-1) MAj * gamma_ijm] + .... =0
% Oi=-delta_i0 K/sqrt(2) - delta(i>=2)*R_i-1 / 4/sqrt(i) - delta(i>=2)*delta(i<=nrm-1)*R_i+1 / 4/sqrt(i+1) 
% 
% Am: K 
for mk=0:nyk-1
ik=0;
for jk=0:nyk-1% j 
AA(1+nrm+(1+mk),1)=AA(1+nrm+(1+mk),1) + GG*QQ* (-1/sqrt(2))* MAm(jk+1)*gammaijm(ik+1,jk+1,mk+1); % (Ok)
end
end
%
% Now Am: Rm 
% dtAm + GG* QQ* [sum(i=0,nyk-1) Oi *[ sum(j=0,nyk-1) MAj * gamma_ijm]  ] +....=0
% Oi=-delta(i=0)*K/sqrt(2) - delta(i>=2)*R_i-1/4/sqrt(i) - delta(i<=nrm-1)*R_i+1/4/sqrt(i+1)
for mk=0:nyk-1
for ik=0:nyk-1; 
%
if ik>=2;
for jk=0:nyk-1% j 
AA(1+nrm+(1+mk),1+(ik-1))=AA(1+nrm+(1+mk),1+(ik-1)) + GG*QQ* ( -1/4/sqrt(ik) )* MAm(jk+1)*gammaijm(ik+1,jk+1,mk+1); % the Ri-1
end
end
%
if ik<=nrm-1; 
for jk=0:nyk-1% j 
AA(1+nrm+(1+mk),1+(ik+1))=AA(1+nrm+(1+mk),1+(ik+1)) + GG*QQ* ( -1/4/sqrt(ik+1) )* MAm(jk+1)*gammaijm(ik+1,jk+1,mk+1); % the Ri+1
end
end% ik<=nrm-1
%
end
end
%

%
%
%%%%%%%%%%%%%%%
% Compute stability
% System: dtX + AA X + CC dxX=0
% Solution : X=Xo*exp(i(kx-wt))
% -iw*X + AA*X + ik*CC*X=0
% dtX = (-AA -i*k*CC)X=MX, eigenvalues l(k)
% w(k)=i*l(k)
MM=AA-ii*kk*CC;
[V,d]=eig(MM);
term=size(d);
nmod=term(1); % can have less modes than nxs
vects=zeros(nxs,nmod);
eigenvals=zeros(nmod,1);
for ixs=1:nmod; 
eigenvals(ixs,1)=d(ixs,ixs); 
end
vects=V;
wr=real(ii*eigenvals);% w(k)=i*l(k)
wi=imag(ii*eigenvals);% w(k)=i*l(k)

% ORDER BY DECREASING SPEED:
%1=Kelvin dry (50 ms-1)
%2=MJO mode (5 ms-1)
%3=Rossby moist mode (-5 ms-1)
%4=Rossby dry mode (-15 ms-1)
if 1==1
cr=wr/kk;
[term,index]=sort(cr,'descend');% sort by decreasing speed
wr=wr(index);
wi=wi(index);
vects(:,:)=vects(:,index);
end

% Notes to recall
%wi=imag(ii*eigenvals);
%cr=wr/kk;
% Solution : X=Xo*exp(i(kx-wt))
%cos(2pi*t/T) has period T
% frequency=2pi/T=wr 
%Cycle per day=1/T=wr/2pi
% wavenumber=2pi/L=k=j*(2pi/40000km)
%wr=wr/ta*oneday;
%wi=wi/ta*oneday;

% Parameters (adim may be wrong !)
%onehour=3600.;% one hour in seconds
%oneday=onehour*24.;% one day in seconds
%onekm=1000.;% one km in meters
%xa=1500.*onekm; %x adim (meters), cf BielloMajda2004 
%dx=625.*onekm/xa;% dx grid step (adim)
%ta=8*onehour; %t adim (seconds)














