% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% plot stability (uses a sq(y) profile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%showgraphsstab=   0 ; % do graphs (if called by another program, must be set to zero)
%
% Parameters
nmodshow=        NaN     ; % index e.g. [1,3] of (fastest) eigenmodes to show in dograph5 (ref=NaN shows all)
npergraph=       4         ; % number of modes to show on each graph   
%
kk=2*pi/LL*(0.1:0.1:10); % values to show (line)
kkp= 2*pi/LL*(0:1:10); % values to overplot (points)
kk=[-fliplr(kk),kk]; kkp=[-fliplr(kkp),kkp];% take also negative values

showgrowth=      0        ; % show growth rate on dograph5

%inix= 20; % where to compute stability (on xg, in 1000km). this is normally set in the ini file
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
xxstab=searchclosest(xg,inix); % where to compute stability
%
% Index of points to show
ikp=kkp*0;
for ik=1:length(kkp)
ikp(ik)=searchclosest(kk,kkp(ik));
end

% Compute
nxs=1+nrm+nyk+nyk;
nkk=length(kk);
if isfinite(nmodshow)==1;
imod1=nmodshow(1); imod2=nmodshow(2);
nmod=imod2-imod1+1;
else 
%[term1,term2,term3]=stabskel(kk(1),ekr,QQ,GG,HH,soref,nyk,nrm);
[term1,term2,term3]=stabskelnew(kk(1),ekr,QQ,GG,HH,squeeze(so(xxstab,:)),nyk,nrm,psim,Hyk);% new is with better parameters
nmod=length(term1); 
imod1=1; imod2=nmod;
end
wr=zeros(nmod,nkk)+NaN; wi=zeros(nmod,nkk)+NaN; vects=zeros(nxs,nmod,nkk)+NaN;

for ik=1:nkk
kki=kk(ik);
%[term1,term2,term3]=stabskel(kki,ekr,QQ,GG,HH,soref,nyk,nrm);
[term1,term2,term3]=stabskelnew(kki,ekr,QQ,GG,HH,squeeze(so(xxstab,:)),nyk,nrm,psim,Hyk);% new is with better parameters
wr(1:nmod,ik)=term1(imod1:imod2);wi(1:nmod,ik)=term2(imod1:imod2); vects(:,1:nmod,ik)=term3(:,imod1:imod2);
end
wcpd=wr/ta*oneday/2/pi;
wgcpd=wi/ta*oneday/2/pi;
kj=kk/xa/2/pi*40000*1000;
cr=wr*0; for ik=1:nkk; cr(:,ik)=wr(:,ik)/kk(ik)*xa/ta;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if showgraphsstab==1
%
% Show
vectshow=abs(vects); % has to show module ! 
figure(iwind); iwind=iwind+1;, clf;
% column1=speed
for imod=1:nmod
if showgrowth==1; subplot(nmod,4,1+(imod-1)*4); else subplot(nmod,3,1+(imod-1)*3); end
plot(kj,cr(imod,:)')
hold on; plot(kj(ikp),cr(imod,ikp)','ko') 
xlabel('k (2pi/40000km)')
title('cr (ms-1)')
xlim([min(kj),max(kj)])
end
% column2=freq
for imod=1:nmod
if showgrowth==1; subplot(nmod,4,2+(imod-1)*4); else subplot(nmod,3,2+(imod-1)*3); end
plot(kj,wcpd(imod,:)')
hold on; plot(kj(ikp),wcpd(imod,ikp)','ko') 
xlabel('k (2pi/40000km)')
title('wr (cpd)')
xlim([min(kj),max(kj)])
end
% column3=eigenvects
ithick=2;
for imod=1:nmod
if showgrowth==1; subplot(nmod,4,3+(imod-1)*4); else subplot(nmod,3,3+(imod-1)*3); end
contourf(kj,1:nxs,squeeze(vectshow(:,imod,:)),'LineStyle','none')
hold on; plot(kj,kj*0+1+0.5,'k-','Linewidth',ithick); % separation line K-Rm
hold on; plot(kj,kj*0+1+nrm+0.5,'k-','Linewidth',ithick); % separation line Rm-Am
hold on; plot(kj,kj*0+1+nrm+nyk+0.5,'k-','Linewidth',ithick); % separation line Am-Zm
hold on; plot([0,0],[1,nxs],'k-','Linewidth',ithick); % vertical line k=0
xlabel('k (2pi/40000km)')
ylabel('K-Rm-Am-Zm')
title('Vects')
xlim([min(kj),max(kj)])
end
% column4=growth (optionnal)
if showgrowth==1;
for imod=1:nmod
subplot(nmod,4,4+(imod-1)*4)
plot(kj,wgcpd(imod,:)')
hold on; plot(kj(ikp),wgcpd(imod,ikp)','ko') 
xlabel('k (2pi/40000km)')
title('wi (cpd)')
xlim([min(kj),max(kj)])
end
end



% Do a k-w plot with all eigenmodes
if 0==1
figure(iwind); iwind=iwind+1;, clf;
ithick=2;
%
for imod=1:6
plot(kj,wcpd(imod,:)','Linewidth',ithick)
hold on; plot(kj(ikp),wcpd(imod,ikp)','ko') 
hold on
end
plot([0,0],[0,0.5],'k-','Linewidth',ithick);
xlabel('k (2pi/40000km)')
ylabel('w (cpd)')
title('k-w')
xlim([min(kj),max(kj)])
ylim([0,0.5])
end







end% graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


