% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Graphi 2: hovmuller y-t quick
% does hovmuller U,o,q,a non filtered, and MJO projection next to it
%
% NB: this code is intented for the truncation M=5 
%   for a truncation M=0 i.e. only m=0 (as in JAS 2013 article) then
%   v (which has there a m=1 component) is not computed. The code may be adapted to get it
% 
% y-t variables 
xshow=           20       ; % strip location for y-t contour (in 1000 km)
timeunits=       1        ; % time units factor (1 is day) for graph
tlim=        [0,2000]     ; % range on t for graphs (adding start time from netcdf)
ts=ncdfgetvar(fileout,'ts'); tlim=tlim+ts(1);
%
doseasog=        0         ; % add seasonal cycle in figures
seathick=       0.5        ; % tihckness for graphs
%
% meridional axis modify
changeyk=        1        ; % change the yk axis to yka
yka=-4:0.5:4; % modification of axis for plotting ! 
yyg=yka*1500/1000; % plotting y axis (in 1000 km, must correspond to yka),
%
% filter
dofilterkw=     1          ; % filter within k,w, for the hovmullers (ref=0) 
kmin=1; kmax=3; wmin=1/90; wmax=1/30; % in cpd, for hovmullers filtering
%
% smooth
dosmooth=   1; % simple smooth for graphs (faster contour)
nysmoo=     1; % for meridional smooth
ntsmoo=     5; % t size of smooth (odd), ref=5
%
% levels
dolevs=          1        ; % put levels on contours, ref=1 
levelsu=     (-1:0.1:1)*2  ; % levels a for contours (m.s-1), dograph2
levelso=     (-1:0.1:1)*0.5  ; % levels a for contours (K), dograph2
levelsq=     (-1:0.1:1)*0.5  ; % levels a for contours (K), dograph2
levelsHa=    (-1:0.1:1)*0.5  ; % levels a for contours (K.d-1), dograph2
levelsv=    (-1:0.1:1)*0.5  ; % levels a for contours (K.d-1), dograph2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs of run: hovmullers
%
% Read Outputs 
'reading 1'
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas');
Zs=ncdfgetvar(fileout,'Zs');
%
% Compute os and us
os=zeros(nx,nyk,nts); us=zeros(nx,nyk,nts);vs=zeros(nx,nyk,nts);
for kts=1:nts
os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim);
us(:,:,kts)=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim);
if doseaso==1
spass=seasowarmpool(ts(kts),yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,1);
Sm=singlecolumnfm(HH*DDA*etas(:,:,kts)-spass,psim,Hyk);
else
Sm=singlecolumnfm(HH*DDA*etas(:,:,kts)-so,psim,Hyk);
end
vs(:,:,kts)=singlecolumnv(Rms(:,:,kts),Sm,psim,dx); 
end
qs=Zs-QQ*os;
%
% Filter kw 
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for ik=1:nyk; 
Zs(:,ik,:)=filterkw(squeeze(Zs(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]); 
etas(:,ik,:)=filterkw(squeeze(etas(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
us(:,ik,:)=filterkw(squeeze(us(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
vs(:,ik,:)=filterkw(squeeze(vs(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
os(:,ik,:)=filterkw(squeeze(os(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
qs(:,ik,:)=filterkw(squeeze(qs(:,ik,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
end
Zs=real(Zs); 
etas=real(etas);
us=real(us);
vs=real(vs);
os=real(os);
qs=real(qs);
end% end filter
%
% arrange arrays
if changeyk==0
yka=yk; yyg=yka*1500/1000;
[ishow,term]=searchclosest(xg,xshow);
us=us(ishow,:,:);
os=os(ishow,:,:);
qs=qs(ishow,:,:);
etas=etas(ishow,:,:);
vs=vs(ishow,:,:);
us=squeeze(us);
os=squeeze(os);
etas=squeeze(etas);
qs=squeeze(qs);
vs=squeeze(vs);
end
etas0=etas; % save for latter
%
if changeyk==1
'changing yk'
[ishow,term]=searchclosest(xg,xshow);
nyka=length(yka);
psima=zeros(nyka,nym); for im=1:nym; psima(:,im)=hermitefunc(im-1,yka); end;%new hermites
usa=zeros(nyka,nts);
osa=zeros(nyka,nts);
qsa=zeros(nyka,nts);
etasa=zeros(nyka,nts);
vsa=zeros(nyka,nts);
for kts=1:nts
pass=singlecolumnfm(reshape(us(ishow,:,kts),1,nyk),psim,Hyk);
usa(:,kts)=singlecolumnf(pass,psima);
pass=singlecolumnfm(reshape(os(ishow,:,kts),1,nyk),psim,Hyk);
osa(:,kts)=singlecolumnf(pass,psima);
pass=singlecolumnfm(reshape(qs(ishow,:,kts),1,nyk),psim,Hyk);
qsa(:,kts)=singlecolumnf(pass,psima);
pass=singlecolumnfm(reshape(etas(ishow,:,kts),1,nyk),psim,Hyk);
etasa(:,kts)=singlecolumnf(pass,psima);
pass=singlecolumnfm(reshape(vs(ishow,:,kts),1,nyk),psim,Hyk);
vsa(:,kts)=singlecolumnf(pass,psima);
end
us=squeeze(usa); usa=0;
os=squeeze(osa); osa=0;
etas=squeeze(etasa); etasa=0;
qs=squeeze(qsa); qsa=0;
vs=squeeze(vsa); vsa=0;
end
%
% smooth 
if dosmooth==1 
'smoothing'
us=calc_smoothn(us,[nysmoo,ntsmoo]); 
os=calc_smoothn(os,[nysmoo,ntsmoo]); 
etas=calc_smoothn(etas,[nysmoo,ntsmoo]); 
qs=calc_smoothn(qs,[nysmoo,ntsmoo]); 
vs=calc_smoothn(vs,[nysmoo,ntsmoo]); 
end
%
% compute warm pool center y position
if doseasog==1; 
ynow=seasowarmpool(ts,yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,2); 
ynow=ynow*ya/1000/1000; % convert to 1000km
end
%
% Hovmuller Contour Physical (on one hermite): u,o,q,Ha (with dims)
figure(iwind); iwind=iwind+1;, clf;
'graph variables'
%colormap('jet')
%colormap('gray')
ngraphs=5;
subplot(1,ngraphs,1)
term=ua*us;
if dolevs; contourf(yyg,ts/timeunits,term',levelsu,'LineStyle','none');
else contourf(yyg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsu(1) levelsu(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits); 
colorbar('location','southoutside');
xlabel('y(1000km)')
ylabel('time (days)')
title('u')
if doseasog==1; hold on; plot(ynow,ts,'k-','Linewidth',seathick); end
%
subplot(1,ngraphs,2)
term=oa*os;
if dolevs; contourf(yyg,ts/timeunits,term',levelso,'LineStyle','none');
else contourf(yyg,ts/timeunits,term','LineStyle','none'); end
caxis([levelso(1) levelso(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('y(1000km)')
%ylabel('t (days)')
title('theta')
if doseasog==1; hold on; plot(ynow,ts,'k-','Linewidth',seathick); end
%
subplot(1,ngraphs,3)
term=oa*qs; 
if dolevs; contourf(yyg,ts/timeunits,term',levelsq,'LineStyle','none');
else contourf(yyg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsq(1) levelsq(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('y(1000km)')
%ylabel('t (days)')
title('q')
if doseasog==1; hold on; plot(ynow,ts,'k-','Linewidth',seathick); end
%axis off;
%
subplot(1,ngraphs,4)
term=HH*DDA*oa/(ta/oneday)*etas; 
if dolevs; contourf(yyg,ts/timeunits,term',levelsHa,'LineStyle','none');
else contourf(yyg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsHa(1) levelsHa(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits);
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('y(1000km)')
%ylabel('t (days)')
title('Ha')
if doseasog==1; hold on; plot(ynow,ts,'k-','Linewidth',seathick); end
%axis off;

subplot(1,ngraphs,5)
term=ua*vs;
if dolevs; contourf(yyg,ts/timeunits,term',levelsv,'LineStyle','none');
else contourf(yyg,ts/timeunits,term','LineStyle','none'); end
caxis([levelsv(1) levelsv(end)]);% VERY IMPORTANT !!!!
ylim(tlim/timeunits); 
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('y(1000km)')
%ylabel('t (days)')
title('v')
if doseasog==1; hold on; plot(ynow,ts,'k-','Linewidth',seathick); end
%axis off;
