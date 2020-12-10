% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% contour x-y of mean fields (averaged over one restart file)
% does not work for M=1 (cannot contour y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
% type to do
dotype=       2            ; % 1=time-mean, 2=time-std
%
% filter before
dofilterkw=     1          ; % filter within k,w, for the hovmullers (ref=0) 
kmin=1; kmax=3; wmin=1/90; wmax=1/30; % in cpd, MJO band
%kmin=-3; kmax=-1; wmin=1/90; wmax=1/30; % in cpd, WEST BAND
%
% Levels
dolevs=          1        ; % put levels on contours, ref=1 
if dotype==1% time-mean levels
levelsu=     (-1:0.1:1)*0.4  ; % levels a for contours (m.s-1), dograph2
levelso=     (0:0.1:1)*0.2  ; % levels a for contours (K), dograph2
levelsq=     (-1:0.1:1)*0.4  ; % levels a for contours (K), dograph2
levelsHa=    (0:0.05:1)*2  ; % levels a for contours (K.d-1), dograph2
levelssq=    (0:0.05:1)*2  ; % levels a for contours (K.d-1), dograph2
levelsv=     (-1:0.1:1)*0.02  ; % levels a for contours (m.s-1), dograph2
levelsv=     (-1:0.1:1)*0.1  ; % levels a for contours (m.s-1), dograph2
end
if dotype==2% time-mean std
levelsu=     (0:0.1:1)*4  ; % levels a for contours (m.s-1), dograph2
levelso=     (0:0.1:1)*1  ; % levels a for contours (K), dograph2
levelsq=     (0:0.1:1)*2  ; % levels a for contours (K), dograph2
levelsHa=    (0:0.05:1)*2  ; % levels a for contours (K.d-1), dograph2
levelssq=    (0:0.05:1)*2  ; % levels a for contours (K.d-1), dograph2
levelsv=     (0:0.1:1)*6  ; % levels a for contours (m.s-1), dograph2
end
if dotype==2% time-mean std, when filtered
if dofilterkw==1
levelsu=     (0:0.1:1)*1  ; % levels a for contours (m.s-1), dograph2
levelso=     (0:0.1:1)*0.1  ; % levels a for contours (K), dograph2
levelsq=     (0:0.1:1)*0.2  ; % levels a for contours (K), dograph2
levelsHa=    (0:0.05:1)*0.3  ; % levels a for contours (K.d-1), dograph2
levelsv=     (0:0.1:1)*0.2  ; % levels a for contours (m.s-1), dograph2
end
end

% modification of axis for plotting 
yka=-4:0.5:4; 
yyg=yka*1500/1000; % plotting y axis (in 1000 km, must correspond to yka),
nyka=length(yka);
psima=zeros(nyka,nym); for im=1:nym; psima(:,im)=hermitefunc(im-1,yka); end;%new hermites
%
ngraphs=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters projection (do not touch)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs of run: hovmullers
%
% READING FILES: THIS HAS TO BE DONE FOR one variable at the time

% Compute time mean/std
mus=zeros(nx,nyka);
mos=zeros(nx,nyka);
mqs=zeros(nx,nyka);
metas=zeros(nx,nyka);
mvs=zeros(nx,nyka);

%%%%%%%%%%%%
% u 
'u'
xcon=zeros(nx,nyka,2); 
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim); end
term2=zeros(nx,nyka,nts);
for kts=1:nts; pass=singlecolumnfm(term(:,:,kts),psim,Hyk); term2(:,:,kts)=singlecolumnf(pass,psima); end
term=term2; term2=0;
xcon=cat(3,xcon,term);
end
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyka;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter
us=xcon; term=0;
if dotype==1; % time mean
for ix=1:nx; for jy=1:nyka
mus(ix,jy)=mean(us(ix,jy,:));
end; end
end
if dotype==2; % time std
for ix=1:nx; for jy=1:nyka
mus(ix,jy)=std(us(ix,jy,:));
end; end
end
us=0;
%%%%%%%%%%%%
% o 
'o'
xcon=zeros(nx,nyka,2); 
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term2=zeros(nx,nyka,nts);
for kts=1:nts; pass=singlecolumnfm(term(:,:,kts),psim,Hyk); term2(:,:,kts)=singlecolumnf(pass,psima); end
term=term2; term2=0;
xcon=cat(3,xcon,term);
end
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyka;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter
os=xcon; term=0;
if dotype==1; % time mean
for ix=1:nx; for jy=1:nyka
mos(ix,jy)=mean(os(ix,jy,:));
end; end
end
if dotype==2; % time std
for ix=1:nx; for jy=1:nyka
mos(ix,jy)=std(os(ix,jy,:));
end; end
end
os=0;
%%%%%%%%%%%%
% q 
'q'
xcon=zeros(nx,nyka,2); 
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
Zs=ncdfgetvar(fileout,'Zs');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term=Zs-QQ*term;
term2=zeros(nx,nyka,nts);
for kts=1:nts; pass=singlecolumnfm(term(:,:,kts),psim,Hyk); term2(:,:,kts)=singlecolumnf(pass,psima); end
term=term2; term2=0;
xcon=cat(3,xcon,term);
end
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyka;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter
qs=xcon; term=0; Zs=0;
if dotype==1; % time mean
for ix=1:nx; for jy=1:nyka
mqs(ix,jy)=mean(qs(ix,jy,:));
end; end
end
if dotype==2; % time std
for ix=1:nx; for jy=1:nyka
mqs(ix,jy)=std(qs(ix,jy,:));
end; end
end
qs=0;
%%%%%%%%%%%%
% Ha
'Ha'
xcon=zeros(nx,nyka,2); 
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
term=ncdfgetvar(fileout,'etas');
term2=zeros(nx,nyka,nts);
for kts=1:nts; pass=singlecolumnfm(term(:,:,kts),psim,Hyk); term2(:,:,kts)=singlecolumnf(pass,psima); end
term=term2; term2=0;
xcon=cat(3,xcon,term);
end
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyka;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter
etas=xcon; term=0; 
if dotype==1; % time mean
for ix=1:nx; for jy=1:nyka
metas(ix,jy)=mean(etas(ix,jy,:));
end; end
end
if dotype==2; % time std
for ix=1:nx; for jy=1:nyka
metas(ix,jy)=std(etas(ix,jy,:));
end; end
end
etas=0;
%%%%%%%%%%%%
% v
'v'
xcon=zeros(nx,nyka,2); 
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas');
term=zeros(nx,nyk,nts);
for kts=1:nts; 
spass=seasowarmpool(ts(kts),yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,1);
Sm=singlecolumnfm(HH*DDA*etas(:,:,kts)-spass,psim,Hyk);
term(:,:,kts)=singlecolumnv(Rms(:,:,kts),Sm,psim,dx); 
end
term2=zeros(nx,nyka,nts);
for kts=1:nts; pass=singlecolumnfm(term(:,:,kts),psim,Hyk); term2(:,:,kts)=singlecolumnf(pass,psima); end
term=term2; term2=0;
xcon=cat(3,xcon,term);
end
if dofilterkw==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin); iwmax=searchclosest(wg,wmax); 
ikmin=searchclosest(kg,kmin); ikmax=searchclosest(kg,kmax); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyka;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin,kmax],[wmin,wmax]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter
vs=xcon; term=0;
if dotype==1; % time mean
for ix=1:nx; for jy=1:nyka
mvs(ix,jy)=mean(vs(ix,jy,:));
end; end
end
if dotype==2; % time std
for ix=1:nx; for jy=1:nyka
mvs(ix,jy)=std(vs(ix,jy,:));
end; end
end
vs=0;
%%%%%%%%%%%%
%

figure(iwind); iwind=iwind+1;, clf;
subplot(1,ngraphs,1)
term=ua*mus;
if dolevs; contourf(xg,yyg,term',levelsu,'LineStyle','none');
else contourf(xg,yyg,term','LineStyle','none'); end
caxis([levelsu(1) levelsu(end)]);% VERY IMPORTANT !!!!
colorbar('location','southoutside');
xlabel('x(1000km)')
ylabel('y(1000km)')
title(' u ')
%
subplot(1,ngraphs,2)
term=oa*mos;
 if dolevs; contourf(xg,yyg,term',levelso,'LineStyle','none');
 else contourf(xg,yyg,term','LineStyle','none'); end
 caxis([levelso(1) levelso(end)]);% VERY IMPORTANT !!!!
%contourf(xg,yyg,term','LineStyle','none');
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('y(1000km)')
title('theta ')
%
subplot(1,ngraphs,3)
term=oa*mqs; 
if dolevs; contourf(xg,yyg,term',levelsq,'LineStyle','none');
else contourf(xg,yyg,term','LineStyle','none'); end
caxis([levelsq(1) levelsq(end)]);% VERY IMPORTANT !!!!
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('y(1000km)')
title('q ')
%axis off;
%
subplot(1,ngraphs,4)
term=HH*DDA*oa/(ta/oneday)*metas; 
 if dolevs; contourf(xg,yyg,term',levelsHa,'LineStyle','none');
 else contourf(xg,yyg,term','LineStyle','none'); end
caxis([levelsHa(1) levelsHa(end)]);% VERY IMPORTANT !!!!
%calc_goodyaxis();
colorbar('location','southoutside');
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
xlabel('x(1000km)')
%ylabel('y(1000km)')
title('Ha')
%axis off;
%
subplot(1,ngraphs,5)
term=ua*mvs;
if dolevs; contourf(xg,yyg,term',levelsv,'LineStyle','none');
else contourf(xg,yyg,term','LineStyle','none'); end
caxis([levelsv(1) levelsv(end)]);% VERY IMPORTANT !!!!
colorbar('location','southoutside');
xlabel('x(1000km)')
ylabel('y(1000km)')
title(' v ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
