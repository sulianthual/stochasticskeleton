% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% regress variables on a signal (data proj) at 180E for example: does covariance or correlation maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%
% variable on which is regressed
xpos=searchclosest(xg,20);% in 1000 km 
yinkm=     0            ; % position y to look in 1000 km
changeypos=        0        ; % change the yk axis for looking ypos (not closest strip but exact position)
varref=    21          ; % 21 for u, 22 for Ha, 23 for q, 24 for theta, otherwise data projections
dofilterkw2=     1          ; % filter within k-w the variable on which is regressed (only for >=21)
kmin2=1; kmax2=3;% in (2pi/40000km)
wmin2=1/90; wmax2=1/30; %
%
%%%%%%%%%%%%
timelag=0; % in days, lag of regressed vs on which regressed (beware it is inverted !)
% timelag=15;
%  timelag=5;
%  timelag=-5;
% timelag=-15;

%%%%%%%%%%%%
% variable that is regressed (doing a x-y map)
varreg=             25               ; % 21 for u, 22 for Ha, 23 for q, 24 for theta, 25 for v otherwise data projections
dofilterkw3=     1          ; % filter within k-w the variable on which is regressed
kmin3=1; kmax3=3;% in (2pi/40000km)
wmin3=1/90; wmax3=1/30; %
%
%% type of y used (for the regressed variable)
changeyk=        1        ; % change the yk axis to yka
yka=-4:0.5:4; % modification of axis for plotting ! 
yyg=yka*1500/1000; % plotting y axis (in 1000 km, must correspond to yka),
%
%%%%%%%%%%%%
dorecalc=   1    ; % recompute
anatype=   3     ; %(3=ref,1 also good) % type of analysis:
%                  1=correlation,2=covariance, 3=regression coeff (<ab>/<a2>)
addconf=    1     ; % add analysis of confidence 95%
%
%%%%%%%%%%%%
% contour
dolevs=1;
if anatype==1; levelsco=(-1:0.1:1)*0.8; end % correlation
if anatype==2; levelsco=(-1:0.1:1)*0.8; end % covariance
if anatype==3; levelsco=(-1:0.1:1)*0.2; end % regression (u others)
if anatype==3; if varreg==24; levelsco=(-1:0.1:1)*0.1;end; end  % levels for regression (u theta)
if anatype==3; if varreg==21; levelsco=(-1:0.1:1)*1;end; end  % levels for regression (u u)
if anatype==3; if varreg==22; levelsco=(-1:0.1:1)*0.2;end; end  % levels for regression (u Ha)
if anatype==3; if varreg==23; levelsco=(-1:0.1:1)*0.2;end; end  % levels for regression (u q)
if anatype==3; if varreg==25; levelsco=(-1:0.1:1)*0.1;end; end  % levels for regression (u v)
if anatype==1; levelsco=(-1:0.1:1)*0.8; end % correlation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get ref signals (on which is regressed)
%
if dorecalc==1

%%%%%%%%%%%%
% data proj
if varref<21
xcon=zeros(nx,2); 
for ilpps=ilpmin:ilpmax
indexrestart=ilpps; run(fileini);% get all infos, notably the restart file
imodshow=varref;
fileproj=strcat(dfolder,'/dataproj_',num2str(indexrestart),'_',num2str(imodshow), '.nc');% file to write
term=ncdfgetvar(fileproj,'xproj');
xcon=cat(2,xcon,term);
end
% remove first dummy
xcon=xcon(:,3:end);
% take only at 180E
xcon=squeeze(xcon(xpos,:));
end

%%%%%%%%%%%%
% other variables
if varref>=21
xcon=zeros(nx,nyk,2); 
for ilpps=ilpmin:ilpmax
indexrestart=ilpps; run(fileini);% get all infos, notably the restart file
% u
if varref==21
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=ua*singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim); end
xcon=cat(3,xcon,term);
end

% Ha
if varref==22
etas=ncdfgetvar(fileout,'etas');
term=HH*DDA*oa/(ta/oneday)*etas; 
xcon=cat(3,xcon,term);
end
%
% q
if varref==23
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
Zs=ncdfgetvar(fileout,'Zs');
os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term=oa*(Zs-QQ*os);
xcon=cat(3,xcon,term);
end
%
% theta
if varref==24
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
Zs=ncdfgetvar(fileout,'Zs');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=oa*singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
xcon=cat(3,xcon,term);
end

end% ilpps
%%%%%%%%%
%
% remove first dummy
xcon=xcon(:,:,3:end);
% Filter kw (only modify at jshow) 
if dofilterkw2==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin2); iwmax=searchclosest(wg,wmax2); 
ikmin=searchclosest(kg,kmin2); ikmax=searchclosest(kg,kmax2); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyk;
pass=filterkw(squeeze(xcon(:,iy,:)),dxg,dtg,[kmin2,kmax2],[wmin2,wmax2]);
xcon(:,iy,:)=real(pass);
end; pass=0;
end% end filter

% after eventual filtering, take only at 180E
% if necessary, compute not on strip but on reconstructed y
if changeypos== 1        
yinadim=yinkm*1000/1500;% adim position
psima=zeros(1,nym); for im=1:nym; psima(:,im)=hermitefunc(im-1,yinadim); end;%new hermites
ntt=size(xcon); ntt=ntt(3);
term=zeros(nx,ntt);
for kt=1:ntt
pass=singlecolumnfm(squeeze(xcon(:,:,kt)),psim,Hyk);
term(:,kt)=singlecolumnf(pass,psima);
end
xcon=term; term=0;
%)
xcon=squeeze(xcon(xpos,:));
size(xcon)
%
else
ypos=searchclosest(yg,yinkm);% in 1000 km 
xcon=squeeze(xcon(xpos,ypos,:))';
end
%
end% >=21
%

end%dorecalc 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get regressed signal

if dorecalc==1
% other variables
xcon1=zeros(nx,nyk,2); 
%%%%%%%%%%
for ilpps=ilpmin:ilpmax
indexrestart=ilpps; run(fileini);% get all infos, notably the restart file
%
% u
if varreg==21
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=ua*singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim); end
xcon1=cat(3,xcon1,term);
end

% Ha
if varreg==22
etas=ncdfgetvar(fileout,'etas');
term=HH*DDA*oa/(ta/oneday)*etas; 
xcon1=cat(3,xcon1,term);
end
%
% q
if varreg==23
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
Zs=ncdfgetvar(fileout,'Zs');
os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term=oa*(Zs-QQ*os);
xcon1=cat(3,xcon1,term);
end
%
% theta
if varreg==24
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
Zs=ncdfgetvar(fileout,'Zs');
term=zeros(nx,nyk,nts);
for kts=1:nts; term(:,:,kts)=oa*singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
xcon1=cat(3,xcon1,term);
end

% v
if varreg==25
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas');
term=zeros(nx,nyk,nts);
for kts=1:nts;
spass=seasowarmpool(ts(kts),yseaso,pseaso,wpint,soref,yk,nx,nyk,LL,xx,1);
Sm=singlecolumnfm(HH*DDA*etas(:,:,kts)-spass,psim,Hyk);
%Sm=singlecolumnfm(HH*DDA*etas(:,:,kts)-so,psim,Hyk);
term(:,:,kts)=ua*singlecolumnv(Rms(:,:,kts),Sm,psim,dx); 
end
xcon1=cat(3,xcon1,term);
end

end% ilpps
%%%%%%%%%
%
% remove first dummy
xcon1=xcon1(:,:,3:end);
% Filter kw (only modify at jshow) 
if dofilterkw3==1
kg=fftkspe(nx,dx)/(2*pi)/xa*40000*1000; % (2pi/40000km)
wg=fftkspe(nts,dt*mts)/(2*pi)/ta*oneday;% (cpd)
dxg=dx*xa/40000/1000;
dtg=dt*mts*ta/oneday;
iwmin=searchclosest(wg,wmin3); iwmax=searchclosest(wg,wmax3); 
ikmin=searchclosest(kg,kmin3); ikmax=searchclosest(kg,kmax3); 
if iwmin>=iwmax; term=iwmax; iwmax=iwmin; iwmin=term; end
if ikmin>=ikmax; term=ikmax; ikmax=ikmin; ikmin=term; end
for iy=1:nyk;
pass=filterkw(squeeze(xcon1(:,iy,:)),dxg,dtg,[kmin3,kmax3],[wmin3,wmax3]);
xcon1(:,iy,:)=real(pass);
end; pass=0;
end% end filter
end% dorecalc
%
if changeyk==1
nyka=length(yka);
psima=zeros(nyka,nym); for im=1:nym; psima(:,im)=hermitefunc(im-1,yka); end;%new hermites
ntt=size(xcon1); ntt=ntt(3);
term=zeros(nx,nyka,ntt);
for kt=1:ntt
pass=singlecolumnfm(squeeze(xcon1(:,:,kt)),psim,Hyk);
term(:,:,kt)=singlecolumnf(pass,psima);
end
xcon1=term; term=0;
else
yka=yk; nyka=nyk;
yyg=yka*1500/1000; % plotting y axis (in 1000 km)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time lag
if timelag~=0
ts=ncdfgetvar(fileout,'ts'); dts=ts(2)-ts(1);
itt=fix(timelag/dts);% lag in index
if itt>0
xcon=xcon(1+itt:end);
xcon1=xcon1(:,:,1:end-itt);
else
xcon=xcon(1:end+itt);
xcon1=xcon1(:,:,1-itt:end);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute covariance map
'compute reg'
comap=zeros(nx,nyka);
confmap=zeros(nx,nyka);
for ix=1:nx
for iy=1:nyka
%
% correlation
if anatype==1
term=xcon; term=detrend(term,'constant');% remove mean
term1=squeeze(xcon1(ix,iy,:));term1=detrend(term1,'constant');% remove mean
nnt=length(term); comap(ix,iy)=sum(term.*term1')/sqrt(sum(term.*term))/sqrt(sum(term1.*term1));% covariance
end

% covariance
if anatype==2
term=xcon; term=detrend(term,'constant');% remove mean
term1=squeeze(xcon1(ix,iy,:));term1=detrend(term1,'constant');% remove mean
comap(ix,iy)=sum(term.*term1')/(nnt-1);% covariance
end
%
% regression coeff
if anatype==3
term=xcon; term=detrend(term,'constant');% remove mean
term1=squeeze(xcon1(ix,iy,:));term1=detrend(term1,'constant');% remove mean
comap(ix,iy)=sum(term.*term1')/sum(term.*term);
%nnt=length(term); comap(ix,iy)=sum(term.*term1')/(nnt-1);% covariance
end
%
% add analysis of confidence interval
if addconf==1
term=xcon; term=detrend(term,'constant');% remove mean
term1=squeeze(xcon1(ix,iy,:));term1=detrend(term1,'constant');% remove mean
%pass=corrcoef([term term1'],'alpha',pass1);
[pass, conf]=corrcoef(term,term1');
%comap(ix,iy)=pass(1,2);
conf=conf(1,2);
if conf<=0.05; confmap(ix,iy)=1; else confmap(ix,iy)=0; end; % 0.05 is 95%conf interval
end
%
end;end;% ix iy
%term=cov(xmat1',xmat2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contour cov map

figure(iwind); iwind=iwind+1;, clf;
term=comap;
if dolevs; [c1,h1]=contourf(xg,yyg,term',levelsco,'LineStyle','none');
caxis([levelsco(1) levelsco(end)]);% VERY IMPORTANT !!!!
else [c1,h1]=contourf(xg,yyg,term','LineStyle','none'); end
colorbar('location','southoutside');
xlabel('x 1000km')
ylabel(' y 1000 km')
if addconf==1
hold on; [c2,h2]=contour(xg,yyg,confmap',0.5,'w-');
end


