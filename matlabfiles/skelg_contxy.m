% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% contour x-y of successive snapshots
% does not work for M=1 (cannot contour y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
% Time (here has to choose the restart file and corresponding time)
%indexrestart= 39 ; timeshow= 77050 ; % restart file and time to show
indexrestart= 2 ; timeshow= 3500 ; % restart file and time to show
%
steps=[0,1,2,3,4,5]; steps=fliplr(steps); % steps to show (at time nstep*dtstep+timeshow)
dtstep=10;% spacing between steps
run(fileini);
%
% Variables
doderived=       0      ; % rather plot some circulation quantities (div, curl....)
dofilterkw=     1          ; % filter within k,w, for the hovmullers (ref=0) 
kmin=1; kmax=3; wmin=1/90; wmax=1/30; % in cpd, for hovmullers filtering
%
% grids
xxrange=[0 40];
xpos=20;% add some straight lines there
dxpos=20; 
changeyk=        1        ; % change the yk axis to yka
yka=-4:0.5:4; % modification of axis for plotting ! 
yyg=yka*1500/1000; % plotting y axis (in 1000 km, must correspond to yka),
yyrange=[-6,6];
%
% graph
dolevs=          1        ; % put levels on contours, ref=1 
levelsu=     (-1:0.1:1)*2  ; % levels a for contours (m.s-1), dograph2
levelsv=     (-1:0.1:1)*0.5  ; % levels a for contours (m.s-1), dograph2
levelso=     (-1:0.1:1)*0.2  ; % levels a for contours (K), dograph2
levelsq=     (-1:0.1:1)*0.5  ; % levels a for contours (K), dograph2
levelsHa=    (-1:0.1:1)*0.5  ; % levels a for contours (K.d-1), dograph2
doquiver=    1    ; % add quiver in all graphs (quiver params within code)
docolorbar=        1          ; % add colorbar
if doderived==1; doquiver=0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Outputs 
['read']
Ks=ncdfgetvar(fileout,'Ks');
Rms=ncdfgetvar(fileout,'Rms');
etas=ncdfgetvar(fileout,'etas');
Zs=ncdfgetvar(fileout,'Zs');
ts=ncdfgetvar(fileout,'ts');
% Compute other variables
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


% Filter kw 
if dofilterkw==1
['filter']
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
osrr=os; usrr=us; vsrr=vs; qsrr=qs; etasrr=etas;
usrr0=us;
vsrr0=vs;
%
% compute derivatives and associated quantities (replace theta and q in graphs)
if doderived==1
['derivatives']
dxus=zeros(nx,nyk,nts); dxvs=zeros(nx,nyk,nts);
dyus=zeros(nx,nyk,nts); dyvs=zeros(nx,nyk,nts);
for kts=1:nts
pass=squeeze(us(:,:,kts));
dxus(:,:,kts)=singlecolumndx(pass,xg);% note there is no need to dimensionalize is xg is dimensional
pass=singlecolumnfm(pass,psim,Hyk); pass=singlecolumndy(pass);
dyus(:,:,kts)=singlecolumnf(pass,psim)*(onekm*1000/ya); % need to dimensionalize and put in (cm.m-1)(1000km-1)
pass=squeeze(vs(:,:,kts));
dxvs(:,:,kts)=singlecolumndx(pass,xg);% note there is no need to dimensionalize is xg is dimensional
pass=singlecolumnfm(pass,psim,Hyk); pass=singlecolumndy(pass);
dyvs(:,:,kts)=singlecolumnf(pass,psim)*(onekm*1000/ya); % need to dimensionalize and put in (cm.m-1)(1000km-1)
end
% finish to dimensionalize
dxus=ua*dxus; dyus=ua*dyus;
dxvs=ua*dxvs; dyvs=ua*dyvs;
% Compute divergence
divuv=dxus+dyvs;
curluv=dxvs-dyus;
% see contributions to div + div+curl + Ha
usrr=dxus; levelsu=     (-1:0.1:1)*0.4  ; % dxu
vsrr=dyvs; levelsv=     (-1:0.1:1)*0.4  ; % dyv
osrr=divuv; levelso=     (-1:0.1:1)*0.4  ; % div
qsrr=curluv; levelsq=     (-1:0.1:1)*1.5  ; % curl
end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs of run: hovmullers
%
figure(iwind); iwind=iwind+1;, clf;
icount=1;
nstep=length(steps);
%
% loop on each step
for istep=1:nstep;
['step=',num2str(istep),'/',num2str(nstep)]
iistep=steps(istep);
ktss=searchclosest(ts,timeshow+dtstep*iistep);

% TAKE AT TIME KTS
os=osrr(:,:,ktss);
us=usrr(:,:,ktss);
vs=vsrr(:,:,ktss);
qs=qsrr(:,:,ktss);
etas=etasrr(:,:,ktss);
us0=usrr0(:,:,ktss);
vs0=vsrr0(:,:,ktss);
% squeeze everything
os=squeeze(os);
qs=squeeze(qs);
etas=squeeze(etas);
us=squeeze(us);
vs=squeeze(vs);
us0=squeeze(us0);
vs0=squeeze(vs0);
%
% change the yk for contouring
etas0=etas;
if changeyk==1
nyka=length(yka);
psima=zeros(nyka,nym); for im=1:nym; psima(:,im)=hermitefunc(im-1,yka); end;%new hermites
us=singlecolumnfm(us,psim,Hyk);
os=singlecolumnfm(os,psim,Hyk);
qs=singlecolumnfm(qs,psim,Hyk);
vs=singlecolumnfm(vs,psim,Hyk);
etas=singlecolumnfm(etas,psim,Hyk);
us0=singlecolumnfm(us0,psim,Hyk);
vs0=singlecolumnfm(vs0,psim,Hyk);
us=singlecolumnf(us,psima);
os=singlecolumnf(os,psima);
qs=singlecolumnf(qs,psima);
vs=singlecolumnf(vs,psima);
etas=singlecolumnf(etas,psima);
us0=singlecolumnf(us0,psima);
vs0=singlecolumnf(vs0,psima);
else
yka=yk;
yyg=yka*1500/1000; % plotting y axis (in 1000 km)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%5
% do graphs
%
iithick=2;
ngraphs=5;
%
%
% QUIVER DEFINITIONS
if doquiver==1
if doderived==1; term1=ua*us0; term2=ua*vs0; else; term1=ua*us; term2=ua*vs; end
%skip1=8; skip2=2; 
skip1=6; skip2=1; scale1=2; scale2=scale1; 
ithick=1;
ixgs=1:skip1:length(xg); iygs=1:skip2:length(yyg); 
scale1=4; scale2=scale1; 
%s;scaale=2; % irrelevant now
end
%
% u
subplot(nstep,ngraphs,icount); icount=icount+1;
if doderived==1; term=us; else; term=ua*us; end
if dolevs; contourf(xg,yyg,term',levelsu,'LineStyle','none');
caxis([levelsu(1) levelsu(end)]);% VERY IMPORTANT !!!!
else contourf(xg,yyg,term','LineStyle','none'); end
%ylim([-6,6]);
if (istep==nstep && docolorbar==1); 
p=get(gca,'position'); % save position
colorbar('location','southoutside'); 
set(gca,'position',p); % restore position
end
%if istep==nstep; xlabel('x(1000km)'); else set(gca,'XTick',[]); end
pass=['t=',num2str(fix(ts(ktss)))];
ylabel(pass)
if istep==1; 
title('u'); 
if doderived==1; title('dxu'); end
end 
%
if doquiver==1 % add quiver
hold on; quiver(xg(ixgs),yyg(iygs),scale1*term1(ixgs,iygs)',scale2*term2(ixgs,iygs)',...
'AutoScale','off','linewidth',ithick,'color','k')
end
%
% v
subplot(nstep,ngraphs,icount); icount=icount+1;
if doderived==1; term=vs; else; term=ua*vs; end
if dolevs; contourf(xg,yyg,term',levelsv,'LineStyle','none');
caxis([levelsv(1) levelsv(end)]);% VERY IMPORTANT !!!!
else contourf(xg,yyg,term','LineStyle','none'); end
ylim(yyrange);
if (istep==nstep && docolorbar==1); 
p=get(gca,'position'); % save position
colorbar('location','southoutside'); 
set(gca,'position',p); % restore position
end
  if istep==nstep; %xlabel('x(1000km)'); 
  else set(gca,'XTick',[]); end; 
set(gca,'YTick',[]);
if istep==1; 
title('v'); 
if doderived==1; title('dyv'); end
end 
hold on; plot([xpos-dxpos,xpos-dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
hold on; plot([xpos+dxpos,xpos+dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
xlim(xxrange);
if doquiver==1 % add quiver
hold on; quiver(xg(ixgs),yyg(iygs),scale1*term1(ixgs,iygs)',scale2*term2(ixgs,iygs)',...
'AutoScale','off','linewidth',ithick,'color','k')
end
%
% o
subplot(nstep,ngraphs,icount); icount=icount+1;
if doderived==1; term=os; else; term=oa*os; end
if dolevs; contourf(xg,yyg,term',levelso,'LineStyle','none');
caxis([levelso(1) levelso(end)]);% VERY IMPORTANT !!!!
else contourf(xg,yyg,term','LineStyle','none'); end
ylim(yyrange);
%calc_goodyaxis();
if (istep==nstep && docolorbar==1); 
p=get(gca,'position'); % save position
colorbar('location','southoutside'); 
set(gca,'position',p); % restore position
end
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
if istep==nstep; %xlabel('x(1000km)'); 
else set(gca,'XTick',[]); end
%ylabel('y(1000km)')
if istep==1; 
title('theta'); 
if doderived==1; title('div'); end
end 
hold on; plot([xpos-dxpos,xpos-dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
hold on; plot([xpos+dxpos,xpos+dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
xlim(xxrange);
if doquiver==1 % add quiver
hold on; quiver(xg(ixgs),yyg(iygs),scale1*term1(ixgs,iygs)',scale2*term2(ixgs,iygs)',...
'AutoScale','off','linewidth',ithick,'color','k')
end
%
% q
subplot(nstep,ngraphs,icount); icount=icount+1;
if doderived==1; term=qs; else; term=oa*qs; end
if dolevs; contourf(xg,yyg,term',levelsq,'LineStyle','none');
caxis([levelsq(1) levelsq(end)]);% VERY IMPORTANT !!!!
else contourf(xg,yyg,term','LineStyle','none'); end
ylim(yyrange);
%calc_goodyaxis();
if (istep==nstep && docolorbar==1); 
p=get(gca,'position'); % save position
colorbar('location','southoutside'); 
set(gca,'position',p); % restore position
end
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
if istep==nstep; %xlabel('x(1000km)'); 
else set(gca,'XTick',[]); end
%ylabel('y(1000km)')
if istep==1; 
title('q'); 
if doderived==1; title('curl'); end
end 
%hold on; plot(xg,xg*0,'k-','linewidth',iithick);
%hold on; plot(xg,xg*0,'k-','linewidth',iithick);
hold on; plot([xpos-dxpos,xpos-dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
hold on; plot([xpos+dxpos,xpos+dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
xlim(xxrange);
if doquiver==1 % add quiver
hold on; quiver(xg(ixgs),yyg(iygs),scale1*term1(ixgs,iygs)',scale2*term2(ixgs,iygs)',...
'AutoScale','off','linewidth',ithick,'color','k')
end
%axis off;
%
% etas
subplot(nstep,ngraphs,icount); icount=icount+1;
term=HH*DDA*oa/(ta/oneday)*etas; 
if dolevs; contourf(xg,yyg,term',levelsHa,'LineStyle','none');
caxis([levelsHa(1) levelsHa(end)]);% VERY IMPORTANT !!!!
else contourf(xg,yyg,term','LineStyle','none'); end
ylim(yyrange);
%calc_goodyaxis();
if (istep==nstep && docolorbar==1); 
p=get(gca,'position'); % save position
colorbar('location','southoutside'); 
set(gca,'position',p); % restore position
end
%set(gca,'XTick',[]); 
set(gca,'YTick',[]);
if istep==nstep; %xlabel('x(1000km)'); 
else set(gca,'XTick',[]); end
%ylabel('y(1000km)')
if istep==1; title('Ha'); end
%hold on; plot(xg,xg*0,'k-','linewidth',iithick);
hold on; plot([xpos-dxpos,xpos-dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
hold on; plot([xpos+dxpos,xpos+dxpos],[min(yyg),max(yyg)]+[-1,1],'k-','linewidth',iithick);
xlim(xxrange);
if doquiver==1 % add quiver
hold on; quiver(xg(ixgs),yyg(iygs),scale1*term1(ixgs,iygs)',scale2*term2(ixgs,iygs)',...
'AutoScale','off','linewidth',ithick,'color','k')
end
%
end% istep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

