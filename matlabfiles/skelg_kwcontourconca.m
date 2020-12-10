% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% Graph: panel of k-w contours, all in one code (with concatenation on time files)
% here shows four variables u, o, q, Ha in one graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
yshow=           0        ; % strip location (yg-1000km) to show (only for skeleton run) (ref=0)
sshowstab =    1   ; % show stability dispersion curves
%
% data treatment
removemean=1; % remove temporal mean (necessary for WP)
donormakw=0; % normalize by global standard dev first (ref=0)
dosmooth=   1; % simple smooth for graphs (faster contour) (ref=1) (done before undersampling)
nxsmoo=     1; % x size of smooth (odd) (ref=1)
ntsmoo=     25; % t size of smooth (odd) (ref=25)
undersamplext=    1 ; % undersample outputs to be faster
nku=1; % k-space between samples 
nwu=11; % w-space between samples 

% graph ranges
dokwlevs=  1; % do levels (is log10) (ref=1)
xxrange=[-5,5]; yyrange=[0,0.1]; %range of graphs (ref= -5,5 and 0., 0.1)
%xxrange=[-10,10]; yyrange=[0,0.2]; 
%
% others
doasym=        0            ; % rather do asymmetric (diff >yshow minus <yshow)
dor0=0; % test specfic (ref=0)
iwindspe=10; 
figure(iwindspe); clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop on each variable
for vvarsshow=1:4
%vvarsshow=     4           ; % 1=u, 2= o, 3= q, 4= Ha, 5=ntaus, is provided here to kw-contour
if vvarsshow==1; kwlevels=-10:0.5:-2; end % panel of u (non norma,10000days), 
if vvarsshow==2; kwlevels=-11:0.5:-3; end% panel of o (non norma,10000days) 
if vvarsshow==3; kwlevels=-11:0.5:-3; end% panel of q (non norma,10000days)
if vvarsshow==4; kwlevels=-11:0.5:-3; end% panel of Ha (non norma,10000days)
% beware (e.g. for u), a 0.5 step is ok, but a 0.25 is very slow
%
[jshow,term]=searchclosest(yg,yshow);% selected location for graphs-print
npanel1=1; npanel2=1;
icount=0;
%for ip2=1:1; % if only DA
for ip2=1:npanel2
for ip1=1:npanel1;% 
%[ip1 npanel1 ip2 npanel2]
icount=icount+1; % for position in graph
% Read variable to concatenate
xcon=((1:nx)*0)'; tcon=[0,0];% first timestep is dummy
for ilpps=ilpmin: ilpmax
indexrestart=ilpps; run(fileini);% get all infos
ts=ncdfgetvar(fileout,'ts'); nts=length(ts);
if vvarsshow==1% u
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); us=zeros(nx,nyk,nts);
for kts=1:nts; us(:,:,kts)=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim); end
if doasym==1; 
term=ua*squeeze( mean(us(:,jshow+1:end,:),2)-mean(us(:,1:jshow-1,:),2) );
else; term=ua*squeeze(us(:,jshow,:)); end
end
if vvarsshow==2% o
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term=oa*squeeze(os(:,jshow,:));
if doasym==1; 
term=oa*squeeze( mean(os(:,jshow+1:end,:),2)-mean(os(:,1:jshow-1,:),2) );
else; term=oa*squeeze(os(:,jshow,:)); end
end
if vvarsshow==3% q
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
Zs=ncdfgetvar(fileout,'Zs'); qs=(Zs-QQ*os); 
if doasym==1; 
term=oa*squeeze( mean(qs(:,jshow+1:end,:),2)-mean(qs(:,1:jshow-1,:),2) );
else; term=oa*squeeze(qs(:,jshow,:)); end
end
if vvarsshow==4% Ha
etas=ncdfgetvar(fileout,'etas'); 
has=HH*DDA*(oa/(ta/oneday))*etas;
if doasym==1; 
term=squeeze( mean(has(:,jshow+1:end,:),2)-mean(has(:,1:jshow-1,:),2) );
else; term=squeeze(has(:,jshow,:)); end
end
if vvarsshow==5% ntau
ntaus=ncdfgetvar(fileout,'ntaus'); 
term=squeeze(ntaus(:,jshow,:)); 
end
xcon=[xcon,term];  
tcon=[tcon,ts]; % concatenate
end
%xcon=xcon(:,2:end); 
tcon=tcon(2:end);
nts=length(tcon);
x=xcon; xcon=0;
ts=tcon; tcon=0;
% remove temporal mean
if removemean==1
for iii=1:nx, x(iii,:)=x(iii,:)-mean(x(iii,:)); end
end
% normalize by global stddev
if donormakw==1
xnorma=std(x(:));
x=x./xnorma;
end
% Compute power
term=zeros(nx,nts);
term=fftshift(fft2(x)); x=0;
term=abs(term).^2/(nx*nts)^2; 
kg=fftkspe(nx,dx)/(2*pi); % adim
wg=fftkspe(nts,dt*mts)/(2*pi); %adim
kg=kg/xa*40000*1000;
wg=wg/ta*oneday;
kg=-kg; % why ? this is no longer due to transpositions '(instead of .') of complex matrix ?
term=log10(term);
if dosmooth==1
term=calc_smoothn(term,[nxsmoo,ntsmoo]);
end
if undersamplext==1% undersample for figure
term=term(1:nku:end,1:nwu:end);
kg=kg(1:nku:end);
wg=wg(1:nwu:end); 
end

%%%%%%%%%%%%%%%%%%%%%
% Graph for one panel
% take the good subplot
subplot(2,2,vvarsshow);
%
if dokwlevs~=1; % do levels
contourf(kg,wg,term','LineStyle','none');
colorbar('eastoutside')
else
contourf(kg,wg,term',kwlevels, 'Linestyle', 'none'); 
caxis([kwlevels(1) kwlevels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
%
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
if 1==1
%set(gcf,'defaulttextinterpreter','latex')
if vvarsshow==1; title('u '); end
if vvarsshow==2; title('theta'); end
if vvarsshow==3; title('q '); end
if vvarsshow==4; title('Ha'); end
end
end
%
%% add lines
ithick=2;
hold on; plot([0,0],yyrange,'k-','Linewidth',ithick);
xlim(xxrange); ylim(yyrange)
%title(truc);

if sshowstab==1
% Additions (after running stability to get disp curves)
showgraphsstab=0; skelg_stabnew;% compute linear stability, without making a figure
% case M=1
if 1==1
for imod=1:4
hold on; 
if imod==1; plot(kj(ikp),wcpd(imod,ikp)','ko'); end
if imod==2; plot(kj(ikp),wcpd(imod,ikp)','ko'); end
if imod==3; plot(kj(ikp),wcpd(imod,ikp)','ko'); end
if imod==4; plot(kj(ikp),wcpd(imod,ikp)','ko'); end
end
end
%
end

% Add horizontal lines
hold on; fsingle=1/10; plot(kg,kg*0+fsingle,'k--','Linewidth',ithick);
hold on; fsingle=1/30; plot(kg,kg*0+fsingle,'k--','Linewidth',ithick);
hold on; fsingle=1/90; plot(kg,kg*0+fsingle,'k--','Linewidth',ithick);
%
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end% on variables









