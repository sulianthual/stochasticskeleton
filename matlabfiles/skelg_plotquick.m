% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% do a quick plot of some timeseries 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ilppsmin=1;ilppsmax=21;% restart files to use, this is set in the skelmain
%
xshow=20; yshow=0; % where to show (in x 1000km and in adim y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(iwind); iwind=iwind+1;, clf;
icount=1;

for vvarsshow=1:4
% some preparation
% things to concatenate over restarts
xcon=[0,0]; tcon=[0,0];
for ilpps=ilpmin:ilpmax      ;% loop index
{'il=',num2str(ilpps),'/',num2str(ilpmax)}
indexrestart=ilpps; run(fileini);
%
[ishow,term]=searchclosest(xg,xshow);[jshow,term]=searchclosest(yg,yshow);
ts=ncdfgetvar(fileout,'ts');
if vvarsshow==1% u
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); us=zeros(nx,nyk,nts);
for kts=1:nts; us(:,:,kts)=singlecolumnu(Ks(:,kts),Rms(:,:,kts),psim); end
term=ua*squeeze(us(ishow,jshow,:));% u
end
if vvarsshow==2% o
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
term=oa*squeeze(os(ishow,jshow,:)); % o
Zs=ncdfgetvar(fileout,'Zs'); term=oa*squeeze(Zs(ishow,jshow,:)-QQ*os(ishow,jshow,:)); % q
end
if vvarsshow==3% q
Ks=ncdfgetvar(fileout,'Ks'); Rms=ncdfgetvar(fileout,'Rms'); os=zeros(nx,nyk,nts);
for kts=1:nts; os(:,:,kts)=singlecolumno(Ks(:,kts),Rms(:,:,kts),psim); end
Zs=ncdfgetvar(fileout,'Zs'); term=oa*squeeze(Zs(ishow,jshow,:)-QQ*os(ishow,jshow,:)); % q
end
if vvarsshow==4% Ha
etas=ncdfgetvar(fileout,'etas'); 
term=HH*DDA*(oa/(ta/oneday))*squeeze(etas(ishow,jshow,:)); 
end
xcon=[xcon,term']; tcon=[tcon,ts]; % concatenate
end% end of loop on restarts
xcon=xcon(3:end); tcon=tcon(3:end);
%
subplot(2,2,icount); icount=icount+1;
plot(tcon/1000,xcon)
xlabel('time (1000 days)');
xlim([min(tcon/1000),max(tcon/1000)])
if vvarsshow==1; title('u'); end
if vvarsshow==2; title('theta'); end
if vvarsshow==3; title('q'); end
if vvarsshow==4; title('Ha'); end
%
end% loop on vvarshow
