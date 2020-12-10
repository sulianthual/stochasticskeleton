% Skeleton model (deterministic or stochastic)
% x-y-t numerical solving and others
% by Sulian Thual
% 
% skelmain: main program for running everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters general 
%
if 0==1
clear all;
%close all;
%colormap('default')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Setup:
%(parameters are in associated fileini.m and simulation files in ./fileini folder)
% loopstart, loopend are the starting and ending index of restart files
%
setup='default'; loopstart=1; loopend=2; % default/test
%setup='HM'; loopstart=1; loopend=41; % with homogeneous
%setup='WP'; loopstart=1; loopend=41; % with Warm Pool
%setup='SC'; loopstart=1; loopend=41; % with Seasonal Cycle Warm Pool
%
% Setups with deterministic skeleton model (Majda Stechmann 2011 JAS)
%setup='HMdet'; loopstart=1; loopend=41; % with homogeneous
setup='WPdet'; loopstart=1; loopend=2;% 41; % with Warm Pool deterministic


doverylongrun=            1                 ; % key for re-running
dopostpro=                0                 ; % key for post-processing
%
fileini=strcat('skelini_',setup); % for calling ini file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs (very long, BEWARE) 
if doverylongrun==1
for indexrestart=loopstart:loopend% loop over restart files 
% Run
run(fileini); % execute ini file to get parameters
if dodefault==1
skelrun_stocha; % run fast skeleton with no keys
else
skelrun; % run skeleton with keys
end
%
end% indexrestart
end%doverylongrun==1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional post-processing (creates sequences of files in data folder)
if dopostpro==1
for indexrestart=loopstart:loopend 
strcat('skelg_vars=',num2str(indexrestart),'/',num2str(loopend))
skelg_vars; % compute files of all variables in spectral coefficients, adim
skelg_varsmjo; % compute files of all variables in spectral coefficients, adim, filtered in mjo band
skelg_varcmjo; % compute (dxu,dyu,dxv,dyv) in spectral coefficients, adim, from vars filtered in mjo band
end
ilpmin=loopstart; ilpmax=loopend; skelg_dataproj; % compute files of data projections 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs from only the ini file
if 0==1
indexrestart=1; run(fileini);
%
showgraphsstab=1; skelg_stabnew; % compute (and show) linear stability
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs FOR ONE RUN, ONE RESTART FILE
if 0==1
% Params
indexrestart=        2      ; % restart file to consider
run(fileini);
%
skelg_hovquick; % hovmuller xt, of variables (u,o,Z,Ha) and of data projection
%skelg_hovquick_lat; % hovmuller y-t (NB: v not working for M=1)
%skelg_contxy; % contour fields at a serie of timesteps (NB: v not working for M=1)
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs FOR ONE RUN, CONCATENATING RESTARTS
if 0==1
% Params
% NB: replace ilpmin and ilpmax with starting/ending restart files used to compute diagnostic
iwind=7; indexrestart=1; run(fileini);% get parameters
%
%ilpmin=loopstart; ilpmax=loopend; skelg_plotquick; % plot single column timeserie
ilpmin=loopstart; ilpmax=loopend; skelg_kwcontourconca; %k-w spectrum (only for M=1)
%ilpmin=loopstart; ilpmax=loopend; skelg_stabprojconcanew;% low-frequency modulation
%ilpmin=loopstart; ilpmax=loopend; skelg_meanstdconca; % contour xy  of time-means or time-stds
%ilpmin=21; ilpmax=39; skelg_regressvar; % regress variables on a 180E signal (doing a x-y map)

%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
