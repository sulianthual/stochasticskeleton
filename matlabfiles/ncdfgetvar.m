function [fout]=ncdfgetvar(filenc,varnc)
% by Sulian Thual
fout=ncread(filenc,varnc);
