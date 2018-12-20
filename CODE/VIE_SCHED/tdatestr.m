% Purpose  
%   Convert [yy,mm,dd,hh,min,secf] to string format.
%   varargin = flag to switch to old output format 
%
% History  
%   2010-04-06   Jing SUN   Created
%	2016-07-11	M. Schartner: nicer output
%   2016-09-06 Matthias Schartner: 	now has 2 possible output formats (new one is default)
%		old: yyyy mm dd hh mm ss.ssssssssss
%		new: yyyy-mm-dd hh:mm:ss



function [datestr] = tdatestr(yy, mm, dd, hh, min, secf, varargin)

dsec = round(secf);
min(dsec==60)=min(dsec==60)+1;
dsec(dsec==60)=0;
hh(min==60)=hh(min==60)+1;
min(min==60)=0;

if nargin == 7 
    [datestr] = sprintf('%4d %02d %02d %02d %02d %14.10f', yy, mm, dd, hh, min, dsec);
else
    [datestr] = sprintf('%4d-%02d-%02d %02d:%02d:%02d', yy, mm, dd, hh, min, dsec);
end
