% ************************************************************************
%   Description:
%   This function reads the OUTLIER file and writes information about those
%   to ini_opt.scan_excl.
%
%   Input:	
%      OOUTLIER filename.
%
%   Output:
%      ini_opt.scan_excl.
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   Jul 2012 by Matthias Madzak
%
%   Revision: 
%   dd mmm yyyy by FIRSTNAME SECONDNAME:
% ************************************************************************
function [scan_excl_field]=readOUT(outfile)
scan_excl_field = [];

fid = fopen(outfile,'r');
  a = 1;
  while ~feof(fid)
    str = fgetl(fid);
    if length(str)>19
        scan_excl_field(a).sta1 = str(1:8);
        scan_excl_field(a).sta2 = str(10:17);
        scan_excl_field(a).mjd = str2num(str(19:end));
        a=a+1;
    end
  end
  fclose(fid);