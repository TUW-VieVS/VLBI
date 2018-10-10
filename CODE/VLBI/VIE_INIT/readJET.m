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
%   11 sep 2013 by Lucia Plank: read jet
% ************************************************************************
function [scan_excl_field2]=readJET(outfile,jamax)

fid=fopen(outfile,'r');
  a=1; ka=1;
  while ~feof(fid)
    str=fgetl(fid);
    if length(str)>19
        scan_all(ka).sta1=str(1:8);
        scan_all(ka).sta2=str(10:17);
        scan_all(ka).mjd=str2num(str(19:36));
        ka=ka+1;
        if abs(str2num(str(38:end)))>jamax
        scan_excl_field(a).sta1=str(1:8);
        scan_excl_field(a).sta2=str(10:17);
        scan_excl_field(a).mjd=str2num(str(19:36));
        scan_excl_field(a).jetang=str2num(str(38:end));
        a=a+1;
        end
    end
  end
  
 
'zufall'

for i=1:(a-1)
     zz=ceil(rand*(ka-i));
     scan_excl_field2(i)=scan_all(zz);
     scan_all(zz)=[];
end 

  
  fclose(fid);