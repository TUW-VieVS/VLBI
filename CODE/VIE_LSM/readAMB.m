% ************************************************************************
%   Description:
%   This function reads the AMB file
% ************************************************************************
function [obsamb]=readAMB(outfile)
obsamb = [];

fid = fopen(outfile,'r');
  a = 1;
  while ~feof(fid)
    str = fgetl(fid);
%     splstr= split(str);
    

        obsamb(a).sta1 = str(1:8);
        obsamb(a).sta2 = str(10:17);
        obsamb(a).mjd = str2num(str(19:36));
        obsamb(a).sou = str(38:45);
        obsamb(a).amb = str2num(str(47:end));
        a=a+1;

  end
  fclose(fid);