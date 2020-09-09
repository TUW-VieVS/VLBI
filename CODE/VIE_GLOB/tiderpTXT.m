% ************************************************************************
%   Description:
%   This function writes estimated tidal ERP terms to a text file
%
%   Input:										
%      globsol              a structure with estimates and all relavant
%                           information about the global adjustment
%      ses_time             order of sessions according to time
%      paths                paths to directories
%
%   Output:                
%      'tiderp_*.txt' TXT file with estimates in
%                           VieVS/OUT/GLOB/_ESTIMATES
%
%
%
%   External calls: 	
%      structure tide is loaded from VieVS/OUT/GLOB             					    											
%       
%   Coded for VieVS: 
%   10 Nov 2011 by Sigrid Boehm
%   
%   Revision:
%   20 Nov 2015 by Sigrid Boehm: text file is written more flexible
%   according to the tides specified in DATA/GLOB/tidalERPvar_list.txt
%
% ************************************************************************


function tiderpTXT(globsol,ses_time,paths)

path_level='../';
load([path_level 'DATA/GLOB/tide'],'tide');
if exist([paths.path_out '_ESTIMATES/' paths.out])~=7
   mkdir([paths.path_out '_ESTIMATES/' paths.out])
end

fid=fopen([paths.path_out '_ESTIMATES/' paths.out '/tidpm_' paths.L2 '.txt'],'wt');
pro = find(tide.gmstpi<2);
np = globsol.tidnum.pro(globsol.tidnum.pro<=pro(end));
nr = globsol.tidnum.ret;
period = tide_per([tide.gmstpi tide.l tide.lp tide.F tide.D tide.OM]);

fprintf(fid,'Solution calculated: %s\n',date);
fprintf(fid,'Date of first session: %c%c%c%c%c%c%c\n',ses_time{1}(1:7));
fprintf(fid,'Date of last session:  %c%c%c%c%c%c%c\n',ses_time{end}(1:7));
fprintf(fid,'Number of sessions in the solution: %1.0f\n',size(globsol.sessions,2));
fprintf(fid,'Number of estimated periods: %2.0f\n',length(np));
fprintf(fid,'********************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'Tide       argument     period     B+   rms       A+   rms        B-   rms       A-   rms\n');
fprintf(fid,'      X  l  l'' F  D OM   [h]     microarcsec    microarcsec     microarcsec    microarcsec\n');
fprintf(fid,'********************************************************************************************** \n');
for j = 1:length(np)
   fprintf(fid,'%4s %2s %2s %2s %2s %2s %2s %7.2f %8.2f %4.2f  %8.2f %4.2f \n',...
           tide.name(np(j),1:4), num2str(tide.gmstpi(np(j))), num2str(tide.l(np(j))),...
           num2str(tide.lp(np(j))), num2str(tide.F(np(j))), num2str(tide.D(np(j))),...
           num2str(tide.OM(np(j))), period(np(j)), globsol.tidBp.val(j),...
           globsol.tidBp.sigma(j), globsol.tidAp.val(j), globsol.tidAp.sigma(j));
end  
npi = length(np);
for j = 1:length(nr)
fprintf(fid,'%4s %2s %2s %2s %2s %2s %2s %7.2f %8.2f %4.2f  %8.2f %4.2f   %8.2f %4.2f  %8.2f %4.2f\n',...
        tide.name(nr(j),1:4), num2str(tide.gmstpi(nr(j))), num2str(tide.l(nr(j))),...
           num2str(tide.lp(nr(j))), num2str(tide.F(nr(j))), num2str(tide.D(nr(j))),...
           num2str(tide.OM(nr(j))), period(nr(j)), globsol.tidBp.val(j+npi),...
           globsol.tidBp.sigma(j+npi), globsol.tidAp.val(j+npi), globsol.tidAp.sigma(j+npi),...
           globsol.tidBm.val(j),globsol.tidBm.sigma(j), globsol.tidAm.val(j), globsol.tidAm.sigma(j));   
end

fclose(fid);

fid=fopen([paths.path_out '_ESTIMATES/' paths.out '/tidut_' paths.L2 '.txt'],'wt');
np = globsol.tidnum.ut1;

fprintf(fid,'Solution calculated: %s\n',date);
fprintf(fid,'Date of first session: %c%c%c%c%c%c%c\n',ses_time{1}(1:7));
fprintf(fid,'Date of last session:  %c%c%c%c%c%c%c\n',ses_time{end}(1:7));
fprintf(fid,'Number of sessions in the solution: %1.0f\n',size(globsol.sessions,2));
fprintf(fid,'Number of estimated periods: %2.0f\n',length(np));
fprintf(fid,'********************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'Tide       argument     period     Us   rms        Uc   rms    \n');
fprintf(fid,'      X  l  l'' F  D OM   [h]       microsec        microsec    \n');
fprintf(fid,'********************************************************************\n');
for j = 1:length(np)
   fprintf(fid,'%4s %2s %2s %2s %2s %2s %2s %7.2f %8.3f %4.3f  %8.3f %4.3f \n',...
           tide.name(np(j),1:4), num2str(tide.gmstpi(np(j))), num2str(tide.l(np(j))),...
           num2str(tide.lp(np(j))), num2str(tide.F(np(j))), num2str(tide.D(np(j))),...
           num2str(tide.OM(np(j))), period(np(j)), globsol.tiduts.val(j),...
           globsol.tiduts.sigma(j), globsol.tidutc.val(j), globsol.tidutc.sigma(j));
end  
fclose(fid);

