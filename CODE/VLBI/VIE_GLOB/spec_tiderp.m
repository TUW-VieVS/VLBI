% ************************************************************************
%   Description:
%   This function reads which tidal terms are to be estimated in ERP
%   from a text file: /DATA/GLOB/tidalERPvar_list.txt
%
%
%   Input:										
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2),
%
%   Output:                
%      parGS               number of tidal terms which should be estimated
%                          are stored in spectid, spectidret
%
%   External calls: 	
%      globind               					    											
%       
%   Coded for VieVS: 
%   13 Nov 2015 by Sigrid Boehm
%
%   Revision: 
% ************************************************************************

function parGS = spec_tiderp(parGS,g,path_level)

load([path_level 'DATA/GLOB/tide'],'tide');


% initialize special tide parameters 'spectid,-ret' for tidal ERP terms
    ret = find(tide.gmstpi>=2); 
    parGS(g.g_tidpm).spectid = (1:size(tide.num))';
    parGS(g.g_tidpm).spectidret = (1:length(ret))';
    parGS(g.g_tidut).spectid = (1:size(tide.num))';
    
     fileID = fopen([path_level 'DATA/GLOB/tidalERPvar_globest.txt']);
     C = textscan(fileID,'%s %d %d %d %d %d %d %d','Delimiter',',','CommentStyle','%');
     fclose(fileID);
     
     t.name = cell2mat(C{1,1});
     t.num = [C{1,8}];
     t.gmstpi = [C{1,2}];
     t.l = [C{1,3}];
     t.lp = [C{1,4}];
     t.F = [C{1,5}];
     t.D = [C{1,6}];
     t.OM = [C{1,7}];
     
     spectid = t.num;
     spectidret=[];
     retest = find(spectid>=ret(1));
     if ~isempty(retest)
     spectidret = spectid(retest) - (ret(1)-1);
     end
     
parGS(g.g_tidpm).spectid = spectid;
parGS(g.g_tidpm).spectidret = spectidret;
parGS(g.g_tidut).spectid = spectid;
