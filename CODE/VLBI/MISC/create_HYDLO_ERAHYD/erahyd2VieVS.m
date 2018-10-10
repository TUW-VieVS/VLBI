% routine to write hydrology data 
% from Strasburg_EOST_Loading_service in mat-structure
%
% It needs to be run without input 
%  and will generate mat-structure from data in current folder
%
%   Coded for VieVS: 
% 09 Aug 2016 by A. Girdiuk
%
%	Revision:
%

list_level=dir(fullfile(['*.erahyd']));
list_stations = {list_level.name}';
ind=0;
load('../../TRF/superstation');
for i = 1:length(list_stations)
    name = list_stations{i,1};
    name1 = name(1:4); name2 = name(6:14);
    for j = 1:length(superstations)   
    erahyd_txt = []; 
        if strcmp(name1,superstations(j).CDP) && strcmp(name2,superstations(j).domes)
            erahyd_txt = load([superstations(j).CDP '_' superstations(j).domes '_NEU.erahyd'],'r');
            ind=ind+1;
            erahyd(ind).ivsname = superstations(j).name;
            erahyd(ind).mjd = erahyd_txt;
        end
   end
   name = []; name1 = []; name2 = [];
end

save erahyd erahyd
