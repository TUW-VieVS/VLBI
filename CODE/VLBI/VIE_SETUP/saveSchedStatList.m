% This function writes the selected stations in scheduling windows to a
% file in VieVS/WORK/STATIONLIST.
% 
% 9.4.2014, Matthias Madzak
%

function saveSchedStatList(hObject,eventdata,handles)

% get chosen stations
stations=get(handles.listbox_vie_sched_stanet, 'String');

if isempty(stations)
    msgbox('Please select a station network first','No station selected','warn');
    return;
end

outPath='../WORK/STATIONLIST/';
if ~exist(outPath, 'dir')
    mkdir(outPath);
end

prompt = {'Input stationlist name:'};
dlg_title = 'Name of station list';
num_lines = 1;
def = {''};
answer = inputdlg(prompt,dlg_title,num_lines,def);

while strcmpi(answer,'')
    prompt = {'Input stationlist name:'};
    dlg_title = 'Name of station list';
    num_lines = 1;
    def = {''};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
end

if ~isempty(answer)
    save([outPath, answer{1}, '.mat'], 'stations');

    msgbox(['Stations saved to WORK/STATIONLIST/', answer{1}, '.mat'], 'Done');
end