% #########################################################################
% #     out_basEop
% #########################################################################
%
% DESCRITPION
% This function handldes writing of EOP files and baseline files.
%
% AUTHOR 
%   Andreas Hellerschmied (2015.01.22)
%
% INPUT
% - hObject
% - handles     : GUI handles structure
%
% OUTPUT
%
% COUPLING
% - eop_out
% - saveProcessList
% - bas_outtxt
% - eop_out_int
%
% CHANGES
% - 2015-08-19, A. Girdiuk: Possibility added to write three differently formatted EOP files
% - 2015-08-20, A. Hellerschmied: Function extraced from vievs2_3.m
% - 2016-09-13, A. Girdiuk: eop_out is changed for intensive sessions usage
% - 2016-09-19, A. Girdiuk: eop_out is changed to get data from chosen subfolder
%
function out_basEop(hObject, handles)

if strcmpi(get(hObject, 'Tag'), 'uipanel_plot_eopOut_write')
    outBut=1;
elseif strcmpi(get(hObject, 'Tag'), 'uipanel_plot_eopOut_writeInt')
    outBut=2;
elseif strcmpi(get(hObject, 'Tag'), 'uipanel_plot_basOut_write')
    outBut=3;
else
    fprintf('NOT DEFINED\nreturning from out_basEOP.m\n');
    return;
end

% This function writes eop data using eop_out.m
% 1. get process list
if get(handles.radiobutton_plot_eopOut_pl_current, 'Value')
    pl4eopOut='../WORK/PROCESSLIST/temporaryProcessListForEopOut.mat';
    saveProcessList(hObject, handles, pl4eopOut);
else
    allInPopupmenu=cellstr(get(handles.popupmenu_plot_eopOut_pl, 'String'));
    pl4eopOut=['../WORK/PROCESSLIST/',...
        allInPopupmenu{get(handles.popupmenu_plot_eopOut_pl, 'Value')}];
end
 
% 2. get subfolder
if get(handles.radiobutton_plot_eopOut_subfolder_current, 'Value')
    curSubfolder=get(handles.edit_run_outDirs_oneSub, 'String');
else
    allInPopupmenu=cellstr(get(handles.popupmenu_plot_eopOut_subfolder, 'String'));
    curSubfolder=allInPopupmenu{get(handles.popupmenu_plot_eopOut_subfolder, 'Value')};
end

% 2_. get process list in curSubfolder
if  get(handles.radiobutton_use_subfolder_data, 'Value')
    allX_Files=dir(['../DATA/LEVEL3/', curSubfolder, '/x_*.mat']);
    sessionnamesShort=(strrep({allX_Files.name}, 'x_', 'yyyy/'));
    sessionnamesShort=strrep(sessionnamesShort, '.mat', '');
    pl4eopOut=char(sessionnamesShort);
end

% 3. get output file
if get(handles.radiobutton_plot_eopOut_outFile_default, 'Value')
    if outBut==3
        eopOutFile=['../OUT/basRep_', curSubfolder, '.txt'];
    else
        eopOutFile=['../OUT/eop_', curSubfolder, '.txt'];
    end
else
    eopOutFile=get(handles.edit_plot_eopOut_outFile, 'String');
end

% 4. #### write files ####

% Write EOP file
if outBut==1
    output_mode(1:3)=0;
	if get(handles.checkbox_plot_eopOut_write_detailed_eop_data,'Value') 
		output_mode(1) = 1;
	end
	if get(handles.checkbox_plot_eopOut_write_sorted_eop_data,'Value')
		output_mode(2) = 1;
	end
	if get(handles.checkbox_plot_eopOut_write_vievs_eop_data,'Value')
		output_mode(3) = 1;
    end
    
	if sum(output_mode) == 0 % No output option is selected (output_mode = [0, 0, 0]) => Error-msg.!
		msgbox('ERROR: Please select at least one of the output options for the EOP file!');
    else
        eop_out(pl4eopOut, curSubfolder, eopOutFile, output_mode);
	end
	
% Write EOP file for intensives
elseif outBut==2
    eop_out(pl4eopOut,curSubfolder,eopOutFile,[],1);    % 1 (the last argument) is flag for intensive
% Write baseline file
elseif outBut==3
%     bas_outtxt(pl4eopOut,curSubfolder,eopOutFile);
    % BLR for longer span:
    limitation=str2double(get(handles.edit_plot_eopOut_basRepOptions_minBasObs,'String')); %  minimum number of observations of baselines per
    makeFig=get(handles.checkbox_plot_eopOut_basRepOptions_simplePlot,'Value');      
	printToCommand=get(handles.checkbox_plot_eopOut_basRepOptions_commandWindowOutput,'Value');
    basOutFname=[];
    if get(handles.checkbox_plot_eopOut_basRepOptions_writeBasOut, 'Value')
        basOutFname=['../OUT/', get(handles.edit_plot_eopOut_basRepOptions_writeBasOutFname, 'String')];
    end
    superstations=[];
    repeatab(curSubfolder,pl4eopOut,limitation,eopOutFile,basOutFname,...
        makeFig,printToCommand,superstations);
    
    fprintf('Baseline length repeatability written to:\n%s\n',...
        eopOutFile);
    if get(handles.checkbox_plot_eopOut_basRepOptions_writeBasOut, 'Value')
        fprintf('Baselines written to:\n%s\n',...
            basOutFname);
    end
    % NICHT MEHR NOTWENDIG: bas=baselineLengthRep(x_, antenna, bas, ip);
    
    
else
    fprintf('ERROR!!\n');
end

% delete temp process list
if get(handles.radiobutton_plot_eopOut_pl_current, 'Value')
    if exist(pl4eopOut)
        delete(pl4eopOut);
    end
end

