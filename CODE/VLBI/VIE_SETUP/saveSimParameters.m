% #########################################################################
% #     saveSimParameters
% #########################################################################
%
% This function saves the GUI parameters for VIE_SIM
%
% AUTHOR 
%   moved in separate function (from vie_setup.m) by A. Hellerschmied
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
% 2016-10-19, A. Hellerschmied: - Added options to define a seed for the random number generator
%                               - Added edit field to specify white noise for satellite observations
% 2016-10-25, A. Hellerschmied: Added checkbox "checkbox_vie_sim_write_vso_file"
%
%

function saveSimParameters(hObject, handles)
% This function saves SIM parameters according to GUI

% what shall happen if the user presses the 'OK' button?
% which parameters does the user want to simulate? -> store
simparam.sswd           = get(handles.checkbox_vie_sim_simSwd,'Value');
simparam.sclk           = get(handles.checkbox_vie_sim_simClock,'Value');
simparam.swn            = get(handles.checkbox_vie_sim_simNoise,'Value');
simparam.ref0           = get(handles.checkbox_vie_sim_setRefClk0,'Value');
simparam.ngs            = get(handles.checkbox_vie_sim_writeNgsFile,'Value');
simparam.write_vso      = get(handles.checkbox_vie_sim_write_vso_file,'Value');
simparam.ssou           = get(handles.checkbox_vie_sim_simSS,'Value');
simparam.rng_use_seed   = get(handles.checkbox_vie_sim_rng_use_seed,'Value');

err = 0;

% store all other parameters - if an entry that is needed is missing ->
% don't continue
if get(handles.radiobutton_vie_sim_fromParameterFile,'Value') == 1
    list_entries = get(handles.listbox_vie_sim_paramFile,'String');
    index_selected = get(handles.listbox_vie_sim_paramFile,'Value');
    simparam.turbfile = list_entries{index_selected};
    if isreal(str2num(get(handles.edit_vie_sim_nDaysSim,'String'))) & str2num(get(handles.edit_vie_sim_nDaysSim,'String')) > 0
        simparam.idays = str2num(get(handles.edit_vie_sim_nDaysSim,'String'));
    else
        err = err + 1;
    end
    if isreal(str2num(get(handles.edit_vie_sim_startInd,'String'))) & str2num(get(handles.edit_vie_sim_startInd,'String')) > 0
        simparam.sind = str2num(get(handles.edit_vie_sim_startInd,'String'));
    else
        err = err + 1;
    end
%     if err == 0
%        save('../DATA/LEVEL4/simparam.mat','simparam');
%     end
elseif get(handles.radiobutton_vie_sim_specifyNow,'Value') == get(handles.radiobutton_vie_sim_specifyNow,'Max')
    if isreal(str2num(get(handles.edit_vie_sim_nDaysSim,'String'))) & str2num(get(handles.edit_vie_sim_nDaysSim,'String')) > 0
        simparam.idays = str2num(get(handles.edit_vie_sim_nDaysSim,'String'));
    else
        % + A. Pany, 29 Aug 2011
        if ~((simparam.sswd == 0) && (simparam.sclk == 0) && (simparam.swn == 0))
            err = err + 1;
        end
        % - A. Pany
    end
    if get(handles.checkbox_vie_sim_simSwd,'Value') == get(handles.checkbox_vie_sim_simSwd,'Max')
        if isreal(str2num(get(handles.edit_vie_sim_tropo_cn,'String'))) & str2num(get(handles.edit_vie_sim_tropo_cn,'String')) > 0
            simparam.Cn    = str2num(get(handles.edit_vie_sim_tropo_cn,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_H,'String'))) & str2num(get(handles.edit_vie_sim_tropo_H,'String')) > 0
            simparam.H     = str2num(get(handles.edit_vie_sim_tropo_H,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_wzd0,'String'))) & str2num(get(handles.edit_vie_sim_tropo_wzd0,'String')) > 0
            simparam.wzd0  = str2num(get(handles.edit_vie_sim_tropo_wzd0,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_vn,'String')))
            simparam.vn    = str2num(get(handles.edit_vie_sim_tropo_vn,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_ve,'String')))
            simparam.ve    = str2num(get(handles.edit_vie_sim_tropo_ve,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_dhseg,'String'))) & str2num(get(handles.edit_vie_sim_tropo_dhseg,'String')) > 0
        simparam.dhseg = str2num(get(handles.edit_vie_sim_tropo_dhseg,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_tropo_dh,'String'))) & str2num(get(handles.edit_vie_sim_tropo_dh,'String')) > 0
            simparam.dh    = str2num(get(handles.edit_vie_sim_tropo_dh,'String'));
        else
            err = err + 1;
        end
    end
    if get(handles.checkbox_vie_sim_simClock,'Value') == get(handles.checkbox_vie_sim_simClock,'Max')
        if isreal(str2num(get(handles.edit_vie_sim_clock_asd,'String'))) & str2num(get(handles.edit_vie_sim_clock_asd,'String')) > 0
            simparam.sy1   = str2num(get(handles.edit_vie_sim_clock_asd,'String'));
        else
            err = err + 1;
        end
        if isreal(str2num(get(handles.edit_vie_sim_clock_asdInt,'String'))) & str2num(get(handles.edit_vie_sim_clock_asdInt,'String')) > 0
            simparam.sy2   = str2num(get(handles.edit_vie_sim_clock_asdInt,'String'));
        else
            err = err + 1;
        end
    end
    if get(handles.checkbox_vie_sim_simNoise,'Value') == get(handles.checkbox_vie_sim_simNoise,'Max')
        % White noise for quasar observations:
        if isreal(str2num(get(handles.edit_vie_sim_noise,'String'))) & str2num(get(handles.edit_vie_sim_noise,'String')) > 0
            simparam.wn    = str2num(get(handles.edit_vie_sim_noise,'String'));
        else
            err = err + 1;
        end
        % White noise for satellite observations:
        [wn_sat, status] = str2num(get(handles.edit_vie_sim_noise_sat,'String'));
        % Check input wn_sat:
        if status
           if wn_sat > 0
               simparam.wn_sat = wn_sat;
           else
               err = err + 1;
           end
        else
            err = err + 1;
        end
    end
    if isreal(str2num(get(handles.edit_vie_sim_startInd,'String'))) & str2num(get(handles.edit_vie_sim_startInd,'String')) > 0
        simparam.sind = str2num(get(handles.edit_vie_sim_startInd,'String'));
    else
        err = err + 1;
    end
    
    simparam.turbfile = [];
%     if err == 0
%         save('../DATA/LEVEL4/simparam.mat','simparam');
%     end
else
    % + A. Pany, 29 Aug 2011
    simparam.turbfile = [];
    if isreal(str2num(get(handles.edit_vie_sim_startInd,'String'))) & str2num(get(handles.edit_vie_sim_startInd,'String')) > 0
        simparam.sind = str2num(get(handles.edit_vie_sim_startInd,'String'));
    else
        err = err + 1;
    end
    if isreal(str2num(get(handles.edit_vie_sim_nDaysSim,'String'))) & str2num(get(handles.edit_vie_sim_nDaysSim,'String')) > 0
        simparam.idays = str2num(get(handles.edit_vie_sim_nDaysSim,'String'));
    else
        err = err + 1;
    end
%     if err == 0
%         save('../DATA/LEVEL4/simparam.mat','simparam');
%     % - A. Pany
%     end
end

% Source structure:
if get(handles.checkbox_vie_sim_simSS,'Value')==1
    list_entries = get(handles.listbox_vie_sim_ss,'String');
    index_selected = get(handles.listbox_vie_sim_ss,'Value');
    simparam.sou_cat = list_entries{index_selected};
    
    save('../DATA/LEVEL4/simparam.mat','simparam');
end

% Random number generator settings:
if simparam.rng_use_seed
    % Check input:
    [wn_sat, status] = str2num(get(handles.edit_vie_sim_rng_seed,'String'));
    if status
        if (wn_sat <= 0) || (mod(wn_sat, 1) ~= 0) 
            error('The seed for the random number generator defined in the VIE_SIM GUI has to be a nonnegative integer! Current seed: %f\n', wn_sat);
        else
            simparam.rng_seed = wn_sat;
        end
    else
        error('The seed for the random number generator defined in the VIE_SIM GUI has to be a nonnegative integer!\n');
    end

end

% Save sim. parameters, if no errors occured:
if err == 0
    save('../DATA/LEVEL4/simparam.mat','simparam');
else
    error('Error while obtaining simulation parameters from the VIE_SIM GUI! ../DATA/LEVEL4/simparam.mat was not saved!')
end





