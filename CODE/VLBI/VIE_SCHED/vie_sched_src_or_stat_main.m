% This function is mainly a wrapper to call all necessary subfunctions of
% vie_sched with source or station based optimization.
% It is necessary for PARFOR looping
%
% CREATED: 21.02.17 - Matthias Schartner
%
% CHANGELOG: 


function [ station, twin, source, obsmode, sched ] = vie_sched_src_or_stat_main( station, twin, source, obsmode, srcat, catpair, PARA, poutfile_main, isched,INFILE_multiSched )
%% preparations
if (PARA.MULTISCHED || PARA.SAVEOUTPUT )
    PARA.fid_body = fopen([poutfile_main sprintf('body_V%03d.txt',isched)],'w');
end

% change parameters for multi-Schedif PARA.MULTISCHED
if PARA.MULTISCHED    
    [ station, twin, source, PARA ] = multiSchedReader( station, twin, source, PARA, INFILE_multiSched, isched );
end

iOpt = 0;
while 1
    %% calculate the schedule
    % remove special station for tag-along mode 
    [station, station_tagalong] = tagalong_splitStations( station, PARA.tagalongname);

    % ##### calculate the schedule #####
    % Select scheduling approach:
    sched = vie_sched_src(station, twin, source, obsmode, srcat, catpair, PARA);                 

    % trim the schedule to selected timespan
    sched = sched_trim( sched, PARA );

    % add special station with tag-along mode 
    if ~isempty(station_tagalong)
        fprintf('scheduling the special station with tag-along mode ...\n');
        sched = vie_sched_tagalong(station, station_tagalong, source, obsmode, sched, PARA);
        station = tagalong_mergeStations( station, station_tagalong );
    end

    % add strong antenna in STAR mode
    if PARA.STARMODE==1
        fprintf('*** fill-in scans for STAR mode ***\n')
        sched = vie_sched_tagalong_star(station, source, obsmode, sched, PARA);
        disp_stat_starmode(source,sched)
    end
    
    %% everything below here is just for optimization conditions
    % check if you have a certain condition
    if isempty(PARA.optimization_condition)
        break
    else
        [ flag_break, sched, source ] = optimizer_main( sched, station, twin, source, obsmode, PARA, iOpt, isched );
        if flag_break
            break
        end
    end
    
end

if PARA.MULTISCHED || PARA.SAVEOUTPUT
    fprintf('finish schedule %d\n',isched);
    fclose(PARA.fid_body);
end

end

