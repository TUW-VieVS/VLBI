% Description:
%   SCHED2010: Next Generation Scheduling for VLBI2010.
% Coded for VieVS:
%   06 Apr 2010 by Jing Sun
% Revision:
%   06 Jan 2014 by Jing Sun:  develop tag-along mode for special station
%                             (/CATALOGS/tagalong.txt)
%   24 Mar 2014 by Jing Sun:  add PARA.FORSI option in param.txt to check
%                             for source structure study
%   24 Mar 2014 by Jing Sun:  use twin.txt to specify the twin/multiple stations
%                             (/CATALOGS/twin.txt)
%   13 Aug 2014 by David Mayer:  use sweight.txt to specify a weight for sources (default is 1)
%                             (/CATALOGS/sweight.txt)
%   Sept 2014 by David Mayer: de-bugging
%   05 Feb 2015 by A. Hellerschmied: Added option to have a sub-directory for
%     the outout (defined via the GUI)
%   04 Mar 2015 by A. Hellerschmied: Observation mode name now defined in
%     param.txt (not hard-coded any more).
%   15 Jul 2015 by A. Hellerschmied: Modifications to add satellite scheduling functions
%   25.May 2016 by M. Schartner: added opportunity for STAR mode
%   15 Jun 2016 by A. Hellerschmied: Minor changes for satellite scheduling
%   11 Jul 2016 by M. Schartner: Changes for station or source based scheduling
%   05 Sep 2016 by M. Schartner: changes for sched_manual
%   05 Okt 2016 by M. Schartner: now saves LEVEL5 in subdir
%   04 Nov 2016 by M. Schartner: now calls checkSkyCoverage.m and checkAndStat.m
%   14 Nov 2016 by M. Schartner: changes for multi scheduling
%   29 Nov 2016 by M. Schartner: changes for deeper integration in vie_batch3.0, now returns process_list
%   30 Nov 2016 by A. Hellerschmied: vie_sched_sat now returns a process_list
%   22 Dec 2016 by A. Hellerschmied: Minor changes for log files (fid_body)
%   21 Feb 2017 by M. Schartner: rearranged all parts for better readability, better multicore calculation sequence, saves statistics for multisched in excel sheet
%   29 Mai 2017 by M. Schartner: bugfix in vie_sched_manual

function [process_list] = vie_sched()
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')


%% ----------------------------- PARAMETERS -----------------------------%
pvievs  = '../';   %%%
pfolder = [pvievs 'DATA/LEVEL5/'];
load([pfolder 'schedparam.mat'], 'PARA');
load([pfolder 'stanetstr.mat'], 'stanetstr');
PARA.pvievs  = pvievs;
PARA.pfolder = pfolder;

% session time
[PARA.startmjd] = tmjd(PARA.year_s, PARA.month_s, PARA.day_s, PARA.hour_s, PARA.minute_s, PARA.second_s);
PARA.endmjd = PARA.startmjd + PARA.duration/24;
[PARA.sessname] = tymdstr(PARA.year_s, PARA.month_s, PARA.day_s);

% input files
pinfile = PARA.pinfile;
INFILE.source   = [pinfile 'source.cat'];     % source positions
INFILE.flux     = [pinfile 'flux.cat'];       % source fluxes
INFILE.antenna  = [pinfile 'antenna.cat'];    % antenna information
INFILE.position = [pinfile 'position.cat'];   % station x,y,z locations
INFILE.equip    = [pinfile 'equip.cat'];      % equipment IDs
INFILE.mask     = [pinfile 'mask.cat'];       % horizon and coordinate masks
INFILE.modes    = [pinfile 'modes.cat'];      % observing modes

INFILE.param    = [pinfile 'param.txt'];      % minor scheduling parameters
INFILE.snrmin   = [pinfile 'snrmin.txt'];     % minimum SNR
INFILE.down     = [pinfile 'down.txt'];       % downtime information
INFILE.psource  = [pinfile 'psource.txt'];    % particular sources
INFILE.tagalong = [pinfile 'tagalong.txt'];   % station for tag-along mode
INFILE.twin     = [pinfile 'twin.txt'];       % twin/multiple station
INFILE.sweight  = [pinfile 'sweight.txt'];    % source weighting, David Mayer, 2014 Aug 13

INFILE.freq     = [pinfile 'freq.cat'];       % frequency sequences
INFILE.rx       = [pinfile 'rx.cat'];         % receiver setups
INFILE.loif     = [pinfile 'loif.cat'];       % station Lo and IF setups
INFILE.rec      = [pinfile 'rec.cat'];        % recording modes
INFILE.hdpos    = [pinfile 'hdpos.cat'];      % head offsets
INFILE.tracks   = [pinfile 'tracks.cat'];     % standard recorded tracks
INFILE.sourcestar = [pinfile 'source_star.cat']; % sources for star mode
INFILE.acceleration = [pinfile 'acceleration.cat']; % catalog for acceleration
INFILE.multiSched = [pinfile 'multiSched.txt']; % multi scheduling parameter file


% create main output folder
if PARA.USE_OUTPUT_SUBDIR
    sub_directroy_main = PARA.OUTPUT_SUBDIR;             % Sub-directory defined in GUI
else
    sub_directroy_main = num2str(round(PARA.year_s));    % Sub-directory = Year
end
if ~isdir([[pvievs 'DATA/SCHED/'],sub_directroy_main])
    mkdir([[pvievs 'DATA/SCHED/'],sub_directroy_main]);
end
poutfile_main = [[pvievs 'DATA/SCHED/'] sub_directroy_main '/'];

% creating files
PARA.fid_header = 1;
PARA.fid_footer = 1;
PARA.fid_body   = 1;
if PARA.MULTISCHED || PARA.SAVEOUTPUT
    PARA.fid_header = fopen([poutfile_main 'header.txt'],'w');
end

% constants
PARA.MAX_FLUXPARA = 30;
PARA.MAX_SEFDPARA = 4;    % (sefd, y, c0, c1)
PARA.MAX_HOR      = 80;   % a maximum of 40 pairs of numbers for horizon masks


% VIE_SCHED start
t = clock;
tstr1 = sprintf('%4d/%02d/%02d', t(1), t(2), t(3));
tstr2 = sprintf('%02d:%02d:%02d', t(4), t(5), round(t(6)));
fprintf('VIE_SCHED starts at %s, %s\n', tstr2, tstr1);

if PARA.MULTISCHED
    fprintf(PARA.fid_header,'VIE_SCHED starts at %s, %s\n', tstr2, tstr1);
    copyfile(INFILE.multiSched,poutfile_main)
end


%% -------------------------- INPUT ---------------------------%
fprintf('---------------------------vie_sched input---------------------------\n');
if PARA.MULTISCHED || PARA.SAVEOUTPUT
    fprintf(PARA.fid_header,'---------------------------vie_sched input---------------------------\n');
    fprintf('busy...\nOutput can be found at: %s\n\n',[poutfile_main 'header.txt'])
end

% read input information and do initial calculation
[station, twin, tagalongname, source, obsmode, srcat, catpair, PARA] = vie_sched_i(stanetstr, INFILE, PARA);
PARA.tagalongname = tagalongname;

% check STAR mode
if PARA.STARMODE == 1
    PARA = checkStarSetup( station, PARA );
end

%% ----------------------- prepare multi Scheduling -------------------%
if PARA.MULTISCHED == 1
    [multisched_start_idx, multisched_end_idx] = multiSchedCounter( INFILE.multiSched, PARA );
    nscheds = multisched_end_idx-multisched_start_idx+1;
else
    multisched_start_idx = 1;
    multisched_end_idx = 1;
    nscheds = 1;
end

%% ----------------------- print optimization condtion -----------------
if ~isempty(PARA.optimization_condition)
    fprintf(PARA.fid_header,'Optimization condition was found:\n');
    fprintf(PARA.fid_header,'    %s\n',PARA.optimization_condition);
end

fprintf(PARA.fid_header,'\n');
if PARA.MULTISCHED || PARA.SAVEOUTPUT
    fclose(PARA.fid_header);
end

%% --------------------------- schedule ----------------------------%
fprintf('--------------------------vie_sched schedule--------------------------\n');
if PARA.MULTISCHED || PARA.SAVEOUTPUT
    fprintf('busy... this might take a while!!!\nLog can be found at: %s\n\n',[poutfile_main 'body*.txt',])
end

if PARA.OPTIMIZATION == 1 || PARA.OPTIMIZATION == 2
%% station or source based scheduling

    if PARA.MULTISCHED
        try
            saveStatAsXLSX([],poutfile_main,PARA,0,0);
        catch
            warning('couldn''t write statistics in summary.xlsx!')
        end
    end

    if PARA.parallel == 0
        % no parallel processing
        for isched = multisched_start_idx:multisched_end_idx
            fprintf('starting schedule number %d (%d total)\n',isched,nscheds);
            [station_, twin_, source_, obsmode_, sched] = vie_sched_src_or_stat_main( station, twin, source, obsmode, srcat, catpair, PARA, poutfile_main, isched, INFILE.multiSched );
            process_list{isched}= vie_sched_o_main( station_, twin_, source_, obsmode_, sched, INFILE, poutfile_main, isched, PARA );
        end

    else
        % parallel processing
        stat = struct('pSlewTime',[],'pScanTime',[],'pConstTime',[],'pIdleTime',[],'totalNumberScans',[],'scansPerH',[],'meanSky',[],'nobs',[],'nscan',[]);
        parfor isched = multisched_start_idx:multisched_end_idx
            fprintf('starting schedule number %d (%d total)\n',isched,nscheds);
            [station_, twin_, source_, obsmode_, sched] = vie_sched_src_or_stat_main( station, twin, source, obsmode, srcat, catpair, PARA, poutfile_main, isched, INFILE.multiSched);
            [process_list{isched},stat(isched)] = vie_sched_o_main( station_, twin_, source_, obsmode_, sched, INFILE, poutfile_main, isched, PARA );
        end
        for isched = multisched_start_idx:multisched_end_idx
            try
                saveStatAsXLSX(stat(isched),poutfile_main,PARA,isched,0);
            catch
                warning('couldn''t write statistics in summary.xlsx!')
            end
        end
    end


    process_list = char(process_list{multisched_start_idx:multisched_end_idx});
elseif PARA.OPTIMIZATION == 3
%% satellite scheduling

    if PARA.SAVEOUTPUT
        PARA.fid_body = fopen([poutfile_main 'body.txt'],'w');
    end
    [sched_data, process_list] = vie_sched_sat(station, source, PARA, INFILE, obsmode);
    if PARA.SAVEOUTPUT
        fclose(PARA.fid_body);
    end

    save([PARA.pfolder 'sched_data.mat'], 'sched_data');

elseif PARA.OPTIMIZATION == 4
%% manual scheduling
    h = sched_manual(station, twin, source, obsmode, PARA);
    uiwait(h);
    try
        load('..\DATA\LEVEL5\sched_manually.mat')
        delete('..\DATA\LEVEL5\sched_manually.mat')
    catch
        error('you have closed sched manually bevore it was finished!')
    end
    
    [ sched ] = sched_trim( sched, PARA );

    process_list = vie_sched_o_main( station, twin, source, obsmode, sched, INFILE, poutfile_main, 1, PARA );

end



fclose('all');
fprintf('--------------------------------------------------------------------\n');
fprintf('VIE_SCHED ends at %s, %s\n', tstr2, tstr1);
