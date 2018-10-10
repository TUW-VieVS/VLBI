% This function is mainly a wrapper to call all necessary subfunctions of
% vie_sched_o. It is necessary for PARFOR looping
%
% CREATED: 21.02.17 - Matthias Schartner
%
% CHANGELOG: 

function [ process_list, stat ] = vie_sched_o_main( station, twin, source, obsmode, sched, INFILE, poutfile_main, isched, PARA )

if PARA.SAVEOUTPUT
    PARA.fid_footer = fopen([poutfile_main sprintf('footer_V%03d.txt',isched)],'w');
end


if PARA.MULTISCHED || PARA.SAVEOUTPUT
    t = clock;
    tstr1 = sprintf('%4d/%02d/%02d', t(1), t(2), t(3));
    tstr2 = sprintf('%02d:%02d:%02d', t(4), t(5), round(t(6)));
    fprintf(PARA.fid_footer,'start: %s, %s\n', tstr2, tstr1);
%     fprintf('--------------------------vie_sched output--------------------------\n');
%     fprintf('vie_sched output...\nOutput can be found at: %s\n\n',[poutfile_main 'footer.txt'])
end

if PARA.MULTISCHED 
% create folder for multi-scheduling
    fprintf(PARA.fid_footer,'*******************************************\n');
    fprintf(PARA.fid_footer,'**          schedule number: %3d         **\n',isched);
    fprintf(PARA.fid_footer,'*******************************************\n');
    poutfile_sub = sprintf('%sv%03d/',poutfile_main,isched);
    if ~isdir(poutfile_sub)
        mkdir(poutfile_sub);
    end
    
else
    
% create standard folders
    poutfile_sub = poutfile_main;
    
end
if ~isdir([poutfile_sub,'/LEVEL5'])
    mkdir([poutfile_sub,'/LEVEL5']);
end

% output directories
PARA.poutfile = poutfile_sub;
OUTFILE.ngs        = PARA.ngs;     % if write output file in NGS format;
OUTFILE.ngsname    = [poutfile_sub PARA.sessname sprintf('VA_V%03d',isched)];
OUTFILE.skd        = PARA.skd;     % if write output file in SKD and VEX format;
OUTFILE.skdname    = [poutfile_sub PARA.sessname '.skd'];
OUTFILE.vexname    = [poutfile_sub PARA.sessname '.vex'];
OUTFILE.skdsum     = PARA.skdsum;  % if write output file in SKD-sum format;
OUTFILE.skdsumname = [poutfile_sub PARA.sessname '-skdsum.txt'];

% ##### write output file #####
vie_sched_o(station, twin, source, obsmode, sched, INFILE, OUTFILE, PARA);

[meanSky,h] = checkSkyCoverage( sched,PARA,station );
stat = checkAndStat( sched,station,source,PARA,meanSky );
if PARA.MULTISCHED && ~PARA.parallel
    try
        saveStatAsXLSX(stat,poutfile_main,PARA,isched,0);
    catch
        warning('couldn''t write statistics in summary.xlsx!')
    end
end

% saves output in LEVEL5 directories
save([PARA.pfolder 'schedparam.mat'], 'PARA');
save([poutfile_sub,'LEVEL5/','schedparam.mat'], 'PARA'); 
save([PARA.pfolder 'station.mat'], 'station');
save([poutfile_sub,'LEVEL5/station.mat'], 'station'); 
save([PARA.pfolder 'sched.mat'], 'sched');
save([poutfile_sub,'LEVEL5/sched.mat'], 'sched'); 
save([PARA.pfolder 'source.mat'], 'source');
save([poutfile_sub,'LEVEL5/source.mat'], 'source'); 
save([PARA.pfolder 'obsmode.mat'], 'obsmode');
save([poutfile_sub,'LEVEL5/obsmode.mat'], 'obsmode'); 
if (PARA.OPTIMIZATION == 1)
    save([PARA.pfolder 'srcat.mat'],   'srcat');
    save([poutfile_sub,'LEVEL5/srcat.mat'], 'srcat'); 
    save([PARA.pfolder 'catpair.mat'], 'catpair');
    save([poutfile_sub,'LEVEL5/catpair.mat'], 'catpair'); 
end

% save path to ngs file in process_list
process_list = OUTFILE.ngsname;

h.Name = sprintf('Sky Coverage V%03d',isched);
print(gcf,sprintf('%sSky Coverage V%03d',PARA.poutfile,isched),'-dpng','-r150')

if PARA.openAnalyser
    sched_analyser();
end

if PARA.MULTISCHED || PARA.SAVEOUTPUT
    fclose(PARA.fid_footer);
end
if PARA.MULTISCHED
    try
        movefile([poutfile_main sprintf('footer_V%03d.txt',isched)],poutfile_sub)
    catch
        warning('couldn''t move file %s in directory %s',sprintf('footer_V%03d.txt',isched),poutfile_sub)
    end
    try
        movefile([poutfile_main sprintf('body_V%03d.txt',isched)],poutfile_sub)
    catch
        warning('couldn''t move file %s in directory %s',sprintf('body_V%03d.txt',isched),poutfile_sub)
    end
end


end

