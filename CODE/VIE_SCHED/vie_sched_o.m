% Purpose  
%   Write output file.
% History  
%   2010-04-06   Jing SUN   Created

function vie_sched_o(station, twin, source, obsmode, sched, INFILE, OUTFILE, PARA)

% write /VieVS/TRF/trf_sched.txt file
fprintf(PARA.fid_footer,'1. write /VieVS/TRF/trf_sched.txt\n');
otrf(station, PARA);

% write /VieVS/CRF/crf_sched.txt file
fprintf(PARA.fid_footer,'2. write /VieVS/CRF/crf_sched.txt\n');
ocrf(source, sched, INFILE, PARA);

filenum = 2;

% write output file in NGS format
if (OUTFILE.ngs == 1)
    filenum = filenum + 1;
    fprintf(PARA.fid_footer,'%1d. write output file in NGS format: %s\n', filenum, OUTFILE.ngsname);
    ongs(source, station, sched, OUTFILE.ngsname, twin);
end

% write output file in SKD-sum format
if (OUTFILE.skdsum == 1)
    filenum = filenum + 1;
    fprintf(PARA.fid_footer,'%1d. write output file in SKD-sum format: %s\n', filenum, OUTFILE.skdsumname);
    oskdsum(source, station, obsmode, sched, OUTFILE.skdsumname, PARA);
end

% write output file in SKD format
if (OUTFILE.skd == 1)
    filenum = filenum + 1;
    fprintf(PARA.fid_footer,'%1d. write output file in SKD format: %s\n', filenum, OUTFILE.skdname);
    oskd(source, station, obsmode, sched, INFILE, OUTFILE.skdname, PARA);
end

% write output file in VEX format
% if (OUTFILE.skd == 1)
%     filenum = filenum + 1;
%     fprintf('%1d. write output file in VEX format: %s\n', filenum, OUTFILE.vexname);
%     oskd2vex(OUTFILE.skdname, OUTFILE.vexname);
% end

