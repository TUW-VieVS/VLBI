% saves statistics in excel file
%
% CREATED: Matthias Schartner
%
% Changelog: 2017-05-29: M. Schartner: Bugfix 

function [  ] = saveStatAsXLSX( stat,poutfile,PARA,line,level3Flag )

fname = [poutfile 'summary.xlsx' ];
sheet = 1;

if line == 0 && level3Flag == 0
    time = clock;
    if exist(fname,'file') == 2
        delete(fname);
    end
    tab = cell(5,14);
    tab{1,1} = ['VieVS multisched summary: ' PARA.sessname];
    tab{2,1} = sprintf('created: %04d-%02d-%02d %02d:%02d',time(1:5));

    tab{5,1} = 'Version';
    tab{5,2} = '# scans';
    tab{5,3} = 'scans/h';
    tab{5,4} = 'min scans/h';
    tab{5,5} = '#obs';
    tab{4,6} = 'sky coverage [1-13]';
    tab{4,8} = 'scan time [%]';
    tab{4,10} = 'slew time [%]';
    tab{4,12} = 'idle time [%]';
    tab{4,14} = 'calibration time [%]';

    for i = 1:5
        tab{5,2*i+4} = 'mean';
        tab{5,2*i+5} = 'std';
    end

    xlswrite(fname,tab,sheet,'A1');

elseif level3Flag == 0
    line = line+5;
    tab = cell(1,14);
    tab{1,1} = sprintf('V%03d',line-5);
    tab{1,2} = stat.nscan;
    tab{1,3} = round(mean(stat.scansPerH),2);
    tab{1,4} = round(min(stat.scansPerH),2);

    tab{1,5} = stat.nobs;

    tab{1,6} = round(mean(stat.meanSky),2);
    tab{1,7} = round(std(stat.meanSky),2);

    tab{1,8} = round(mean(stat.pScanTime),2);
    tab{1,9} = round(std(stat.pScanTime),2);

    tab{1,10} = round(mean(stat.pSlewTime),2);
    tab{1,11} = round(std(stat.pSlewTime),2);

    tab{1,12} = round(mean(stat.pIdleTime),2);
    tab{1,13} = round(std(stat.pIdleTime),2);

    tab{1,14} = round(mean(stat.pConstTime),2);
    tab{1,15} = round(std(stat.pConstTime),2);

    xlswrite(fname,tab,sheet,sprintf('A%d',line));

elseif line == 0 && level3Flag == 1
    tab = cell(2,33);
    tab{2,1} = '#simulations';
    tab{2,2} = '#stations';
    tab{2,3} = '#sources';
    tab{1,4} = 'station coordinates';
    tab{1,6} = 'source coordinates';

    tab{1,8} = 'xpol';
    tab{1,10} = 'ypol';
    tab{1,12} = 'dut1';
    tab{1,14} = 'nutdx';
    tab{1,16} = 'nutdy';

    tab{1,18} = 'station coordinate x';
    tab{1,20} = 'station coordinate y';
    tab{1,22} = 'station coordinate z';

    tab{1,24} = 'source ra';
    tab{1,26} = 'source de';
    tab{1,28} = 'source ra in NNR';
    tab{1,30} = 'source de in NNR';
    tab{1,32} = 'source ra not in NNR';
    tab{1,34} = 'source de not in NNR';
    for i = 1:16
        tab{2,2*i+2} = 'mean std estimates';
        tab{2,2*i+3} = 'mean sigma';
    end

    xlswrite(fname,tab,sheet,'P4');

elseif level3Flag == 1
    line = line+5;
    tab = cell(1,32);

    if isnan(stat.coor3d)
        tab{1,1} = '-';
        tab{1,2} = '-';
        tab{1,15} = '-';
        tab{1,16} = '-';
        tab{1,17} = '-';
        tab{1,18} = '-';
        tab{1,19} = '-';
        tab{1,21} = '-';
    else
        tab{1,1} = round(stat.coor3d,4);
        tab{1,2} = round(stat.coor3dm,4);

        tab{1,15} = round(stat.coorx,4);
        tab{1,16} = round(stat.coorxm,4);
        tab{1,17} = round(stat.coory,4);
        tab{1,18} = round(stat.coorym,4);
        tab{1,19} = round(stat.coorz,4);
        tab{1,20} = round(stat.coorzm,4);
    end

    if isnan(stat.soude)
        tab{1,3} = '-';
        tab{1,4} = '-';

        tab{1,21} = '-';
        tab{1,22} = '-';
        tab{1,23} = '-';
        tab{1,24} = '-';
        tab{1,25} = '-';
        tab{1,26} = '-';
        tab{1,27} = '-';
        tab{1,28} = '-';
        tab{1,29} = '-';
        tab{1,30} = '-';
        tab{1,31} = '-';
        tab{1,32} = '-';
    else
        tab{1,3} = round(stat.sou2d,4);
        tab{1,4} = round(stat.sou2dm,4);

        tab{1,21} = round(stat.soura,4);
        tab{1,22} = round(stat.souram,4);
        tab{1,23} = round(stat.soude,4);
        tab{1,24} = round(stat.soudem,4);
        tab{1,25} = round(stat.sourainNNR,4);
        tab{1,26} = round(stat.sourainNNRm,4);
        tab{1,27} = round(stat.soudeinNNR,4);
        tab{1,28} = round(stat.soudeinNNRm,4);
        tab{1,29} = round(stat.souranotinNNR,4);
        tab{1,30} = round(stat.souranotinNNRm,4);
        tab{1,31} = round(stat.soudenotinNNR,4);
        tab{1,32} = round(stat.soudenotinNNRm,4);
    end

    if isnan(stat.xpol)
        tab{1,5} = '-';
        tab{1,6} = '-';
    else
        tab{1,5} = round(stat.xpol,6);
        tab{1,6} = round(stat.xpolm,6);
    end

    if isnan(stat.ypol)
        tab{1,7} = '-';
        tab{1,8} = '-';
    else
        tab{1,7} = round(stat.ypol,6);
        tab{1,8} = round(stat.ypolm,6);
    end

    if isnan(stat.dut1)
        tab{1,9} = '-';
        tab{1,10} = '-';
    else
        tab{1,9} = round(stat.dut1,6);
        tab{1,10} = round(stat.dut1m,6);
    end

    if isnan(stat.nutdx)
        tab{1,11} = '-';
        tab{1,12} = '-';
    else
        tab{1,11} = round(stat.nutdx,6);
        tab{1,12} = round(stat.nutdxm,6);
    end

    if isnan(stat.nutdy)
        tab{1,13} = '-';
        tab{1,14} = '-';
    else
        tab{1,13} = round(stat.nutdy,6);
        tab{1,14} = round(stat.nutdym,6);
    end

    tab = {stat.isim,stat.ista,stat.isrc,tab{:}};
    xlswrite(fname,tab,sheet,sprintf('P%d',line));

end


end
