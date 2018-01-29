% Purpose  
%   Read source information from the STAR-source catalog file.
% History  
%   2016-05-15 M. Schartner & L. Plank: created
%   2016-09-06 Matthias Schartner: bugfix


function [source] = isource_star(INFILE, PARA, source1)

% read source file
fprintf('2.1* read star source file: %s \n', INFILE.sourcestar);
fid = fopen(INFILE.sourcestar, 'r');
if (fid < 0)
    error('    no star source catalog file %s !\n', INFILE.sourcestar);
end

srcnum1 = length(source1);
srcnumBeginning = srcnum1; 

while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength > 9) & (line(1) == ' '))
        % exclude the repeated source information
        ifre = 0;
        for isrc = 1 : srcnum1
            if strcmp(source1(isrc).name, line(2:9))
               source1(isrc).star = 1;
               ifre = 1;
               break;
            end
        end
        if ifre == 1
            continue;
        end
        srcnum1 = srcnum1 + 1;
        source1(srcnum1).name = line(2:9);
        if (line(11) == '$')
            source1(srcnum1).commoname = '        ';
        else
            source1(srcnum1).commoname = line(11:18);
        end
        [tmprd, count, errmsg, nextindex] = sscanf(line(20:linelength), '%f', 3);
        source1(srcnum1).rah = tmprd(1);   
        source1(srcnum1).ram = tmprd(2);   
        source1(srcnum1).ras = tmprd(3);
        source1(srcnum1).ra = (source1(srcnum1).rah + source1(srcnum1).ram/60 + source1(srcnum1).ras/3600) * 15 * pi / 180; 
        index = 20 + nextindex - 1;
        [il] = findstr(line(index:linelength), '-');
        if (length(il) > 0)
            sign = -1;
            source1(srcnum1).sign = '-';
        else
            sign = +1;
            source1(srcnum1).sign = '+';
        end
        [tmprd, count, errmsg, nextindex] = sscanf(line(index:linelength), '%f', 3);
        source1(srcnum1).ded = abs(tmprd(1)); 
        source1(srcnum1).dem = tmprd(2); 
        source1(srcnum1).des = tmprd(3);
        source1(srcnum1).de = sign * (source1(srcnum1).ded + source1(srcnum1).dem/60 + source1(srcnum1).des/3600) * pi /180;
        index = index + nextindex - 1;
        source1(srcnum1).info = line(index:linelength);
        source1(srcnum1).fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = 0.0;
        source1(srcnum1).ifp = 0;
        source1(srcnum1).pst = 0.0;
        source1(srcnum1).pdt = 0.0;
        source1(srcnum1).weight = 1.0;      % David Mayer, 2014 Aug 13
        source1(srcnum1).star = 1;
    end
end
fclose(fid);
fprintf('  number of sources in %s is %d\n', INFILE.source, srcnum1);

if srcnumBeginning ~= length(source1)
    % read file flux
    fprintf('2.2* read flux file for STAR sources: %s \n', INFILE.flux);
    [source1_Star] = iflux(INFILE.flux, source1(srcnumBeginning+1:end), PARA);

    % check source flux information
    srcnum = 0;
    for isrc = 1 : length(source1_Star)
        ifuse = 1;
        for ib = 1 : PARA.MAX_BANDNUM   % check all bands
            if (source1_Star(isrc).fluxpara(ib,PARA.MAX_FLUXPARA) < 0.5)
                ifuse = 0;
            end
        end
        if (ifuse == 1)
            srcnum = srcnum + 1;
            source_Star(srcnum) = source1_Star(isrc);
        end
    end
    if exist('source_Star','var')
        source = [source1(1:srcnumBeginning) source_Star];
    else
        source = source1(1:srcnumBeginning);
    end
else
    source = source1;
end


