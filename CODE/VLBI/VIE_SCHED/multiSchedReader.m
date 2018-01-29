% Purpose  
%   overwrites parameters from multiSched.txt    
% History  
%	2016-11-14 Matthias Scharner: created

function [ station, twin, source, PARA ] = multiSchedReader( station, twin, source, PARA, filename, currentLoopCount )

fid = fopen(filename);
lnumber = 0;
loopCount = 0;

while ~feof(fid)
    line = fgetl(fid);
    lnumber = lnumber +1;
    
    % check for comment
    if strcmp(line(1),'%')
        continue;
    end
    
    % check for new schedule
    if strcmp(line(1),'#')
        loopCount = loopCount + 1;
    end
    
    % if you fund your current schedule
    if loopCount == currentLoopCount
        t = clock;
        tstr1 = sprintf('%4d/%02d/%02d', t(1), t(2), t(3));
        tstr2 = sprintf('%02d:%02d:%02d', t(4), t(5), round(t(6)));
        fprintf(PARA.fid_body,'start: %s, %s\n', tstr2, tstr1);
        fprintf(PARA.fid_body,'*******************************************\n');
        fprintf(PARA.fid_body,'**   Now starting schedule number: %3d   **\n',loopCount);
        fprintf(PARA.fid_body,'*******************************************\n');
        fprintf(PARA.fid_body,'changes:\n');
        while ~feof(fid)
            line = fgetl(fid);
            lnumber = lnumber +1;
            
            % check for comment
            if strcmp(line(1),'%')
                continue;
            end
            
            % check for new schedule
            if strcmp(line(1),'#')
                break;
            end
            
            % try to evaluate text as MATLAB expression  
            try
                eval(line)
                fprintf(PARA.fid_body,['    ' line '\n']);
            catch ME
                warning([ME.message,'\n',line,'\n','In multiSched.txt at %d'],lnumber)
            end
        end
        break
    end
    
end

fclose(fid);
end

