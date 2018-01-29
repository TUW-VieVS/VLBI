% read information from freq.cat

% created 22-07-2015 by Lucia Plank

function obsmode=ofreq(obsmode,filename)

% open freq catalog file
fid = fopen(filename, 'r');


% get freq information
    frewind(fid);
    entryfinished=0;
    while ~feof(fid)  && entryfinished==0;
        line = fgetl(fid);
        if length(line>0)
        if (line(1) ~= '*' && line(1) ~= '-' && line(1) ~= '&' && line(1)~= '!')
            freqname=textscan(line,'%s %s %s %s');
            if strcmp(freqname{1},obsmode.freqname)
                obsmode.sx=freqname{2};
                obsmode.rxname = freqname{4};
                i=1;
                while entryfinished==0
                line=fgetl(fid);
                if strcmp(line(1),'*')
                    entryfinished=1;
%                 elseif strcmp(line(1),'*')
                else    
                    if strcmp(line(1),'-') 
                        obsmode.freq(i,:)=textscan(line(3:end),'%s %s %s %s %s %s %s');
                        i=i+1;
                    end
                end
                end
                
            end
        end
        end

    end


% close file
fclose(fid);