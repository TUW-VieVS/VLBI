% read information from tracks.cat

% created 22-07-2015 by Lucia Plank
% revised 05-09-2016 by Lucia Plank: fanout added
% revised 11-9-2017 by Lucia McCallum: correct reading of empty lines

function [tracks,fanout]=otrack(filename, trackmode)

% open rec catalog file
fid = fopen(filename, 'r');

% get tracks information
    frewind(fid);
    entryfinished=0;
    while ~feof(fid)  && entryfinished==0;
        line = fgetl(fid)
        if length(line>4)
        if (line(1) ~= '*' && line(2) ~= '-' && line(1) ~= '&' && line(1)~= '!')
            trackslin=textscan(line,'%s %s %s');
            if strcmp(trackslin{1},trackmode)
                fanout=cell2mat(trackslin{2});
                i=1;
                while entryfinished==0
                line=fgetl(fid)
                if line==-1 
                    entryfinished=1;
                elseif strcmp(line,'')
                    entryfinished=1;
                elseif strcmp(line(1),'*')
                    entryfinished=1;
                else    
                    if strcmp(line(2),'-') 
                        tracks(i,:)=textscan(line(3:end),'%s %s');
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