% read information from rec.cat

% created 22-07-2015 by Lucia Plank

% CHANGES:
% 2015-07-23: A. Hellerschmied & D. Mayer: minor revision (problem with lines in rec.cat with more than 4 columns)

function station=orec(station, filename, obsmode)

% open rec catalog file
fid = fopen(filename, 'r');

%length recname
lcmp=length(obsmode.recname)+1;

% get rec information
    frewind(fid);
    entryfinished=0;
    while ~feof(fid)  && entryfinished==0;
        line = fgetl(fid);
        linelength = length(line);
        %line(2:lcmp+2)
%         obsmode.recname
        if (linelength == lcmp) && (line(1) ~= '*' && strcmp(line(2:lcmp),obsmode.recname))
            entryfinished=0; i=1;
            while entryfinished==0
                line = fgetl(fid);
                if strcmp(line(1:4),'****')
                    entryfinished=1;
                elseif strcmp(line(1),'*')
                else    
                    if strcmp(line(2),'-') 
                        rec(i,:)=textscan(line(3:end),'%s %s %s %s %*[^\n]');
                        i=i+1;
                    end
                end
            end
            if entryfinished==1
              for is=1:length(station)
                  for ind = 1 : length(rec)
                      if strcmp(deblank(station(is).name), rec{ind}{1})
                         break; 
                      end
                  end
                  
                  
                  station(is).trackmode=rec{ind,3};                  
              end
            end     
        end
    end


% close file
fclose(fid);