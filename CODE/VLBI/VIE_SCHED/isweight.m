% Purpose  
%   Read source weight information from the sweight.txt catalog file.
% History  
%   2014-06-20   David MAYER   Created
%   2014-09-09


function [source] = isweight(filename, source)

fid = fopen(filename, 'r');
if (fid < 0)
    fprintf('    no sweight.txt file !\n');
    return;
end

% read sweight.txt file
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if isempty(line)
        continue;
    end
    if ~strcmp(line(1), ' ') 
        continue;
    end
    strcmp1 = sscanf(line(2:9), '%s', 1);
    for i=1:length(source)
        if strcmp(deblank(source(i).name),deblank(strcmp1))
            [w] = sscanf(line(10:linelength), '%f');
            source(i).weight = w;
        end
    end
end

% close file
fclose(fid);


