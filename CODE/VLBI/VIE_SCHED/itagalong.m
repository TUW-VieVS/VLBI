% Purpose  
%   Read tagalong.txt file.
% History  
%   2013-12-30   Jing SUN   Created
%   


function [tagalongname] = itagalong(filename,PARA)

% open tagalong.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf(PARA.fid_header,'No tagalong.txt file !\n');
    tagalongname = '';
    return;
end

% read the tag-along information
num = 0;
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ~strcmp(line(1), ' ') 
        continue;
    end
    num = num + 1;
    tagalongname(num,:) = sscanf(line(2:linelength), '%s', 1);
    fprintf(PARA.fid_header,'There is tag-along mode for station : %s', tagalongname(num,:));
    fprintf(PARA.fid_header,'\n');
end

% close file
fclose(fid);