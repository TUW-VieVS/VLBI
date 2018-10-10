% Purpose  
%   counts schedules in  multiSched.txt        
% History  
%	2016-11-14 Matthias Schartner: created
%   2017-02-23 Matthias Schartner: possibillity to calculate only parts of
%                                  the multisched list

function [ start_index, end_index ] = multiSchedCounter( filename,PARA )

boolIdxStart = false;
start_index = 0;
boolIdxEnd = false;
end_index = 0;

fid = fopen(filename);
loopCount = 0;

while ~feof(fid)
    line = fgetl(fid);    
    % check for comment
    if strcmp(line(1),'+')
        val = strsplit(line);
        if strcmp(val{2},'start_index')
            start_index = str2double(val{3});
            boolIdxStart = true;
        end
        if strcmp(val{2},'end_index')
            end_index = str2double(val{3});
            boolIdxEnd = true;
        end
    end
    if strcmp(line(1),'#')
        loopCount = loopCount + 1;
    end
end

if ~(boolIdxStart && boolIdxEnd)
    error('couldn''t find start-index and/or end-index in multiSched.txt');
end

if start_index <1
    start_index = 1;
end

if end_index == Inf
    end_index = loopCount;
end
if end_index > loopCount
    end_index = loopCount;
end
fprintf(PARA.fid_header,'\nMULTISCHED session:\n');
fprintf(PARA.fid_header,'    A total of %d schedules were found\n',loopCount);
fprintf(PARA.fid_header,'    processing schedule %d to %d \n\n',start_index,end_index);
fprintf('\nMULTISCHED session:\n');
fprintf('    A total of %d schedules were found\n',loopCount);
fprintf('    processing schedule %d to %d\n\n',start_index,end_index);

end

