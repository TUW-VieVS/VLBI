% ONLY READ-AS-IS FUNCTION! NO SMART MERGING!!
% BETTER (MAYBE): readSnx.m
% This function reads in SINEX files and provides the data of that file in
% a struct.
% It works as follows: In the beginning there should be defined the formats
% of (if possible) all different blocks. Below is a loop going through all 
% lines of the file and reading in line-by-line the data according to the
% format specified.
% So there is an easy handling (and adding of new blocks).
% --> the column.
%
% e.g. block 'SITE/ID' at the beginning
%  would give an output:
%  block(1).name ... 'SITE/ID'
%  block(1).description ... 'CODE PT DOMES____ T STATION_DESCRIPTION___ APPROX_LON_ APPROX_LAT_ APP_H__'
%  block(1).data{1} ... {7225, 7298, 7380, 7224, ...}
%  block(1).data{2} ... {'A', 'A', 'A', 'A', ...}
%  block(1).data{3} ... {'40408S002', '40424S007', ...}
%  block(1).name .....
%
% - if block is empty (only +blockname -blockname): nothing is read in
% - mind the ' ' in the beginning of each format description


function block=sinex_reader(sinexFile)

if ~exist(sinexFile, 'file')
    fprintf('sinex file does not exist\nreturning\n');
    return;
else
    fid=fopen(sinexFile);
end

% 1. Format definintions of blocks
blockNames={'FILE/REFERENCE',...
    'SITE/ID', ...
    'FILE/COMMENT', ...
    'INPUT/ACKNOWLEDGMENTS', ...
    'SOLUTION/DISCONTINUITY'};
blockFormats={' %17s %60s',...
    ' %4f %2s %9s %1s %22s %3f %2f %4f %3f %2f %4f %7f',...
    ' %s',...
    ' %3s %75s',...
    ' %4f %2s %4f %1s %2f:%3f:%5f %2f:%3f:%5f %50s'};

% 2. Preallocating block struct
block=struct;
blockIndex=1;

% 3. Read SINEX lines
while ~feof(fid)
    % read one line
    curLine=[fgetl(fid), ' '];
    
    % if we have a '+' at the beginning of the line-> we have a new block
    if strcmp(curLine(1), '+')
        % try to find the block name in the definitions
        logBlockFound=~cellfun(@isempty, strfind(blockNames, deblank(curLine(2:end))));
        
        % if we "know" the block (ie it is defined above)
        if sum(logBlockFound)>0
            % get index of block in 'blockNames'
            indBlockNames=find(logBlockFound);
            
            
            
            % get line after +blockname line
            curLine=[fgetl(fid), ' '];

            % if we have actually data (and not '-blockname', which is end of block)
            if ~strcmp(curLine(1), '-')
                % write name and description to struct
                block(blockIndex).name=blockNames{indBlockNames};
                
                % if we have description line (ie first character is not ' ')
                if ~strcmp(curLine(1), ' ')
                    % write description line to block
                    block(blockIndex).description=curLine;
                    
                    % and read following line
                    curLine=fgetl(fid);
                end
                
                % read in all lines that start with a blank
                rowInd=1;
                while strcmp(curLine(1), ' ')

                    % read in data according to file format
                    dataFromRow=textscan(curLine, blockFormats{indBlockNames}, 'delimiter', '\n');

                    % put data into blockstruct
                    for k=1:size(dataFromRow,2)
                        block(blockIndex).data{k}{rowInd}=dataFromRow{k};  
                    end

                    % increase row index for all data rows per block
                    rowInd=rowInd+1;

                    % read next line
                    curLine=[fgetl(fid), ' '];
                end

                % increase block index by one (next block should be written to
                % next block
                blockIndex=blockIndex+1;
            end
        else
            fprintf('block ''%s'' was not found in description\nand therefore not read in\n', deblank(curLine(2:end)));
        end
    end
end

fclose(fid);

        
    


