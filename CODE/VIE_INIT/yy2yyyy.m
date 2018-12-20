function outYrNum = yy2yyyy(inYY)

% This function converts a two-digit year to the full number (matlab:
% double)
%
% INPUT can be char: '15' or double [15], can be cellstr or vector
% OUTPUT will be double: [2015] (same size as input)
%
% 79 to 99 -> 1979 to 1999
% 00 to 78 -> 2000 to 2078
%
% coded: 04.11.2015 by Matthias Madzak
%

if ~isnumeric(inYY)
    inYY=cellstr(inYY);
    outYrNum=zeros(size(inYY));
    for k=1:numel(inYY)
        outYrNum(k)=str2double(inYY{k});
    end
else
    outYrNum=inYY;
end

larger79=outYrNum>=79;
outYrNum(larger79)=outYrNum(larger79)+1900;
outYrNum(~larger79)=outYrNum(~larger79)+2000;