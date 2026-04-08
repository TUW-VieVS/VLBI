% ************************************************************************
%   Description:
%   Read the file VieVS\DATA\GLOB\TRF\VELOC\TIES\ and make a structure of
%   it. It containts names of antennas, where the velocity is the same.
%
%
%   Input:
%      file_velties      name of the file with path
%
%         FORMAT of the file: (station names (8character) and 2 spaces between the names and the line has to end with "\")
%         ...
%         PIETOWN   VLA       VLA-N8    \
%         WETTZELL  TIGOWTZL   \
%         YLOW7296  YELLOWKN     \
%         ...
%
%   Output:                
%      velties           structure with names of antennas, where the
%                        velocity should be equal
%
%   Coded for VieVS: 
%   21 Jun 2011 by Hana Spicakova
%
%   Revision: 
%%


function velties = read_velocity_ties(file_velties)

velties=[];

fid=fopen(file_velties,'r');

if fid~=-1
    j=0;

    while ~feof(fid);
        j=j+1;
        k=1;
        str=fgetl(fid);

        velties(j).aname = split(str);
        clear str    
    end

    fclose(fid);
end