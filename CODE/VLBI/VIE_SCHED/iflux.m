% #########################################################################
% #     iflux
% #########################################################################
%
% DESCRIPTION
%    Read source flux information from the flux catalog file.
%
% CREATED  
%   - 2010-04-06   Jing SUN
%
% REFERENCES
%
%
% COUPLING
%   - TLE_propagation
%
%
% INPUT
%  
% - filename            - filename of flux catalog (flux.cat)
% - source              - source data structure
% - PARA                - Global Parameter structure.
%
% OUTPUT
% - source              - source data structure, flux parameters added
%
% CHANGES:
%   - 2015-02-25: A. Hellerschmied: Modifications to speed up the loading
%       process significantly. Waitbar added.
%   - 2016-06-06: M. Schartner: Modifications to speed up the loading
%       process extremely, bugfix


function [source] = iflux(filename, source, PARA)

% open flux catalog file
fid = fopen(filename, 'r');
if (fid < 0)
    error('    no flux catalog file %s !\n', filename);
end

% Read data flux catalog file
fid = fopen(filename, 'r');
C = textscan(fid,'%s %s %s %[^\n]','commentstyle','*','Delimiter',' ','multipledelimsasone',1);
fclose(fid);
src_Name = C{1};
XorS = C{2};
BorM = C{3};
fluxInfo = C{4};

srcnum = length(source);
X = strcmp(XorS,'X');
S = strcmp(XorS,'S');

for isrc = 1:srcnum
    idx1 = strcmp(src_Name,strtrim(source(isrc).name)) | strcmp(src_Name,strtrim(source(isrc).commoname));
    
    % find matching Band
    for iband = 1 : PARA.MAX_BANDNUM
        if(PARA.BAND(iband)=='S')
            idx2 = S;
        elseif(PARA.BAND(iband)=='X')
            idx2 = X;
        end
        idx = idx1&idx2;
        sidx = sum(idx);
        idx = find(idx);
        
        for i = 1:sidx
            [ftmp,np] = sscanf(fluxInfo{idx(i)}, '%f');  
            if strcmp(BorM{idx(i)},'B')     
                source(isrc).fluxpara(iband,1:np) = ftmp(1:np);
                source(isrc).fluxpara(iband,PARA.MAX_FLUXPARA) = 1.0;

            elseif strcmp(BorM{idx(i)},'M')    
                if (source(isrc).fluxpara(iband,1) < 1.0d-3)
                    source(isrc).fluxpara(iband,1:6) = ftmp(1:6,1);
                elseif (source(isrc).fluxpara(iband,7) < 1.0d-3)
                    source(isrc).fluxpara(iband,7:12) = ftmp(1:6,1);
                elseif (source(isrc).fluxpara(iband,13) < 1.0d-3)
                    source(isrc).fluxpara(iband,13:18) = ftmp(1:6,1);
                end 
                source(isrc).fluxpara(iband,PARA.MAX_FLUXPARA) = 2.0;
            end
        end
    end
end