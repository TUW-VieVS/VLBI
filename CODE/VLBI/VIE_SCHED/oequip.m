% Purpose  
%   Write $STATIONS (list of stations in this experiment).
%          the T lines : station datea acquisition terminal information
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-10 Matthias SCHARTNER 	rewrote function

function oequip(staname, staid, filename, fid_skd, PARA)

% open equip catalog file
fid = fopen(filename, 'r');
C = textscan(fid,'%s %s %s %s %s %s %s %s %s %[^\n]','delimiter',' ','multipledelimsasone',1,'commentstyle','*');
stanameCat = C{1};
staIDCat = C{2};
band1Cat = C{6};
band2Cat = C{8};
frewind(fid);
Cline = textscan(fid,'%*s %[^\n]','commentstyle','*');
fclose(fid);

thisBand = cell(PARA.MAX_BANDNUM,1);
for ib = 1 : PARA.MAX_BANDNUM
    thisBand{ib} = PARA.BAND(ib);
end

for i = 1:size(staname,1);
    thisStaname = staname(i,1:8);
    thisId = staid{i};
    bool1 = strcmp(stanameCat,strtrim(thisStaname));
    bool3 = false(size(bool1));
    for id = 1 : PARA.MAX_BANDNUM
        bool3 = bool3 | strcmp(band1Cat,thisBand{id});
        bool3 = bool3 | strcmp(band2Cat,thisBand{id});
    end
    bool = bool1 & bool3;
    if sum(bool)==1
        if ~strcmp(strtrim(thisId),staIDCat{bool})
            warning('Station ID for %s in equip.cat does not match with ID in antenna.cat',thisStaname)
        end
        fprintf(fid_skd, 'T %s\n', Cline{1}{bool});
    elseif sum(bool)==0
        warning('There is no Station %s in equip.cat',thisStaname)
        fprintf(fid_skd, 'T %s %s *WARNING: no match in equip.cat! Check catalogs\n',thisId,thisStaname);
    else
        bool2 = strcmp(staIDCat,strtrim(thisId));
        bool = bool2 & bool;
        if sum(bool)==1
            fprintf(fid_skd, 'T %s\n', Cline{1}{bool});
        elseif sum(bool)>1
            boolidFirst = find(bool,1);
            warning('More than one station with name %s and ID %s was found in the catalogs! Check catalogs\n',thisStaname,thisId)
            fprintf(fid_skd, 'T %s\n', Cline{1}{boolidFirst});
        else
            warning('There is no Station %s with ID %s in equip.cat',thisStaname,thisId)
            fprintf(fid_skd, 'T %s %s\n',thisId,thisStaname);
        end
    end
end


