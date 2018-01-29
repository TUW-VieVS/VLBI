% Purpose  
%   Get the distribution of sources [1/2/4].
% History  
%   2010-04-06   Jing SUN   Created
%   


function [cat, srcat, catpair] = isrcseq(source, PARA)

% generate the spatial grid 'cat'
gridnum = 0;
dg_de = 20;
% level 1
ide = 90; 
dg_ra = 360/2;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 2
ide = ide - dg_de; 
dg_ra = 360/6;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 3
ide = ide - dg_de; 
dg_ra = 360/8;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 4
ide = ide - dg_de; 
dg_ra = 360/10;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 5
ide = ide - dg_de; 
dg_ra = 360/12;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 6
ide = ide - dg_de; 
dg_ra = 360/10;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end 
% level 7
ide = ide - dg_de; 
dg_ra = 360/8;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end
% level 8
ide = ide - dg_de; 
dg_ra = 360/6;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end
% level 9
ide = ide - dg_de; 
dg_ra = 360/2;
for ira = 0 : dg_ra : (360-dg_ra)
    gridnum = gridnum + 1;
    cat(gridnum).demax = (ide)*pi/180;
    cat(gridnum).demin = (ide - dg_de)*pi/180;
    cat(gridnum).ramin = (ira)*pi/180;
    cat(gridnum).ramax = (ira + dg_ra)*pi/180;
    cat(gridnum).demean = (cat(gridnum).demin+cat(gridnum).demax)/2;
    cat(gridnum).ramean = (cat(gridnum).ramin+cat(gridnum).ramax)/2;
end

% calculate the 'srcat'
srcnum = length(source);
for igrid = 1 : gridnum
    srcat(igrid).num = 0;
    srcat(igrid).src(1:srcnum) = 0;
end
for isrc = 1 : srcnum
    ra = source(isrc).ra;
    de = source(isrc).de;
    for igrid = 1 : gridnum
        if ((ra>=cat(igrid).ramin)&(ra<cat(igrid).ramax)&(de>=cat(igrid).demin)&(de<cat(igrid).demax))
            srcat(igrid).num = srcat(igrid).num + 1;
            srcat(igrid).src(srcat(igrid).num) = isrc;
            break;
        end
    end
end

% analyze the subnet 'catpair' [1/2/4]
if (PARA.SRCNUM == 1)
    [catpair] = isub1grid(cat);
elseif(PARA.SRCNUM == 2)
    [catpair] = isub2grid(cat);
elseif (PARA.SRCNUM == 4)
    [catpair] = isub4grid(cat);
end


