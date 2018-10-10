% This function checks all optimization conditions and returns a list of
% good and bad sources. 
%
% CREATED: 21.02.17 - Matthias Schartner
%
% CHANGELOG: 
function [ flag, goodSources, badSources ] = optimizer_get_good_sources( sched, con )

allScans = [sched.scan];
allsrcid = [allScans.srcid];
[unique_Srcid,~,idx] = unique(allsrcid);
allnsta  = [allScans.nsta];

for i = 1:length(allScans)
    allScans(i).endmjd = max([allScans(i).sta.endmjd]);
    allScans(i).nobs = (allScans(i).nsta*(allScans(i).nsta-1))/2;
end

ntime = zeros(size(unique_Srcid));
nscans = zeros(size(unique_Srcid));
nsta = zeros(size(unique_Srcid));
nbl = zeros(size(unique_Srcid));

for i = 1:length(unique_Srcid)
    thisSrc = unique_Srcid(i);
    
    thisScans = allScans(idx==i);
    
    minMJD = min([thisScans.startmjd]);
    maxMJD = max([thisScans.endmjd]);
    ntime(i) = (maxMJD-minMJD)*24;
    
    nscans(i) = length(thisScans);
    
    this_nsta = [thisScans.nsta];
    nbl(i) = sum( (this_nsta.*(this_nsta-1))/2 );
    
    allSta = [thisScans.sta];
    nsta(i) = length(unique([allSta.staid]));
end


allTypes = {'nbl','nscan','nsta','ntime'};

split = strsplit(con);
split = split(2:2:end);

con = struct('type',{},'value',{},'combination',{});
counter = 1;
for i=1:(length(split)-2)/3+1
    con(1).type(i) = split(counter); 
    counter = counter+1;
    con(1).value(i) =  split(counter); 
    counter = counter+1;
    if counter < length(split)
        con(1).combination(i) = split(counter); 
        counter = counter+1;
    end
end

bool = false(length(unique_Srcid),length(con.type));
for i = 1:length(con.type)
    
    idx = find(strcmp(con.type(i),allTypes));
    value = str2double(con.value(i));
    
    if idx == 1
        bool(:,i) = nbl>value;
    elseif idx == 2
        bool(:,i) = nscans>value;
    elseif idx == 3
        bool(:,i) = nsta>value;
    elseif idx == 4
        bool(:,i) = ntime>value;
    end
end

boolfinal = bool(:,1);
for i = 2:length(con.type)
    if strcmp(con.combination(i-1),'AND')
        boolfinal = boolfinal & bool(:,i);
    else
        boolfinal = boolfinal | bool(:,i);
    end
end


goodSources = unique_Srcid(boolfinal);
badSources = unique_Srcid(~boolfinal);
if all(boolfinal)
    flag = 1;
else
    flag = 0;
end

end

