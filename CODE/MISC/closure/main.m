format compact;
clc, clear, close all;

%% generate list of closure delays from CONT17 XA sessions:
XA = dir(fullfile('Data','14*XA'))

%for i = 1:length(XA) % <-- run over all XA sessions
for i = 1:1 % <-- run only one XA session
    fprintf('Read: %s \n',XA(i).name);
    scans = getScanList(fullfile(XA(i).folder,XA(i).name));
    fprintf('Find closures: %s \n',XA(i).name);
    scans2closures(XA(i).name, scans);
end

%% Read input of closure delays: 
sources = [];
triangles = [];
closure_delay = [];
time = [];
closures_XA = dir(fullfile('Closures','closures_*XA.txt'));
for i = 1:length(closures_XA)
    fid = fopen(fullfile(closures_XA(i).folder, closures_XA(i).name));
    C = textscan(fid,'%s %s %s %s %f %s');
    fclose(fid);
    sources = [sources; C{1}];
    triangles_u = [C{2} C{3} C{4}];
    triangles = [triangles; join(triangles_u)];
    closure_delay = [closure_delay; C{5}];
    time_str = C{6};
    time = [time; datetime(time_str,'InputFormat','yyyy.MM.dd_HH:mm:SS')];    
end

all_sources = unique(sources);
all_triangles = unique(triangles);

%% randomly pick one source and plot all closure delays
src = sources{1};

bool = strcmp(src, sources);
plot(time(bool), closure_delay(bool)*1e12)
title(sprintf('%s', src))
ylabel('\tau_{ABC} [ps]')

