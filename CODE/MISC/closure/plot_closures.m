% Matthias Schartner, Hana Krasna

% process_list=['2020/20OCT29PI [vgosDB]'
%               '2020/20OCT29VE [vgosDB]'];
% typ = 1; % 1 - baseline delay, 2 - geocentric delay


function plot_closures(process_list,typ)
close all

pth = '../OUT/CLOSURES/';


%% Read input of closure delays: 
sources = [];
triangles = [];
closure_delay = [];
time = [];
%closures_XA = dir(fullfile('Closures','closures_*XA.txt'));


for i = 1:size(process_list,1)
%     fid = fopen(fullfile(closures_XA(i).folder, closures_XA(i).name));

    if typ == 1
        fid = fopen([pth 'VALUES/closures_B_' process_list(i,6:14) '.txt']);
    else
        fid = fopen([pth 'VALUES/closures_G_' process_list(i,6:14) '.txt']);
    end
    C = textscan(fid,'%s %s %s %s %f %f %f %f %f %f %s','CommentStyle','%');
    fclose(fid);
    
    
    
    sources=[];
    triangles_u = [];
    triangles = [];
    closure_delay = [];
    sigma_closure_delay = [];
    mjd = [];
    time_str = [];
    time = []; 
    delayFlag = [];
    
    sources = [sources; C{1}];
    triangles_u = [C{2} C{3} C{4}];
    triangles = [triangles; join(triangles_u)];
    closure_delay = [closure_delay; C{5}];
    sigma_closure_delay = [sigma_closure_delay; C{6}];
    delayFlag = [delayFlag; C{7} C{8} C{9}];       
    mjd = [mjd; C{10}];
    time_str = [time_str; C{11}];
    time = [time; datetime(time_str,'InputFormat','yyyy.MM.dd_HH:mm:SS')];    
    
    
    


    all_sources = unique(sources);
    all_triangles = unique(triangles);

    %% plot closure delays
    tria=['all stations'];
    
    src = 'all sources';
    bool = logical(1:length(mjd)); % all sources
    
    
    
    % src = sources{1}; % pick a source and plot all closure delays
    % src = '0749+540'
    % src = '1300+580'
    % src = '1851+488'
    % bool = strcmp(src, sources);


    % tria = ['ONSA13SW RAEGYEB WETTZ13S'];
 %    tria = ['ONSA13NE RAEGYEB WETTZ13S'];
  %   tria = ['ISHIOKA ONSA13SW WETTZ13S'];
     
     
   %  bool = strcmp(tria, triangles);
    
    % bool = strcmp(src, sources) & strcmp(tria, triangles); % specific source and specific triangle

    time_bool=time(bool);
    mjd_bool=mjd(bool);
    closure_delay_bool = closure_delay(bool);
    sigma_closure_delay_bool = sigma_closure_delay(bool);
    delayFlag_bool = delayFlag(bool,:);



    % check delay flag
    dflimit = 0; % set delay flag limit
    r=[];
    [r,~,~] = find(delayFlag_bool>dflimit); % find the bad
    id = setdiff(1:length(delayFlag_bool),unique(r));
    if ~isempty(unique(r)) & size(delayFlag_bool,1)>3
            tm = time_bool(id);
            mjdplot = mjd_bool(id);
            clodel = closure_delay_bool(id);
            sigclodel = sigma_closure_delay_bool(id);
    else
            tm = time_bool;
            mjdplot = mjd_bool;
            clodel = closure_delay_bool;
            sigclodel = sigma_closure_delay_bool;
    end


    [uni_mjd_bool, id_uni_mjd_bool]=unique(mjdplot);
    uni_time_bool = tm(id_uni_mjd_bool);


    figure(i)
    errorbar(mjdplot, clodel*1e12, sigclodel*1e12,'.')
    
    xticks(uni_mjd_bool)
    xticklabels(datestr(uni_time_bool))
    xtickangle(45)

    title([sprintf('%s', src) ' ' sprintf('%s', tria)])
    ylabel('\tau_{ABC} [ps]')

    if typ == 1
        print( '-dpdf' ,'-r800',[pth 'PLOTS/' src '_' tria '_B_' process_list(i,6:14)])
    else
        print( '-dpdf' ,'-r800',[pth 'PLOTS/' src '_' tria '_G_' process_list(i,6:14)])

    end

 end