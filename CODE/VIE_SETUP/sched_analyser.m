% 2016-07-50, M. Schartner: creator
% 2016-10-19, M. Schartner: bugfix, now loads source instead of source5
% 2017-05-29, M. Schartner: bugfix

function varargout = sched_analyser(varargin)
% SCHED_ANALYSER MATLAB code for sched_analyser.fig
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sched_analyser_OpeningFcn, ...
                   'gui_OutputFcn',  @sched_analyser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before sched_analyser is made visible.
function sched_analyser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sched_analyser (see VARARGIN)

% general parameters
handles.slider_endmjd.SliderStep = [0.05 0.05];
handles.slider_startmjd.SliderStep = [0.05 0.05];
handles.plotdown = 1;
handles.show_baseline.Value = 1;
handles.show_station_id.Value = 1;

% Choose default command line output for sched_analyser
handles.output = hObject;

%% load the Data
if isempty(varargin)
    load('../DATA/LEVEL5/sched.mat');
    handles.sched = sched;
    load('../DATA/LEVEL5/source.mat');
    handles.source = source;
    load('../DATA/LEVEL5/station.mat');
    handles.station = station;
    handles.figure1.Name='VieVS SCHED Analyser: ../DATA/LEVEL5';
    handles.input = 1;
    handles.pushbutton_load_star_catalog.Visible = 'off';
    handles.pushbutton_browse_for_star_catalog.Visible = 'off';
elseif strcmp(varargin{1},'file')
    fid = fopen(varargin{2});

    flagSource = false;
    flagStat = false;
    flagSched = false;
    while ~feof(fid)
        fline = fgetl(fid);

        if strcmp(strtrim(fline),'$SOURCES')
            source = struct();
            counter = 1;
            flagSource = true;
            while 1
                fline = fgetl(fid);
                if strcmp(fline(1),'$')
                    break
                end
                if strcmp(fline(1),'*')
                    continue
                end
                C = textscan(fline,'%s %s %f %f %f %f %f %f %[^\n]');
                source(counter).name = C{1}{1};
                if strcmp(C{2},'$')
                    source(counter).commoname = '        ';
                else
                    source(counter).commoname = C{2}{1};
                    source(counter).name = C{2}{1};
                end
                source(counter).ra = (C{3} + C{4}/60 + C{5}/3600) * 15 * pi / 180;
                source(counter).de = (C{6} + C{7}/60 + C{8}/3600) * pi /180;
                source(counter).info = C{9}{1};
                source(counter).star = 0;
                counter = counter+1;
            end
        end

        if strcmp(strtrim(fline),'$STATIONS')
            flagStat = true;
            station = struct();
            counter = 1;
            while 1
                fline = fgetl(fid);
                if strcmp(fline(1),'$')
                    break
                end
                if strcmp(fline(1),'*')
                    continue
                end

                if strcmp(fline(1),'A')
                    C = textscan(fline,'%s %s %s %s %f %f %f %f %f %f %f %f %f %s %s');
                    station(counter).name = C{3}{1};
                    station(counter).id = C{2}{1};
                    counter = counter+1;
                end


                if strcmp(fline(1),'P')
                    C = textscan(fline,'%s %s %s %f %f %f %f %f %f %[^\n]');
                    for i = 1:counter
                        if strcmp(station(i).name,C{3})
                            [station(i).llh(2),station(i).llh(1),station(i).llh(3)] = xyz2ell([C{4} C{5} C{6}]);
                            station(i).po = C{2}{1};
                            break;
                        end
                    end
                end
            end

            if length({station.id})~=length(unique({station.id}))
                error('No unique id for stations')
            end
        end

        if strcmp(strtrim(fline),'$SKED')
            flagSched = true;
            sched = struct();
            counter = 1;
            while 1
                fline = fgetl(fid);
                if strcmp(fline(1),'$')
                    break
                end
                if strcmp(fline(1),'*')
                    continue
                end
                C = textscan(fline,'%s %d %s %s %s %d %s %d %s %s %[^\n]');
                srcname = C{1};
                date = C{5};
%                 maxdur = C{6};
                stat = C{10};
                ids = stat{1}(1:2:end);
                info = C{11};
                info = strsplit(strtrim(info{1}));
                sched(counter).nscan = 1;

                scan = struct();
                yy = str2double(date{1}(1:2));
                d1 = datenum(sprintf('01-Jan-20%02d 00:00:00',yy));

                ddd = str2double(date{1}(3:5));
                hh = str2double(date{1}(6:7));
                mm = str2double(date{1}(8:9));
                ss = str2double(date{1}(10:11));
                d2 = d1+ddd-1+hh/24+mm/1440+ss/86400;
                startTime = datetime(d2,'ConvertFrom','datenum');
                dVector = datevec(startTime);

                scan.startmjd = date2mjd(dVector);
                scan.srcid = find(strcmp(srcname,{source.name}));
                ra = source(scan.srcid).ra;

                de = source(scan.srcid).de;
                scan.nsta = length(ids);

                sta = struct();
                for i = 1:scan.nsta
                    id = ids(i);
                    sta(i).staid = find(strcmp(id,{station.id}));
                    sta(i).startmjd = scan.startmjd;
                    lat = station(sta(i).staid).llh(2);
                    lon = station(sta(i).staid).llh(1);

                    [az, el, ha, dc] = zazel_s(scan.startmjd, lon, lat, ra, de);
                    sta(i).az = az;
                    sta(i).el = el;
                    sta(i).ha = ha;
                    sta(i).dc = dc;

                    sta(i).duration = str2double(info(end-scan.nsta+i));
                    sta(i).endmjd = scan.startmjd + sta(i).duration/86400;
                end
                scan.sta = sta;
                sched(counter).scan = scan;

                counter = counter+1;
            end
        end

    end
    handles.source = source;
    handles.station = station;
    handles.sched = sched;
    fclose(fid);
    handles.figure1.Name=sprintf('VieVS SCHED Analyser: %s',varargin{2});
    handles.input = 2;
elseif strcmp(varargin{1},'folder')
    fname = [varargin{2} '/sched.mat'];
    if exist(fname, 'file') == 2
        load(fname);
        handles.sched = sched;
    else
        error('cannot find file "sched.mat"!')
    end

    fname = [varargin{2} '/source.mat'];
    if exist(fname, 'file') == 2
        load(fname);
        handles.source = source;
    else
        error('cannot find file "source.mat"!')
    end

    fname = [varargin{2} '/station.mat'];
    if exist(fname, 'file') == 2
        load(fname);
        handles.station = station;
    else
        error('cannot find file "station.mat"!')
    end
    handles.figure1.Name=sprintf('VieVS SCHED Analyser: %s',varargin{2});
    handles.pushbutton_load_star_catalog.Visible = 'off';
    handles.pushbutton_browse_for_star_catalog.Visible = 'off';
    handles.input = 3;
end

% change name to commonname if it exists
commonnames = {handles.source.commoname};
boolcommon = ~strcmp(commonnames,'        ');
idx = find(boolcommon);
for i = 1:length(idx)
    handles.source(idx(i)).name = commonnames{idx(i)};
end

% save defining flag
defining = {handles.source.info};
defining = strfind(defining,'ICRF2 def');
defining = ~cellfun(@isempty,defining);
for i = 1:length(handles.source)
    handles.source(i).defining = defining(i);
end

% save scans individually to improve performance
allScans = [handles.sched.scan];
handles.scans = allScans;
for i = 1:length(handles.scans)
    thisScan = handles.source(handles.scans(i).srcid);
    handles.scans(i).endmjd = max([handles.scans(i).sta.endmjd]);
    handles.scans(i).nobs = (handles.scans(i).nsta*(handles.scans(i).nsta-1))/2;
    handles.scans(i).srcname = [' ' thisScan.name];
    handles.scans(i).star = thisScan.star;
    handles.scans(i).defining = thisScan.defining;
end
handles.select{2} = true(1,length(handles.scans));

% save observation individually to improve performance
allStanum = [allScans.nsta];
allStart = [allScans.startmjd];
allSource = [allScans.srcid];
allnSta = [allScans.nsta];
allnObs = (allnSta.*(allnSta-1))/2;
allStar = [handles.source(allSource).star];
alldefining = [handles.source(allSource).defining];
handles.obs = [allScans.sta];
c = 1;
for i = 1:length(allStanum)
    for j = 1:allStanum(i)
        handles.obs(c).startmjd=allStart(i);
        handles.obs(c).srcid = allSource(i);
        handles.obs(c).star = allStar(i);
        handles.obs(c).defining = alldefining(i);
        handles.obs(c).nsta = [' ' num2str(allnSta(i))];
        handles.obs(c).nobs = [' ' num2str(allnObs(i))];
        handles.obs(c).srcname = [' ' handles.source(allSource(i)).name];
        c = c+1;
    end
end
handles.select{1} = true(1,length(handles.obs));

% general settings
handles.startsched = handles.sched(1).scan.startmjd;
handles.startmjd = handles.sched(1).scan.startmjd;
if length([handles.sched(end).scan])==1
    end1 = max([handles.sched(end).scan.sta.endmjd]);
    handles.endsched = end1;
end
if length([handles.sched(end).scan])==2
    end1 = max([handles.sched(end).scan(1).sta.endmjd]);
    end2 = max([handles.sched(end).scan(2).sta.endmjd]);
    handles.endsched = max([end1,end2]);
end
if length([handles.sched(end).scan])==3
    end1 = max([handles.sched(end).scan(1).sta.endmjd]);
    end2 = max([handles.sched(end).scan(2).sta.endmjd]);
    end3 = max([handles.sched(end).scan(3).sta.endmjd]);
    handles.endsched = max([end1,end2,end3]);
end
if length([handles.sched(end).scan])==4
    end1 = max([handles.sched(end).scan(1).sta.endmjd]);
    end2 = max([handles.sched(end).scan(2).sta.endmjd]);
    end3 = max([handles.sched(end).scan(3).sta.endmjd]);
    end4 = max([handles.sched(end).scan(4).sta.endmjd]);
    handles.endsched = max([end1,end2,end3,end4]);
end

handles.endmjd = handles.endsched;

[year, month, day, hour, minute, second] = tymdhms(handles.startmjd);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.txt_start.String = sprintf('start: mjd: %.4f date: %s',handles.startmjd,datestr);

[year, month, day, hour, minute, second] = tymdhms(handles.endsched);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.txt_end.String = sprintf('end: mjd: %.4f date: %s',handles.endmjd,datestr);

handles.highlight{1} = false(size(handles.obs));
handles.highlight{2} = false(size(handles.scans));

%% plot station network
axes(handles.axes_overview_station)
axesm('MapProjection','robinson','Grid','on','MeridianLabel','off',...
    'ParallelLabel','off','PLineLocation',30,'MLineLocation',60,...
    'Frame', 'on','FontSize',8);
tightmap();
axis off
coast = load('coast');
geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor',[1 1 1],'EdgeColor',[.3 .3 .3]);
llh = [handles.station.llh];
lon = llh(1:3:end)*180/pi;
lat = llh(2:3:end)*180/pi;
id = {handles.station.po};
handles.plot_sta = geoshow(lat,lon,'DisplayType','point','marker','o','MarkerEdgeColor','k','MarkerFaceColor','g');
for i = 1:length(id)
    sta_id(i)=textm(lat(i),lon(i),[' ' id{i}]);
end
handles.plot_sta_id = sta_id;
counter = 1;
for i = 1:length(id)
    for j = i+1:length(id)
        sta_bas(counter)=geoshow([lat(i) lat(j)],[lon(i) lon(j)],'LineStyle','--','LineWidth',1,'Color','r');
        counter = counter+1;
    end
end
handles.plot_sta_bas = sta_bas;
zoom on

%% plot source Overview
axes(handles.axes_source_overview)
axesm('MapProjection','mollweid','Grid','on','MeridianLabel','off',...
    'ParallelLabel','off','PLineLocation',30,'MLineLocation',60,...
    'Frame', 'on','FontSize',8,'FFaceColor',[1 1 1]);
tightmap();
axis off

% create plots
p1 = line();
set(p1,'XData',[],'YData',[]);
p1Star = line();
set(p1Star,'XData',[],'YData',[]);
handles.plot_source_overview = [p1,p1Star];

pdef = line();
set(pdef,'XData',[],'YData',[]);
handles.plot_source_overview_defining = pdef;

p2 = line();
set(p2,'XData',[],'YData',[]);
p2Star = line();
set(p2Star,'XData',[],'YData',[]);
handles.plot_source_overview_highlight = [p2,p2Star];

handles.plot_source_overview_txt_nr = 0;
handles.plot_source_overview_txt = text();

%% plot skyplot
name = {handles.station.name};
handles.popup_station.String = name;
axes(handles.axes_sky_coverage)
polar([0:0.1:2*pi 2*pi],90*ones(1,64),'k');
hold on
h = polar([0:0.1:2*pi 2*pi],30*ones(1,64));
set(h,'Color',[0.8725    0.8725    0.8725])
h = polar([0:0.1:2*pi 2*pi],60*ones(1,64));
set(h,'Color',[0.8725    0.8725    0.8725])
h = polar([0:0.1:2*pi 2*pi],60*ones(1,64));
set(h,'Color',[0.8725    0.8725    0.8725])
hold off
h = findall(gca,'type','line');
delete(h(end));
h = findall(gca,'type','text');
delete(h([13,14]));
view([90 -90])

% create empty polar plots
pdef = line();
set(pdef,'XData',[],'YData',[])
set(pdef,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','Linestyle','none');
handles.plot_sky_coverage_defining = pdef;

p1 = line();
set(p1,'XData',[],'YData',[])
set(p1,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','Linestyle','none');
p1Star = line();
set(p1Star,'XData',[],'YData',[])
set(p1Star,'marker','p','MarkerEdgeColor','k','markerSize',10,'MarkerFaceColor','r','Linestyle','none')
handles.plot_sky_coverage = [p1,p1Star];

p2 = line();
set(p2,'XData',[],'YData',[])
set(p2,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','y','Linestyle','none')
p2Star = line();
set(p2Star,'XData',[],'YData',[])
set(p2Star,'marker','p','MarkerEdgeColor','k','markerSize',10,'MarkerFaceColor','y','Linestyle','none')
handles.plot_sky_coverage_highlight = [p2,p2Star];

handles.plot_sky_coverage_txt_nr = 0;
handles.plot_sky_coverage_txt = text();

%% update stats
handles = updateStats(hObject, eventdata, handles);
%% update plots
handles = updatePlots(hObject, eventdata, handles);
%% update handles structure
guidata(hObject, handles);


function handles = updateStats(hObject, eventdata, handles)
%% Fill source stat table
allSource = [handles.scans.srcid];
allScans = [handles.scans];
defining = [handles.source.defining];

names = {handles.source(allSource).name};
nsta = [allScans.nsta];
baselines = (nsta.*(nsta-1))/2;

[names_unique,ia,ic] = unique(names);
srcid_unique = allSource(ia);
defining = defining(srcid_unique);
star = [handles.source(srcid_unique).star];
star = star==1;
n = length(ia);
nscans = histcounts(ic,n);
nbaselines = zeros(n,1);
for i = 1:n
    nbaselines(i)=sum(baselines(ic==i));
end

handles.source_list.Data(1:length(names_unique),1)=names_unique;
handles.source_list.Data(1:length(names_unique),2)=num2cell(nscans);
handles.source_list.Data(1:length(names_unique),3)=num2cell(nbaselines);
handles.source_list.Data(1:length(names_unique),4)=num2cell(false);
handles.source_list.Data(1:length(names_unique),5)=num2cell(defining);
handles.source_list.Data(1:length(names_unique),6)=num2cell(star);

%% Fill station stat table
names = {handles.station.name};
n = length(names);

handles.uitable_stat_station.Data(1:n,1)=names;
staids = [handles.obs.staid];
scans = cell(n,1);
source = cell(n,1);
obs = cell(n,1);
STAR = cell(n,1);
ICRF2 = cell(n,1);
for i = 1:n
    bool = staids==i;
    scans(i)=num2cell(sum(bool));
    source{i} = sprintf('%.1f [%%]',(sum(~[handles.obs(bool).star]&~[handles.obs(bool).defining])/scans{i})*100);
    ICRF2{i} = sprintf('%.1f [%%]',(sum([handles.obs(bool).defining])/scans{i})*100);
    STAR{i} = sprintf('%.1f [%%]',(sum([handles.obs(bool).star])/scans{i})*100);
end
handles.uitable_stat_station.Data(1:n,2)=scans;
handles.uitable_stat_station.Data(1:n,3)=source;
handles.uitable_stat_station.Data(1:n,4)=STAR;
handles.uitable_stat_station.Data(1:n,5)=ICRF2;

scansta = [handles.scans.nsta];
for i = 2:n
    handles.uitable_stat_station_nscan.Data(i-1,1)={sprintf('Number of %2d-scan subnets:',i)};
    handles.uitable_stat_station_nscan.Data(i-1,2)={sum(scansta==i)};
end
handles.uitable_stat_station_nscan.ColumnWidth = {150,30};
%% Fill baseline table
name = {handles.station.po};
n = length(name);
handles.uitable_stat_baseline.Data = cell(n);
handles.uitable_stat_baseline.ColumnName = {name{:} 'SUM'};
handles.uitable_stat_baseline.RowName = name;
handles.uitable_stat_baseline.ColumnWidth = {35};

stationPerScan = {allScans.sta};
counter = zeros(n);
for i = 1:length(stationPerScan)
    cr = nchoosek([stationPerScan{i}.staid],2);
    ind = sub2ind([n,n],cr(:,1),cr(:,2));
    counter(ind)=counter(ind)+1;
end
counter = counter+counter';
counter = triu(counter);
s1 = sum(counter,1);
s2 = sum(counter,2);
s = s1'+s2;
counter = [counter s];
handles.uitable_stat_baseline.Data = counter;
handles.text_stat_baseline_sum.String = sprintf('Total Nr. of Baselines: %d',sum(s)/2);

%% statistics general
handles.text_stat_general_startmjd.String = sprintf('%.5f',handles.startsched);
[year, month, day, hour, minute, second] = tymdhms(handles.startsched);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.text_stat_general_startdate.String = sprintf('%s',datestr);

handles.text_stat_general_endmjd.String = sprintf('%.5f',handles.endsched);
[year, month, day, hour, minute, second] = tymdhms(handles.endsched);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.text_stat_general_enddate.String = sprintf('%s',datestr);

handles.text_stat_general_totalTime.String = sprintf('%.3f [h]',(handles.endsched-handles.startsched)*24);

handles.text_stat_general_nsta.String = sprintf('%d',length(handles.station));
handles.text_stat_general_nscans.String = sprintf('%d',length(handles.scans));
handles.text_stat_general_nbaselines.String = sprintf('%d',sum(s)/2);
handles.text_stat_general_nSource.String = sprintf('%d',length(handles.source));
sched = handles.sched;
num = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        thisScan = sched(isched).scan(iscan);
        idsource = sched(isched).scan(iscan).srcid;
        for sta1 = 1:thisScan.nsta-1;
            id1 = sched(isched).scan(iscan).sta(sta1).staid;
            for sta2 = sta1+1:thisScan.nsta
                id2 = sched(isched).scan(iscan).sta(sta2).staid;
                num = num+1;
                blsrc(num,1) = id1;
                blsrc(num,2) = id2;
                blsrc(num,3)  = idsource;
            end
        end
    end
end
tmp = blsrc(:,1:2);
blsrc(:,1:2) = sort(tmp,2);
[C,ia,ic] = unique(blsrc,'rows');

xtmp1 = length(ic);
xtmp2 = 0;
for i=1:max(ic)
    xtmp2 = xtmp2 + numel(ic(ic==i))^2;
end
num = length(ia);
M = mode(ic);
most = C(M,:);
sta1 = most(1);
sta2 = most(2);
src  = most(3);
maxva = numel(ic(ic==M));
minva = 10000;
for i=1:max(ic)
    tmp = numel(ic(ic==i));
    if tmp==1
        minva = 1;
        break
    end
    if tmp<minva
        minva=tmp;
    end
end

handles.text_stat_general_dBaseline.String = sprintf('Average number of obs. per baseline...');
handles.text_stat_general_dBaseline2.String = sprintf('...per source(normalized by up-time) =   %5.1f\n', (xtmp1/num));
handles.text_stat_general_dBaseline_min.String = sprintf('%6d', minva);
handles.text_stat_general_dBaseline_max.String = sprintf('%6d', maxva);
handles.text_stat_general_dBaseline_bline.String = sprintf('(Baseline %s-%s on %s)', handles.station(sta1).po, handles.station(sta2).po, strtrim(handles.source(src).name));
handles.text_stat_general_dBaseline_rms.String = sprintf('%.1f', sqrt(xtmp2/num));

if any([handles.obs.star])
    handles.pushbutton_disp_star_stats.Visible = 'on';
end

function handles = updatePlots(hObject, eventdata, handles)
%% update sky coverage
axes(handles.axes_sky_coverage)

id = get(handles.popup_station,'Value');
bool1 = [handles.obs.staid]==id;
bool2 = handles.startmjd<=[handles.obs.endmjd] & handles.endmjd>=[handles.obs.startmjd];

boolStar = [handles.obs.star];
boolDef = [handles.obs.defining];

bool = bool1 & bool2 & handles.select{1} & ~boolStar & ~boolDef;
boolDef = bool1 & bool2 & handles.select{1} & boolDef;
boolStar = bool1 & bool2 & handles.select{1} & boolStar & ~boolDef;

az = [handles.obs(bool).az];
zd = 90-[handles.obs(bool).el]*180/pi;

azDef = [handles.obs(boolDef).az];
zdDef = 90-[handles.obs(boolDef).el]*180/pi;

azStar = [handles.obs(boolStar).az];
zdStar = 90-[handles.obs(boolStar).el]*180/pi;

set(handles.plot_sky_coverage_defining,'XData',zdDef.*cos(azDef),'YData',zdDef.*sin(azDef))
set(handles.plot_sky_coverage(1),'XData',zd.*cos(az),'YData',zd.*sin(az))
set(handles.plot_sky_coverage(2),'XData',zdStar.*cos(azStar),'YData',zdStar.*sin(azStar))

%% update text
boolAll = bool | boolStar | boolDef;
switch handles.plot_sky_coverage_txt_nr
    case 0
        delete(handles.plot_sky_coverage_txt);
        handles.plot_sky_coverage_txt = text();
    case 1
        delete(handles.plot_sky_coverage_txt);
        txt = {handles.obs(boolAll).srcname};
        az = [handles.obs(boolAll).az];
        zd = 90-[handles.obs(boolAll).el]*180/pi;
        t = text(zd.*cos(az),zd.*sin(az),txt);
        handles.plot_sky_coverage_txt = t;
    case 2
        delete(handles.plot_sky_coverage_txt);
        txt = {handles.obs(boolAll).nsta};
        az = [handles.obs(boolAll).az];
        zd = 90-[handles.obs(boolAll).el]*180/pi;
        t = text(zd.*cos(az),zd.*sin(az),txt);
        handles.plot_sky_coverage_txt = t;
    case 3
        delete(handles.plot_sky_coverage_txt);
        txt = {handles.obs(boolAll).nobs};
        az = [handles.obs(boolAll).az];
        zd = 90-[handles.obs(boolAll).el]*180/pi;
        t = text(zd.*cos(az),zd.*sin(az),txt);
        handles.plot_sky_coverage_txt = t;
end

%% update highlight{1}
boolHighlight = (bool | boolDef) & handles.highlight{1};
boolHighlightStar = boolStar & handles.highlight{1};
az = [handles.obs(boolHighlight).az];
zd = 90-[handles.obs(boolHighlight).el]*180/pi;
azStar = [handles.obs(boolHighlightStar).az];
zdStar = 90-[handles.obs(boolHighlightStar).el]*180/pi;

set(handles.plot_sky_coverage_highlight(1),'XData',zd.*cos(az),'YData',zd.*sin(az))
set(handles.plot_sky_coverage_highlight(2),'XData',zdStar.*cos(azStar),'YData',zdStar.*sin(azStar))


%% Legende
l = legend([handles.plot_sky_coverage handles.plot_sky_coverage_defining handles.plot_sky_coverage_highlight],{'source','STAR-source','ICRF2 defining source','highlight source','highlight STAR source'});
set(l,'position',[.05 .15 .3 .1])


%% source_overview_plot
if handles.plotdown == 2;
    bool2 = handles.startmjd<=[handles.scans.endmjd] & handles.endmjd>=[handles.scans.startmjd];
    boolStar = [handles.scans.star];
    boolDef = [handles.scans.defining];

    bool =  handles.select{2} & bool2 & ~boolStar & ~boolDef;
    boolDef =  handles.select{2} & bool2 & boolDef;
    boolStar =  handles.select{2} & bool2 & boolStar & ~boolDef;
    boolHighlight = (bool | boolDef) & handles.highlight{2};
    boolHighlightStar = boolStar & handles.highlight{2};

    axes(handles.axes_source_overview);
    delete(handles.plot_source_overview);
    delete(handles.plot_source_overview_highlight);
    delete(handles.plot_source_overview_defining);
    delete(handles.plot_source_overview_txt);

    p1 = line();
    set(p1,'XData',[],'YData',[]);
    p1Star = line();
    set(p1Star,'XData',[],'YData',[]);
    handles.plot_source_overview = [p1,p1Star];

    pdef = line();
    set(pdef,'XData',[],'YData',[]);
    handles.plot_source_overview_defining = pdef;

    p2 = line();
    set(p2,'XData',[],'YData',[]);
    p2Star = line();
    set(p2Star,'XData',[],'YData',[]);
    handles.plot_source_overview_highlight = [p2,p2Star];

    % source
    srcid = [handles.scans(bool).srcid];
    srcid = unique(srcid);
    ra = [handles.source(srcid).ra]*180/pi;
    de = [handles.source(srcid).de]*180/pi;
    g = geoshow(de,ra,'DisplayType','Point','marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','Linestyle','none');
    if ~isempty(g)
        handles.plot_source_overview(1) = g;
    else
        set(handles.plot_source_overview(1),'XData',[],'YData',[]);
    end

    % STAR-source
    srcid = [handles.scans(boolStar).srcid];
    srcid = unique(srcid);
    ra = [handles.source(srcid).ra]*180/pi;
    de = [handles.source(srcid).de]*180/pi;
    g = geoshow(de,ra,'DisplayType','Point','marker','p','markerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r','Linestyle','none');
    if ~isempty(g)
        handles.plot_source_overview(2) = g;
    else
        set(handles.plot_source_overview(2),'XData',[],'YData',[]);
    end

    % defining
    srcid = [handles.scans(boolDef).srcid];
    srcid = unique(srcid);
    ra = [handles.source(srcid).ra]*180/pi;
    de = [handles.source(srcid).de]*180/pi;
    g = geoshow(de,ra,'DisplayType','Point','marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','Linestyle','none');
    if ~isempty(g)
        handles.plot_source_overview_defining(1) = g;
    else
        set(handles.plot_source_overview_defining(1),'XData',[],'YData',[]);
    end

    % highlight source
    srcid = [handles.scans(boolHighlight).srcid];
    srcid = unique(srcid);
    ra = [handles.source(srcid).ra]*180/pi;
    de = [handles.source(srcid).de]*180/pi;
    g = geoshow(de,ra,'DisplayType','Point','marker','o','MarkerEdgeColor','k','MarkerFaceColor','y','Linestyle','none');
    if ~isempty(g)
        handles.plot_source_overview_highlight(1) = g;
    else
        set(handles.plot_source_overview_highlight(1),'XData',[],'YData',[]);
    end

    % highlight STAR-source
    srcid = [handles.scans(boolHighlightStar).srcid];
    srcid = unique(srcid);
    ra = [handles.source(srcid).ra]*180/pi;
    de = [handles.source(srcid).de]*180/pi;
    g = geoshow(de,ra,'DisplayType','Point','marker','p','markerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Linestyle','none');
    if ~isempty(g)
        handles.plot_source_overview_highlight(2) = g;
    else
        set(handles.plot_source_overview_highlight(2),'XData',[],'YData',[]);
    end

    % text
    boolAll = bool | boolStar | boolDef;
    switch handles.plot_source_overview_txt_nr
        case 0
            delete(handles.plot_source_overview_txt);
            handles.plot_source_overview_txt = text();
        case 1
            delete(handles.plot_source_overview_txt);
            thisScans = handles.scans(boolAll);
            thisScansSrcId = [thisScans.srcid];
            [thisScansSrcId] = unique(thisScansSrcId);
            txt = {handles.source(thisScansSrcId).name};
            ra = [handles.source(thisScansSrcId).ra]*180/pi;
            de = [handles.source(thisScansSrcId).de]*180/pi;
            t = textm(de,ra,txt);
            handles.plot_source_overview_txt = t;
        case 2
            delete(handles.plot_source_overview_txt);
            thisScans = handles.scans(boolAll);
            thisScansSrcId = [thisScans.srcid];
            [thisScansSrcId,ic,ia] = unique(thisScansSrcId);
            txt = cell(length(thisScansSrcId),1);
            for i = 1:length(ic)
                idx = find(ia==i);
                txt{i}=sprintf('-%d',length(idx));
            end
            ra = [handles.source(thisScansSrcId).ra]*180/pi;
            de = [handles.source(thisScansSrcId).de]*180/pi;
            t = textm(de,ra,txt);
            handles.plot_source_overview_txt = t;
        case 3
            delete(handles.plot_source_overview_txt);
            thisScans = handles.scans(boolAll);
            thisScansSrcId = [thisScans.srcid];
            [thisScansSrcId,ic,ia] = unique(thisScansSrcId);
            txt = cell(length(thisScansSrcId),1);
            for i = 1:length(ic)
                idx = find(ia==i);
                txt{i}='';
                for j = 1:length(idx)
                    txt{i}=sprintf('%s,%d',txt{i},thisScans(idx(j)).nsta);
                end
                txt{i}(1)=' ';
            end
            ra = [handles.source(thisScansSrcId).ra]*180/pi;
            de = [handles.source(thisScansSrcId).de]*180/pi;
            t = textm(de,ra,txt);
            handles.plot_source_overview_txt = t;
        case 4
            delete(handles.plot_source_overview_txt);
            thisScans = handles.scans(boolAll);
            thisScansSrcId = [thisScans.srcid];
            [thisScansSrcId,ic,ia] = unique(thisScansSrcId);
            txt = cell(length(thisScansSrcId),1);
            for i = 1:length(ic)
                idx = find(ia==i);
                txt{i}='';
                for j = 1:length(idx)
                    txt{i}=sprintf('%s,%d',txt{i},thisScans(idx(j)).nobs);
                end
                txt{i}(1)=' ';
            end
            ra = [handles.source(thisScansSrcId).ra]*180/pi;
            de = [handles.source(thisScansSrcId).de]*180/pi;
            t = textm(de,ra,txt);
            handles.plot_source_overview_txt = t;
    end

end
%% --- Outputs from this function are returned to the command line.
function varargout = sched_analyser_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on selection change in popup_station.
function popup_station_Callback(hObject, eventdata, handles)

% hObject    handle to popup_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popup_station contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_station

%% --- Executes during object creation, after setting all properties.
function popup_station_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on slider movement.
function slider_startmjd_Callback(hObject, eventdata, handles)
% hObject    handle to slider_startmjd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(handles.slider_startmjd,'Value');
mjd = handles.startsched+(handles.endsched-handles.startsched)*value;
if ~get(handles.button_fixed_length,'Value');
    if mjd>handles.endmjd
        mjd = handles.endmjd;
        value = (mjd-handles.startsched)/(handles.endsched-handles.startsched);
        hObject.Value = value;
    end
    handles.startmjd = mjd;
else
    handles.startmjd = mjd;
    minutes = str2double(get(handles.edittxt_minutes,'String'));
    handles.endmjd = min([handles.startmjd+minutes/1440 handles.endsched]);
    value = (handles.endmjd-handles.startsched)/(handles.endsched-handles.startsched);
    handles.slider_endmjd.Value = value;
    [year, month, day, hour, minute, second] = tymdhms(handles.endmjd);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    handles.txt_end.String = sprintf('end: mjd: %.4f date: %s',handles.endmjd,datestr);
end

[year, month, day, hour, minute, second] = tymdhms(handles.startmjd);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.txt_start.String = sprintf('start: mjd: %.4f date: %s',handles.startmjd,datestr);

handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%% --- Executes during object creation, after setting all properties.
function slider_startmjd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_startmjd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% --- Executes on slider movement.
function slider_endmjd_Callback(hObject, eventdata, handles)
% hObject    handle to slider_endmjd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
mjd = handles.startsched+(handles.endsched-handles.startsched)*value;

if ~get(handles.button_fixed_length,'Value');
    if mjd<handles.startmjd
        mjd = handles.startmjd;
        value = (mjd-handles.startsched)/(handles.endsched-handles.startsched);
        hObject.Value = value;
    end
    handles.endmjd = mjd;
else
    handles.endmjd = mjd;
    minutes = str2double(get(handles.edittxt_minutes,'String'));
    handles.startmjd = max([handles.endmjd-minutes/1440 handles.startsched]);
    value = (handles.startmjd-handles.startsched)/(handles.endsched-handles.startsched);
    handles.slider_startmjd.Value = value;
    [year, month, day, hour, minute, second] = tymdhms(handles.startmjd);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    handles.txt_start.String = sprintf('start: mjd: %.4f date: %s',handles.startmjd,datestr);
end

[year, month, day, hour, minute, second] = tymdhms(handles.endmjd);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.txt_end.String = sprintf('end: mjd: %.4f date: %s',handles.endmjd,datestr);

handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%% --- Executes during object creation, after setting all properties.
function slider_endmjd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_endmjd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Value = 1;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% --- Executes on button press in button_fixed_length.
function button_fixed_length_Callback(hObject, eventdata, handles)
% hObject    handle to button_fixed_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if value
    handles.txt_min.Enable = 'on';
    handles.edittxt_minutes.Enable = 'on';
%     handles.txt_end.String = '';
    minutes = str2double(get(handles.edittxt_minutes,'String'));
    wholeTime = (handles.endsched-handles.startsched)*1440;
    oneStep = minutes/wholeTime/2;
    handles.slider_endmjd.SliderStep = [oneStep oneStep];
    handles.slider_startmjd.SliderStep = [oneStep oneStep];
    handles.endmjd = min([handles.startmjd+minutes/1440 handles.endsched]);
    value = (handles.endmjd-handles.startsched)/(handles.endsched-handles.startsched);
    handles.slider_endmjd.Value = value;
    [year, month, day, hour, minute, second] = tymdhms(handles.endmjd);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    handles.txt_end.String = sprintf('end: mjd: %.4f date: %s',handles.endmjd,datestr);
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.txt_min.Enable = 'off';
    handles.edittxt_minutes.Enable = 'off';
    handles.slider_endmjd.SliderStep = [0.05 0.05];
    handles.slider_startmjd.SliderStep = [0.05 0.05];
end
% Hint: get(hObject,'Value') returns toggle state of button_fixed_length

%% --- Executes on editing the text in edittxt_minutes
function edittxt_minutes_Callback(hObject, eventdata, handles)
% hObject    handle to edittxt_minutes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject,'String'));
if ~isnan(value) && value>0
    wholeTime = (handles.endsched-handles.startsched)*1440;
    oneStep = value/wholeTime/2;
    handles.slider_endmjd.SliderStep = [oneStep oneStep];
    handles.slider_startmjd.SliderStep = [oneStep oneStep];
    slider_startmjd_Callback(hObject, eventdata, handles)
else
    warning('numeric input required!')
end
% Hints: get(hObject,'String') returns contents of edittxt_minutes as text
%        str2double(get(hObject,'String')) returns contents of edittxt_minutes as a double


%% --- Executes during object creation, after setting all properties.
function edittxt_minutes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittxt_minutes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in show_baseline.
function show_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to show_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    for i = 1:length(handles.plot_sta_bas)
        handles.plot_sta_bas(i).Visible = 'off';
    end
else
    for i = 1:length(handles.plot_sta_bas)
        handles.plot_sta_bas(i).Visible = 'on';
    end
end
% Hint: get(hObject,'Value') returns toggle state of show_baseline


%% --- Executes on button press in show_station_id.
function show_station_id_Callback(hObject, eventdata, handles)
% hObject    handle to show_station_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    for i = 1:length(handles.plot_sta_id)
        handles.plot_sta_id(i).Visible = 'off';
    end
else
    for i = 1:length(handles.plot_sta_id)
        handles.plot_sta_id(i).Visible = 'on';
    end
end
% Hint: get(hObject,'Value') returns toggle state of show_station_id


%% --- Executes on button press in button_sky_cov_src_name.
function button_sky_cov_src_name_Callback(hObject, eventdata, handles)
% hObject    handle to button_sky_cov_src_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_sky_coverage_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_sky_coverage_txt_nr=1;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.radiobutton_nsta,'Value',0);
    set(handles.radiobutton_nobs,'Value',0);
    guidata(hObject, handles);
end

%% --- Executes on button press in radiobutton_nsta.
function radiobutton_nsta_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_nsta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_sky_coverage_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_sky_coverage_txt_nr=2;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.button_sky_cov_src_name,'Value',0);
    set(handles.radiobutton_nobs,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_nsta


%% --- Executes on button press in radiobutton_nobs.
function radiobutton_nobs_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_nobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_sky_coverage_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_sky_coverage_txt_nr=3;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.button_sky_cov_src_name,'Value',0);
    set(handles.radiobutton_nsta,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_nobs

%% --- Executes on button press in source_overview_name.
function source_overview_name_Callback(hObject, eventdata, handles)
% hObject    handle to source_overview_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = get(hObject,'Value');
if ~value
    handles.plot_source_overview_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_source_overview_txt_nr=1;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.radiobutton_source_overview_scans,'Value',0);
    set(handles.source_overview_obs,'Value',0);
    set(handles.radiobutton_source_overview_nsta,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of source_overview_name


%% --- Executes on button press in radiobutton_source_overview_scans.
function radiobutton_source_overview_scans_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_source_overview_scans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_source_overview_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_source_overview_txt_nr=2;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.source_overview_name,'Value',0);
    set(handles.source_overview_obs,'Value',0);
    set(handles.radiobutton_source_overview_nsta,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_source_overview_scans

%% --- Executes on button press in radiobutton_source_overview_nsta.
function radiobutton_source_overview_nsta_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_source_overview_nsta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_source_overview_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_source_overview_txt_nr=3;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.radiobutton_source_overview_scans,'Value',0);
    set(handles.source_overview_name,'Value',0);
    set(handles.source_overview_obs,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton_source_overview_nsta

%% --- Executes on button press in source_overview_obs.
function source_overview_obs_Callback(hObject, eventdata, handles)
% hObject    handle to source_overview_obs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
if ~value
    handles.plot_source_overview_txt_nr=0;
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
else
    handles.plot_source_overview_txt_nr=4;
    handles = updatePlots(hObject, eventdata, handles);
    set(handles.radiobutton_source_overview_scans,'Value',0);
    set(handles.source_overview_name,'Value',0);
    set(handles.radiobutton_source_overview_nsta,'Value',0);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of source_overview_obs

%% --- Executes on selection change in pupup_sortby.
function pupup_sortby_Callback(hObject, eventdata, handles)
% hObject    handle to pupup_sortby (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
value = contents{get(hObject,'Value')};
if strcmp(value,'name')
    d = handles.source_list.Data(:,1);
    [~,idx] = sort(d);
    handles.source_list.Data = handles.source_list.Data(idx,:);
elseif strcmp(value,'scans')
    d = handles.source_list.Data(:,2);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.source_list.Data = handles.source_list.Data(idx,:);
elseif strcmp(value,'obs')
    d = handles.source_list.Data(:,3);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.source_list.Data = handles.source_list.Data(idx,:);
elseif strcmp(value,'highlight')
    d = handles.source_list.Data(:,4);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.source_list.Data = handles.source_list.Data(idx,:);
elseif strcmp(value,'ICRF2 defining')
    d = handles.source_list.Data(:,5);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.source_list.Data = handles.source_list.Data(idx,:);
elseif strcmp(value,'STAR-source')
    d = handles.source_list.Data(:,6);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.source_list.Data = handles.source_list.Data(idx,:);
end




% Hints: contents = cellstr(get(hObject,'String')) returns pupup_sortby contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pupup_sortby


%% --- Executes during object creation, after setting all properties.
function pupup_sortby_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupup_sortby (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pushbutton_highlight.
function pushbutton_highlight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_highlight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.source_list.Data(:,4);
h = cell2mat(h);
hid = find(h);
hname = handles.source_list.Data(hid,1);

srcid = [handles.obs.srcid];
names = {handles.source(srcid).name};
boolHighlight = false(size(srcid));

srcid2 = [handles.scans.srcid];
names2 = {handles.source(srcid2).name};
boolHighlightSrc = false(size(handles.scans));
for i = 1:length(hname)
    highlightID = strcmp(hname(i),names);
    boolHighlight = boolHighlight | highlightID;
    highlightID2 = strcmp(hname(i),names2);
    boolHighlightSrc = boolHighlightSrc | highlightID2;
end
handles.highlight{1} = boolHighlight;
handles.highlight{2} = boolHighlightSrc;
handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_select_plot_down.
function popupmenu_select_plot_down_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_plot_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
value = contents{get(hObject,'Value')};

if strcmp(value,'plot: station overview')
    handles.uipanel_station_overview.Visible = 'on';
    handles.uipanel_source_overview.Visible = 'off';
    handles.plotdown = 1;
elseif strcmp(value,'plot: source overview')
    handles.uipanel_station_overview.Visible = 'off';
    handles.uipanel_source_overview.Visible = 'on';
    handles.plotdown = 2;
    handles = updatePlots(hObject, eventdata, handles);
end
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_select_plot_down contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_select_plot_down


%% --- Executes during object creation, after setting all properties.
function popupmenu_select_plot_down_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_select_plot_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_open_skd_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open_skd_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_timespan_Callback(hObject, eventdata, handles)
% hObject    handle to menu_timespan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'start: "yyyy-mm-dd hh:mm:ss"','end: "yyyy-mm-dd hh:mm:ss"'};
dlg_title = 'timespan';
num_lines = [1 30];

[year, month, day, hour, minute, second] = tymdhms(handles.startsched);
[datestr1] = tdatestr(year, month, day, hour, minute, second);

[year, month, day, hour, minute, second] = tymdhms(handles.endsched);
[datestr2] = tdatestr(year, month, day, hour, minute, second);

defaultans = {datestr1,datestr2};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
if ~isempty(answer)
    startT = textscan(answer{1},'%f-%f-%f %f:%f:%f');
    startT  = date2mjd(cell2mat(startT));
    endT = textscan(answer{2},'%f-%f-%f %f:%f:%f');
    endT = date2mjd(cell2mat(endT));

    handles.startmjd = startT;
    handles.endmjd = endT;

    value = max([0,(handles.startmjd-handles.startsched)/(handles.endsched-handles.startsched)]);
    handles.slider_startmjd.Value = value;
    handles.txt_start.String = sprintf('start: mjd: %.4f date: %s',handles.startmjd,answer{1});

    value = min([1 (handles.endmjd-handles.startsched)/(handles.endsched-handles.startsched)]);
    handles.slider_endmjd.Value = value;
    handles.txt_end.String = sprintf('end: mjd: %.4f date: %s',handles.endmjd,answer{2});

    set(handles.button_fixed_length,'Value',0);
    handles.txt_min.Enable = 'off';
    handles.edittxt_minutes.Enable = 'off';

    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);
end

%% --- Executes on selection change in popupmenu_statistics.
function popupmenu_statistics_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
val = contents{get(hObject,'Value')};
if strcmp(val,'statistic: general informations')
    handles.uipanel_stats_source.Visible = 'off';
    handles.uipanel_stat_station.Visible = 'off';
    handles.uipanel_stat_baseline.Visible = 'off';
    handles.uipanel_stat_general.Visible = 'on';
elseif strcmp(val,'statistic: source based')
    handles.uipanel_stats_source.Visible = 'on';
    handles.uipanel_stat_station.Visible = 'off';
    handles.uipanel_stat_baseline.Visible = 'off';
    handles.uipanel_stat_general.Visible = 'off';
elseif strcmp(val,'statistic: baseline based')
    handles.uipanel_stats_source.Visible = 'off';
    handles.uipanel_stat_station.Visible = 'off';
    handles.uipanel_stat_baseline.Visible = 'on';
    handles.uipanel_stat_general.Visible = 'off';
elseif strcmp(val,'statistic: station based')
    handles.uipanel_stat_station.Visible = 'on';
    handles.uipanel_stats_source.Visible = 'off';
    handles.uipanel_stat_baseline.Visible = 'off';
    handles.uipanel_stat_general.Visible = 'off';
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_statistics contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_statistics


%% --- Executes during object creation, after setting all properties.
function popupmenu_statistics_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pushbutton_sky_coverage_select_plot.
function pushbutton_sky_coverage_select_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sky_coverage_select_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Selection,ok] = listdlg('PromptString','Select plots: (STRG+Click for multiple)','ListSize',[190 100],...
                'OKString','plot','InitialValue',[1 2 3 4 5],...
                'ListString',{'source','STAR-source','ICRF2 defining source','highlight source','highlight STAR source'});
if ok
    sum = false(1,length(handles.obs));
    highB = handles.highlight{1};
    star = [handles.obs.star];
    highStar = highB & star;
    high = highB & ~star;
    def = [handles.obs.defining];
    src = ~high&~star&~highStar&~def;

    if any(Selection==1)
        sum = sum | src;
    end
    if any(Selection==2)
        sum = sum | star;
    end
    if any(Selection==3)
        sum = sum | def;
    end
    if any(Selection==4)
        sum = sum | high;
    end
    if any(Selection==5)
        sum = sum | highStar;
    end
    handles.select{1} = sum;
end
handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);


%% --- Executes on button press in pushbutton_source_overview_select_plot.
function pushbutton_source_overview_select_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_source_overview_select_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Selection,ok] = listdlg('PromptString','Select plots: (STRG+Click for multiple)','ListSize',[190 100],...
                'OKString','plot','InitialValue',[1 2 3 4 5],...
                'ListString',{'source','STAR-source','ICRF2 defining source','highlight source','highlight STAR source'});
if ok
    sum = false(1,length(handles.scans));
    highB = handles.highlight{2};
    star = [handles.scans.star];
    highStar = highB & star;
    high = highB & ~star;
    def = [handles.scans.defining];
    src = ~high&~star&~highStar&~def;
    if any(Selection==1)
        sum = sum | src;
    end
    if any(Selection==2)
        sum = sum | star;
    end
    if any(Selection==3)
        sum = sum | def;
    end
    if any(Selection==4)
        sum = sum | high;
    end
    if any(Selection==5)
        sum = sum | highStar;
    end
    handles.select{2} = sum;
end
handles = updatePlots(hObject, eventdata, handles);
guidata(hObject, handles);


%% --- Executes on selection change in popupmenu_stat_station_sort.
function popupmenu_stat_station_sort_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_stat_station_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
value = contents{get(hObject,'Value')};
if strcmp(value,'name')
    d = handles.uitable_stat_station.Data(:,1);
    [~,idx] = sort(d);
    handles.uitable_stat_station.Data = handles.uitable_stat_station.Data(idx,:);
elseif strcmp(value,'scans')
    d = handles.uitable_stat_station.Data(:,2);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.uitable_stat_station.Data = handles.uitable_stat_station.Data(idx,:);
elseif strcmp(value,'source')
    d = handles.uitable_stat_station.Data(:,3);
    d = strcat(d,{' '});
    C = textscan([d{:}],'%f %s');
    d = C{1};
    [~,idx] = sort(d,'descend');
    handles.uitable_stat_station.Data = handles.uitable_stat_station.Data(idx,:);
elseif strcmp(value,'STAR')
    d = handles.uitable_stat_station.Data(:,4);
    d = strcat(d,{' '});
    C = textscan([d{:}],'%f %s');
    d = C{1};
    [~,idx] = sort(d,'descend');
    handles.uitable_stat_station.Data = handles.uitable_stat_station.Data(idx,:);
elseif strcmp(value,'ICRF2 defining')
    d = handles.uitable_stat_station.Data(:,5);
    d = strcat(d,{' '});
    C = textscan([d{:}],'%f %s');
    d = C{1};
    [~,idx] = sort(d,'descend');
    handles.uitable_stat_station.Data = handles.uitable_stat_station.Data(idx,:);
end



%% --- Executes during object creation, after setting all properties.
function popupmenu_stat_station_sort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stat_station_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --------------------------------------------------------------------
function menu_browse_Callback(hObject, eventdata, handles)
% hObject    handle to menu_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('../DATA/SCHED/*.skd','Selekt skd file');
if ~isequal(FileName,0)
    sched_analyser('file',[PathName,FileName])
end

%% --------------------------------------------------------------------
function menu_level5_Callback(hObject, eventdata, handles)
% hObject    handle to menu_level5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir('../DATA/');
if folder_name ~=0
    sched_analyser('folder',folder_name);
end

%% --- Executes on button press in pushbutton_load_star_catalog.
function pushbutton_load_star_catalog_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to pushbutton_load_star_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(varargin)
    [FileName,PathName,FilterIndex] = uigetfile('../CATALOGS/*.*','Selekt STAR catalog file');
else
    FileName = 'source_star.cat';
    PathName = '../CATALOGS/';
end
if ~isequal(FileName,0)
    fid = fopen([PathName FileName]);
    C = textscan(fid,'%s %s %*[^\n]','CommentStyle','*');
    fclose(fid);
    for i = 1:length(C{1})
        if strcmp(C{2},'$')
            starNames(i) = C{1}(i);
        else
            starNames(i) = C{2}(i);
        end
    end
    sourceNames = {handles.source.name};

    for i = 1:length(starNames)
        bool = strcmp(starNames(i),sourceNames);
        if any(bool)
            handles.source(bool).star = 1;
        end
    end

    starId = find([handles.source.star]);
    bool = false(1,length(handles.obs));
    for i = 1:length(starId)
        bool = bool | [handles.obs.srcid]==starId(i);
    end
    for i = 1:length(handles.obs)
        if bool(i)
            handles.obs(i).star = 1;
        end
    end

    bool = false(1,length(handles.scans));
    for i = 1:length(starId)
        bool = bool | [handles.scans.srcid]==starId(i);
    end
    for i = 1:length(handles.scans)
        if bool(i)
            handles.scans(i).star = 1;
        end
    end
    handles.pushbutton_load_star_catalog.Visible = 'off';
    handles.pushbutton_browse_for_star_catalog.Visible = 'off';
    handles = updateStats(hObject, eventdata, handles);
    handles = updatePlots(hObject, eventdata, handles);
    guidata(hObject, handles);

end


% --- Executes on button press in pushbutton_browse_for_star_catalog.
function pushbutton_browse_for_star_catalog_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse_for_star_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton_load_star_catalog_Callback(hObject, eventdata, handles, 1)


% --- Executes on button press in pushbutton_disp_star_stats.
function pushbutton_disp_star_stats_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_disp_star_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp_stat_starmode( handles.source,handles.sched )
