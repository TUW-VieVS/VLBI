% called by L3plot.m

function fig = plot_scan_meteo(antenna, scan)


nSta = numel(antenna);
nScan = numel(scan);

antNames = cell(nSta,1);
gpt3 = cell(nSta,1);
for k = 1:nSta
    antNames{k} = strtrim(antenna(k).name);
    gpt3{k}.p = antenna(k).gpt3.p;
    gpt3{k}.T = antenna(k).gpt3.T;
    %gpt3{k}.Tm = antenna(k).gpt3.Tm - 273.15;
    gpt3{k}.e = antenna(k).gpt3.e;
end

%% Colors

colors1 = [55,126,184;
          255,127,0;
          77,175,74;
          247,129,191;
          166,86,40;
          152,78,163;
          153,153,153;
          228,26,28;
          222,222,0] ./ 255;

colors2 = ([55,126,184;
          255,127,0;
          77,175,74;
          247,129,191;
          166,86,40;
          152,78,163;
          153,153,153;
          228,26,28;
          222,222,0] +...
          [50 50 50
          0 50 50
          50 50 50
          0 50 50
          50 50 50
          50 50 50
          50 50 50
          10 50 50
          10 10 50])./ 255;

colors=[colors1;colors2]; % 18 ant max


%% Collect scan meteorological data

mjdData  = cell(nSta,1);
tempData = cell(nSta,1);
presData = cell(nSta,1);
eData    = cell(nSta,1);

for i = 1:nScan

    mjd = double(scan(i).mjd);

    for k = 1:nSta

        if k <= numel(scan(i).stat)

            st = scan(i).stat(k);

            if ~isempty(st.temp)
                mjdData{k}(end+1,1)  = mjd;
                tempData{k}(end+1,1) = double(st.temp);
                presData{k}(end+1,1) = double(st.pres);
                eData{k}(end+1,1)    = double(st.e);
            end

        end

    end

end

minmjd=min(cell2mat(mjdData));
maxmjd=max(cell2mat(mjdData));

%% ============================================================
%% FIGURE: Temperature / Pressure / Humidity e
%% ============================================================

fig = figure('Color',[1 1 1], ...
             'Position',[50 100 1400 900]);

sgtitle('VLBI Scan Meteorological Data', ...
        'FontWeight','bold');

%% Temperature
t = tiledlayout(3,1);

ax1 = nexttile;
%subplot(3,1,1);


for k = 1:nSta
    if ~isempty(mjdData{k})
        plot(mjdData{k}, tempData{k}, ...
             '.-', ...
             'Color',colors(k,:), ...
             'LineWidth',1.0);
        %hold on
        %plot(mjdData{k}(1), gpt3{k}.Tm, 'o','MarkerSize',5,'MarkerEdgeColor',colors(k,:),'MarkerFaceColor',[1 1 1]); % gpt3 mean temperature weighted with the water vapor in degrees Kelvin 
        hold on
        plot(mjdData{k}(1), gpt3{k}.T, 'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:)); % gpt3 temperature 
        hold on
    end
end
hold off
title('Temperature');
%xlabel('MJD');
ylabel('Temperature [°C]');

styleAxes(minmjd, maxmjd);

%% Pressure
ax2 = nexttile;
%subplot(3,1,2);


allPres = [];

for k = 1:nSta

    if ~isempty(mjdData{k})

        plot(mjdData{k}, presData{k}, ...
             '.-', ...
             'Color',colors(k,:), ...
             'LineWidth',1.0);
        hold on
        plot(mjdData{k}(1), gpt3{k}.p, 'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:)); % gpt3 p 

        allPres = [allPres; presData{k}(:)];
        hold on
    end

end
hold off
ymin = min(allPres);
ymax = max(allPres);

yrange = ymax - ymin;

padding = max(1.0,0.3*yrange);

ylim([ymin-padding ymax+padding])

title('Pressure');
%xlabel('MJD');
ylabel('Pressure [hPa]');

styleAxes(minmjd, maxmjd);

%% Humidity / e
ax3 = nexttile;
%subplot(3,1,3);


for k = 1:nSta
    if ~isempty(mjdData{k})
        p1{k}=plot(mjdData{k}, eData{k}, ...
             '.-', ...
             'Color',colors(k,:), ...
             'LineWidth',1.0);
        hold on                
        plot(mjdData{k}(1), gpt3{k}.e, 'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(k,:)); % gpt3 temperature 
        hold on
    end
end
hold off
title('Water vapor pressure e');
xlabel('MJD');
ylabel('e [hPa]');

lgd = legend([p1{:}],antNames);
lgd.Layout.Tile = 'east';  
styleAxes(minmjd, maxmjd);

%% Save

% Achsen-Toolbar vor Export ausblenden
axs = findall(gcf,'Type','axes');

for a = 1:numel(axs)
    if isprop(axs(a),'Toolbar')
        axs(a).Toolbar.Visible = 'off';
    end
end

% Aktuelle Figur speichern
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');


%% Local function

function styleAxes(minmjd, maxmjd)

    xtickformat('%.1f')
    ax = gca;

    ax.Color  = [1 1 1];
    ax.XColor = [0.2 0.2 0.2];
    ax.YColor = [0.2 0.2 0.2];

    ax.XAxis.Exponent = 0;
    ax.XAxis.Limits=[minmjd-0.01 maxmjd+0.01];

    box on;
    grid on;

end

end