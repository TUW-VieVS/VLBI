% called by L3plot.m

function fig = plot_scan_cable(antenna, scan)


nSta = numel(antenna);
nScan = numel(scan);

antNames = cell(nSta,1);
for k = 1:nSta
    antNames{k} = strtrim(antenna(k).name);
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
cableData = cell(nSta,1);
cableCalData = cell(nSta,1);
cableCDMSData = cell(nSta,1);

for i = 1:nScan

    mjd = double(scan(i).mjd);

    for k = 1:nSta

        if k <= numel(scan(i).stat)

            st = scan(i).stat(k);

            if ~isempty(st.cab)
                mjdData{k}(end+1,1)  = mjd;
                cableData{k}(end+1,1) = double(st.cab);
                if isfield(st,'cab_cablecal')
                    cableCalData{k}(end+1,1) = double(st.cab_cablecal);
                else
                    cableCalData{k}(end+1,1) =NaN;
                end

                if isfield(st,'cab_CDMS')
                    cableCDMSData{k}(end+1,1) = double(st.cab_CDMS);
                else
                    cableCDMSData{k}(end+1,1) =NaN;
                end
            end
        end
    end
end

minmjd=min(cell2mat(mjdData));
maxmjd=max(cell2mat(mjdData));

%% ============================================================
%% FIGURE: Cable calibration
%% ============================================================

fig = figure('Color',[1 1 1], ...
             'Position',[50 00 1400 900]);

sgtitle('VLBI Scan Meteorological Data', ...
        'FontWeight','bold');

%% Temperature
tiledlayout("vertical");

for k = 1:nSta
    nexttile;

    if ~isempty(mjdData{k})
        plot(mjdData{k}, cableData{k}, ...
             'o', ...
             'MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
        hold on
        plot(mjdData{k}, cableCalData{k}, ...
             '.', ...
             'Color',colors(2,:), ...
             'LineWidth',1.0);
        hold on
        plot(mjdData{k}, cableCDMSData{k}, ...
             '.', ...
             'Color',colors(3,:), ...
             'LineWidth',1.0);

        styleAxes(minmjd, maxmjd);
        title(antNames(k))
        ylabel('[ns]');
    end
end



lgd = legend([{'VieVS applied'}; {'Cal-Cable'}; {'CDMS'}]);
lgd.Layout.Tile = 'east';  


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
    ytickformat('%.2f')
    ax = gca;

    ax.Color  = [1 1 1];
    ax.XColor = [0.2 0.2 0.2];
    ax.YColor = [0.2 0.2 0.2];

    ax.XAxis.Exponent = 0;
    ax.XAxis.Limits=[minmjd-0.01 maxmjd+0.01];

    set(gca,'XTickLabel',[;;;])

    box on;
    grid on;

end


end