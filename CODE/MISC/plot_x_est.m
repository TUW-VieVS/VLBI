% called by L3plot.m

function [fig1, fig2, fig3, fig4, fig5] = plot_x_est(x,opt_,sources)

nSta = numel(x.antenna);

antNames = cell(nSta,1);
for k = 1:nSta
    antNames{k} = strtrim(x.antenna(k).name);
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



%% MJD span from ZWD
allZwdMjd = [];
for i = 1:nSta
    allZwdMjd = [allZwdMjd; double(x.zwd(i).mjd(:))];
end

firstMJD = min(allZwdMjd);
lastMJD  = max(allZwdMjd);

%% ============================================================
% FIGURE 1: Atmosphere / ZWD
%% ============================================================

fig1 = figure('Color',[1 1 1], ...
              'Position',[50 100 1400 600]);



for i = 1:nSta

    z = x.zwd(i);

    errorbar(double(z.mjd), ...
             z.val, ...
             z.mx, ...
             '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)
    hold on
    p1{i}=plot(double(z.mjd), ...
             z.val, ...
             'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(i,:),'LineWidth',0.1);
    hold on
end
hold off
yline(0,'--');

title(sprintf('Zenith Wet Delay | MJD %.3f - %.3f', ...
      firstMJD,lastMJD), ...
      'FontWeight','bold');

xlabel('MJD');
ylabel('ZWD [cm]');
           
xtickformat('%.2f')
ax = gca;
ax.XAxis.Exponent = 0;

legend([p1{:}],antNames,'Location','eastoutside');
styleAxes();

%exportgraphics(fig1,[ses '_zwd.png'],'Resolution',800);
% exportgraphics(fig1,[ses '_zwd.pdf'], ...
               % 'ContentType','vector', ...
               % 'BackgroundColor','none');
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');

%% ============================================================
% FIGURE 2: EOP + Nutation
%% ============================================================

fig2 = figure('Color',[1 1 1], ...
              'Position',[50 100 1400 800]);

sgtitle(sprintf('EOP and Nutation | nobs=%d   nscans=%d   WRMS=%.4f cm', ...
        x.nobs,x.nscans,x.wrms), ...
        'FontWeight','bold');

%% x-pole
subplot(2,2,1);
yline(0,'--');
hold on
errorbar(double(x.xpol.mjd), ...
         x.xpol.val, ...
         x.xpol.mx, ...
         '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)
hold on
plot(double(x.xpol.mjd), ...
         x.xpol.val, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'LineWidth',0.1)

hold off
xlim([x.xpol.mjd(1)-0.05  x.xpol.mjd(end)+0.05])

title('x pole');
xlabel('MJD');
ylabel('x-pole [mas]');
xtickformat('%.2f')
ax = gca;
ax.XAxis.Exponent = 0;

styleAxes();

%% y-pole
subplot(2,2,2);

yline(0,'--');
hold on
errorbar(double(x.ypol.mjd), ...
         x.ypol.val, ...
         x.ypol.mx, ...
         '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)

hold on
plot(double(x.ypol.mjd), ...
         x.ypol.val, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'LineWidth',0.1)

hold off
xlim([x.ypol.mjd(1)-0.05  x.ypol.mjd(end)+0.05])

title('y pole');
xlabel('MJD');
ylabel('y-pole [mas]');
xtickformat('%.2f')
ax = gca;
ax.XAxis.Exponent = 0;

styleAxes();

%% dUT1
subplot(2,2,3);
yline(0,'--');
hold on
errorbar(double(x.dut1.mjd), ...
         x.dut1.val, ...
         x.dut1.mx, ...
         '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)
hold on
plot(double(x.dut1.mjd), ...
         x.dut1.val, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'LineWidth',0.1)
hold on
yline(0,'--');
hold off
xlim([x.dut1.mjd(1)-0.05  x.dut1.mjd(end)+0.05])



title('UT1-UTC');
xlabel('MJD');
ylabel('UT1-UTC [ms]');
xtickformat('%.2f')
ax = gca;
ax.XAxis.Exponent = 0;

styleAxes();

%% Nutation dX / dY
subplot(2,2,4);

% Kleiner horizontaler Offset nur für die Darstellung:
% 0.08 MJD = ca. 1 h 55 min
% Dadurch liegen dX und dY optisch nicht direkt übereinander.
nutOffset = 0.005;

yline(0,'--');
hold on
errorbar(double(x.nutdx.mjd) - nutOffset, ...
         x.nutdx.val, ...
         x.nutdx.mx, ...
         '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)
hold on
errorbar(double(x.nutdy.mjd) + nutOffset, ...
         x.nutdy.val, ...
         x.nutdy.mx, ...
         '-.','color',[0.6 0.6 0.6],'LineWidth',0.2)
hold on

plot(double(x.nutdx.mjd)- nutOffset, ...
         x.nutdx.val, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:),'LineWidth',0.1)
hold on
plot(double(x.nutdy.mjd)+ nutOffset, ...
         x.nutdy.val, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(3,:),'LineWidth',0.1)
hold off

xlim([x.nutdx.mjd(1)-0.05  x.nutdx.mjd(end)+0.05])

title('Nutation offsets dX / dY');
xlabel('MJD');
ylabel('[mas]');
legend('','','','dX','dY','Location','best');

xtickformat('%.2f')
ax = gca;
ax.XAxis.Exponent = 0;
styleAxes();

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');

%% ============================================================
% FIGURE 3: Station coordinate estimates
%% ============================================================

antNNR = [opt_.stat.nnr_inc];
idantnoNNR = find(antNNR==0);

idx = ~cellfun(@isempty, {x.coorx.val});
nStae=sum(idx);
ordstat=[1:nStae];

fig3 = figure('Color',[1 1 1], ...
              'Position',[50 100 1400 700]);

sgtitle('Station coordinate estimates', ...
        'FontWeight','bold');

coordNames  = {'coorx','coory','coorz'};
coordLabels = {'\DeltaX [cm]','\DeltaY [cm]','\DeltaZ [cm]'};

for p = 1:3

    subplot(1,3,p);

    vals = zeros(nStae,1);
    errs = zeros(nStae,1);

    for k = 1:nStae
        vals(k) = x.(coordNames{p})(k).val;
        errs(k) = x.(coordNames{p})(k).mx;
    end

    yline(0,'--');
    hold on;
    errorbar(ordstat, ...
             vals, ...
             errs, ...
             '.','color',[0.6 0.6 0.6],'LineWidth',0.2)
    hold on
    plot(ordstat, ...
             vals, ...
             'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'LineWidth',0.1)
    if nStae>0
        hold on
        plot(ordstat(idantnoNNR), ...
                 vals(idantnoNNR), ...
                 'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(8,:),'LineWidth',0.1)
    end
    hold off

    set(gca, ...
        'XTick',1:nSta, ...
        'XTickLabel',antNames, ...
        'XTickLabelRotation',45);

    ylabel(coordLabels{p});
    title(coordLabels{p});
    xlim([0 nSta+1])

    styleAxes();

end



set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');


%% ============================================================
% FIGURE: Source coordinate estimates
%% ============================================================


idinNNR = cellfun(@(xx) isequal(xx,1), {x.soura.inNNR});

idx = ~cellfun(@isempty, {x.soura.val});
nSou=sum(idx);
ordsou=[1:nSou];



if opt_.pw_sou
    souNames = {x.soura.name};
    for i=1:nSou
        idsou=strcmp(souNames(i),{sources.q.name});
        RAapr(i) = [sources.q(idsou).ra2000];
        DEapr(i) = [sources.q(idsou).de2000];
    end
else
    RAapr = [sources.q.ra2000];
    DEapr = [sources.q.de2000];
    souNames = {sources.q.name};
end



fig4 = figure('Color',[1 1 1], ...
              'Position',[50 100 1400 700]);

sgtitle('Source coordinate estimates', ...
        'FontWeight','bold');

coordNames  = {'soura','soude'};
coordLabels = {'\DeltaRA* [mas]','\DeltaDe [mas]'};

for p = 1:2

    subplot(2,1,p);

    vals = zeros(nSou,1);
    errs = zeros(nSou,1);

    for k = 1:nSou
        vals(k) = x.(coordNames{p})(k).val(1);
        errs(k) = x.(coordNames{p})(k).mx(1);
    end

    yline(0,'--');

    if opt_.est_sourceNNR | opt_.pw_sou
        if strcmp('soura',coordNames{p}) & ~isempty(vals)
            vals = vals.*cos(DEapr)'.*15;
            errs = errs.*cos(DEapr)'.*15;
        end
    
        hold on;
        errorbar(DEapr, ...
                 vals, ...
                 errs, ...
                 '.','color',[0.6 0.6 0.6],'LineWidth',0.2)
        hold on
        plot(DEapr, ...
                 vals, ...
                 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'LineWidth',0.1)
        if opt_.est_sourceNNR
            hold on
            plot(DEapr(idinNNR), ...
                     vals(idinNNR), ...
                     'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(8,:),'LineWidth',0.1)
            legend('','','NNR','wo NNR')
        end
        hold off
    end
    set(gca, ...
        'XTick',[-pi/2:pi/6:pi/2], ...
        'XTickLabel',[-90:30:90])
    

    xlabel('De')
    ylabel(coordLabels{p});
    xlim([-pi/2 pi/2]);

    styleAxes();

end
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');


%% ============================================================
% FIGURE 5: Atmosphere / TGR
%% ============================================================

fig5 = figure('Color',[1 1 1], ...
              'Position',[50 100 1400 600]);

for g=1:2
    subplot(2,1,g)
    for i = 1:nSta
    
        if g==1
            z = x.egr(i);
        else
            z = x.ngr(i);
        end
        errorbar(double(z.mjd), ...
                 z.val, ...
                 z.mx, ...
                 '-.','color',[0.6 0.6 0.6],'LineWidth',0.2);
        hold on
        p1{i}=plot(double(z.mjd), ...
                 z.val, ...
                 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',colors(i,:),'LineWidth',0.1);
        hold on
    end
    hold on
    yline(0,'--');
    hold off
    
    if g==1
        ylabel('TRGe [cm]');
        title(sprintf('Tropospheric Gradient Delay | MJD %.3f - %.3f', ...
          firstMJD,lastMJD), ...
          'FontWeight','bold');
    else
        ylabel('TRGn [cm]');
        xlabel('MJD');
    end
               
    xtickformat('%.2f')
    ax = gca;
    ax.XAxis.Exponent = 0;
    
    legend([p1{:}],antNames,'Location','eastoutside');
    styleAxes();
end

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperOrientation','landscape');


%% ============================================================
%% Local function
%% ============================================================

function styleAxes()

    ytickformat('%.2f')
    ax = gca;

    ax.Color  = [1 1 1];
    ax.XColor = [0.2 0.2 0.2];
    ax.YColor = [0.2 0.2 0.2];
    
    box on;
    grid on;

end


end