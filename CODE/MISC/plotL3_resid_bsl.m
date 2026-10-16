% Hana Krasna, 2026-09-24
% Creates .pdf file with post-fit residual plots from the main solution per
% baselines

% sol = 'first', 'main'

function plotL3_resid_bsl(L3,ses,sol)

pdf_filename = ['../OUT/' ses '_' sol '_postfitresid.pdf'];
txt_filename = ['../OUT/' ses '_' sol '_resid.txt'];

load(['../DATA/LEVEL3/' L3 '/res_' ses '.mat' ]);

mm2ps = 1/299792458*10^9;

% residuals per baseline
baselines = double(res.baselineOfObs);
if strcmp(sol,'main')
    x  = double(res.mainVal(:));
    mx = double(res.sigma_residuals_aposteriori(:));
else
    x  = double(res.firstVal(:));
    mx = double(res.sigma_from_fringe_fitting(:))./mm2ps .*10^12 .*10; %cm
end
mjd  = double(res.mjd(:));
sourceNames = cellstr(strtrim(res.allSourceNames));
statNames   = cellstr(strtrim(res.allStatNames));


baselineID_all = strings(size(baselines,1),1);


for i = 1:size(baselines,1)

    s1 = baselines(i,1);
    s2 = baselines(i,2);

    if s1 <= numel(statNames) && s2 <= numel(statNames)
        name1 = statNames{s1};
        name2 = statNames{s2};

        pair = sort({name1,name2});
        baselineID_all(i) = string(pair{1}) + "-" + string(pair{2});
    end
end

uniqueBase = unique(baselineID_all);

baseName = {};
baseWRMS = [];
baseNobs = [];
baseNproblem = [];

% Create a figure that is invisible to speed up processing
fig = figure('Visible', 'off');
sfont=10;


mjd_datetimeALL = datetime(mjd, 'ConvertFrom', 'mjd');

for i = 1:numel(uniqueBase)

    idx = baselineID_all == uniqueBase(i);

    baseName{end+1,1} = char(uniqueBase(i));
    baseNobs(end+1,1) = sum(idx);
   

    resbsl=[]; mresbsl=[]; mjdresbsl=[];
    resbsl = x(idx).*10; %mm
    mresbsl = mx(idx).*10; %mm
    mjdresbsl = mjd(idx);

    P=1./(mresbsl.^2);
    wrms1 = sqrt(((resbsl.*resbsl)'*P)/sum(P)); 

    % 1. Convert MJD to MATLAB datetime
    % MJD is the number of days since November 17, 1858
    mjd_datetime = datetime(mjdresbsl, 'ConvertFrom', 'mjd');
    xla = string(datetime(mjd_datetime(1),'Format','MMMM d, yyyy'));



    % Clear current figure for a new page
    clf(fig);
    
    yyaxis left
    %errorbar(mjd_datetime, resbsl, mresbsl,'.','color',[0.7 0.7 0.7])
    %hold on
    plot(mjd_datetime, resbsl,'o','Color',[55./255,  126./255, 184./255],'Markersize',3,'MarkerEdgeColor','k','MarkerFaceColor',[55./255,  126./255, 184./255],'LineWidth',0.1)
    hold off
    ylabel('residuals (mm)')
    
    % Get the current tick locations
    ticksy = get(gca, 'YTick');
    % Convert those numbers to strings
    la = string(num2str(ticksy', '%.0f'));
    % Apply the labels
    yticklabels(la)

    title([char(uniqueBase(i))])

text(0.05, 0.98, ['obs: ', num2str(length(resbsl)) ], ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'FontWeight', 'normal', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none');
text(0.05, 0.93, ['wrms: ', num2str(wrms1,'%6.2f'), ' mm (' num2str(wrms1*mm2ps,'%5.1f') ' ps)' ], ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, ...
    'FontWeight', 'normal', ...
    'BackgroundColor', 'none', ...
    'EdgeColor', 'none');


    
    % Get the current limits of the left axis
    yl = ylim; 

    % Create the right axis scaling
    yyaxis right
    % Set the right axis limits by multiplying the left limits by 1e12
    ylim([yl(1)*mm2ps, yl(2)*mm2ps]); 
    ylabel('residuals (ps)');

    % Get the current tick locations
    ticksy = get(gca, 'YTick');
    % Convert those numbers to strings
    la = string(num2str(ticksy', '%.0f'));
    % Apply the labels
    yticklabels(la)

    % Use datetick or xtickformat to make the date display readable
    xtickformat('HH:mm'); 
    xlabel(xla)
    xlim([mjd_datetimeALL(1)-0.01, mjd_datetimeALL(end)+0.01])
    % If the data span is very small, you might prefer:
    % xtickformat('HH:mm:ss'); 
    set(gca,'FontSize',sfont)

    if i == 1
        exportgraphics(fig, pdf_filename, 'ContentType', 'vector');
    else
        exportgraphics(fig, pdf_filename, 'ContentType', 'vector', 'Append', true);
    end


end

close(fig);


fid=fopen(txt_filename,'wt');
fprintf(fid,'#  baseline          mjd                 date-time               source         resid (mm)    sigma (mm) \n');
for i=1:length(mjd)
    parts = strsplit(baselineID_all(i), '-');
    fprintf(fid,'%-8s %-8s    %17.12f    %s     %s    %10.2f    %10.2f  \n',...
             parts{1}, parts{2}, mjd(i), mjd_datetimeALL(i), string(sourceNames(res.source(i))), x(i)*10, mx(i)*10 );
end
fclose(fid);
