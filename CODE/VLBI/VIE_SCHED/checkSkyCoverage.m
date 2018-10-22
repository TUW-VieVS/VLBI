% #########################################################################
% #     checkSkyCoverage
% #########################################################################
%
% DESCRIPTION
%   This function calculates the mean sky coverage in a 0.5 hour interval per station and
%   plots this information ta a figure.
%
%
% CREATED
%   2016     M. Schartner
%
% REFERENCES
%
%
% COUPLING
%   - zcoverage
%
%
% INPUT
%   -
%
%   - sched                         : sched data structure
%   - station                       : station data structure
%   - PARA                          : Global scheduling parameter strucutre
%
%
% OUTPUT
%   - meanSky                       : mean sky coverage
%   - f                             : figure handle
%
% CHANGES:
% - 2016-12-22, A. Hellerschmied: - Header added
%                                 - Function now works correctly for sessions shorter than 0.5 hours
% - 2017-05-28, M. Schartner: -remove dependency from statistics and machine learning toolbox 

function [ meanSky,f ] = checkSkyCoverage( sched,PARA,station )

allScans = [sched.scan];

allStanum = [allScans.nsta];
allStart = [allScans.startmjd];
obs = [allScans.sta];
c = 1;
for i = 1:length(allStanum)
    for j = 1:allStanum(i)
        obs(c).startmjd=allStart(i);
        c = c+1;
    end
end

nsta = max([obs.staid]);
startTime = PARA.startmjd;
endTime =   PARA.endmjd;
steps = round((endTime-startTime)*24/.5);
if steps == 0
    steps = 1;
end

storage = nan(ceil(steps/12)*12,nsta);

for i= 1:nsta
    thisObs = obs([obs.staid] == i);
    for j = 1:steps
        startStep = startTime + (j-1)*30/1440;
        endStep =   startTime + (j)*30/1440;
        obsInTime = thisObs([thisObs.startmjd]>=startStep & [thisObs.startmjd]<endStep);
        cat = zeros(1,length(obsInTime));
        for k = 1:length(obsInTime)
            cat(k) = zcoverage(obsInTime(k).az, obsInTime(k).el);
        end
        storage(j,i) = length(unique(cat));
    end
end

for i = 1:length(station)
    meanSky(i) = mean(storage(~isnan(storage(:,i)),i));
end

if PARA.MULTISCHED || ~PARA.disp_sky_coverage
    f = figure('Units','normalized','Position',[0.05 .25 .9 .5],'visible','off');
else
    f = figure('Units','normalized','Position',[0.05 .25 .9 .5]);
end
plots = ceil(steps/12);
h = (0.85-(plots-1)*0.05)/plots;

for i = 1:plots
    subplot('Position',[0.025 1-i*0.05-(i)*h .8 h])
    bar((i-1)*6+(0:.5:5.5),storage((i-1)*12+[1:12],:))
    set(gca,'XTick',(i-1)*6+(0:.5:5.5))
    set(gca,'YTick',1:2:13)
    xlim([(i-1)*6-.5 (i-1)*6+6]);
    ylim([0 13])
    set(gca,'Ygrid','on')
end
meanmean = mean(meanSky);
meanvar = var(meanSky);

meanSkyStr = strcat(' (mean:',{' '},strsplit(num2str(meanSky)),')');

legend(strcat({station.name},meanSkyStr),'Location',[0.85 .5 .125 .45]);
if nsta>7
    colormap(colorcube(nsta))
else
    colormap(jet(nsta))
end
annotation('textbox',[0.025 0.025 0.8 0.03],...
    'String',{'hours since start'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off','FontWeight','bold','EdgeColor','none');

annotation('textbox',[0.025 0.975 0.8 0.03],...
    'String',{'Number of sky coverage - 60min intervall - 30min overlap'},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off','FontWeight','bold','EdgeColor','none','FontSize',12);

annotation('textbox',[0.85 .45 .125 .05],...
    'String',{sprintf('total mean: %.3f',meanmean)},...
    'HorizontalAlignment','left',...
    'FitBoxToText','off','FontWeight','bold','EdgeColor','none','FontSize',12);

annotation('textbox',[0.85 .4 .125 .05],...
    'String',{sprintf('var: %.3f',meanvar)},...
    'HorizontalAlignment','left',...
    'FitBoxToText','off','FontWeight','normal','EdgeColor','none','FontSize',12);

annotation('textbox',[0.85 .35 .125 .05],...
    'String',{sprintf('std: %.3f',sqrt(meanvar))},...
    'HorizontalAlignment','left',...
    'FitBoxToText','off','FontWeight','normal','EdgeColor','none','FontSize',12);

[minVal,idx] = min(meanSky);
annotation('textbox',[0.85 .30 .125 .05],...
    'String',{sprintf('min: %.3f (%s)',minVal,station(idx).name)},...
    'HorizontalAlignment','left',...
    'FitBoxToText','off','FontWeight','normal','EdgeColor','none','FontSize',12);

[maxVal,idx] = max(meanSky);
annotation('textbox',[0.85 .25 .125 .05],...
    'String',{sprintf('max: %.3f (%s)',maxVal,station(idx).name)},...
    'HorizontalAlignment','left',...
    'FitBoxToText','off','FontWeight','normal','EdgeColor','none','FontSize',12);


end
