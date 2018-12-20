function AxesLimitsChanged( hobj, evt, plink )
%LinkPropChangedFcn Links 2 compatible properties

curSession=get(plink.handles.popupmenu_plot_residuals_session, 'Value');

SessionStartTimeMJD =  plink.handles.data.plot.res(curSession).mjd(1);

xTicksLabel = datetime(get(plink.handles.axes_plot_residuals,'XTick')/24+SessionStartTimeMJD,'ConvertFrom','ModifiedJulianDate');
hhmm = datestr(xTicksLabel,'hh:MM');
date = unique(datestr(xTicksLabel,'dd.mm.yyyy'),'rows');
set(plink.handles.axes_plot_residuals, 'XTickLabel',hhmm);

dateLab = '';
for i=1:size(date,1)
    if i == 1
        dateLab = [dateLab date(i,:)];
    else
        dateLab = [dateLab ' - ' date(i,:)];
    end
end

set(plink.handles.axes_plot_residuals.XLabel, 'String', dateLab);
end
