
% Barbara Tumler, Hana Krasna
% 2026-06-30

% d - name of LEVEL3 dir
% ses - name of session
% plotScanRes - 1/0 - create outPdf
% plotX - 1/0 - create outPdfEst
% tableout - 1/0 - save tables as csv

function L3plots(d,ses,plotScanRes,plotX,tableout)
close all

lev = 'LEVEL3';

scanFile = ['../DATA/' lev '/' d  '/' ses '_scan.mat'];
antFile  = ['../DATA/' lev '/' d  '/' ses '_antenna.mat'];
souFile  = ['../DATA/' lev '/' d  '/' ses '_sources.mat'];
resFile  = ['../DATA/' lev '/' d  '/res_' ses '.mat'];
optFile  = ['../DATA/' lev '/' d  '/opt_' ses '.mat'];
xFile    = ['../DATA/' lev '/' d  '/x_' ses '.mat'];


%% Load data

S=load(scanFile);
scan=S.scan;
S=load(antFile);
antenna=S.antenna;
S=load(souFile);
sources=S.sources;
S=load(resFile);
res=S.res;
S=load(optFile);
opt_=S.opt_;


if plotScanRes
    outPdf = [ses '_plots.pdf'];
    
    fig1 = plot_scan_meteo(antenna, scan);
    exportFigureToPDF(fig1,outPdf,false)
    
    fig2 = plot_scan_cable(antenna, scan);
    exportFigureToPDF(fig2,outPdf,true)
    
    % Plot wrms from res
    [wrmsfig1, wrmsfig2, wrmsfig3, wrmsfig4, wrmsfig5, wrmsfig6] = plot_res_wrms(res, opt_,tableout);
    exportFigureToPDF(wrmsfig1,outPdf,true)
    exportFigureToPDF(wrmsfig2,outPdf,true)
    exportFigureToPDF(wrmsfig3,outPdf,true)
    exportFigureToPDF(wrmsfig4,outPdf,true)
    exportFigureToPDF(wrmsfig5,outPdf,true)
    exportFigureToPDF(wrmsfig6,outPdf,true)
    fprintf('\nSaved combined PDF:\n')
    fprintf('  %s\n',outPdf)
end


if plotX
    outPdfEst = [ses '_plots_est.pdf'];
    S=load(xFile);
    x_=S.x_;

    % Plot estimates from x_
    [xfig1, xfig2, xfig3, xfig4, xfig5] = plot_x_est(x_,opt_,sources);
    exportFigureToPDF(xfig1,outPdfEst,false)
    exportFigureToPDF(xfig2,outPdfEst,true)
    exportFigureToPDF(xfig3,outPdfEst,true)
    exportFigureToPDF(xfig4,outPdfEst,true)
    exportFigureToPDF(xfig5,outPdfEst,true)
    fprintf('\nSaved combined PDF:\n')
    fprintf('  %s\n',outPdfEst)
end

function exportFigureToPDF(fig,outPdf,appendFlag)
    figure(fig);
    axs = findall(fig,'Type','axes');
    for k = 1:numel(axs)
        if isprop(axs(k),'Toolbar')
            axs(k).Toolbar.Visible = 'off';
        end
    end
    drawnow;
    if appendFlag
        exportgraphics(fig,['../OUT/' outPdf],'ContentType','vector','Append',true);
    else
        exportgraphics(fig,['../OUT/' outPdf],'ContentType','vector');
    end
end

end