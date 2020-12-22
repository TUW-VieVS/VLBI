% #########################################################################
% #     plotSessionAnalysisToAxes
% #########################################################################
%
% DESCRITPION
% This function is used to handle all the session analysis tools provided by 
% the VieVS GUI. 
%
% AUTHOR 
%   Matthias Madzak (?)
%
% INPUT
%   handles      structure from the GUI (also containing data, e.g. residuals)
%
% OUTPUT
%   handles      structure from the GUI (also containing data, e.g. residuals)
% 
% COUPLING:
%   - repeatab.m
%
%
% CHANGES
%   2015-16-12, A. Hellerschmied: Exctact code from "vievs2_3.m" file and create this file.
%   2016-08-31, A. Girdiuk: bug-fix and repeatability is plotted correctly
%   2016-09-21, A. Girdiuk: warning messages are added in case of absence of station coordinates estimates
%   2017-03-03, A. Hellerschmied: Removed absolute positioning of the plot window.
%   2017-05-29, A. Hellerschmied: The correlation matrix is now plotted correctly again 
%   2017-06-14, A. Hellerschmied: Bug-fix: indexing error for plotting the correlation matirx corrected
%   2017-06-20, A. Hellerschmied: Solved problem with indices when plotting the correlation matrix


function handles = plotSessionAnalysisToAxes(handles)
% This function plots the session analysis as selected in the interface.

% make axes empty
cla(handles.axes_plot_sessionAnalysis)
colorbar('off')

set(handles.axes_plot_sessionAnalysis, 'Projection', 'orthographic', 'Ydir', 'normal');
axis(handles.axes_plot_sessionAnalysis, 'auto')

% get index of current selected session
curSelSesInd=get(handles.popupmenu_plot_sessionAnalysis_session, 'Value');

% plot others (2nd, 3rd, 4th) session if chosen (and available)
others2plot=[get(handles.checkbox_plot_sessionAnalysis_add2, 'Value');...
    get(handles.checkbox_plot_sessionAnalysis_add3, 'Value');...
    get(handles.checkbox_plot_sessionAnalysis_add4, 'Value')];

% get indices of those which are to plot ([1,0,1] means 2nd and 4th to plot
indMore2plot=find(others2plot);
        
% define colros /sizes/...
myColors={'r', 'b', 'g'};
myCirclesizes=[38,22,6];
mySymbols={'^', 'x', 's'};
        
% set the axes of session analysis to current axis
set(handles.figure_vievs2,'CurrentAxes',handles.axes_plot_sessionAnalysis)

if get(handles.radiobutton_plot_sessionAnalysis_network, 'Value')
    % plot station network
    antenna=handles.data.plot.sessionAnalysis.antennaFiles(1).antenna(curSelSesInd).antenna;
    
    [phi,lam]=xyz2ell([[antenna.x]', [antenna.y]', [antenna.z]']);
    lam=lam*180/pi;
    phi=phi*180/pi;
    
    handles.axes_plot_sessionAnalysis=...
    axesm('MapProjection','giso','MeridianLabel','on','ParallelLabel','on');
    gridm on
    framem on
    % check if coast.mat is there, if not use coastlines.mat (from R2020b)
    if exist('coast.mat','file')            
      load coast
      plotm(lat,long, 'Color',[0.6,0.6,0.6]);
    else 
      load coastlines
      plotm(coastlat,coastlon, 'Color',[0.6,0.6,0.6]);
    end
    hold on
    scatterm(phi,lam,70,'k','filled')
    textm(phi+2,lam+2,{antenna.name},'FontSize',10)
    
    % if more than the main are to be plotted
    if ~isempty(indMore2plot)

        for iMore=1:length(indMore2plot)
            curInd=indMore2plot(iMore)+1; % 2|3|4
            
            % get current selected session index
            curSelSesIndMore=get(eval(['handles.popupmenu_plot_sessionAnalysis_session', num2str(curInd)]), 'Value');
            
            antenna=handles.data.plot.sessionAnalysis.antennaFiles(curInd).antenna(curSelSesIndMore).antenna;

            [phi,lam]=xyz2ell([[antenna.x]', [antenna.y]', [antenna.z]']);
            lam=lam*180/pi;
            phi=phi*180/pi;

            scatterm(phi,lam,myCirclesizes(curInd-1),myColors{curInd-1},'filled')
            textm(phi+2,lam+2,{antenna.name},'FontSize',10)
        end
    end
    hold off
    
elseif get(handles.radiobutton_plot_sessionAnalysis_baselLeRep, 'Value')
    % plot baseline length repeatability            
    % get sessions with estimated coordinates
    allCoorx={handles.data.plot.sessionAnalysis.x_files(1).x_.coorx};
    allCoory={handles.data.plot.sessionAnalysis.x_files(1).x_.coory};
    allCoorz={handles.data.plot.sessionAnalysis.x_files(1).x_.coorz};
    coordsAvail=false;
    for ind_x_=1:length(allCoorx)
        if ~isempty([allCoorx{1,ind_x_}.col])
            coordsAvail(ind_x_)=true;
        else
            fprintf('You need to estimate the station coordinates for %s session !\n',...
                handles.data.plot.sessionAnalysis.sessionnamesShort{1}.list(ind_x_,6:end));
        end
    end
    if sum(coordsAvail)~=length(allCoorx)
        fprintf('In this analysis, %i session(s) will be used for calculation instead of %i session(s) in the absence of station coordinates estimates !\n',sum(coordsAvail),length(allCoorx));
    end
    
    limitation=10;
    outfile=[]; % no file: []
    rigFormErr=0;
    printToCommand=0;
    process_list = handles.data.plot.sessionAnalysis.sessionnamesShort{1}.list;
    [blr,wblr,bl,blnames]=repeatab([], process_list(coordsAvail,:), limitation, outfile, [],  0, printToCommand, [], rigFormErr,...
        { {handles.data.plot.sessionAnalysis.x_files(1).x_(coordsAvail)},...
        {handles.data.plot.sessionAnalysis.antennaFiles(1).antenna(coordsAvail).antenna},... % CHeck here!
        {handles.data.plot.sessionAnalysis.atpaFiles(1).atpa(coordsAvail).atpa},...
        {handles.data.plot.sessionAnalysis.optFiles(1).opt(coordsAvail).opt} });
    nanVals=isnan(blr);
    blr=blr(~nanVals);
    wblr=wblr(~nanVals);
    bl=bl(~nanVals);
    blnames=blnames(~nanVals);

    if isempty(bl)
        fprintf('Non-empty bas_out.m demands, please change folder in the main panel (black)!\n');
    end
    
    plot(bl/1000,blr*100, 'o', 'color', 'k', 'linewidth',2)
    hold on

    % quadratic curve (LS)
    A=[(bl(:)/1000).^2, bl(:)/1000,ones(length(bl),1)];
    est=A\(blr(:)*100);
    x_BLR=min(bl)/1000:100:max(bl)/1000;
    plot(x_BLR, est(1)*x_BLR.^2+est(2)*x_BLR+est(3), '-', 'color', 'k', 'linewidth', 2)

    % if baseline names should be plotted
    if get(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Value')
        text(bl/1000,blr*100,blnames);
    end

    % if more than the main are to be plotted
    if ~isempty(indMore2plot)
        
        for iMore=1:length(indMore2plot)
            curInd=indMore2plot(iMore)+1; % 2|3|4
            
            allCoorx={handles.data.plot.sessionAnalysis.x_files(curInd).x_.coorx};
            allCoory={handles.data.plot.sessionAnalysis.x_files(curInd).x_.coory};
            allCoorz={handles.data.plot.sessionAnalysis.x_files(curInd).x_.coorz};
            
            coordsAvail=false;
            flag_speaking=0;
            for ind_x_=1:length(allCoorx)
                if ~isempty([allCoorx{1,ind_x_}.col])
                    coordsAvail(ind_x_)=true;
                else
                    if ~flag_speaking
                        flag_speaking=1;
                        switch indMore2plot(iMore)
                            case 1
                                color = 'red';
                            case 2
                                color = 'blue';
                            case 3
                                color = 'green';
                        end
                        fprintf('Additional plot panel %i marked with %s speaking:\n',indMore2plot(iMore),color);
                    end
                    fprintf('You need to estimate the station coordinates for %s session !\n',...
                        handles.data.plot.sessionAnalysis.sessionnamesShort{curInd}.list(ind_x_,6:end));
                end
            end
            if sum(coordsAvail)~=length(allCoorx)
                fprintf('In this analysis, %i session(s) will be used for calculation instead of %i session(s) in the absence of station coordinates estimates !\n',sum(coordsAvail),length(allCoorx));
            end

            limitation=10;
            outfile=[]; % no file: []
            process_list = handles.data.plot.sessionAnalysis.sessionnamesShort{curInd}.list;
            [blr,wblr,bl,blnames]=repeatab([], process_list(coordsAvail,:), limitation, outfile, [], ...
                0, printToCommand, [], rigFormErr,...
                { {handles.data.plot.sessionAnalysis.x_files(curInd).x_(coordsAvail)},...
                {handles.data.plot.sessionAnalysis.antennaFiles(curInd).antenna(coordsAvail).antenna},...
                {handles.data.plot.sessionAnalysis.atpaFiles(curInd).atpa(coordsAvail).atpa},...
                {handles.data.plot.sessionAnalysis.optFiles(curInd).opt(coordsAvail).opt} });
            nanVals=isnan(blr);
            blr=blr(~nanVals);
            wblr=wblr(~nanVals);
            bl=bl(~nanVals);
            blnames=blnames(~nanVals);

            plot(bl/1000,blr*100, ...
                [myColors{curInd-1},mySymbols{curInd-1}], 'linewidth',2)

            % quadratic curve (LS)
            A=[(bl(:)/1000).^2, bl(:)/1000,ones(length(bl),1)];
            est=A\(blr(:)*100);
            x_BLR=min(bl)/1000:100:max(bl)/1000;
            plot(x_BLR,est(1)*x_BLR.^2+est(2)*x_BLR+est(3), '-', 'color',...
                myColors{curInd-1}, 'linewidth', 2)

            % if baseline names should be plotted
            if get(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Value')
                text(bl/1000,blr*100,blnames);
            end
        end
    end
    
    hold off
    
elseif get(handles.radiobutton_plot_sessionAnalysis_corMatrix, 'Value')
    % plot correlation matrix
    
    % get first and last parameter to be plotted
    allFirstParInMenu   = get(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'String');
    allLastParInMenu    = get(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'String');
    firstpar            = allFirstParInMenu{get(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Value')};
    lastpar             = allLastParInMenu{get(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value')};
    
    % get required x_ and atpa_ file
    x_=handles.data.plot.sessionAnalysis.x_files(1).x_(curSelSesInd);
    atpa_=handles.data.plot.sessionAnalysis.atpaFiles(1).atpa(curSelSesInd).atpa;
    
    % get column indices
    if size(x_.(firstpar),2)==1
        colfirst=x_.(firstpar).col(1);
    else
        partials        = {x_.(firstpar).col};
        ind_partials    = find(~cellfun(@isempty,partials));
        colfirst        = x_.(firstpar)(ind_partials(1)).col(1);
    end
    if size(x_.(lastpar),2)==1
        collast = x_.(lastpar).col(end);
    else
        collast = x_.(lastpar)(end).col;
        if isempty(collast) % If there are no estimates for the "last" station, e.g. ref. station for clock estimates
            collast = x_.(lastpar)(end-1).col;
            if isempty(collast) % If there are no estimates for the "last" station, e.g. ref. station for clock estimates
                collast = x_.(lastpar)(end-2).col;
            end
        end
    end
    collast = collast(end);
    
    % get parameters which may be plotted
    posspar = {'pwclk'; 'rqclk'; 'zwd'; 'ngr'; 'egr'; 'xpol'; 'ypol'; 'dut1'; 'nutdx'; 'nutdy'; 'soura'; 'soude'; 'coorx'; 'coory'; 'coorz'; 'sat_pos1'; 'sat_pos2'; 'sat_pos3'; 'scale'; 'bdclko'}; 
    
    % get parameters which were extimated for current session
    parest=cell(length(posspar),1);
    for iPossiblePar=1:length(posspar)
        if isfield(x_,posspar{iPossiblePar})
            if isfield(x_.(posspar{iPossiblePar}), 'val')
                parest{iPossiblePar}=posspar{iPossiblePar};
            end
        end
    end
       
    % delete cell entries of parameters which were not estimated in session
    parest(cellfun(@isempty, parest))=[];
    ind1 = find(strcmp(parest,firstpar));
    ind2 = find(strcmp(parest,lastpar));
    parint = parest(ind1:ind2);
    nint(1) = 0;
    
    for k = 2:length(parint)+1
          nint(k) = nint(k-1)+length([x_.(parint{k-1}).col]);
    end
    nint = nint(1:end-1)+0.5;

    % Extract the Qxx matrix of your interest
    qxx = inv(atpa_.mat);
    qxxred = qxx(colfirst:collast, colfirst:collast);

    % Calculate the correlation matrix
    varii = diag(qxxred);
    rxy = qxxred./sqrt(varii*varii');

    % Plot the correlation matrix
    imagesc(rxy,[-1 1]);
    colorbar(handles.axes_plot_sessionAnalysis);
    axis image;
    set(handles.axes_plot_sessionAnalysis,'YTickLabel',parint(logical([~(nint(2:end)-nint(1:end-1)==0),1])));
    set(handles.axes_plot_sessionAnalysis,'YTick',nint(logical([~(nint(2:end)-nint(1:end-1)==0),1])));
    set(handles.axes_plot_sessionAnalysis,'YGrid','on');
    set(handles.axes_plot_sessionAnalysis,'XTick',nint(logical([~(nint(2:end)-nint(1:end-1)==0),1])));
    set(handles.axes_plot_sessionAnalysis,'XTickLabel',[]);
    set(handles.axes_plot_sessionAnalysis,'XGrid','on');
    set(handles.axes_plot_sessionAnalysis,'LineWidth',1.5);
end
