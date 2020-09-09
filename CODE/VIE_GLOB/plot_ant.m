% ************************************************************************
%   Description:
%   This function plots a matlab figure which shows the location on the world map
%   of the antennas, whose coordinates/velocities are estimated in the global adjustment
%
%   Input:										
%       X0                  station coordinates
%       Y0
%       Z0
%       excidant            indices of excluded stations from the NNT/NNR
%                           condition
%       refname             names of the antennas
%
%
%   Output:                
%      figure(2)
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%%



function plot_ant(X0,Y0,Z0,excidant,refname)

[phi,lam,h]=xyz2ell([X0;Y0;Z0]'); %rad rad m
lam=lam/pi*180;
phi=phi/pi*180;


figure (2)
    landareas = shaperead('landareas.shp','UseGeoCoords',true);
    geoshow(landareas,'FaceColor',[.8 .8 .8],'EdgeColor',[.6 .6 .6]);
    hold on
    scatter(lam,phi,'o','SizeData',25,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1])
    hold on
    scatter(lam(excidant),phi(excidant),'o','SizeData',25,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
    hold on
    text(lam,phi,refname,'HorizontalAlignment','left','VerticalAlignment','Top','FontSize',8)
    hold on
    hold on
    title({'Stations in global adjustment';'blue: included in NNT/NNR condition; red: excluded from NNT/NNR condition'})
    axis equal;
    set(gca,'XLim',[-180 180],'YLim',[-90 90])
    set(get(gca,'XLabel'),'String','longitude [°]');
    set(get(gca,'YLabel'),'String','latitude [°]');
    orient landscape
    hold off



