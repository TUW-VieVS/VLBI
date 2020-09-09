% ************************************************************************
%   Description:
%   This function plots a matlab figure which shows the location on the sky
%   of the sources, whose coordinates are estimated in the global adjustment
%
%   Input:										
%       RA_all              source coordinates [rad,rad]
%       De_all
%       excidsouc           indices of excluded sources from the NNR
%                           condition
%
%
%   Output:                
%      figure(4)
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%%

function plot_sou(RA_all,De_all,excidsouc)

figure(4)
    scatter(RA_all,De_all,'o','SizeData',25,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 1])
    hold on
    scatter(RA_all(excidsouc),De_all(excidsouc),'o','SizeData',25,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
    hold on
    title({'Radio sources in global adjustment';'blue: included in NNR condition; red: excluded from NNR condition'})
    set(gca,'XLim',[0 2*pi],'YLim',[-pi/2 pi/2])
    set(get(gca,'XLabel'),'String','RA [h]');
    set(gca,'XTick',[0;pi/2;pi;1.5*pi;2*pi])
    set(gca,'XTickLabel',[0;6;12;18;24])
    set(get(gca,'YLabel'),'String','De [°]');
    set(gca,'YTick',[-pi/2;-pi/4;0;pi/4;pi/2])
    set(gca,'YTickLabel',[-90;-45;0;45;90])
    orient landscape
    hold off





