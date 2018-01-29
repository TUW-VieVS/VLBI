% ************************************************************************
%   Description:
%   This function plots a matlab figure which shows stations in sessions
%   included in the global adjustment
%
%   Input:										
%      antname_plot        names of the antennas
%      antactiv_plot       matrix [num. of antennas, num. of sessions]
%                          1/0 if the antenna was included in the session
%                          or not
%   Output:                
%      figure(1)
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%%


function plot_antactiv(antname_plot,antactiv_plot)
   
figure(1)
    spy(antactiv_plot,8,'k');
    axis fill
    xlabel('Sessions')
    set(gca,'YTick',1:size(antname_plot,1))
    set(gca,'YTickLabel',antname_plot)
    set(gca,'xlim', [0 size(antactiv_plot,2)+1])
    title('Activity of the antennas')
    orient landscape
