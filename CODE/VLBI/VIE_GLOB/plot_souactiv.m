% ************************************************************************
%   Description:
%   This function plots a matlab figure which shows sources in sessions
%   included in the global adjustment
%
%   Input:
%      qrefname            names of the sources (sorted by order of sessions in processlist)
%      souactiv            matrix [num. of sources, num. of sessions]
%                          1/0 if the source was included in the session
%                          or not (sorted by order of sessions in
%                          processlist)
%      mjd_all             mjd of each session
%
%   Output:										
%      souname_plot        names of the sources (sorted for the plot)
%      souactiv_plot       matrix [num. of sources, num. of sessions]
%                          1/0 if the source was included in the session
%                          or not
%      figure(3)
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   09 Jun 2011 by Hana Spicakova   the sorting of sources for plotting is
%           done in this function. It was removed from the main program vie_glob.
%%


function [souname_plot,souactiv_plot]=plot_souactiv(qrefname,souactiv,mjd_all)

lns = size(qrefname,1); % final number of sources

souactiv(lns+1,:)=mjd_all;
souactiv=sortrows(souactiv',lns+1);
souactiv(:,lns+1)=[];

[souactiv_plot,id_sou]=sortrows(souactiv',-[1:size(souactiv,1)]);
souname_plot=qrefname(id_sou,:);


figure(3)
    spy(souactiv_plot,8,'k');
    axis fill
    xlabel('Sessions')
    set(gca,'YTick',1:size(souname_plot,1))
    set(gca,'YTickLabel',souname_plot)
    set(gca,'xlim', [0 size(souactiv_plot,2)+1])
    title('Observed sources')
    orient landscape

