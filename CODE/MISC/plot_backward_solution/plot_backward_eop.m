% Function for plotting of the EOP estimates from the back solution.
% Input: eop_DIRIN.txt file created by function backward_solution.m


% Created for VieVS by Hana Spicakova
% 01 Aug 2011
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory
% WORK, script changed to a function; 
% plot_backward_eop(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%

function plot_backward_eop(DIRIN, DIROUT)

path_outglob = '../OUT/GLOB/';

EOPs = ['xpol'; 'ypol';'dut1'; 'dX  '; 'dY  '];

%%
formatEOP = '%4c %f %f %f %s';
[param, val(:,1),val(:,2),val(:,3), ses]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '.txt'], formatEOP, 'commentstyle','matlab');
%val = mjd, val, std

neop=size(EOPs,1);

xmin = 45650;
xmax = 58849;
sfont=14;

for i=1:neop
   eop=EOPs(i,:);
   id=strcmp(cellstr(eop),cellstr(param));
   if sum(id)>0
       data=val(id,:);
       
       [sval(:,1) idsr]=sort(data(:,1));
       sval(:,2)=data(idsr,2);
       sval(:,3)=data(idsr,3);
              
       figure(i)
       errorbar(sval(:,1),sval(:,2),sval(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       title(eop)
       if strcmp(eop,'dut1')
           ylabel(['[ms]'],'FontSize',sfont);
%            ylim([-0.2   0.2])
       else
           ylabel('[mas]','FontSize',sfont);
%            ylim([-3   3])
       end

       % delete blank spaces at the end: 'dX  ' --> 'dX'
       steop=char(cellstr(eop));

       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;
       
%         hgsave([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '_' steop]) %.fig
%         orient landscape      
%         print('-depsc' ,'-r400',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '_' steop]); %.eps
      print('-dpng' ,'-r500',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '_' steop]);

%       close all
       clear data  sval idsr 
   end
end



