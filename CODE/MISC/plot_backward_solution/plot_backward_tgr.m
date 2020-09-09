% Function for plotting of the antenna coordinate estimates from the back solution.
% Input: ant_DIRIN.txt file created by function backward_solution.m


% Created for VieVS by Hana Spicakova
% 01 Aug 2011
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory
% WORK, script changed to a function; 
% plot_backward_tgr(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%

function plot_backward_tgr(DIRIN, DIROUT)

path_outglob = '../OUT/GLOB/';


%%
formatTGR = '%3c %8c %f %f %f %s';
[param, anam, val(:,1),val(:,2),val(:,3),ses ]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/tgr_' DIRIN '.txt'], formatTGR, 'commentstyle','matlab');
%val = mjd, val, std




%% for which stations should the plot be done?
% all stations
for i=1:size(anam,1)
   if i==1
       anames=anam(i,:);
   else
       if ~strcmp(cellstr(anames),cellstr(anam(i,:)))
           anames=[anames; anam(i,:)];
       end
   end
end

% anames='MEDICINA'


%%

sfont=14;
xmin = 45650;
xmax = 58849;

nant=size(anames,1);

for i=1:nant
   aname=anames(i,:);
   id=strcmp(cellstr(anam),cellstr(aname));
   if sum(id)>0

       data=val(id,:);
       par=param(id,:);

       idngr=strcmp(cellstr('ngr'),cellstr(par));
       idegr=strcmp(cellstr('egr'),cellstr(par));

       dataNGR=data(idngr,:);
       dataEGR=data(idegr,:);

       [svalNGR(:,1) idsr]=sort(dataNGR(:,1));
       svalNGR(:,2)=dataNGR(idsr,2);
       svalNGR(:,3)=dataNGR(idsr,3);

       [svalEGR(:,1) idsr]=sort(dataEGR(:,1));
       svalEGR(:,2)=dataEGR(idsr,2);
       svalEGR(:,3)=dataEGR(idsr,3);


       figure(i)
       subplot(2,1,1)
       errorbar(svalNGR(:,1),svalNGR(:,2),svalNGR(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       title(aname)
       ylabel('\Deltangr [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;

       subplot(2,1,2)
       errorbar(svalEGR(:,1),svalEGR(:,2),svalEGR(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       ylabel('\Deltaegr [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;

       
       % ' ' --> '_' for file name
       idblank=find(aname(1:max(find(aname(:)~=' ')))==' ');
       aname(idblank)='_';
       % delete blank spaces at the end of station name: 'KOKEE   ' --> 'KOKEE'
       stname=char(cellstr(aname));

%        hgsave([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/tgr_' DIRIN '_' stname]) %.fig
%        orient landscape
%        print('-depsc' ,'-r400',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/tgr_' DIRIN '_' stname]); %.eps
       print('-dpng' ,'-r500',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/tgr_' DIRIN '_' stname]);

%        close all
       clear data dataNGR dataEGR svalNGR svalEGR idsr
   end
end



