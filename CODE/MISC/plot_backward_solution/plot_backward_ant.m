% Function for plotting of the antenna coordinate estimates from the back solution.
% Input: ant_DIRIN.txt file created by function backward_solution.m


% Created for VieVS by Hana Spicakova
% 01 Aug 2011
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory
% WORK, script changed to a function; 
% plot_backward_ant(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%

function plot_backward_ant(DIRIN, DIROUT)

path_outglob = '../OUT/GLOB/';


%%
formatANT = '%5c %8c %f %f %f %s';
[param, anam, val(:,1),val(:,2),val(:,3), ses]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '.txt'], formatANT, 'commentstyle','matlab');
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

%  anames='VERAMZSW'

%%

fidCoord = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/averCoord_' DIRIN '.txt'],'wt');
 fprintf(fidCoord,'%% station          dx            dy           dz  [m]  \n');
 formatCoord = '%c%c%c%c%c%c%c%c      %10.4f   %10.4f   %10.4f  \n';


nant=size(anames,1);

       
xmin = 45650;
xmax = 58849;
sfont=14;

for i=1:nant
   aname=anames(i,:);
   id=strcmp(cellstr(anam),cellstr(aname));
   if sum(id)>0
       data=val(id,:);
       par=param(id,:);

       idx=strcmp(cellstr('ant_x'),cellstr(par));
       idy=strcmp(cellstr('ant_y'),cellstr(par));
       idz=strcmp(cellstr('ant_z'),cellstr(par));

       dataX=data(idx,:);
       dataY=data(idy,:);
       dataZ=data(idz,:);

       
       [svalX(:,1) idsr]=sort(dataX(:,1));
       svalX(:,2)=dataX(idsr,2);
       svalX(:,3)=dataX(idsr,3);

       [svalY(:,1) idsr]=sort(dataY(:,1));
       svalY(:,2)=dataY(idsr,2);
       svalY(:,3)=dataY(idsr,3);

       [svalZ(:,1) idsr]=sort(dataZ(:,1));
       svalZ(:,2)=dataZ(idsr,2);
       svalZ(:,3)=dataZ(idsr,3);
       
       
       figure(i)
       subplot(3,1,1)
       errorbar(svalX(:,1),svalX(:,2),svalX(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',12);
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       title(aname)
       ylabel('\Deltax [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;
       
       subplot(3,1,2)
       errorbar(svalY(:,1),svalY(:,2),svalY(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',12);
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       ylabel('\Deltay [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;
       
       subplot(3,1,3)
       errorbar(svalZ(:,1),svalZ(:,2),svalZ(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',12);
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       ylabel('\Deltaz [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;
      
       % ' ' --> '_' for file name
       idblank=find(aname(1:max(find(aname(:)~=' ')))==' ');
       aname(idblank)='_';
       % delete blank spaces at the end of station name: 'KOKEE   ' --> 'KOKEE'
       stname=char(cellstr(aname));


%         hgsave([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '_' stname]) %.fig
%         orient landscape      
%         print('-depsc' ,'-r400',['path_outglob BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '_' stname]); %.eps
          print('-dpng' ,'-r500',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '_' stname]);

%        close all
      
      
       % COMPUTE AVERAGE POSITION
       averESTm=[mean(svalX(:,2))/100 mean(svalY(:,2))/100   mean(svalZ(:,2))/100];
      
       fprintf(fidCoord,formatCoord, aname,averESTm);
     
      
       clear data dataX dataY dataZ svalX svalY svalZ idsr aname averESTm
   end
end

 fclose(fidCoord);

