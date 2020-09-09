% Function for plotting of the source coordinate estimates from the back solution.
% Input: sou_DIRIN.txt file created by function backward_solution.m


% Created for VieVS by Hana Spicakova
% 01 Aug 2011
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory
% WORK, script changed to a function; 
% plot_backward_sou(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%

function plot_backward_sou(DIRIN, DIROUT)

path_outglob = '../OUT/GLOB/';

%%
formatSOU = '%3c %8c %f %f %f %s';
[param, sou, val(:,1),val(:,2),val(:,3), ses]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/sou_' DIRIN '.txt'], formatSOU, 'commentstyle','matlab');
%val = mjd, val, std


%% for which sources should the plot be done?
% all sources
for i=1:size(sou,1)
   if i==1
       sources=sou(i,:);
   else
       if ~strcmp(cellstr(sources),cellstr(sou(i,:)))
           sources=[sources; sou(i,:)];
       end
   end
end

% sources='0014+813';

%%

sfont=14;
xmin = 45650;
xmax = 58849;

nsou=size(sources,1);

for i=1:nsou
   source=sources(i,:);
   id=strcmp(cellstr(sou),cellstr(source));
   if sum(id)>0
       data=val(id,:);
       par=param(id,:);

       idra=strcmp(cellstr('sra'),cellstr(par));
       idde=strcmp(cellstr('sde'),cellstr(par));

       dataRA=data(idra,:);
       dataDE=data(idde,:);

       
       [svalRA(:,1) idsr]=sort(dataRA(:,1));
       svalRA(:,2)=dataRA(idsr,2);
       svalRA(:,3)=dataRA(idsr,3);

       [svalDE(:,1) idsr]=sort(dataDE(:,1));
       svalDE(:,2)=dataDE(idsr,2);
       svalDE(:,3)=dataDE(idsr,3);

       
       figure(i)
       subplot(2,1,1)
       errorbar(svalRA(:,1),svalRA(:,2),svalRA(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       title(source)
       ylabel('\DeltaRA [mas]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;

       subplot(2,1,2)
       errorbar(svalDE(:,1),svalDE(:,2),svalDE(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       ylabel('\DeltaDe [mas]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;

       
       % delete blank spaces at the end of station name: 'KOKEE   ' --> 'KOKEE'
       fnsou=[ char(cellstr(source(1:4))) char(cellstr(source(6:8))) ];

%        hgsave([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '_' stname]) %.fig
%        orient landscape      
%        print('-depsc' ,'-r400',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/sou_' DIRIN '_' fnsou]); %.eps
       print('-dpng' ,'-r500',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/sou_' DIRIN '_' fnsou]);

%        close all
       clear data dataRA dataDE svalRA svalDE idsr source
   end
end



