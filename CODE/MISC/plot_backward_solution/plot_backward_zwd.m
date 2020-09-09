% Function for plotting of the antenna coordinate estimates from the back solution.
% Input: ant_DIRIN.txt file created by function backward_solution.m


% Created for VieVS by Hana Spicakova
% 01 Aug 2011
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory
% WORK, script changed to a function; 
% plot_backward_zwd(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%

function plot_backward_zwd(DIRIN, DIROUT)

path_outglob = '../OUT/GLOB/';

%%

formatZWD = '%3c %8c %f %f %f %s';
[param, anam, val(:,1),val(:,2),val(:,3), ses]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '.txt'], formatZWD, 'commentstyle','matlab');
%val = mjd, val, std


%% for which stations should the plot be done?
for i=1:size(anam,1)
   if i==1
       anames=anam(i,:);
   else
       if ~strcmp(cellstr(anames),cellstr(anam(i,:)))
           anames=[anames; anam(i,:)];
       end
   end
end

%  anames=['FORTLEZA' ];

%%

sfont=14;
xmin = 45650;
xmax = 58849;

for i=1:size(anames,1)
   aname=anames(i,:);
   id=strcmp(cellstr(anam),cellstr(aname));
   
   if sum(id)>0
       data=val(id,:);

       [sval(:,1) idsr]=sort(data(:,1));
       sval(:,2)=data(idsr,2);
       sval(:,3)=data(idsr,3);

       
       % ' ' --> '_' for file name
       idblank=find(aname(1:max(find(aname(:)~=' ')))==' ');
       aname(idblank)='_';
       % delete blank spaces at the end of station name: 'KOKEE   ' --> 'KOKEE'
       stname=char(cellstr(aname));

%        % write separate zwd.txt
%        fid=fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '_' stname '.txt'],'wt');
%        fprintf(fid,'%c%c%c%c%c%c%c%c \n', aname);
%        fprintf(fid,'%s \n', 'mjd      zwd[cm]     formal error [cm]');
%        fprintf(fid,'%10.4f    %7.2f    %7.2f \n', sval');     
%        fclose(fid);
       
       
       figure(i)
       subplot(2,1,1)
       errorbar(sval(:,1),sval(:,2),sval(:,3),'.','Color',[.4 .4 .4 ],'MarkerEdgeColor','k');
%        xlim([xmin xmax]) 
%        set(gca,'XTick',[51544-(365.25*10) 51544 51544+(365.25*10) 51544+(365.25*20)])
%        set(gca,'XTickLabel',[1990;2000;2010;2020],'FontSize',sfont)
       title(aname)
       ylabel('\Deltazwd [cm]','FontSize',sfont)
       set(gca,'LineWidth',1,'FontSize',sfont);%,'fontweight','bold')  ;

%        hgsave([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '_' stname]) %.fig
%        orient landscape
%        print('-depsc' ,'-r400',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '_' stname]); %.eps
       print('-dpng' ,'-r500',[path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '_' stname]);

%        close all
       clear data sval idsr
   end
end





