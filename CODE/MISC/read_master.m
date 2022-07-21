
function [sesname, sescode] = read_master(yr)

fileID = fopen(['../DATA/MASTER/master' num2str(yr,'%02.0f') '.txt']);
C = textscan(fileID,'%s %s %s %s %f %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','|','CommentStyle','-----','HeaderLines',10);
fclose(fileID);

sescode = C{3};
sescode(end) = []; % last line

sesnam = char(join([C{4} C{13}],''));
id7 = isspace(sesnam(:,7));
sesnam(id7,7) = '_';
sesnam(:,8:end)=[];
sesnam(end,:)=[]; % last line

sy(1:size(sesnam,1),1) = yr;
sesname = string([num2str(sy,'%02.0f') sesnam]);

if exist(['../DATA/MASTER/master' num2str(yr,'%02.0f') '-vgos.txt'])
    fileID = fopen(['../DATA/MASTER/master' num2str(yr,'%02.0f') '-vgos.txt']);
    V = textscan(fileID,'%s %s %s %s %f %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','|','CommentStyle','-----','HeaderLines',10);
    fclose(fileID);

    vsescode = V{3};
    vsescode(end) = []; % last line

    vsesnam = char(join([V{4} V{13}],''));
    vid7 = isspace(vsesnam(:,7));
    vsesnam(vid7,7) = '_';
    vsesnam(:,8:end)=[];
    vsesnam(end,:)=[]; % last line

    vsy(1:size(vsesnam,1),1) = yr;
    vsesname = string([num2str(vsy,'%02.0f') vsesnam]);
    
    
    sesname = [sesname; vsesname];
    sescode = [sescode; vsescode];
end



