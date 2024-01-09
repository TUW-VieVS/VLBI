% ## Master file format version 2.0

function [sesname, sescode, sestype, sesdate, sesstat] = read_master_v2(yr)

%   SESSION      DATE     SESSION    DOY TIME   DUR                         STATIONS                        SKED CORR  STATUS  DBC  SUBM DEL
%     TYPE     yyyymmdd     CODE     ddd hh:mm  h:mm                                                                  yyyymmdd CODE      days

fileID = fopen(['../DATA/MASTER/master' num2str(yr,'%04.0f') '.txt']);
C = textscan(fileID,'|%s %s %s %f %s %s   %s    %s %s %s %s %s %s', 'Delimiter','|','CommentStyle','-----','HeaderLines',9);
fclose(fileID);

sestype = deblank(C{1});
sesdate = deblank(C{2});
sescode = deblank(C{3});
sesstat = deblank(C{7});

sesname = char(join([C{2} deblank(C{3})]));
id7 = isspace(sesname(:,9));
sesname(id7,9) = '-';





