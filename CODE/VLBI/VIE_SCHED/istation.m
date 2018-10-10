% Purpose  
%   Read station information.
% History  
%   2010-04-06   Jing SUN   Created
%   2016-06-11	M. Schartner: new field acc1 and acc2 in station struct, now reads acceleration.cat 


function [station] = istation(staname, INFILE, PARA)

% station number
stanum = size(staname, 1);

% initialize station information struct
for ista = 1 : stanum
    station(ista).name(1:8) = staname(ista,1:8); 
    % antenna
    station(ista).id(1)     = ' ';
    station(ista).axis(1:4) = ' ';
    station(ista).offset    = 0.0;
    station(ista).rate1     = 0.0;
    station(ista).acc1      = 0.0;
    station(ista).c1        = 0.0;
    station(ista).lim11     = 0.0;
    station(ista).lim12     = 0.0;
    station(ista).rate2     = 0.0;
    station(ista).acc2      = 0.0;
    station(ista).c2        = 0.0;
    station(ista).lim21     = 0.0;
    station(ista).lim22     = 0.0;
    station(ista).diam      = 0.0;
    station(ista).po(1:2)   = ' ';
    station(ista).eq(1:3)   = ' ';
    station(ista).ms(1:2)   = ' ';
    station(ista).aznp      = 0.0;
    station(ista).azn1      = 0.0;
    station(ista).azn2      = 0.0;
    station(ista).azc1      = 0.0;
    station(ista).azc2      = 0.0;
    station(ista).azw1      = 0.0;
    station(ista).azw2      = 0.0;
    % position
    station(ista).xyz(1:3)  = 0.0; 
    station(ista).llh(1:3)  = 0.0; 
    % equip
    station(ista).sefdpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = 0.0; 
    station(ista).equip(1:8) = ' ';
    % mask
    station(ista).hmasknum              = 0;
    station(ista).hmask(1:PARA.MAX_HOR) = 0.0;   
    % minimum SNR 
    for i=1:PARA.MAX_BANDNUM
        station(ista).minsnr(i) = PARA.MIN_SNR(i);
    end
    % downtime 
    station(ista).downum       = 0;
    station(ista).downstart(1) = 0.0;
    station(ista).downend(1)   = 0.0;
end

% read file 'antenna.cat'
fprintf(PARA.fid_header,'1.1 read antenna file: %s \n', INFILE.antenna);
[station] = iantenna(INFILE.antenna, station, PARA);

% read file 'position.cat'
fprintf(PARA.fid_header,'1.2 read position file: %s \n', INFILE.position);
[station] = iposition(INFILE.position, station, PARA);

% read file 'equip.cat'
fprintf(PARA.fid_header,'1.3 read equip file: %s \n', INFILE.equip);
[station] = iequip(INFILE.equip, station, PARA);

% read file 'mask.cat'
fprintf(PARA.fid_header,'1.4 read mask file: %s \n', INFILE.mask);
[station] = imask(INFILE.mask, station, PARA);

% read snrmin.txt
fprintf(PARA.fid_header,'1.5 read snrmin file: %s \n', INFILE.snrmin);
[station] = isnrmin(INFILE.snrmin, station, PARA);

% read down.txt
fprintf(PARA.fid_header,'1.6 read down file: %s \n', INFILE.down);
[station] = idown(INFILE.down, station, PARA);

% read acceleration.txt
fprintf(PARA.fid_header,'1.7 read acceleration file: %s \n', INFILE.acceleration);
[station] = iacceleration(INFILE.acceleration, station, PARA);

% check 1-letter code
% letters = 'A':1:'Z';
% id = {station.id};
% for i = 1:length(id)-1
%     bool = strcmp(id{i},{id{1+i:end}});
%     if any(bool)
%         for j = 1:length(letters)
%             if ~any(strcmp(letters(j),{id{1+i:end}}))
%                 station(i).id = letters(j);
%                 warning('same 1-letter code for %s and %s! New 1-letter code for %s is now %s',...
%                     station(i).name,station(find(bool)+i).name,station(i).name,letters(j));
%                 break
%             end
%         end
%         
%     end
% end


