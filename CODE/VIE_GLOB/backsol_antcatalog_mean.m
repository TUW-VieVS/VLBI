
% Creates a catalog (mean position) for stations which were reduced in the global
% adjustment.
% Input: Session-wise absolute coordinates
% Hana Krasna, 2021-11-08



function backsol_antcatalog_mean(path_outglob,DIROUT,DIRIN)

%clear all

fileID = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_bcksol_sessionwise_' DIRIN '.txt']);
C = textscan(fileID,'%s %f %f %f %f %f %f %f %f    %s   %f %f %f %f %f %f %f %f %f ','Headerlines',6);
fclose(fileID);




fidCoord = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_bcksol_mean_' DIRIN '.txt'],'wt');
curDate=clock;          % current date and time
fprintf(fidCoord, 'Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
fprintf(fidCoord,'station              X [m]              Y [m]            Z [m]            vx [m/y]    vy [m/y]    vz [m/y]   epoch   start  end\n');
formatCoord = '%8s      %15.4f   %15.4f   %15.4f   %10.4f   %10.4f   %10.4f     %5.0f  %5.0f  %5.0f    %1.0f  \n';



% list of the stations
statlist = unique(C{1});


for i = 1:length(statlist)

    ids =[];
    ids = strcmp(statlist(i),C{1});
    
    X_ = C{2}(ids);
    Y_ = C{3}(ids);
    Z_ = C{4}(ids);
    vX_ = C{14}(ids);
    vY_ = C{15}(ids);
    vZ_ = C{16}(ids);
    
    ep_ = C{8}(ids);
    estart_ = C{18}(ids);
    eend_ = C{19}(ids);
    
    
    int1 = unique(estart_); 
    % how many intervals?
    for j = 1:size(int1,1)
        
        idbr = [];
        idbr = find(int1(j)==estart_);
        
        X = mean(X_(idbr));
        Y = mean(Y_(idbr));
        Z = mean(Z_(idbr));
        vX = mean(vX_(idbr));
        vY = mean(vY_(idbr));
        vZ = mean(vZ_(idbr));

        ep = mean(ep_(idbr)); % mean epoch over the included sessions
        
        estart = estart_(idbr(1));
        eend = eend_(idbr(1));
        
        
        aname=[statlist{i} blanks(8-size(statlist{i},2))];
        fprintf(fidCoord,formatCoord, aname, X, Y, Z, vX, vY, vZ, ep, estart, eend, 0 );
  
        
    end
end


 fclose(fidCoord);


