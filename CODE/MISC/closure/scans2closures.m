% Matthias Schartner, Hana Krasna

function [] = scans2closures(pthOUT, name, scans, typ)

if typ == 1
    fid = fopen(fullfile([pthOUT 'VALUES'],['closures_B_' name '.txt']),'w');
else
    fid = fopen(fullfile([pthOUT 'VALUES'],['closures_G_' name '.txt']),'w');
end

fprintf(fid,'%% source sta1     sta2     sta3     closure delay [s]  fe [s] delay flags   mjd       time\n');


c = 0;
c_max = 0;
for i_scan = 1 : length(scans)
    obs = [scans(i_scan).obs];
    stas = unique([[obs.sta1] [obs.sta2]]);
    
    triangles = nchoosek(stas,3);
    for i_triangle = 1:size(triangles,1)
        triangle = triangles(i_triangle,:);
        bool_ab = strcmp([obs.sta1], triangle{1}) & strcmp([obs.sta2], triangle{2});
        bool_ac = strcmp([obs.sta1], triangle{1}) & strcmp([obs.sta2], triangle{3});
        bool_bc = strcmp([obs.sta1], triangle{2}) & strcmp([obs.sta2], triangle{3});
        bc_order=1;
        if sum(bool_bc)==0
            bool_bc = strcmp([obs.sta2], triangle{2}) & strcmp([obs.sta1], triangle{3}); % change order of bc stations
            bc_order = 0;
        end
        if any(bool_ab) && any(bool_ac) && any(bool_bc)
            ab = obs(bool_ab);
            ac = obs(bool_ac);
            bc = obs(bool_bc);
            if typ ==1
                if bc_order == 1
                    tau_abc = ab.delay + bc.delay - ac.delay + (bc.rate * ab.delay_theoretical);
                    sigma_tau_abc = sqrt(ab.sigmaDelay^2 + bc.sigmaDelay^2 + ac.sigmaDelay^2 + (ab.delay_theoretical * bc.sigmaRate)^2 );
                else
                    tau_abc = ab.delay - bc.delay - ac.delay - (bc.rate * ac.delay_theoretical);     
                    sigma_tau_abc = sqrt(ab.sigmaDelay^2 + bc.sigmaDelay^2 + ac.sigmaDelay^2 + (ac.delay_theoretical * bc.sigmaRate)^2 );                     
                end                
            else
                if bc_order == 1
                    tau_abc = ab.delay + bc.delay - ac.delay;  
                else
                    tau_abc = ab.delay - bc.delay - ac.delay;  
                end                
                sigma_tau_abc = sqrt(ab.sigmaDelay^2 + bc.sigmaDelay^2 + ac.sigmaDelay^2 );
            end
            fprintf(fid,'%8s %8s %8s %8s %12.5e %.5e    %1.0f %1.0f %1.0f %15.7f %s\n', ab.src{1}, triangle{1}, triangle{2}, triangle{3}, tau_abc, sigma_tau_abc, ab.delayFlag, bc.delayFlag, ac.delayFlag, ab.mjd, datestr(ab.time,'yyyy.mm.dd_HH:MM:SS'));
            c = c+1;
        end
    end
    c_max = c_max + size(triangles,1);
end
fprintf('    %d of %d closures found\n', c, c_max);
fclose(fid);

end

