function [] = scans2closures(name, scans)

fid = fopen(fullfile('Closures',['closures_' name '.txt']),'w');
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
        if any(bool_ab) && any(bool_ac) && any(bool_bc)
            ab = obs(bool_ab);
            ac = obs(bool_ac);
            bc = obs(bool_bc);
            tau_abc = ab.delay + bc.delay - ac.delay + (bc.rate * ab.delay_theoretical);
            fprintf(fid,'%8s %8s %8s %8s %.5e %s\n', ab.src{1}, triangle{1}, triangle{2}, triangle{3}, tau_abc, datestr(ab.time,'yyyy.mm.dd_HH:MM:SS'));
            c = c+1;
        end
    end
    c_max = c_max + size(triangles,1);
end
fprintf('    %d of %d closures found\n', c, c_max);
fclose(fid);

end

