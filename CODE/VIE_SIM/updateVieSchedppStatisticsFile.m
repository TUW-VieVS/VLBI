function [] = updateVieSchedppStatisticsFile(file, t_mean_sig, t_rep)

fid = fopen(file);
if(fid == -1)
    error('Cannot open file: %s', file);
end
header = strsplit(fgetl(fid),',');
header(end) = [];
frewind(fid);
data = dlmread(file,',',1,0);
data(:,end) = [];
fclose(fid);

[~,idx] = sort(data(:,1));
data = data(idx,:);

simHeader = t_mean_sig.Properties.VariableNames(1:end);
simHeader = {'n_sim','dUT1_[mus]','x_pol_[muas]','y_pol_[muas]','x_nut_[muas]','y_nut_[muas]', 'average_3d_coordinates_[mm]', simHeader{9:end}};

mean_sig = t_mean_sig.Variables;
average = mean(mean_sig(:,9:end),2);
mean_sig = [mean_sig(:,[1 4:8]), average, mean_sig(:,9:end)];

rep = t_rep.Variables;
average = mean(rep(:,9:end),2);
rep = [rep(:,[1 4:8]), average, rep(:,9:end)];

for i=1:length(simHeader)
    header{end+1} = ['vievs_sim_mean_formal_error_' simHeader{i}];
end
for i=1:length(simHeader)
    header{end+1} = ['vievs_sim_repeatability_' simHeader{i}];
end
data = [data, mean_sig, rep];



outfile = [file(1:end-4) '_sim.csv'];
fid = fopen(outfile, 'wt');
fprintf(fid, '%s,',header{:});
fprintf(fid, '\n');
fclose(fid);
dlmwrite(outfile,data,'delimiter',',','-append');


end

