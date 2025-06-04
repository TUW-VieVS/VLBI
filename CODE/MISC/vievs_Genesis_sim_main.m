%% Options
vievs_root_str = '/data/USERS/hwolf/VieVS/VLBI';
%scripts_path_str = '/home/auto/auto_proc/scripts/';
%addpath(scripts_path_str);

% Process List:
pathProcessList = "/data/USERS/hwolf/VieVS/VLBI/WORK/PROCESSLIST";
process_list = load(pathProcessList + "/GEN_i60_Nc.mat");
process_list = process_list.process_list;

% Parameter File: 
pathParaFile = "/data/USERS/hwolf/VieVS/VLBI/WORK/PARAMETERS";
parameter = load(pathParaFile + "/GEN_i60.mat");
parameter = parameter.parameter;
%parameter.lsmopt.level1OutDir = strcat('GENESIS_2yr_wn10_v', num2str(i));

%runp
runp = load(vievs_root_str + "/WORK/runp_GEN.mat");
runp = runp.runp;
runp.init = 1;
runp.mod = 1;
runp.lsm = 1;
runp.sim = 1;
runp.parallel = 1; 
runp.nCores{1,1} = '3';
save(vievs_root_str + "/WORK/runp.mat",'runp')

pathSimFile = '/data/USERS/hwolf/VieVS/VLBI/DATA/LEVEL4';
simparam = load(pathSimFile + "/simparam_GEN.mat");

simparam = simparam.simparam;
simparam.idays = 1;
simparam.wn = 10;
simparam.wn_sat = 10; % adjust
save(pathSimFile + "/simparam.mat",'simparam')

% ##### Session analysis #####
for i = 1:1
	fprintf(1, '\n\n########################################################\n');
	fprintf(1, '##### Simulation %d #####\n', i);
	fprintf(1, '########################################################\n');
    %parameter.lsmopt.level1OutDir = strcat('GENESIS_2yr_wn10_v', num2str(i));
    name = strcat('GEN_i60_Nc_v', num2str(i));
    parameter.lsmopt.level1OutDir = name;
    runp.init_path = name; 
    runp.mod_path = name;
    runp.lsm_path = name;
    runp.glob_path = name; 
    save([vievs_root_str,'/WORK/guiparameter.mat'], 'parameter');
    save([vievs_root_str,'/WORK/process_list.mat'], 'process_list');
    save(vievs_root_str + "/WORK/runp.mat",'runp');
	
    vie_batch

    % ##### Delete temp. folders and data ##### 
    level0_path_str = [vievs_root_str, 'DATA/LEVEL0/', name, '/'];
    level1_path_str = [vievs_root_str, 'DATA/LEVEL1/', name, '/'];
    level3_path_str = [vievs_root_str, 'DATA/LEVEL3/', name, '/'];
    if isfolder(level0_path_str)
        rmdir(level0_path_str, 's');
    end
    if isfolder(level1_path_str)
        rmdir(level1_path_str, 's');
    end
    if isfolder(level3_path_str)
        delete([level3_path_str,'*_sources.mat']);
        delete([level3_path_str,'*_antenna.mat']);
        delete([level3_path_str,'*_scan.mat']);
        delete([level3_path_str,'atpa_*']);
        delete([level3_path_str,'res_*']);
    end
end


