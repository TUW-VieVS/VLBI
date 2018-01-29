% Purpose  
%   Read param.txt file.
%
% CREATED  
%   2010-04-06   Jing SUN
%
% CHANGES
%   2015-02-11, A. Hellerschmied: removed parameters PARA.EXPER and
%     PARA.DESCRIPTION (now via GUI input)
%   2015-03-04, A. Hellerschmied: Parameter OBSMODE_NAME added.
%   2015-03-19, D. Mayer: Parameter SCHEDULER added.
%   2015-07-08, A. Hellerschmied: Parameters for VieVS satellite scheduling added
%   2015-11-19, A. Hellerschmied: Parameters for VieVS auto satellite scheduling added
%   2016-05-20, M. Schartner: STAR mode added
%   2016-06-16, A. Hellerschmied: Parameters for VieVS auto satellite scheduling added (CALC_SKYCOV_SAT, SKYDT_SAT, MIN_SRCRP_SAT, CENTER_TRACK_EPOCHS)
%   2016-06-10, M. Schartner: band information is now in GUI
%   2016-10-19, M. Schartner: additional parameters added
%   2016-12-22, A. Hellerschmied: PARA.VEX_PARA_FILEPATH and PARA.SAT_FREQ_FILEPATH removed (path from vie_sched GUI is used instead)

function [PARA] = iparam(filename, stanum, PARA)

% open param.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    error('    no param.txt file !\n');
end

% read the scheduling parameters
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength < 21) | (~strcmp(line(1), ' ')))   %%%%%
        continue;
    end
    [paramname, count, errmsg, nextindex1] = sscanf(line(2:linelength), '%s', 1);
    index = 2 + nextindex1 - 1;
    [paramvalue, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%f', 1); 
    switch (paramname)
        case 'PARA.RATE1A'
            PARA.RATE1A  = paramvalue;
        case 'PARA.RATE2A'
            PARA.RATE2A = paramvalue;
        case 'PARA.MARGEL1'
            PARA.MARGEL1 = paramvalue;
        case 'PARA.MARGEL2'    
            PARA.MARGEL2 = paramvalue;
        case 'PARA.MIN_SRCRP'
            PARA.MIN_SRCRP = paramvalue;
        case 'PARA.SOURCE'
            PARA.SOURCE = paramvalue;
        case 'PARA.TAPETM'
            PARA.TAPETM = paramvalue;
        case 'PARA.IDLE'
            PARA.IDLE = paramvalue;
        case 'PARA.CALIBRATION'
            PARA.CALIBRATION = paramvalue;
        case 'PARA.MAXSLEWTIME'
            PARA.MAXSLEWTIME = paramvalue;
        case 'PARA.MAX_WAIT'
            PARA.MAX_WAIT = paramvalue;
        case 'PARA.CORSYNCH'
            PARA.CORSYNCH = paramvalue;
        case 'PARA.MAX_SCAN'
            PARA.MAX_SCAN = paramvalue;
        case 'PARA.MIN_SCAN'
            PARA.MIN_SCAN = paramvalue;
        case 'PARA.FILLINMODE'
            PARA.FILLINMODE = round(paramvalue);
        case 'PARA.FILLENDT'
            PARA.FILLENDT = paramvalue;
        case 'PARA.SCREEN'
            PARA.SCREEN = round(paramvalue);
        case 'PARA.MIN_STANUM'
            PARA.MIN_STANUM = paramvalue;
        case 'PARA.MIN_STASCAN'
            PARA.MIN_STASCAN = paramvalue;
        case 'PARA.MIN_STANUM_FI'
            PARA.MIN_STANUM_FI = paramvalue;
        case 'PARA.SUBNETTING'
            PARA.SUBNETTING = paramvalue;
        case 'PARA.SKYDT'
            PARA.SKYDT = paramvalue;
        case 'PARA.MIN_SRC2ANG'
            PARA.MIN_SRC2ANG = paramvalue * pi / 180.0;
        case 'PARA.SORTNUM'
            PARA.SORTNUM = paramvalue;
        case 'PARA.FORSI'
            PARA.FORSI = round(paramvalue);
        case 'PARA.UPSTA'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1); 
            PARA.UPSTA = paramstr;
        case 'PARA.DOWNSTA'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1); 
            PARA.DOWNSTA = paramstr;
        case 'PARA.SRCFRINGE'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1); 
            PARA.SRCFRINGE = paramstr;
        case 'PARA.STARMODE'
            PARA.STARMODE = paramvalue;
        case 'PARA.STRONGANT'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1); 
            PARA.STRONGANT = '        ';
            PARA.STRONGANT(1:length(paramstr))=paramstr;
        case 'PARA.CADENCE'
            PARA.CADENCE = paramvalue;
        case 'PARA.SCANDURA'
            PARA.SCANDURA = paramvalue;
        case 'PARA.TRACKSMODE'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.TRACKSMODE = paramstr;
        case 'PARA.CORRELATOR'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.CORRELATOR = paramstr;
        case 'PARA.TRACKSNAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.TRACKSNAME = paramstr;
        case 'PARA.OBSMODE_NAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.OBSMODE_NAME = paramstr;
        case 'PARA.SCHEDULER'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.SCHEDULER = paramstr;
		case 'PARA.WEIGHT_NUMBER_OF_OBS'	
			PARA.WEIGHT_NUMBER_OF_OBS = paramvalue;
		case 'PARA.WEIGHT_SKY_COVERAGE'	
			PARA.WEIGHT_SKY_COVERAGE = paramvalue;
		case 'PARA.WEIGHT_SCAN_END_TIME'	
			PARA.WEIGHT_SCAN_END_TIME = paramvalue;

		% #### Parameters specific for VieVS satellite scheduling ####     
        case 'PARA.TRF_FILEPATH'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.TRF_FILEPATH = paramstr;
        case 'PARA.TRF_FILENAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.TRF_FILENAME = paramstr;
        case 'PARA.EOP_FILENAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.EOP_FILENAME = paramstr;
        case 'PARA.EOP_FILEPATH'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.EOP_FILEPATH = paramstr;
        case 'PARA.REPOS_INT'
            PARA.REPOS_INT = paramvalue;
        case 'PARA.DELTA_REPOS_INT'
            PARA.DELTA_REPOS_INT = paramvalue;
        case 'PARA.MIN_SAT_SCAN_TIME'
            PARA.MIN_SAT_SCAN_TIME = paramvalue;
        case 'PARA.SAT_FREQ_FILENAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.SAT_FREQ_FILENAME = paramstr;
        case 'PARA.VEX_PARA_FILENAME'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.VEX_PARA_FILENAME = paramstr;
        case 'PARA.CHECK_SCAN_INT'
            PARA.CHECK_SCAN_INT = paramvalue;
        case 'PARA.PREOB_NEW_SOURCE'
            PARA.PREOB_NEW_SOURCE = paramvalue;
        case 'PARA.SAT_BLOCK_DUR_MIN'
            PARA.SAT_BLOCK_DUR_MIN = paramvalue;
        case 'PARA.Q_BLOCK_DUR_MIN'
            PARA.Q_BLOCK_DUR_MIN = paramvalue;
        case 'PARA.SAT_SCAN_DUR_SEC'
            PARA.SAT_SCAN_DUR_SEC = paramvalue;
        case 'PARA.INIT_SCAN_TYPE'
            index = 2 + nextindex1 - 1;
            [paramstr, count, errmsg, nextindex2] = sscanf(line(index:linelength), '%s', 1);
            PARA.INIT_SCAN_TYPE = paramstr;
        case 'PARA.CALC_SKYCOV_SAT'
            PARA.CALC_SKYCOV_SAT = paramvalue;
        case 'PARA.SKYDT_SAT'
            PARA.SKYDT_SAT = paramvalue;
        case 'PARA.MIN_SRCRP_SAT'
            PARA.MIN_SRCRP_SAT = paramvalue; 
        case 'PARA.CENTER_TRACK_EPOCHS'
            PARA.CENTER_TRACK_EPOCHS = paramvalue;
    end 
end

% close param.txt file
fclose(fid);


