% Function to read additional TRF data using snx data.
%
% INPUT
%  ns_codes			main loop file of creating superstation file (mk_superstatFile.m)
%  datafile   		path of TRF
%  datafilename		name of TRF as it will display at structure
%  break0  			for preallocating: the empty struct for a TRF break
% 
% OUTPUT
%  ns_codes(updated)
%
% CREATED
%  11.05.2016 A. Girdiuk

 	function [ns_codes] = trf_by_snx_reader(ns_codes,datafile,datafilename,break0)

	% read sinex data
    [blocknames blockdata]=readSnx(datafile,'blocks2read', {'SITE/ID',...
    'SOLUTION/EPOCHS', 'SOLUTION/ESTIMATE'});

    for kSta=1:length(blockdata{3}.stat)
    % find current station in ns_codes
        foundStatLog=strcmpi(blockdata{1}(kSta).domes, {ns_codes.domes});
    
        if sum(foundStatLog)==0
        % try to find in first 8 characters of desription (might work for
        % VTRF2014/ITRF2015 as the first 8 chars look like ivs name)
            foundStatLog=strcmpi(blockdata{1}(kSta).descr(1:8), {ns_codes.name});
            if sum(foundStatLog)==0
                fprintf('Station %s (%s) not found in ns_codes!\n',...
                    blockdata{1}(kSta).descr,datafilename);
                keyboard;
            end
        end
    
        if sum(foundStatLog)>0
            foundInd=find(foundStatLog);
            for iCorrectStat=1:length(foundInd)
            
                curIStat=foundInd(iCorrectStat);

            % for all breaks
                for breakInd=1:length(blockdata{3}.stat(kSta).break)

                % preallocate break (first break only)
                    if breakInd==1
                        ns_codes(curIStat).(datafilename).break(1)=break0;
                    end
                    if isempty(blockdata{2}(kSta).break)
                        fprintf(' Station %s has no epoch! (but coordinates might be there)\n',...
                        	blockdata{1}(kSta).descr);
                    else
                    % make (by hand) the start of 2, 3,.. break equal to
                    % end of previous break (that there are no "empty"
                    % epochs)
                        if breakInd>1
                            ns_codes(curIStat).(datafilename).break(breakInd).start=ns_codes(curIStat).(datafilename).break(breakInd-1).end;
                        else
                            ns_codes(curIStat).(datafilename).break(breakInd).start = blockdata{2}(kSta).break(breakInd).start;
                        end
                    
                        ns_codes(curIStat).(datafilename).break(breakInd).end = blockdata{2}(kSta).break(breakInd).end;
                    end      
                
                    ns_codes(curIStat).(datafilename).break(breakInd).epoch = blockdata{3}.stat(kSta).break(breakInd).epoch;
                    ns_codes(curIStat).(datafilename).break(breakInd).x = blockdata{3}.stat(kSta).break(breakInd).x;
                    ns_codes(curIStat).(datafilename).break(breakInd).y = blockdata{3}.stat(kSta).break(breakInd).y;
                    ns_codes(curIStat).(datafilename).break(breakInd).z = blockdata{3}.stat(kSta).break(breakInd).z;
                    ns_codes(curIStat).(datafilename).break(breakInd).x_sigma = blockdata{3}.stat(kSta).break(breakInd).x_sigma;
                    ns_codes(curIStat).(datafilename).break(breakInd).y_sigma = blockdata{3}.stat(kSta).break(breakInd).y_sigma;
                    ns_codes(curIStat).(datafilename).break(breakInd).z_sigma = blockdata{3}.stat(kSta).break(breakInd).z_sigma;

                    ns_codes(curIStat).(datafilename).break(breakInd).vx = blockdata{3}.stat(kSta).break(breakInd).vx;
                    ns_codes(curIStat).(datafilename).break(breakInd).vy = blockdata{3}.stat(kSta).break(breakInd).vy;
                    ns_codes(curIStat).(datafilename).break(breakInd).vz = blockdata{3}.stat(kSta).break(breakInd).vz;
                    ns_codes(curIStat).(datafilename).break(breakInd).vx_sigma = blockdata{3}.stat(kSta).break(breakInd).vx_sigma;
                    ns_codes(curIStat).(datafilename).break(breakInd).vy_sigma = blockdata{3}.stat(kSta).break(breakInd).vy_sigma;
                    ns_codes(curIStat).(datafilename).break(breakInd).vz_sigma = blockdata{3}.stat(kSta).break(breakInd).vz_sigma;
                end
            end % for all found (strcmpi) stations
        end
    end
