
% Function to read TRF data from SSC textfiles. Taken from mk_superstatFile.m (where
% this function is now called)
%
% INPUT
%  fileID  filename of the textfile
%  year    double containing the year of the epoch of the chosen TRF (is usually given
%          in the header of the textfile)
%          eg 2000 if itrf2005
%             2005 if itrf2008
%             2000 if vtrf2008
%  ns_codes main loop file of creating superstation file (mk_superstatFile.m)
%  trf     string: name of the trf (this will be the fieldname in the superstation file)
%          e.g. 'itrf2008'
%  break0  for preallocating: the empty struct for a TRF break
% 
% OUTPUT
%  ns_codes(updated)
%
% CREATED
%  09.05.2016 A. Girdiuk
%
% CHANGES
%  10.05.2016 A. Girdiuk: domes, site_name and id are not duplicated in subfields anymore

function [ns_codes] = ITRF_VLBI_SSC_reader(fileID,year,ns_codes,trf,break0)

    fid = fopen(fileID);

    epoch = modjuldat(year,1,1,0,0);
    if strcmp(trf,'vtrf2008')
        for j = 1: 5
            line = [fgetl(fid),'                                               '];
        end
    elseif strcmp(trf,'ivsTrf2014b')
        for j = 1: 4
            line = [fgetl(fid),'                                               '];
        end
    else
        for j = 1: 7
            line = [fgetl(fid),'                                               '];
        end
	end
    while ~feof(fid)
        line = [fgetl(fid),'                                                                                                                '];
    
        % find station index
        iStat=find(strcmp({ns_codes.domes}, line(1:9)));

        % get break number if exist and define if we have just one
        % break in total== there is empty at line(96:96) or more
        a = {char(line(96:96))} ;
        if strcmp(trf,'vtrf2008') || strcmp(trf,'ivsTrf2014b')
            a = {char(line(102:102))} ;
        end
        
        if strcmp(a, ' ');
            breakInd=1;
            justOneBreak=1; % 1 == no break number is given --> use 0 and 99999 as epochs
        else
            breakInd=str2double(a);
            justOneBreak=0;
        end
    
        % write coordinates to all ns_codes entries where to domes was found
        for iCorrectStat=1:size(iStat,2)
            curIStat=iStat(iCorrectStat);
        
            % preallocate break (first break only)
            if breakInd==1
                ns_codes(curIStat).(trf).break(1)=break0;
            end
        
            % only do once -> when first time go through one station
%            if breakInd==1
 %               ns_codes(curIStat).(trf).domes_nr = (line(1:9));
  %              ns_codes(curIStat).(trf).site_name = (line(11:25));
   %             ns_codes(curIStat).(trf).id = str2double(line(33:36));
    %        end

            % if there is no number given -> start=0, end=9999
            if justOneBreak==1
                startEpoch=0;
                endEpoch=99999;
            else %-> there are more breaks
                % get start date
                
                if strcmp(trf,'vtrf2008') || strcmp(trf,'ivsTrf2014b')
                   % get start date
                    startYr=str2double(line(104:105));
                    startDoy=str2double(line(107:109));
                    startSecod=str2double(line(111:115));

                    % get end date
                    endYr=str2double(line(117:118));
                    endDoy=str2double(line(120:122));
                    endSecod=str2double(line(124:128)); 
                else
                    startYr=str2double(line(98:99));
                    startDoy=str2double(line(101:103));
                    startSecod=str2double(line(105:109));
                
                    % get end date
                    endYr=str2double(line(111:112));
                    endDoy=str2double(line(114:116));
                    endSecod=str2double(line(118:121));
                end
                
                % get start epoch
                if startDoy==0
                    startEpoch=0;
                else
                    % convert yy:doy:secod to mjd
                    % make yy -> yyyy
                    if startYr>79
                        startYr=1900+startYr;
                    else
                        startYr=2000+startYr;
                    end
                    startEpoch=doy2jd(startYr,startDoy)-2400000.5;
                end

                % get end epoch
                if endDoy==0
                    endEpoch=99999;
                else
                    % convert yy:doy:secod to mjd
                    % yy -> yyyy
                    if endYr>79
                        endYr=1900+endYr;
                    else
                        endYr=2000+endYr;
                    end
                    try
                    endEpoch=doy2jd(endYr,endDoy)-2400000.5;
                    catch
                        keyboard;
                    end
                end
            end
            ns_codes(curIStat).(trf).break(breakInd).start = startEpoch;
            ns_codes(curIStat).(trf).break(breakInd).epoch = epoch;
            ns_codes(curIStat).(trf).break(breakInd).end = endEpoch;
            
            % if we are in 2nd, 3rd... break: write start epoch to end epoch of
            % previous break (ITRF2008 file does not include the correct end
            % dates of the first breaks (the end of the line is left empty even 
            % if there is a second break))
            if strcmp(trf,'itrf2008') && breakInd>1
                ns_codes(curIStat).(trf).break(breakInd-1).end=startEpoch;
            end

            ns_codes(curIStat).(trf).break(breakInd).x = str2double(line(38:49));
            ns_codes(curIStat).(trf).break(breakInd).y = str2double(line(51:62));
            ns_codes(curIStat).(trf).break(breakInd).z = str2double(line(64:75));
            
            if strcmp(trf,'vtrf2008') || strcmp(trf,'ivsTrf2014b')
                ns_codes(curIStat).(trf).break(breakInd).x_sigma = str2double(line(77:83));
                ns_codes(curIStat).(trf).break(breakInd).y_sigma = str2double(line(85:91));
                ns_codes(curIStat).(trf).break(breakInd).z_sigma = str2double(line(93:99));
            else
                ns_codes(curIStat).(trf).break(breakInd).x_sigma = str2double(line(77:81));
                ns_codes(curIStat).(trf).break(breakInd).y_sigma = str2double(line(83:87));
                ns_codes(curIStat).(trf).break(breakInd).z_sigma = str2double(line(89:93));
            end
        end
    
        % get the second line (with velocities)
        line = [fgetl(fid),'                                                                                                                '];

        % do the same loop again for the second row (velocities)
        for iCorrectStat=1:size(iStat,2)
        
            % write it to station I want
            curIStat=iStat(iCorrectStat);
            if strcmp(trf,'ivsTrf2014b')
                ns_codes(curIStat).(trf).break(breakInd).vx = str2double(line(43:49));
                ns_codes(curIStat).(trf).break(breakInd).vy = str2double(line(56:62));
                ns_codes(curIStat).(trf).break(breakInd).vz = str2double(line(69:75));
            else
                ns_codes(curIStat).(trf).break(breakInd).vx = str2double(line(44:49));
                ns_codes(curIStat).(trf).break(breakInd).vy = str2double(line(57:62));
                ns_codes(curIStat).(trf).break(breakInd).vz = str2double(line(70:75));
            end
            if strcmp(trf,'vtrf2008') || strcmp(trf,'ivsTrf2014b')
                ns_codes(curIStat).(trf).break(breakInd).x_sigma = str2double(line(77:83));
                ns_codes(curIStat).(trf).break(breakInd).y_sigma = str2double(line(85:91));
                ns_codes(curIStat).(trf).break(breakInd).z_sigma = str2double(line(93:99));
            else
                ns_codes(curIStat).(trf).break(breakInd).x_sigma = str2double(line(77:81));
                ns_codes(curIStat).(trf).break(breakInd).y_sigma = str2double(line(83:87));
                ns_codes(curIStat).(trf).break(breakInd).z_sigma = str2double(line(89:93));
            end
        end
    
    end

    fclose(fid);
