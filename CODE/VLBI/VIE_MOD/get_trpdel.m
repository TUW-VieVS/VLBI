% ************************************************************************
%   Description:
%   Finds & writes external tropospheric delays from the file into scan.
% 
%   References: 
%
%   Input:										
% 
%   Output:
% 
%   External calls: 	   
%       
%   Coded for VieVS: 
%   22 May 2012 by Lucia Plank
%
%   Revision: 
%   16 Dec 2014 by Daniel Landskron: partial Gn and Ge deleted (not used)
%   01 Sep 2015 by Daniel Landskron: curLine gets more conditions
%   03 Nov 2015 by Armin Hofmeister: Add handling of cases where 0, 1 or more matching
%          trp entries are found.
%          Avoid setting of input pressure and temperature to [] in case no matching trp entry is
%          found, otherwise VieVS would crash. Display a message in case no matching trp entry is
%          found.
%          Add use of only the first matching trp entry in case of multiple matching trp entries
%          and display a message about the use of the first matching trp entry.
%   27 Jan 2016 by Armin Hofmeister: Increase allowed difference in seconds between observation from
%          VieVS and from trp-file to up to 0.1 seconds.
%   
% ************************************************************************
function [scan] = get_trpdel(trpdata,scan,isc,ist,antenna,trpFileFoundLog,sourceNames)

%        % define if vmf1 is needed because no trp data is found in trp
%        % file (for current station in current scan)
%        vmf1NeededLog=0;

% but first: if we have trp file -> try to get data
if trpFileFoundLog==1
    
    % get correct line in trp file
    tmjd = []; year = []; month = []; day = []; hour = []; sec = [];
    tmjd = modjuldat(scan(isc).tim(1),scan(isc).tim(2),scan(isc).tim(3),scan(isc).tim(4),scan(isc).tim(5),scan(isc).tim(6));
    [year, month, day, hour, minu, sec] = mjd2date(tmjd);
    
    % note: The last condition here is because there are some few NGS-files where two scans of one
    %       station take place at exactly the same time, which is wrong due to interpolational
    %       problems, therefore it is needed to check also for matching source names.
    % note: The check of the seconds being different by up to 0.1 is needed since the trp-file
    %       reports only 1 decimal place for the seconds and rounding effects may lead to a
    %       difference compared to the VieVS seconds of up to 0.1 seconds.
    curLine = trpdata{3} == year   &   trpdata{4} == month   &   trpdata{5} == day   &   trpdata{6} == hour   &   trpdata{7} == minu   &   abs(sec - trpdata{8}) <= 0.1   &   strcmpi(trpdata{9}, antenna(ist).name)   &   strcmpi(trpdata{1},sourceNames(isc));
    
%                         if sum(curLine)~=1
%                             fprintf('No (or more than one) trp line found for scan %1.0f (station %s)\nTry to find offset of 0.5 -> some intensives...', isc, deblank(antenna(ist).name));
%                             % IF we are here, something might be wrong with trp files - usually it's something with mjd2date (because of rounding)
%                             curLine=scanYrs(isc)==trpdata{3} & scanMonths(isc)==trpdata{4} & scanDays(isc)==trpdata{5} & scanHrs(isc)==trpdata{6} & scanMinutes(isc)==trpdata{7} & abs(scanSecs(isc)-trpdata{8})==0.5 & strcmpi(trpdata{9}, antenna(ist).name);
%                             if sum(curLine)==1;
%                                 fprintf(' FOUND!\n');
%                             else
%                                 keyboard;
%                                 fprintf(' NOT FOUND!\nUsing ECMWFh*VMF1h instead\n\n');
%                                 vmf1NeededLog=1;
%                             end
%                         end
%
%            %return if finished
%            
%                         if vmf1NeededLog==0 % if there is exactly one line (what we want)
%                             if vahabsRaytracingFiles==1
%                                 scan(isc).stat(ist).trop =trpdata{14}(curLine);
%                                 vmf1NeededLog=1;
%                                 % the rest does not exist (load vmf1 later)
%                             else
%                                 scan(isc).stat(ist).trop  = trpdata{14}(curLine);   % slant path delay [sec]   
%                                 scan(isc).stat(ist).mfw   = trpdata{15}(curLine);   % wet mapping function
%                                 scan(isc).stat(ist).temp  = trpdata{13}(curLine);   % temperature [°C]
%                                 scan(isc).stat(ist).pres  = trpdata{12}(curLine);   % pressure [hPa]
%                             end
%                         end
            
    % determine the number of matching trp entries
    nr_matches = sum(curLine);


    % test how many matching trp entries have been found and initiate according processes

    % no matching entry
    if nr_matches == 0
        % display message
        fprintf('No trp entry has been found for scan %s for station %s! ',num2str(isc),strtrim(antenna(ist).name));

        % No data assignments necessary.
        % Note: Already existing values of pressure and temperature are preserved as they are
        %       needed by VieVS.
        % Backup of missing trp data by mapping function and according message about
        % backup will be handled in vie_mod.

    % 1 matching entry as it should always be
    elseif nr_matches == 1
        % assign the data
        scan(isc).stat(ist).trop  = trpdata{14}(curLine);   % slant path delay [sec]   
        scan(isc).stat(ist).mfw   = trpdata{15}(curLine);   % wet mapping function
        scan(isc).stat(ist).temp  = trpdata{13}(curLine);   % temperature [°C]
        scan(isc).stat(ist).pres  = trpdata{12}(curLine);   % pressure [hPa]

    % more than 1 matching entry
    elseif nr_matches > 1
        % display message
        fprintf('More than 1 matching trp entry has been found for scan %s for station %s! Only the first matching entry is considered!\n',num2str(isc),strtrim(antenna(ist).name));

        % assign the data of all matching entries
        scan(isc).stat(ist).trop  = trpdata{14}(curLine);   % slant path delay [sec]   
        scan(isc).stat(ist).mfw   = trpdata{15}(curLine);   % wet mapping function
        scan(isc).stat(ist).temp  = trpdata{13}(curLine);   % temperature [°C]
        scan(isc).stat(ist).pres  = trpdata{12}(curLine);   % pressure [hPa]

        % reduce the data to the first matching trp entry
        scan(isc).stat(ist).trop  = scan(isc).stat(ist).trop(1);   % slant path delay [sec]   
        scan(isc).stat(ist).mfw   = scan(isc).stat(ist).mfw(1);   % wet mapping function
        scan(isc).stat(ist).temp  = scan(isc).stat(ist).temp(1);   % temperature [°C]
        scan(isc).stat(ist).pres  = scan(isc).stat(ist).pres(1);   % pressure [hPa]

    end % test of different number of matches

end % if we have trp file