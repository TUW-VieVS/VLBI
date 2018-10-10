% Function to read additional TRF data introduced by user from textfiles and mat-file.
%
% INPUT
%  ns_codes			main loop file of creating superstation file (mk_superstatFile.m)
%  userOwnTrfData   filename of TRF data
% 
% OUTPUT
%  ns_codes(updated)
%
% CREATED
%  11.05.2016 A. Girdiuk
% 07 Jun 2017, H. Krasna: bug fixed


    % format / stationname x y z vx vy vz epoch start end
    % for all lines of the vieTrf file
    
    function [ns_codes] = ownTRFread(ns_codes,userOwnTrfData)

    
    for iLine=1:length(userOwnTrfData{1})

        % find index of current station (try also _ instead of ' '  and
        % vice-versa)
%         iStat=strcmp({ns_codes.name}, userOwnTrfData{1}{iLine})|...
%             strcmp({ns_codes.name}, strrep(userOwnTrfData{1}{iLine}, ' ', '_'))|...
%             strcmp({ns_codes.name}, strrep(userOwnTrfData{1}{iLine}, '_', ' '));
        
        iStat=strcmp({ns_codes.name}, userOwnTrfData{1}(iLine,:))|...
            strcmp({ns_codes.name}, strrep(userOwnTrfData{1}(iLine,:), ' ', '_'))|...
            strcmp({ns_codes.name}, strrep(userOwnTrfData{1}(iLine,:), '_', ' '));
        
        
        % if not found -> append to end of ns_codes
        if sum(iStat)==0
            indInNsCodes=size(ns_codes,2)+1;
        else
            indInNsCodes=find(iStat);
        end

        % for all indices where the name was found (always? only one!)
        for iStatFound=1:length(indInNsCodes)
            curInd=indInNsCodes(iStatFound);

            % +++ get break number +++
            % if we append station to ns_codes, it's also a new break
            if curInd>size(ns_codes,2)
                curBreak=1;
            % if there exists ns_codes(x).vieTrf (means there is already a
            % break
            elseif isfield(ns_codes(curInd), 'userOwnTrf')
                if ~isempty(ns_codes(curInd).userOwnTrf)
                    curBreak=size(ns_codes(curInd).userOwnTrf.break,2)+1;
                else
                    curBreak=1;
                end
            else
                curBreak=1;
            end
            % --- get break number ---

            % write coordinates and velocities to ns_codes
            if ~isempty(userOwnTrfData{9})
                ns_codes(curInd).userOwnTrf.break(curBreak).start=userOwnTrfData{9}(iLine);
            end
            if ~isempty(userOwnTrfData{10})
                ns_codes(curInd).userOwnTrf.break(curBreak).end=userOwnTrfData{10}(iLine);
            end
            if ~isempty(userOwnTrfData{5})
                ns_codes(curInd).userOwnTrf.break(curBreak).vx=userOwnTrfData{5}(iLine);
            end
            if ~isempty(userOwnTrfData{6})
                ns_codes(curInd).userOwnTrf.break(curBreak).vy=userOwnTrfData{6}(iLine);
            end
            if ~isempty(userOwnTrfData{7})
                ns_codes(curInd).userOwnTrf.break(curBreak).vz=userOwnTrfData{7}(iLine);
            end

            ns_codes(curInd).userOwnTrf.break(curBreak).epoch=userOwnTrfData{8}(iLine);
            ns_codes(curInd).userOwnTrf.break(curBreak).x=userOwnTrfData{2}(iLine);
            ns_codes(curInd).userOwnTrf.break(curBreak).y=userOwnTrfData{3}(iLine);
            ns_codes(curInd).userOwnTrf.break(curBreak).z=userOwnTrfData{4}(iLine);
                        
%            ns_codes(curInd).name=''; % needed for cellfun strfind over
%            .name ???????????????????
 %           ns_codes(curInd).code='';

        end    
    end
