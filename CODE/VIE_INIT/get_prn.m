% ************************************************************************
%   Description:
%	Finds the corresponding prn and fso name based on the name and the id 
%   of a satellite in the NGS file.
%
%   Input:										
%     name             name of the satellite (from NGS file)
%     id               id of the satellite (from NGS file)
% 
%   Output:
%     prn              prn of the satellite (for sp3 file)
%     fso			   name/number for satellite used in fso file
%					   nn: GPS
%                      1nn: GLONASS
%                      2nn: Galileo
%                      3nn: SBAS
%                      4nn: BeiDou
%                      5nn: QZSS
% 
%   External calls: 	
%       
%       
%   Coded for VieVS: 
%   12 Jan 2025 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [prn, fso] = get_prn(name, id)
    
    Galileo = {'GSAT0101','37846U','11060A','E11';
    'GSAT0102','37847U','11060B','E12';
    'GSAT0103','38857U','12055A','E19';
    'GSAT0104','38858U','12055B','E20';
    'GSAT0201','40128U','14050A','E18';
    'GSAT0202','40129U','14050B','E14';
    'GSAT0203','40544U','15017A','E26';
    'GSAT0204','40545U','15017B','E22';
    'GSAT0205','40889U','15045A','E24';
    'GSAT0206','40890U','15045B','E30';
    'GSAT0209','41174U','15079A','E09';
    'GSAT0208','41175U','15079B','E08';
    'GSAT0211','41549U','16030A','E02';
    'GSAT0210','41550U','16030B','E01';
    'GSAT0207','41859U','16069A','E07';
    'GSAT0212','41860U','16069B','E03';
    'GSAT0213','41861U','16069C','E04';
    'GSAT0214','41862U','16069D','E05';
    'GSAT0215','43055U','17079A','E21';
    'GSAT0216','43056U','17079B','E25';
    'GSAT0217','43057U','17079C','E27';
    'GSAT0218','43058U','17079D','E31';
    'GSAT0221','43564U','18060A','E15';
    'GSAT0222','43565U','18060B','E33';
    'GSAT0219','43566U','18060C','E36';
    'GSAT0220','43567U','18060D','E13';
    'GSAT0223','49809U','21116A','E34';
    'GSAT0224','49810U','21116B','E10';
    'GSAT0225','59598U','24079A','E29';
    'GSAT0227','59600U','24079C','E06';
    'GSAT0232','61182U','24167A','E16';
    'GSAT0226','61183U','24167B','E23'};
    Galileo = cell2table(Galileo, 'VariableNames', {'Name','NORAD','ID','PRN'});

    Sentinel = {'SENTI-6A', '46984C', '20086A', 'L40'};
    Sentinel = cell2table(Sentinel, 'VariableNames', {'Name','NORAD','ID','PRN'});

    Lageos = {'LAGEOS-1', '08820U', '76039A', 'L51'};
    Lageos = cell2table(Lageos, 'VariableNames', {'Name','NORAD','ID','PRN'});

    Genesis = {'GEN-01', '11111U', '22001A', 'L01'};
    Genesis = cell2table(Genesis, 'VariableNames', {'Name','NORAD','ID','PRN'});

    if startsWith(name, "GSAT", 'IgnoreCase', true)
            prn = checkEntry(Galileo, false, id, name(1:8));
            % Galileo: PRN z.B. 'E11' -> number 11, FSO = 200 + 11
            num = sscanf(prn, '%*1c%d'); 
            if isempty(num)
                fso = "";
            else
                fso = num + 200;
            end
    
        elseif startsWith(name, "SENT", 'IgnoreCase', true) || startsWith(name, "SENTI", 'IgnoreCase', true)
            prn = checkEntry(Sentinel, true, id, name);
            fso = "";  
        elseif startsWith(name, "LAGEOS", 'IgnoreCase', true)
            prn = checkEntry(Lageos, true, id, name);
            fso = "";
        elseif startsWith(name, "GEN", 'IgnoreCase', true)
            prn = checkEntry(Genesis, true, id, name);
            fso = "";
        else
            prn = "";
            fso = "";
     end

end


function prnOut = checkEntry(tbl, exactNameMatch, id, name)
    if exactNameMatch
        maskName = strcmp(tbl.Name, name);
    else
        maskName = startsWith(tbl.Name, name, 'IgnoreCase', true);
    end
    maskID = strcmp(tbl.NORAD, id) | strcmp(tbl.ID, id);

    if ~any(maskName) && ~any(maskID)
        error("Keine Übereinstimmung in der Liste. Name oder ID nicht gefunden.");
    end
    maskBoth = maskName & maskID;
    if ~any(maskBoth)
        error("Name und ID stimmen nicht überein (Liste vorhanden, aber kein passender Eintrag).");
    end
    prnOut = string(tbl.PRN{find(maskBoth,1)});
end
