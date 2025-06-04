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

[numRowsGalileo, ~] = size(Galileo);

    if contains(name, 'GSAT') 
        idx1 = find(ismember(Galileo(:,1),name(1:8)));
        idx2 = find(ismember(Galileo(:,2),id));
        if isempty(idx1) || isempty(idx2)
            fprintf('%-10s %-10s %-10s %-10s\n', 'Name', 'NORAD', 'ID', 'PRN');
            for i = 1:numRowsGalileo
                fprintf('%-10s %-10s %-10s %-10s\n', Galileo{i, :});
            end
            error(' *** There is no satellite found which matches with an entry in Galileo satellite list (see above). Please check the name and id of the satellite used for scheduling.')
        elseif idx1 ~= idx2
            fprintf('%-10s %-10s %-10s %-10s\n', 'Name', 'NORAD', 'ID', 'PRN');
            for i = 1:numRowsGalileo
                fprintf('%-10s %-10s %-10s %-10s\n', Galileo{i, :});
            end
            error('The name of the satellite and the id of this satellite do not match with an entry in the Galileo satellite list (see above). Please check the name and id of the satellite used for scheduling.')
        else
            prn = char(string(Galileo(idx1,4)));
            fso = str2double(prn(2:end))+200;
        end
	else
		prn = '';
		fso = '';
    end
end