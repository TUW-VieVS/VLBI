%
% Computation of ionospheric contribution for dual band VLBI observations
% for the higher frequency
% I.e., ionospheric contribution for X-band in case of S/X observations
%
% Hana Krasna, 2025 Aug 10


function [ionoX, sigma_ionoX, qflag_ion,vs,vx] = iono_contribution_dualB(out_struct,ctau_s,ctau_x,csigma_tau_s,csigma_tau_x,parameter)
fileoutput=false;

tau_s=[ctau_s{:}]' .*1e9; % ns
tau_x=[ctau_x{:}]' .*1e9; % ns
sigma_tau_s=[csigma_tau_s{:}]' .*1e9; % ns
sigma_tau_x=[csigma_tau_x{:}]' .*1e9; % ns

qflag_ion = zeros(length(tau_x), 1);

[vx, sigma_vx, qflag_vx]= effective_freq(out_struct,'bX',parameter); % GHz
[vs, sigma_vs, qflag_vs]= effective_freq(out_struct,'bS',parameter); % GHz



if fileoutput
    if length(vx) ~= length(tau_x)
            fid = fopen(['sessions_different_length_ObsEdit.txt'],'a');
            fprintf(fid,'%s   %s \n', out_struct.head.Session.val, parameter.session_name);
            fclose(fid);
            %return
    end
end



% Calculate ionospheric CONTRIBUTION for taux
% contribution = -correction
ionoX = -((vs.^2) ./ (vs.^2 - vx.^2) .* (tau_s - tau_x)); % [ns]  

% Calculate partial derivatives with respect to each variable
partial_tau_s = -(vs.^2) ./ (vs.^2 - vx.^2); %(-)
partial_tau_x = +(vs.^2) ./ (vs.^2 - vx.^2);

partial_vs = (tau_s - tau_x) .* (2*vs.*vx.^2) ./ ((vs.^2-vx.^2).^2);
partial_vx = -(tau_s - tau_x) .* (2*vs.^2.*vx) ./ ((vs.^2-vx.^2).^2);



% % Propagation of uncertainty 
% sigma_ionoX = sqrt((partial_tau_x .* sigma_tau_x).^2 + ...
%                  (partial_tau_s .* sigma_tau_s).^2 + ...
%                  (partial_vx .* sigma_vx).^2 + ...
%                  (partial_vs .* sigma_vs).^2); % [ns]

sigma_ionoX = sqrt((partial_tau_x .* sigma_tau_x).^2 + ...
                 (partial_tau_s .* sigma_tau_s).^2 + ...
                 (sqrt(2)*1e-3)^2 + ...            % take constant 2ps²
                 (sqrt(2)*1e-3)^2); % [ns]


%% plots
% if isfield (out_struct.head, 'ExpName')
%       exper=out_struct.head.ExpName.val';
% else
%     exper=out_struct.head.Session.val';
% end
% f1=figure;
% plot( (partial_tau_x .* sigma_tau_x).*10^3,'.')
% hold on
% plot((partial_tau_s .* sigma_tau_s).*10^3,'.')
% hold on
% plot((partial_vx .* sigma_vx).*10^3,'.')
% hold on
% plot((partial_vs .* sigma_vs).*10^3,'.')
% hold off
% legend('taux','taus','vx','vs')
% ylabel('[ps]')
% ylabel(['partials x sigma [ps]'])
% title([exper ])
% print(f1,'-dpdf' ,'-r500',[exper '_partials_x_sigma']);
% 
% f2=figure;
% plot( sigma_ionoX.*10^3,'.')
% ylabel('[ps]')
% ylabel(['sigma iono [ps]'])
% title([exper ])
% print(f2,'-dpdf' ,'-r500',[exper '_sigmaIono']);
%%



% flags  
qflag_ion(isnan(ionoX) | ionoX == 0) = -1;
qflag_ion(isnan(sigma_ionoX) | sigma_ionoX == 0) = -1;
qflag_ion(qflag_vx == -1) = -1;
qflag_ion(qflag_vs == -1) = -1;

% set NaN values to zero, they are already flaged
ionoX(isnan(ionoX)) = 0;
sigma_ionoX(isnan(sigma_ionoX)) = 0;

sum_qflag_ion = sum(qflag_ion);

if sum_qflag_ion == 0
    fprintf(' - Same ionospheric delay flag (= 0) used for all scans! \n')
end

fprintf('\t Iono corr: directly calculated in VieVS \n')

end


