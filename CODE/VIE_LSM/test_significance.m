% ************************************************************************
%   Description:
%   Function to test certain estimates for significance
%
%
%   Input:                                                                              
%   x_ structure        (contains the estimates)
%   opt_ structure      (contains the source names)
%   significance_sigma  (level of significance test)
% 
%   Coded for VieVS: 
%   23 Sep 2016 by David Mayer
%
%   Revisions: 
%    - 2016-10-10, A. Girdiuk: adjusted for sources estimates number more than one
%    - 2016-10-10, A. Hellerschmied: Check, if "opt_.source" exists.
%    - 2017-04-04, D. Mayer: added information about datum
%    - 2017-09-11, D. Mayer: compatibility with source- and stationwise parameterization 
% ************************************************************************


function test_significance(x_, opt_, significance_sigma)
% if ~isempty(x_.xpol.col)
%     for i = 1: length(x_.xpol.col)
%         if abs(x_.xpol.val(i)) > significance_sigma*x_.xpol.mx(i)
%             fprintf(1,'xpol estimate is significant on %2.0f sigma level\n',significance_sigma);
%             break
%         end
%     end
% end
% 
% if ~isempty(x_.ypol.col)
%     for i = 1: length(x_.ypol.col)
%         if abs(x_.ypol.val(i)) > significance_sigma*x_.ypol.mx(i)
%             fprintf(1,'ypol estimate is significant on %2.0f sigma level\n',significance_sigma);
%             break
%         end
%     end
% end
% 
% if ~isempty(x_.dut1.col)
%     for i = 1: length(x_.dut1.col)
%         if abs(x_.dut1.val(i)) > significance_sigma*x_.dut1.mx(i)
%             fprintf(1,'dut1 estimate is significant on %2.0f sigma level\n',significance_sigma);
%             break
%         end
%     end
% end
% 
% if ~isempty(x_.nutdx.col)
%     for i = 1: length(x_.nutdx.col)
%         if abs(x_.nutdx.val(i)) > significance_sigma*x_.nutdx.mx(i)
%             fprintf(1,'nutdx estimate is significant on %2.0f sigma level\n',significance_sigma);
%             break
%         end
%     end
% end
% 
% if ~isempty(x_.nutdy.col)
%     for i = 1: length(x_.nutdy.col)
%         if abs(x_.nutdy.val(i)) > significance_sigma*x_.nutdy.mx(i)
%             fprintf(1,'nutdy estimate is significant on %2.0f sigma level\n',significance_sigma);
%             break
%         end
%     end
% end
try
    for i = 1:length(x_.antenna)
        if opt_.stat(i).nnt_inc && opt_.stat(i).nnr_inc
            nnr_nnt_string = '(in NNR/NNT)';
        elseif opt_.stat(i).nnt_inc && ~opt_.stat(i).nnr_inc
            nnr_nnt_string = '(in NNT but not in NNR)';
        elseif ~opt_.stat(i).nnt_inc && opt_.stat(i).nnr_inc
            nnr_nnt_string = '(not in NNT but in NNR)';
        elseif ~opt_.stat(i).nnt_inc && ~opt_.stat(i).nnr_inc
            nnr_nnt_string = '(not in NNR/NNT)';
        end

        if length(x_.coorx(i).col) == 1 && ~isempty(x_.coorx(i).col) && abs(x_.coorx(i).val) > significance_sigma*x_.coorx(i).mx
            fprintf(1,'X coordinate estimate (%5.2f +- %5.2f cm) of %s %s is significant on %2.0f sigma level\n', x_.coorx(i).val, x_.coorx(i).mx, x_.antenna(i).name, nnr_nnt_string, significance_sigma);
        end
        if length(x_.coory(i).col) == 1 && ~isempty(x_.coory(i).col) && abs(x_.coory(i).val) > significance_sigma*x_.coory(i).mx
            fprintf(1,'Y coordinate estimate (%5.2f +- %5.2f cm) of %s %s is significant on %2.0f sigma level\n', x_.coory(i).val, x_.coory(i).mx, x_.antenna(i).name, nnr_nnt_string, significance_sigma);
        end
        if length(x_.coorz(i).col) == 1 && ~isempty(x_.coorz(i).col) && abs(x_.coorz(i).val) > significance_sigma*x_.coorz(i).mx
            fprintf(1,'Z coordinate estimate (%5.2f +- %5.2f cm) of %s %s is significant on %2.0f sigma level\n', x_.coorz(i).val, x_.coorz(i).mx, x_.antenna(i).name, nnr_nnt_string, significance_sigma);
        end
    end
    if isfield(opt_, 'source')
        for i = 1:length({opt_.source.name}) 
            if length(x_.soude(i).val) == 1 && ~isempty(x_.soude(i).val) && any(abs(x_.soude(i).val) > significance_sigma*x_.soude(i).mx)
                for i_est = 1:length(x_.soude(i).val)
                    if abs(x_.soude(i).val(i_est)) > significance_sigma*x_.soude(i).mx(i_est)
                        fprintf(1,'DEC coordinate estimate (%6.3f +- %6.3f mas) of %s is significant on %2.0f sigma level\n', x_.soude(i).val(i_est), x_.soude(i).mx(i_est), opt_.source(i).IERSname, significance_sigma);
                    end
                end
            end
            if length(x_.soura(i).val) == 1 && ~isempty(x_.soura(i).val) && any(abs(x_.soura(i).val) > significance_sigma*x_.soura(i).mx)
                for i_est = 1: length(x_.soura(i).val)
                    if abs(x_.soura(i).val(i_est)) > significance_sigma*x_.soura(i).mx(i_est)
                        fprintf(1,'RA coordinate estimate (%6.3f +- %6.3f ms) of %s is significant on %2.0f sigma level\n', x_.soura(i).val(i_est), x_.soura(i).mx(i_est), opt_.source(i).IERSname, significance_sigma);
                    end
                end
            end
        end
    end
    BCOsl = 3;
    BCOcount = 0;
    if isfield(x_, 'bdclko')
        for i=1:length(x_.bdclko)
           if abs(x_.bdclko(i).val)>BCOsl*x_.bdclko(i).mx
               BCOcount = BCOcount + 1;
               fprintf(1,'BCO estimate (%6.3f +- %6.3f cm) between %s and %s is significant on %2.0f sigma level\n', x_.bdclko(i).val, x_.bdclko(i).mx, x_.bdclko(i).namest1, x_.bdclko(i).namest2, BCOsl);                
           end
        end
        if BCOcount>0
            fprintf(1,'----------\n', opt_.fixed_clock);   
            fprintf(1,'BCOs block for OPT file with CLOCK REFERENCE %s:\n', opt_.fixed_clock);   
            fprintf(1,'+BASELINE-DEPENDENT CLOCK OFFSET\n'); 
            for i=1:length(x_.bdclko)
                if abs(x_.bdclko(i).val)>BCOsl*x_.bdclko(i).mx         
                    fprintf(1,'%s  %s \n',  x_.bdclko(i).namest1, x_.bdclko(i).namest2);                
                end
            end
            fprintf(1,'-BASELINE-DEPENDENT CLOCK OFFSET\n'); 
        end
    end
catch
	fprintf(1,'Significance test failed.');
end
fprintf(1,'\n');

