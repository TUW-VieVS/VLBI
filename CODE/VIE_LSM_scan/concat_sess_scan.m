%**************************************************************************
%   Description:
%    Function to create from a design matrix per scan a its corresponding 
%    design matrix for the whole session
%
%   Input:
%     'scan'                 structure array      (for info. /DOC/scan.doc)
%     'iscan'                (1,1)                 Number of the scan that is being analysed
%     'A'                    structure array       design matrix of real observations
%     'opt'                  structure array      (for info. /DOC/opt.doc)
%     't'                    (1,ntim)              minutes from midnight to the time of the scans
%     't_first'              (1,na)                First observation of antenna
%     'na'                   (1,1)                 Number of antennas
%     'ns'                   (1,1)                 Number of sources
%     'num_inter_....'       matrix                Number of estimation intervals
%
%   Output:
%     'A_session'           structure array        Design matrix of real observations
%
%   Coded for VieVS: 
%   11 oct 2012 by Claudia Tierno Ros
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
%**************************************************************************


function  [A_session]=concat_sess_scan(scan,iscan,A,opt,t,t_first,na,ns...
        ,num_inter_xpol,num_inter_ypol,num_inter_dut1,num_inter_nutdx...
        ,num_inter_nutdy,num_inter_egr,num_inter_ngr,num_inter_clk,num_inter_zwd,per_source,num_inter_sou,t_first_sou)
    
    
     
   % Make 0-matrixes with their corresponding size
    A_session(6).sm= zeros(scan(iscan).nobs,num_inter_xpol+1);
    A_session(7).sm= zeros(scan(iscan).nobs,num_inter_ypol+1);
    A_session(8).sm= zeros(scan(iscan).nobs,num_inter_dut1+1);
    A_session(9).sm= zeros(scan(iscan).nobs,num_inter_nutdx+1);
    A_session(10).sm= zeros(scan(iscan).nobs,num_inter_nutdy+1);

    A_session(4).sm= zeros(scan(iscan).nobs,sum(num_inter_ngr)+na);
    A_session(5).sm= zeros(scan(iscan).nobs,sum(num_inter_egr)+na);
    
   
    A_session(1).sm= zeros(scan(iscan).nobs,sum(num_inter_clk)+na);
    A_session(3).sm=zeros(scan(iscan).nobs,sum(num_inter_zwd)+na);
    
    if opt.pw_sou == 1
    A_session(11).sm= zeros(scan(iscan).nobs,sum(num_inter_sou)+ns);
    A_session(12).sm= zeros(scan(iscan).nobs,sum(num_inter_sou)+ns);
    end
    
    if opt.est_sourceNNR == 1
    A_session(11).sm= zeros(scan(iscan).nobs,sum(num_inter_sou)+ns);
    A_session(12).sm= zeros(scan(iscan).nobs,sum(num_inter_sou)+ns);
    end
    
  
    % Number of column correspoding to the observation interval
    kxpol=(floor(t(iscan)/opt.xpol.int)-floor(t(1)/opt.xpol.int)+1);  
    kypol=(floor(t(iscan)/opt.ypol.int)-floor(t(1)/opt.ypol.int)+1);
    kdut1=(floor(t(iscan)/opt.dut1.int)-floor(t(1)/opt.dut1.int)+1);
    knutdx=(floor(t(iscan)/opt.nutdx.int)-floor(t(1)/opt.nutdx.int)+1);
    knutdy=(floor(t(iscan)/opt.nutdy.int)-floor(t(1)/opt.nutdy.int)+1);
  
  for j=1:scan(iscan).nobs
             
           

        % Xpol
        %kxpol=(floor(t(iscan)/opt.xpol.int)-floor(t(1)/opt.xpol.int)+1);  % Number of column correspoding to the observation interval
        A_session(6).sm(j,kxpol)=A(6).sm(j,1);                            % Put value in its corresponding column
        A_session(6).sm(j,kxpol+1)=A(6).sm(j,2);                          % Put value in its corresponding column

        % Ypol
        %kypol=(floor(t(iscan)/opt.ypol.int)-floor(t(1)/opt.ypol.int)+1);  % Number of column correspoding to the observation interval
        A_session(7).sm(j,kypol)=A(7).sm(j,1);                            % Put value in its corresponding column
        A_session(7).sm(j,kypol+1)=A(7).sm(j,2);                          % Put value in its corresponding column

        % dut1
        %kdut1=(floor(t(iscan)/opt.dut1.int)-floor(t(1)/opt.dut1.int)+1);  % Number of column correspoding to the observation interval
        A_session(8).sm(j,kdut1)=A(8).sm(j,1);                            % Put value in its corresponding column
        A_session(8).sm(j,kdut1+1)=A(8).sm(j,2);                          % Put value in its corresponding column

        % nutdx
        %knutdx=(floor(t(iscan)/opt.nutdx.int)-floor(t(1)/opt.nutdx.int)+1); % Number of column correspoding to the observation interval
        A_session(9).sm(j,knutdx)=A(9).sm(j,1);                             % Put value in its corresponding column
        A_session(9).sm(j,knutdx+1)=A(9).sm(j,2);                           % Put value in its corresponding column

        % nutdy
        %knutdy=(floor(t(iscan)/opt.nutdy.int)-floor(t(1)/opt.nutdy.int)+1); % Number of column correspoding to the observation interval
        A_session(10).sm(j,knutdy)=A(10).sm(j,1);                           % Put value in its corresponding column
        A_session(10).sm(j,knutdy+1)=A(10).sm(j,2);                         % Put value in its corresponding column

        
        
        % rate and quadratic term of clock
        A_session(2).sm=A(2).sm; %A_session(2).sm(j,:)=A(2).sm(j,:); 
     
        % Xstat
        A_session(13).sm=A(13).sm;

        % Ystat
        A_session(14).sm=A(14).sm;
 
        % Zstat
        A_session(15).sm=A(15).sm;
  
     
        XX=scan(iscan).obs(j).i1; % number of first station
        YY=scan(iscan).obs(j).i2; %number of second station
        
               
        
        % first station
        if XX==1
            
            % ngr
            ncol=(floor(t(iscan)/opt.int_ngr)-floor(t_first(XX)/opt.int_ngr)+1);  % Number of column correspoding to the observation interval
            A_session(4).sm(j,ncol)=A(4).sm(j,XX*2-1);                            % Put value in its corresponding column
            A_session(4).sm(j,ncol+1)=A(4).sm(j,XX*2);                            % Put value in its corresponding column
            
            % egr
            ncol=(floor(t(iscan)/opt.int_egr)-floor(t_first(XX)/opt.int_egr)+1);  % Number of column correspoding to the observation interval
            A_session(5).sm(j,ncol)=A(5).sm(j,XX*2-1);                            % Put value in its corresponding column
            A_session(5).sm(j,ncol+1)=A(5).sm(j,XX*2);                            % Put value in its corresponding column
            
            % zwd
            ncol=(floor(t(iscan)/opt.int_zwd)-floor(t_first(XX)/opt.int_zwd)+1);  % Number of column correspoding to the observation interval
            A_session(3).sm(j,ncol)=A(3).sm(j,XX*2-1);                            % Put value in its corresponding column
            A_session(3).sm(j,ncol+1)=A(3).sm(j,XX*2);                            % Put value in its corresponding column
            
            % pwl clock offsets
            ncol=(floor(t(iscan)/opt.int_clk)-floor(t_first(XX)/opt.int_clk)+1);  % Number of column correspoding to the observation interval
            A_session(1).sm(j,ncol)=A(1).sm(j,XX*2-1);                            % Put value in its corresponding column
            A_session(1).sm(j,ncol+1)=A(1).sm(j,XX*2);                            % Put value in its corresponding column
        
        else
            
             % ngr
             ni1=(floor(t(iscan)/opt.int_ngr)-floor(t_first(XX)/opt.int_ngr)+1);  % Observation interval
             ncol=sum(num_inter_ngr(1:XX-1))+ni1+(XX-1);                          % Number of column correspoding to the observation interval
             A_session(4).sm(j,ncol)=A(4).sm(j,XX*2-1);                           % Put value in its corresponding column
             A_session(4).sm(j,ncol+1)=A(4).sm(j,XX*2);                           % Put value in its corresponding column
             
             % egr
             ni1=(floor(t(iscan)/opt.int_egr)-floor(t_first(XX)/opt.int_egr)+1);  % Observation interval
             ncol=sum(num_inter_egr(1:XX-1))+ni1+(XX-1);                          % Number of column correspoding to the observation interval
             A_session(5).sm(j,ncol)=A(5).sm(j,XX*2-1);                           % Put value in its corresponding column
             A_session(5).sm(j,ncol+1)=A(5).sm(j,XX*2);                           % Put value in its corresponding column
             
             % zwd
             ni1=(floor(t(iscan)/opt.int_zwd)-floor(t_first(XX)/opt.int_zwd)+1);  % Observation interval
             ncol=sum(num_inter_zwd(1:XX-1))+ni1+(XX-1);                          % Number of column correspoding to the observation interval
             A_session(3).sm(j,ncol)=A(3).sm(j,XX*2-1);                           % Put value in its corresponding column
             A_session(3).sm(j,ncol+1)=A(3).sm(j,XX*2);                           % Put value in its corresponding column
             
             % pwl clock offsets
             ni1=(floor(t(iscan)/opt.int_clk)-floor(t_first(XX)/opt.int_clk)+1);  % Observation interval
             ncol=sum(num_inter_clk(1:XX-1))+ni1+(XX-1);                          % Number of column correspoding to the observation interval
             A_session(1).sm(j,ncol)=A(1).sm(j,XX*2-1);                           % Put value in its corresponding column
             A_session(1).sm(j,ncol+1)=A(1).sm(j,XX*2);                           % Put value in its corresponding column
        
        end
        

        % second station       
        if YY==1
            
             % ngr
             ncol2=(floor(t(iscan)/opt.int_ngr)-floor(t_first(YY)/opt.int_ngr)+1);  % Number of column correspoding to the observation interval
             A_session(4).sm(j,ncol2)=A(4).sm(j,YY*2-1);                            % Put value in its corresponding column
             A_session(4).sm(j,ncol2+1)=A(4).sm(j,YY*2);                            % Put value in its corresponding column
             
             % egr
             ncol2=(floor(t(iscan)/opt.int_egr)-floor(t_first(YY)/opt.int_egr)+1);  % Number of column correspoding to the observation interval
             A_session(5).sm(j,ncol2)=A(5).sm(j,YY*2-1);                            % Put value in its corresponding column
             A_session(5).sm(j,ncol2+1)=A(5).sm(j,YY*2);                            % Put value in its corresponding column
             
             % zwd
             ncol2=(floor(t(iscan)/opt.int_zwd)-floor(t_first(YY)/opt.int_zwd)+1);  % Number of column correspoding to the observation interval
             A_session(3).sm(j,ncol2)=A(3).sm(j,YY*2-1);                            % Put value in its corresponding column
             A_session(3).sm(j,ncol2+1)=A(3).sm(j,YY*2);                            % Put value in its corresponding column
             
             % pwl clock offsets
             ncol2=(floor(t(iscan)/opt.int_clk)-floor(t_first(YY)/opt.int_clk)+1);  % Number of column correspoding to the observation interval
             A_session(1).sm(j,ncol2)=A(1).sm(j,YY*2-1);                            % Put value in its corresponding column
             A_session(1).sm(j,ncol2+1)=A(1).sm(j,YY*2);                            % Put value in its corresponding column
             
        else
            
            % ngr
             ni2=(floor(t(iscan)/opt.int_ngr)-floor(t_first(YY)/opt.int_ngr)+1);  % Observation interval
             ncol2=sum(num_inter_ngr(1:YY-1))+ni2+(YY-1);                         % Number of column correspoding to the observation interval
             A_session(4).sm(j,ncol2)=A(4).sm(j,YY*2-1);                          % Put value in its corresponding column
             A_session(4).sm(j,ncol2+1)=A(4).sm(j,YY*2);                          % Put value in its corresponding column
             
            % egr
             ni2=(floor(t(iscan)/opt.int_egr)-floor(t_first(YY)/opt.int_egr)+1);  % Observation interval
             ncol2=sum(num_inter_egr(1:YY-1))+ni2+(YY-1);                         % Number of column correspoding to the observation interval
             A_session(5).sm(j,ncol2)=A(5).sm(j,YY*2-1);                          % Put value in its corresponding column
             A_session(5).sm(j,ncol2+1)=A(5).sm(j,YY*2);                          % Put value in its corresponding column
             
             % zwd
             ni2=(floor(t(iscan)/opt.int_zwd)-floor(t_first(YY)/opt.int_zwd)+1);  % Observation interval
             ncol2=sum(num_inter_zwd(1:YY-1))+ni2+(YY-1);                         % Number of column correspoding to the observation interval
             A_session(3).sm(j,ncol2)=A(3).sm(j,YY*2-1);                          % Put value in its corresponding column
             A_session(3).sm(j,ncol2+1)=A(3).sm(j,YY*2);                          % Put value in its corresponding column
             
             % pwl clock offsets
             ni2=(floor(t(iscan)/opt.int_clk)-floor(t_first(YY)/opt.int_clk)+1);  % Observation interval
             ncol2=sum(num_inter_clk(1:YY-1))+ni2+(YY-1);                         % Number of column correspoding to the observation interval
             A_session(1).sm(j,ncol2)=A(1).sm(j,YY*2-1);                          % Put value in its corresponding column
             A_session(1).sm(j,ncol2+1)=A(1).sm(j,YY*2);                          % Put value in its corresponding column
             
        end
        


        if opt.pw_sou == 1 % if pwl offsets are estimated
             
            ni1=(floor(t(iscan)/opt.sour_int_rade)-floor(t_first_sou(per_source.iso)/opt.sour_int_rade)+1); % Observation interval
             
            % Number of column correspoding to the observation interval
            if per_source.iso==1;
                 ncol=ni1;
            else
                 ncol=sum(num_inter_sou(1:per_source.iso-1))+ni1+((per_source.iso)-1);
            end
             
             % Put value in its corresponding column
             A_session(11).sm(j,ncol)=A(11).sm(j,1);
             A_session(11).sm(j,ncol+1)=A(11).sm(j,2);
             A_session(12).sm(j,ncol)=A(12).sm(j,1);
             A_session(12).sm(j,ncol+1)=A(12).sm(j,2);
             

        end
    
         
        if opt.est_sourceNNR==1
            ncol = scan(iscan).iso;
            A_session(11).sm=zeros(size(A(11).sm,1),ns);
            A_session(12).sm=zeros(size(A(12).sm,1),ns);

            A_session(11).sm(:,ncol)=A(11).sm';
            A_session(12).sm(:,ncol)=A(12).sm';
           
            
            
            
        end
        
        
 end