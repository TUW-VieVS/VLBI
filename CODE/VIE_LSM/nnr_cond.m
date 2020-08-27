% ************************************************************************
%   Description:
%   function to form NNR condition equations for datum definition of CRF
%
%   Reference: 
%
%   Input:	
%       'ns'         (1,1)                      number of sources
%       'ra'         (1,ns)                     apriori RA of all sources in the session
%       'de'         (1,ns)                     apriori De of all sources in the session
%       'opt'        structure array            (for info. /DOC/opt.doc)
%       'sum_dj'     (1,number of models)       total number of estimates for each included model
%       'N'          (sum_dj(end),sum_dj(end))  datum free normal equation matrix         
%
%   Output:
%       'N' (sum_dj(end)+cond,sum_dj(end)+cond) normal equation matrix with NNR condition equations        
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   18 Jun 2014 by Hana Krasna
%
%   Revision: 
%   24 Sep 2014 by Hana Krasna; bug fixed occuring if all other parameters
%   have to be fixed to their a priori values
% ************************************************************************
function [N] = nnr_cond(ns,ra,de,opt,sum_dj,N)
             Bra = []; Bde = []; B = [];
    
                                      
	for isou = 1 : ns
        
        %--------------------
        bra(3,1) = 0; bde = bra;
        %--------------------
        bra(1,1) = tan(de(isou))*cos(ra(isou)); 
        bra(2,1) = tan(de(isou))*sin(ra(isou)); 
        bra(3,1) = -1; 
        %--------------------
        bde(1,1) = -sin(ra(isou)); 
        bde(2,1) =  cos(ra(isou)); 
        bde(3,1) = 0; 
        %--------------------
  
        if opt.source(isou).nnr_inc == 0;
           clear bra, bra(3,1) = 0;
           clear bde, bde = bra;
        end

        
        Bra = horzcat(Bra,bra); Bde = horzcat(Bde,bde);
        
        clear bra bde
    end    
   B1 = horzcat(Bra,Bde); 
   
   B1(~any(B1,2),:) = [];
   if sum_dj(11) ~= 0
       B2(size(B1,1),sum_dj(11)) =  0; 
       if size(N,2)-sum_dj(13) ~= 0
           B3(size(B1,1),size(N,2)-sum_dj(13)) =  0; 
       else
           B3=[];
       end
       B = horzcat(B2,B1,B3);
   elseif sum_dj(11) == 0 && size(N,2)-sum_dj(13)~=0
       B3(size(B1,1),size(N,2)-sum_dj(13)) =  0;     
       B = horzcat(B1,B3);
   else
       B=B1;
   end

   
   K(size(B,1),size(B,1)) = 0; 
   N2 = horzcat(vertcat(N,B),vertcat(B',K)); 
   clear N
   N = N2; 

end     
