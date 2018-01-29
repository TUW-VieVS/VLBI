function [RQ,DRQDRA,DRQDDE] = sourcevec(decl,ra)
 
% vectorised
% Input:    declination       [:,1]
%           right ascension   [:,1]
%
% Output:   RQ                source vector [:,3] barycentric/ICRS
%           DRQDRA            partial derivative after ra [:,3]
%           DRQDDE            partial derivative after de [:,3]
%
      SID = sin(decl);             
      COD = cos(decl);             
      SIR = sin(ra);
      COR = cos(ra);

      RQ(:,1) = COD.* COR;           % cos(de) cos(ra)
      RQ(:,2) = COD.* SIR;           % cos(de) sin(ra)
      RQ(:,3) = SID;                 %     sin(de) 
     
      DRQDRA(:,1) = -COD.* SIR;
      DRQDRA(:,2) =  COD.* COR;
      DRQDRA(:,3) =  0;
      
      DRQDDE(:,1) = -SID.* COR;
      DRQDDE(:,2) = -SID.* SIR;
      DRQDDE(:,3) =  COD;
end
      