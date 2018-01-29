function [vmf1h,vmf1w] = vmf1_ht (ah,aw,dmjd,dlat,ht,zd)
      
% C     !!! This is the version with height correction !!!  
% C     !!! It has to be used with the grid !!!    
% C
% C     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
% C     Reference: Boehm, J., B. Werl, H. Schuh (2006), 
% C     Troposphere mapping functions for GPS and very long baseline interferometry 
% C     from European Centre for Medium-Range Weather Forecasts operational analysis data,
% C     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
% C
% C     Please mind that the coefficients in this paper are wrong. The corrected version of
% C     the paper can be found at: 
% C     http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf
% C
% C     input data
% C     ----------
% C     ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
% C     aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)  
% C     dmjd: modified julian date
% C     dlat: ellipsoidal latitude in radians
% C     ht:   ellipsoidal height in meter
% C     zd:   zenith distance in radians
% C
% C     output data
% C     -----------
% C     vmf1h: hydrostatic mapping function
% C     vmf1w: wet mapping function
% C
% C     Johannes Boehm, 2005 October 2
% C     Rev 2011 July 21: latitude -> ellipsoidal latitude
% C
% 
%       implicit double precision (a-h,o-z)
% 
%       pi = 3.14159265359d0
%       
% C     reference day is 28 January
% C     this is taken from Niell (1996) to be consistent

      doy = dmjd  - 44239.d0 + 1 - 28;
      
      bh = 0.0029;
      c0h = 0.062;
      if (dlat<0)      %   ! southern hemisphere
          phh  = pi;
          c11h = 0.007;
          c10h = 0.002;
      else             %   ! northern hemisphere
          phh  = 0.d0;
          c11h = 0.005;
          c10h = 0.001;
      end
          
      ch = c0h + ((cos(doy/365.25d0*2.d0*pi + phh)+1.d0)*c11h/2.d0 ... 
          + c10h)*(1.d0-cos(dlat));

      sine   = sin(pi/2.d0 - zd);
      beta   = bh/( sine + ch  );
      gamma  = ah/( sine + beta);
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)));
      vmf1h   = topcon/(sine+gamma);

% C  height correction for hydrotatic part [Niell, 1996]     
      a_ht = 2.53d-5;
      b_ht = 5.49d-3;
      c_ht = 1.14d-3;
      hs_km  = ht/1000.d0;
      beta         = b_ht/( sine + c_ht);
      gamma        = a_ht/( sine + beta);
      topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)));
      ht_corr_coef = 1.d0/sine - topcon/(sine + gamma);
      ht_corr      = ht_corr_coef * hs_km;
      vmf1h        = vmf1h + ht_corr;
      
      bw = 0.00146;
      cw = 0.04391;
      beta   = bw/( sine + cw );
      gamma  = aw/( sine + beta);
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));
      vmf1w   = topcon/(sine+gamma);
      

