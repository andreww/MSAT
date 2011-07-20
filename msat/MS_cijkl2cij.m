% MS_CIJKL2CIJ - Convert from elastic tensor to Voigt elasticity matrix
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Converts between a 3x3x3x3 tensor representation and a 6x6 Voigt 
%     representation of anisotropic elasticity. 
%
%  % [C] = MS_cijkl2cij(CC)
%
% Usage: 
%     CC must be a rank 4 array with size (3,3,3,3). The
%     returned array C will be rank 2 with size (6,6). The elasticity
%     matrix representation is used as input for almost all MSAT functions.
%
% Notes:
%     Do not use this function for the elastic compliance as additional
%     terms are needed in this case.
%
% See also: MS_cij2cijkl

% (C) James Wookey and Andrew Walker, 2011

function [C] = MS_cijkl2cij(CC)
 C = zeros(6,6);
 for im=1:3
   for jm=1:3
      for km=1:3
         for lm=1:3
            if ( CC(im,jm,km,lm) ~= 0.0) 
               [iv,jv]=ijkl2ij_local(im,jm,km,lm) ;
               C(iv,jv) = CC(im,jm,km,lm);
            end
         end
      end
   end
 end
return

function [iv,jv] = ijkl2ij_local(ii,jj,kk,ll)   
   if (ii==1 & jj==1) 
      iv=1;
   end
   if (ii==1 & jj==2) 
      iv=6;
   end
   if (ii==1 & jj==3) 
      iv=5;
   end
   if (ii==2 & jj==1) 
      iv=6;
   end

   if (ii==2 & jj==2) 
      iv=2;
   end
   if (ii==2 & jj==3) 
      iv=4;
   end
   
   if (ii==3 & jj==1) 
      iv=5;
   end
   if (ii==3 & jj==2) 
      iv=4;
   end
   if (ii==3 & jj==3) 
      iv=3;
   end
   
   if (kk==1 & ll==1) 
      jv=1;
   end
   if (kk==1 & ll==2) 
      jv=6;
   end
   if (kk==1 & ll==3) 
      jv=5;
   end

   if (kk==2 & ll==1) 
      jv=6;
   end
   if (kk==2 & ll==2) 
      jv=2;
   end
   if (kk==2 & ll==3) 
      jv=4;
   end

   if (kk==3 & ll==1) 
      jv=5;
   end
   if (kk==3 & ll==2) 
      jv=4;
   end
   if (kk==3 & ll==3) 
      jv=3;
   end

return
