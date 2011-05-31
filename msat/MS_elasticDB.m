% Database of various elastic constants. 
% These are reference by a string uid, output is stiffness matrix (GPa) and density (kg/m3)

function [C,rh] = MS_elasticDB(uid)

   switch lower(uid)
      case {'olivine','ol'}
      % single crystal olivine at 0GPa from ab initio (Abrams 97)
         C = [320.5  68.1  71.6   0.0   0.0   0.0 ; ...
               68.1 196.5  76.8   0.0   0.0   0.0 ; ...
               71.6  76.8 233.5   0.0   0.0   0.0 ; ...
                0.0   0.0   0.0  64.0   0.0   0.0 ; ...
                0.0   0.0   0.0   0.0  77.0   0.0 ; ...
                0.0   0.0   0.0   0.0   0.0  78.7 ];
         rh = 3355 ;   
      case {'albite','alb'}
      % single crystal ablite at 0GPa from ab initio (Brown et al PCM 2006)
         C = [ 69.9  34.0  30.8   5.1   -2.4   -0.9 ; ...
               34.0 183.5   5.5  -3.9   -7.7   -5.8 ; ...
               30.8   5.5 179.5  -8.7    7.1   -9.8 ; ...
                5.1  -3.9  -8.7  24.9   -2.4   -7.2 ; ...
               -2.4  -7.7  7.1   -2.4   26.8    0.5 ; ...
               -0.9  -5.8  -9.8  -7.2    0.5   33.5 ];
         rh = 2623 ;                
      otherwise
         error('MS_elasticDB: Unknown identifier') ;
   end      
return






 
                             
                             
