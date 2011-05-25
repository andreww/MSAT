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
      otherwise
         error('MS_elasticDB: Unknown identifier') ;
   end      
return