% MS_iso - generate elastic matrix from isotropic velocities.
%
%   [C]=MS_iso(vp,vs,rh)
%
%   Inputs: 
%       rh  : Density (kg/m3) 
%       vp  : p-wave velocity (km/s) 
%       vs  : s-wave velocity (km/s)
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)

function [C]=MS_iso(vp,vs,rh)

   vp=vp*1e3;
   vs=vs*1e3;
         
   C=zeros(6,6) ;
   C(3,3) = vp*vp ;
   C(6,6) = vs*vs ;
   
   C=MS_expand(C,'iso') ;

%  convert to GPa
   C = C.*rh./1e9 ;

return
