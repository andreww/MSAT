%   MS_VTI - generate elastic constants for a vertically transverse isotropic
%             medium from Thomsen parameters. Symmetry is in the 3-axis. 
%
%   [C]=MS_VTI(vp,vs,rh,eps,gam,del)
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s 
%       vs  : km/s 
%       eps, gam, del : Dimensionless
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
% See also: MS_iso, MS_elasticDB

function [C]=MS_VTI(vp,vs,rh,eps,gam,del)

%  convert to m/s
   vp=vp*1e3;
   vs=vs*1e3;

   C=zeros(6,6) ;
   C(3,3) = vp*vp ;
   C(4,4) = vs*vs ;
   C(6,6) = C(4,4)*(2.0*gam +1.0) ;
   C(1,1) = C(3,3)*(2.0*eps +1.0) ;
   
   btm = 2.0*C(4,4) ;
   term = C(3,3) - C(4,4) ;
   ctm = C(4,4)*C(4,4) - (2.0*del*C(3,3)*term + term*term) ;
   dsrmt = (btm*btm - 4.0*ctm) ;
   
	if dsrmt < 0
		error('MS:VTI:HiVS',...
		   'S-velocity too high or delta too negative for Thomsen routine.') ;
	end
   
   C(1,3) = -btm/2.0 + sqrt(dsrmt)/2.0 ;
   
   C(1,2) = C(1,1) - 2.0*C(6,6) ; 
   C(2,3) = C(1,3) ;
   C(5,5) = C(4,4) ;
   C(2,2) = C(1,1) ;

   % make symmetrical
   for i=1:6
      for j=i:6
         C(j,i) = C(i,j) ;
      end
   end

%  convert to GPa
   C = C.*rh./1e9 ;

return
