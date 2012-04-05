% MS_TI - generate elastic constants for a vertically transverse isotropic
%             medium from specified parameter sets. Symmetry is in the 3-axis. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%  [C]=MS_TI( list_of_parameters , parameter_set_string )
%
%  where parameter_set_string defines the set which precede it: 
%
%-------------------------------------------------------------------------------
%   'thomsen' (default)
%-------------------------------------------------------------------------------
%
%   [C]=MS_TI(vp,vs,rh,eps,gam,del) -or-
%   [C]=MS_TI(vp,vs,rh,eps,gam,del,'thomsen')
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s (vertical) 
%       vs  : km/s (vertical)
%       eps, gam, del : Dimensionless
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Given Thomsen (1986) parameters for a weakly anisotropic VTI medium, 
%   return the elasticity matrix
%
%-------------------------------------------------------------------------------
%   'panning'
%-------------------------------------------------------------------------------
%
%   [C]=MS_TI(vp,vs,rh,xi,phi,eta,'panning')
%
%   Inputs: 
%       rh  : Density (kg/m2) 
%       vp  : km/s (isotropic average)
%       vs  : km/s (isotropic average)
%       xi, phi, eta : Dimensionless anisotropy parameters
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from average Vp and Vs,
%   and anisotropic parameters xi, phi and eta (see, e.g., Panning and 
%   Romanowicz, 2006)
%
%-------------------------------------------------------------------------------
%   'love'
%-------------------------------------------------------------------------------
%
%   [C]=MS_TI(A,C,L,N,F,'love')
%
%   Inputs: 
%        A,C,L,N,F : Love parameters (GPa)
%
%   Output:
%        C : Stiffness tensor (6x6 notation, GPa)
%
%   Calculates the elastic tensor for a VTI medium from Love (1927) parameters.
%
%-------------------------------------------------------------------------------
% References: 
%     Thomsen, L. (1986) "Weak elastic anisotropy" Geophysics 
%         vol.51 pp.1954-1966
%
%     Mark Panning and Barbara Romanowicz (2006) A three-dimensional radially 
%        anisotropic model of shear velocity in the whole mantle. Geophysical  
%        Journal International v167, 361–379. 
%        doi: 10.1111/j.1365-246X.2006.03100.x
%
%     Love, A.E.H., (1927). A Treatise on the Theory of Elasticity, 
%        Cambridge Univ. Press, Cambridge.
%
% See also: MS_iso, MS_elasticDB

% Copyright (c) 2011-2012, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [C]=MS_TI(varargin)
   
% see if a parameter set type was selected. 
if ischar(varargin{nargin})
   pset = varargin{nargin} ;
   ncheck = nargin-1;
else
   pset = 'thomsen' ;
   ncheck = nargin ;
end   

% check that all other inputs are scalars
for icheck=1:ncheck
   if ~isscalar(varargin{icheck})
      error('MS:TI:BadScalarInput',...
      'MS_TI requires scalar inputs') ;
   end   
end

%  call the appropriate routine.

switch lower(pset)
%-------------------------------------------------------------------------------
   case {'thomsen', 'thom'}
      if length(varargin)~=7 & length(varargin)~=6 % need to check 
         error('MS:TI:ThomsenWrongArgs', ...
         'Thomsen (1986) requires 6 input parameters.') ;
      end
      vp = varargin{1} ; vs = varargin{2} ; rh = varargin{3} ; 
      eps = varargin{4} ; gam = varargin{5} ; del = varargin{6} ;
      [C]=MS_thomsen(vp,vs,rh,eps,gam,del) ;
%-------------------------------------------------------------------------------
   case {'panning'}
      if length(varargin)~=7 % need to check 
         error('MS:TI:PanningWrongArgs', ...
         'Panning (2006) requires 6 input parameters.') ;
      end
      vp = varargin{1} ; vs = varargin{2} ; rh = varargin{3} ; 
      xi = varargin{4} ; phi = varargin{5} ; eta = varargin{6} ;
      [C]=MS_panning(vp,vs,rh,xi,phi,eta) ;
%-------------------------------------------------------------------------------
   case {'love'}
      if length(varargin)~=6 % need to check 
         error('MS:TI:LoveWrongArgs', ...
         'Love (1926) requires 5 input parameters.') ;
      end
      A = varargin{1} ; C_love = varargin{2} ; L = varargin{3} ; 
      N = varargin{4} ; F = varargin{5} ; 
      [C]=MS_love(A,C_love,L,N,F) ;
%-------------------------------------------------------------------------------
   otherwise
      error('MS:TI:UnknownParSet', ...
         'Specified parameter set is not supported.') ;
%-------------------------------------------------------------------------------
end % of switch
   
% check resulting matrix.
MS_checkC(C) ;
   
end

function [CC]=MS_love(A,C,L,N,F)
   CC= [ A      A-2.*N F  0  0  0  ; ...
         A-2.*N A      F  0  0  0  ; ...
         F      F      C  0  0  0  ; ...
         0      0      0  L  0  0  ; ...
         0      0      0  0  L  0  ; ...
         0      0      0  0  0  N ] ;
end

function [C]=MS_panning(vp,vs,rh,xi,phi,eta)

%  convert to m/s
   vp=vp*1e3;
   vs=vs*1e3;

%  note: xi and phi are defined oppositely; i.e.: 
%  xi = vsv^2/vsh^2 and phi = vph^2/vpv^2
   vsv = sqrt((3.*vs.^2)./(2+xi)) ;
   vsh = sqrt(xi.*vsv.^2) ;
   
   vph = sqrt((5.*vp.^2)./(4+phi)) ;
   vpv = sqrt(phi.*vph.^2) ;
   
   C11 = vph.^2.*rh ; % A
   C33 = vpv.^2.*rh ; % C
   C44 = vsv.^2.*rh ; % L
   C66 = vsh.^2.*rh ; % N
   
   C12 = C11-2.*C66 ;
   C13 = eta.*(C11-2.*C44) ; % F
   
   C22 = C11 ;
   C23 = C13 ;
   C55 = C44 ;
   
   C = [C11 C12 C13  0   0   0  ; ...
        C12 C22 C13  0   0   0  ; ...
        C13 C13 C33  0   0   0  ; ...
         0   0   0  C44  0   0  ; ...
         0   0   0   0  C55  0  ; ...
         0   0   0   0   0  C66 ] ;
   
   C = C./1e9 ; % convert to GPa

end

function [C]=MS_thomsen(vp,vs,rh,eps,gam,del)

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
		error('MS:TI:ThomsenHiVS',...
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

end
