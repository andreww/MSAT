% MS_VTI - generate elastic constants for a vertically transverse isotropic
%             medium from Thomsen parameters. Symmetry is in the 3-axis. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given Thomsen (1986) parameters for a weakly anisotropic 
% VTI medium, return the elasticity matrx
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
% References: 
%     Thomsen, L. (1986) "Weak elastic anisotropy" Geophysics 
%         vol.51 pp.1954-1966
%
% See also: MS_iso, MS_elasticDB

% Copyright (c) 2011, James Wookey and Andrew Walker
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
