% MS_VTI2 - generate elastic constants for a vertically transverse isotropic
%             medium from xi, phi, eta parameters. Symmetry is in the 3-axis. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
%   [C]=MS_VTI2(vp,vs,rh,xi,phi,eta)
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
%   Calculates the elastic tensor for a hexagonal (aka VTI) medium from
%   average Vp and Vs, and anisotropic parameters xi, phi and eta. 
%
%   See, for example, Panning and Romanowicz (2006) for a definition of the 
%   parameters. 
%
%   ** Note: this is wrapper for MS_TI, provided for back-compatibility. **
%
%   References:
%     Mark Panning and Barbara Romanowicz (2006) A three-dimensional radially 
%        anisotropic model of shear velocity in the whole mantle. Geophysical  
%        Journal International v167, 361–379. 
%        doi: 10.1111/j.1365-246X.2006.03100.x
%
% See also: MS_TI, MS_VTI, MS_iso, MS_elasticDB

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

function [C]=MS_VTI2(vp,vs,rh,xi,phi,eta)

   C=MS_TI(vp,vs,rh,xi,phi,eta,'panning') ;

return

%function [C]=MS_VTI2(vp,vs,rh,xi,phi,eta)
%
%%  convert to m/s
%   vp=vp*1e3;
%   vs=vs*1e3;
%
%%  note: xi and phi are defined oppositely; i.e.: 
%%  xi = vsv^2/vsh^2 and phi = vph^2/vpv^2
%   vsv = sqrt((3.*vs.^2)./(2+xi)) ;
%   vsh = sqrt(xi.*vsv.^2) ;
%   
%   vph = sqrt((5.*vp.^2)./(4+phi)) ;
%   vpv = sqrt(phi.*vph.^2) ;
%   
%   C11 = vph.^2.*rh ; % A
%   C33 = vpv.^2.*rh ; % C
%   C44 = vsv.^2.*rh ; % L
%   C66 = vsh.^2.*rh ; % N
%   
%   C12 = C11-2.*C66 ;
%   C13 = eta.*(C11-2.*C44) ; % F
%   
%   C22 = C11 ;
%   C23 = C13 ;
%   C55 = C44 ;
%   
%   C = [C11 C12 C13  0   0   0  ; ...
%        C12 C22 C13  0   0   0  ; ...
%        C13 C13 C33  0   0   0  ; ...
%         0   0   0  C44  0   0  ; ...
%         0   0   0   0  C55  0  ; ...
%         0   0   0   0   0  C66 ] ;
%   
%   C = C./1e9 ; % convert to GPa
%
%return
