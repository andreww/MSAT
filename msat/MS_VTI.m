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
%   ** Note: this is wrapper for MS_TI, provided for back-compatibility. **
%
% References: 
%     Thomsen, L. (1986) "Weak elastic anisotropy" Geophysics 
%         vol.51 pp.1954-1966
%
% See also: MS_TI, MS_iso, MS_elasticDB

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

function [C]=MS_VTI(vp,vs,rh,eps,gam,del)
%  use MS_TI routine. 
   [C]=MS_TI(vp,vs,rh,eps,gam,del,'thomsen') ;
return
