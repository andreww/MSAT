% MS_UNWIND_PM_90 - Unwind an angle until it is between -90 and 90 degrees
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Utility function to find the equivelent angle between -90 and 90 degrees
% assuming 180 defgree periodicity. 
%
%  % [angle] = unwind_pm_90(angle_in)
%
% Usage: 
%     for a given value of angle_in return a value greater than 
%     -90 and less than or equal to +90 by repetedly adding or subtracting
%     180. Angle can be a scalar or a vector.

% Copyright (c) 2011, 2014 James Wookey and Andrew Walker
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

%===============================================================================
function [angle] = MS_unwind_pm_90(angle_in)
%===============================================================================
%  ** check input
      if (~isvector(angle_in) | ~isnumeric(angle_in))
         error('MS:unwind_pm_90:BadInput',...
            'Input is required to be a numeric scalar/vector') ;
      end
      
      angle = angle_in ;
   
%     shift to -180 -> 180   
      angle = angle - 180.*fix(angle./180) ;

%     refine: too small
      ind = find(angle<=-90);
      if ~isempty(ind), angle(ind) = angle(ind) + 180 ; , end    

%     or too large      
      ind = find(angle>90); 
      if ~isempty(ind), angle(ind) = angle(ind) - 180 ; , end 
   
return
%===============================================================================

% See also: MS_UNWIND_0_180, MS_UNWIND_360, MS_UNWIND_PM_180

