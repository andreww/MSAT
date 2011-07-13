% MS_UNWIND_PM_90 - Unwind an angle until it is between 0 and 360 degrees
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Utility function to find the equivelent angle between 0 and 360 degrees. 
%
%  % [angle] = unwind_pm_90(angle_in)
%
% Usage: 
%     for a given value of angle_in return a value greater than 
%     -90 and less than or equal to 900 by repetedly adding or subtracting
%     180. Angle can be a scalar or a vector.

% (C) James Wookey and Andrew Walker, 2011

%===============================================================================
function [angle] = MS_unwind_pm_90(angle_in)
%===============================================================================
%  ** check input
      if ~isvector(angle_in)
         error('MS:unwind_pm_90:BadInput',...
            'Input is required to be a scalar/vector') ;
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

