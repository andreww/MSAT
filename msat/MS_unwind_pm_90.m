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
%
% See also: MS_UNWIND_0_180, MS_UNWIND_360, MS_UNWIND_PM_180

%
%  Created by James Wookey on 2007-05-08.
%
function [angle] = MS_unwind_pm_90(angle_in)
%===============================================================================
   angle = angle_in ;
   for i=1:1000 % should never reach this
      
%  ** check for completion
      ind = find(angle <= -90.0 | angle > 90.0) ;
      if isempty(ind), return, end % done

%     else refine: too small
      ind = find(angle<=-90);
      if ~isempty(ind), angle(ind) = angle(ind) + 180 ; , end    

%     and too large      
      ind = find(angle>90); 
      if ~isempty(ind), angle(ind) = angle(ind) - 180 ; , end 
   end

%  if we got to here something is probably wrong:
error('MS:UNWIND:toomanyits', ...
    'UNWIND_PM_90 completed 1000 iterations without angle reaching +/-90 deg')   
   
return
%===============================================================================
