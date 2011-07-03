%  UNWIND_PM_90 - Unwind an angle until it is between -90 and 90 degrees
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [angle] = MS_unwind_pm_90(angle_in)
%
%  This (obviously) assumes a 180 degree periodicity. Angle can be a scalar
%  or a vector. 
%
% See also: 

% (C) James Wookey and Andrew Walker, 2011

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
error('UNWIND_PM_90 completed 1000 iterations without angle reaching +/-90 deg')   
   
return
%===============================================================================
