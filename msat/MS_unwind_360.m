%  UNWIND_360 - Unwind an angle until it is between 0 and 360 degrees
%
%  [angle] = unwind_pm_360(angle_in)
%
%  Angle can be a scalar or a vector. 
%
%  Created by James Wookey on 2007-05-08.
%
function [angle] = unwind_360(angle_in)
%===============================================================================
   angle = angle_in ;
   for i=1:1000 % should never reach this
      
%  ** check for completion
      ind = find(angle < 0.0 | angle >= 360.0) ;
      if isempty(ind), return, end % done

%     else refine: too small
      ind = find(angle < 0);
      if ~isempty(ind), angle(ind) = angle(ind) + 360 ; , end    

%     and too large      
      ind = find(angle>=360); 
      if ~isempty(ind), angle(ind) = angle(ind) - 360 ; , end 
   end

%  if we got to here something is probably wrong:
error('UNWIND_360 completed 1000 iterations without angle reaching 0-360 deg')   
   
return
%===============================================================================
