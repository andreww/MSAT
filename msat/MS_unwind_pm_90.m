%  MS_UNWIND_PM_90 - Unwind an angle until it is between -90 and 90 degrees
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [angle] = MS_unwind_pm_90(angle_in)
%
% This assumes a 180 degree periodicity. Angle can be a scalar or a vector. 
%
% See also: 

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
