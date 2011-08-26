% MS_ROTM - CREATE A CARTESIAN ROTATION MATRIX
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [ M ] = MS_rotM( alp, bet, gam )
%
%  Inputs:
%     alp - clockwise about 1-axis (looking at origin, == yaw)
%     bet - clockwise about 2-axis (looking at origin, == -dip)
%     gam - clockwise about 3-axis (looking at origin, == azimuth)
%           angles are in ** degrees **, they should be scalars. 
%
%  Output: is M, a 3x3 rotation matrix
%
%     [M] = MS_rotM(...,'order',V)                    
%         order: V is a 3 element vector containing the order in which to apply the 3 rotations.
%                Default is [1 2 3].
%
%  Notes: 
%    The rotations are applied in order, ie: alpha, then beta then gamma (by default).
%

% (C) James Wookey and Andrew Walker, 2011

function [ M ] = MS_rotM( alp, bet, gam, varargin )


orderV = [1 2 3] ;

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'order' 
         orderV = varargin{iarg+1} ;
         if (length(find(orderV==1))+length(find(orderV==2))+length(find(orderV==3)))~=3
            error('MS:ROTM:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;   
         end
         if (length(orderV)~=3)
            error('MS:ROTM:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;  
         end                  
         iarg = iarg + 2 ;
      otherwise 
         error('MS:ROTM:UnknownOption',...
            ['Unknown option: ' varargin{iarg}]) ;   
   end   
end

a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

%  Make rotation matrices

R(1,:,:) = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R(2,:,:) = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R(3,:,:) = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

M =  squeeze(R(orderV(3),:,:)) * ...
         squeeze(R(orderV(2),:,:)) * ...
            squeeze(R(orderV(1),:,:));

return

