% MS_ROT3 - Elasticity matrix rotation.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotates an elasticity matrix around the three axes. 
%
% Usage: 
%     [CR] = MS_rot3(C,alp,bet,gam)                    
%         C: input 6x6 elasticity matrix 
%         alp: clockwise rotation about 1-axis, looking at origin
%         bet: clockwise rotation about 2-axis, looking at origin
%         gam: clockwise rotation about 3-axis, looking at origin
%
%     [CR] = MS_rot3(...,'order',V)                    
%         order: V is a 3 element vector containing the order in which to apply the 3 rotations.
%                Default is [1 2 3].
%
% Notes:
%     Angles are given in degrees and correspond to yaw, -dip and aximuth,
%     respectvly. The rotations are applied in order, ie: alpha, then beta
%     then gamma (by default).

% (C) James Wookey and Andrew Walker, 2011.
% 
%

function [CR] = MS_rot3(C,alp,bet,gam,varargin)

try 
   MS_checkC(C) ;
catch ME
   error('MS:ROT3BadInputMatrix', ...
      ['Invalid input elasticity matrix: ' ME.message]) ;
end

orderV = [1 2 3] ;

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'order' 
         orderV = varargin{iarg+1} ;
         if (length(find(orderV==1))+length(find(orderV==2))+length(find(orderV==3)))~=3
            error('MS:ROT3:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;   
         end
         if (length(orderV)~=3)
            error('MS:ROT3:BadOrder',...
            'Order vector must be 3-element, containing 1,2 and 3') ;  
         end                  
         iarg = iarg + 2 ;
      otherwise 
         error('MS:ROT3:UnknownOption',...
            ['Unknown option: ' varargin{iarg}]) ;   
   end   
end

a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

RM = zeros(3,3,3) ;

R(1,:,:) = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R(2,:,:) = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R(3,:,:) = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  squeeze(R(orderV(3),:,:)) * ...
         squeeze(R(orderV(2),:,:)) * ...
            squeeze(R(orderV(1),:,:));

%  Delegate the rotation
CR = MS_rotR(C, RR) ;


return

