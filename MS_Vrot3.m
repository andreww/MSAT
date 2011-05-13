% Script to rotate a ROW 3-vector by 3 angles:
%
% alp = clockwise about 1-axis (looking at origin, == yaw)
% bet = clockwise about 2-axis (looking at origin, == -dip)
% gam = clockwise about 3-axis (looking at origin, == azimuth)
%
% [VR] = CIJ_rot3(V,alp,bet,gam)
%
%  Inputs: the vector, 3 angles (degrees)
%
%  Output is the rotated vector                 
%
%  NOTE: The rotations are applied in order, ie: alpha, then beta then gamma
%
%

function [VR] = V_rot3(V,alp,bet,gam)

%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;

VR = RR * V';
VR = VR';
 
return

