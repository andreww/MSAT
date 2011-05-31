% Script to rotate a set of elastic constants (Cij) by 3 angles:
%
% alp = clockwise about 1-axis (looking at origin, == yaw)
% bet = clockwise about 2-axis (looking at origin, == -dip)
% gam = clockwise about 3-axis (looking at origin, == azimuth)
%
% [CR] = CIJ_rot3(C,alp,bet,gam)
%
%  Inputs: the Cij elastic stiffness tensor, 3 angles (degrees)
%
%  Output is the rotated Cij elastic stiffness tensor                  
%
%  NOTE: The rotations are applied in order, ie: alpha, then beta then gamma
%
%

function [CR] = MS_rot3C4(C,alp,bet,gam)


%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;

CR = MS_rotRC4(C,RR) ;

return

