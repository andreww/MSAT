% MS_ROT3 - Elasticity matrix rotation.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Rotates an elasticity matrix around the three axies. 
%
%  % [CR] = CIJ_rot3(C,alp,bet,gam)
%
% Usage: 
%     [CR] = CIJ_rot3(C,alp,bet,gam)                    
%         C: input 6x6 elasticity matrix 
%         alp: clockwise rotation about 1-axis, looking at origin
%         bet: clockwise rotation about 2-axis, looking at origin
%         gam: clockwise rotation about 3-axis, looking at origin
%
% Notes:
%     Angles are given in degrees and correspond to yaw, -dip and aximuth,
%     respectvly. The rotations are applied in order, ie: alpha, then beta
%     then gamma

% (C) James Wookey and Andrew Walker, 2011.
% 
%

function [CR] = MS_rot3(C,alp,bet,gam)

%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;

%  Delegate the rotation
CR = MS_rotR(C, RR) ;


return

