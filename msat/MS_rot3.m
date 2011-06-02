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

function [CR] = MS_rot3(C,alp,bet,gam)

try 
   MS_checkC(C) ;
catch ME
   error('MS:ROT3BadInputMatrix', ...
      ['Invalid input elasticity matrix: ' ME.message]) ;
end

%  Make 6x6 rotation matrices
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

% These are from Bower Chapter 3, section 3.2. Note that that e2 and e3 matrices
% in that section contain a number of errors; these have been corrected. 

c = cos(a) ; s = sin(a) ;
R1 = [ 1      0      0      0      0      0 ; ...
       0     c^2    s^2   2*c*s    0      0 ; ...
       0     s^2    c^2  -2*c*s    0      0 ; ...
       0    -c*s    c*s  c^2-s^2   0      0 ; ...
       0      0      0      0      c     -s ; ...
       0      0      0      0      s      c ] 

c = cos(b) ; s = sin(b) ;
R2 = [c^2     0     s^2     0   -2*c*s    0 ; ...
       0      1      0      0      0      0 ; ...
      s^2     0     c^2     0    2*c*s    0 ; ...
       0      0      0      c      0      s ; ...
      c*s     0    -c*s     0   c^2-s^2   0 ; ...
       0      0      0     -s      0      c ] ;

c = cos(g) ; s = sin(g) ;
R3 = [c^2    s^2     0      0      0    2*c*s ; ...
      s^2    c^2     0      0      0   -2*c*s ; ...
       0      0      1      0      0      0 ; ...
       0      0      0      c      s      0 ; ...
       0      0      0     -s      c      0 ; ...
     -c*s    c*s     0      0      0   c^2-s^2 ] ;
              
RR =  R3 * R2 * R1 ;


CR = RR * C * transpose(RR) ;


return

