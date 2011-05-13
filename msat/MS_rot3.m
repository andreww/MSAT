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

function [CR] = CIJ_rot3(C,alp,bet,gam)

[CC] = cij2cijkl(C) ;

%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;
 
 
% rotate the elastic contants
for M=1:3
 for N=1:3
  for R=1:3
   for S=1:3
    CSUM = 0.0 ;
    for I=1:3
     for J=1:3
      for K=1:3
       for L=1:3
        AA = RR(M,I)*RR(N,J)*RR(R,K)*RR(S,L) ;
        CSUM = CSUM + AA * CC(I,J,K,L) ;
       end
      end
     end
    end
    if (abs(CSUM) < 10.0) 
     CSUM=0.0 ;
    end
    CCR(M,N,R,S) = CSUM ;
   end
  end
 end
end

[CR] = cijkl2cij(CCR) ;

return

