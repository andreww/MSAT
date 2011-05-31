% Script to rotate a set of elastic constants (Cij) by a rotation matrix
% 
% [CR] = CIJ_rot3(C,R)
%
%  Inputs: the Cij elastic stiffness tensor, 3x3 rotation matrix
%
%  Output is the rotated Cij elastic stiffness tensor                  
%

function [CR] = MS_rotR(C,R)


RR =  R;
 
[CC] = cij2cijkl(C) ;
 
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
    CCR(M,N,R,S) = CSUM ;
   end
  end
 end
end

[CR] = cijkl2cij(CCR) ;



return

