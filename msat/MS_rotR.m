% Script to rotate a set of elastic constants (Cij) by a rotation matrix
% 
% [CR] = CIJ_rot3(C,R)
%
%  Inputs: the Cij elastic stiffness tensor, 3x3 rotation matrix
%
%  Output is the rotated Cij elastic stiffness tensor                  
%

function [CR] = MS_rotR(C,R)

% form the K matrix (based on Bowers 'Applied Mechanics of Solids' chapter 3)
K = [R(1,1).^2      R(1,2).^2      R(1,3).^2             2.*R(1,2).*R(1,3)             2.*R(1,3).*R(1,1)              2.*R(1,1).*R(1,2) ; ...
     R(2,1).^2      R(2,2).^2      R(2,3).^2             2.*R(2,2).*R(2,3)             2.*R(2,3).*R(2,1)              2.*R(2,1).*R(2,2) ; ...
     R(3,1).^2      R(3,2).^2      R(3,3).^2             2.*R(3,2).*R(3,3)             2.*R(3,3).*R(3,1)              2.*R(3,1).*R(3,2) ; ...
R(2,1).*R(3,1) R(2,2).*R(3,2) R(2,3).*R(3,3) R(2,2).*R(3,3)+R(2,3).*R(3,2) R(2,3).*R(3,1)+R(2,1).*R(3,3)  R(2,1).*R(3,2)+R(2,2).*R(3,1) ; ...
R(3,1).*R(1,1) R(3,2).*R(1,2) R(3,3).*R(1,3) R(3,2).*R(1,3)+R(3,3).*R(1,2) R(3,3).*R(1,1)+R(3,1).*R(1,3)  R(3,1).*R(1,2)+R(3,2).*R(1,1) ; ...
R(1,1).*R(2,1) R(1,2).*R(2,2) R(1,3).*R(2,3) R(1,2).*R(2,3)+R(1,3).*R(2,2) R(1,3).*R(2,1)+R(1,1).*R(2,3)  R(1,1).*R(2,2)+R(1,2).*R(2,1)] ;
      
% apply the rotation
CR = K * C * transpose(K) ;

return

