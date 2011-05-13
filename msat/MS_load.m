% [C,rho] = CIJ_load(str)
function [C,rho] = CIJ_load(str) ;

C = zeros(6,6) ;

a=load(str) ;

[nec ndum] = size(a) ;

for i=1:nec-1
   C(a(i,1),a(i,2)) = a(i,3) ;
end

% if NEC = 2 (i.e. isotropic, need to fill out matrix)
if (nec-1 == 2)
      C(1,1) = C(3,3) ; C(2,2) = C(3,3) ;
      C(5,5) = C(4,4) ; C(6,6) = C(4,4) ;
      C(1,2) = (C(3,3)-2.*C(4,4)) ;
      C(1,3) = C(1,2) ; C(2,3) = C(1,2) ;   
end



% make symmetrical
for i=1:6
   for j=i:6
      C(j,i) = C(i,j) ;
   end
end
   
rho = a(nec,3) ;
