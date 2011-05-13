% CIJ_list(C,rho)

function CIJ_list(C,rho);

for i=1:6
   for j=i:6
      
      fprintf('%1i %1i %e\n',i,j,C(i,j))
      
   end
end   
fprintf('%1i %1i %f\n',7,7,rho)
