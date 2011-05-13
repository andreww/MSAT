% [C] = cijkl2cij(CC)

function [C] = cijkl2cij(CC)
 for im=1:3
   for jm=1:3
      for km=1:3
         for lm=1:3
            if ( CC(im,jm,km,lm) ~= 0.0) 
               [iv,jv]=ijkl2ij_local(im,jm,km,lm) ;
               C(iv,jv) = CC(im,jm,km,lm);
            end
         end
      end
   end
 end
return

function [iv,jv] = ijkl2ij_local(ii,jj,kk,ll)   
	if (ii==1 & jj==1) 
      iv=1;
   end
   if (ii==1 & jj==2) 
      iv=6;
   end
   if (ii==1 & jj==3) 
      iv=5;
   end
   if (ii==2 & jj==1) 
      iv=6;
   end

   if (ii==2 & jj==2) 
      iv=2;
   end
   if (ii==2 & jj==3) 
      iv=4;
   end
   
   if (ii==3 & jj==1) 
      iv=5;
   end
   if (ii==3 & jj==2) 
      iv=4;
   end
   if (ii==3 & jj==3) 
      iv=3;
   end
   
   if (kk==1 & ll==1) 
      jv=1;
   end
   if (kk==1 & ll==2) 
      jv=6;
   end
   if (kk==1 & ll==3) 
      jv=5;
   end

   if (kk==2 & ll==1) 
      jv=6;
   end
   if (kk==2 & ll==2) 
      jv=2;
   end
   if (kk==2 & ll==3) 
      jv=4;
   end

   if (kk==3 & ll==1) 
      jv=5;
   end
   if (kk==3 & ll==2) 
      jv=4;
   end
   if (kk==3 & ll==3) 
      jv=3;
   end

return
