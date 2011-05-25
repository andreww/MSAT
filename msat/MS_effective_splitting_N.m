% EFFECT_SPLITTING_N - N-layer effective splitting operator calculation from
% Silver and Savage (GJI, v119 pp 949-963, 1994)  
%
% [fast_eff,tlag_eff]=effective_splitting_N(f,spol,fast,tlag)
%  
%  Fast and tlag can be scalars or vectors but must all be the same length  
%  spol and f must be scalars. CAUTION! The method gives unstable results when
%  the source polarisation is near the effecive fast direction. This would,
%  however probably be seen as a null result anyway.

%===============================================================================
function [fast_eff,tlag_eff]=MS_effective_splitting(f,spol,fast,tlag)
%===============================================================================
%  check inputs are the same size   
   if ~isequal(length(fast),length(tlag))
      error('Input fast and tlag vectors must be of the same length')            
   end
   
%  check for just one layer
   if length(fast)==1
      fast_eff = fast ;
      tlag_eff = tlag ;
      return
   end                  

%**first unwind the fast directions
   fast = unwind_pm_90(fast) ;
   
%**deal with the situation where all fast directions are the same
   if length(find(fast==fast(1))) == length(fast)
      fast_eff = fast(1) ;
      tlag_eff = sum(tlag) ;
   else   

%**process   
      w = 2.0 .* pi .* f ;
      th = w .* tlag ./ 2. ;
      al = 2.*(fast-spol) .* pi/180.0 ;
            
      S = prod(cos(th)) ;
      Cc = S.*sum((tan(th).*cos(al))) ;
      Cs = S.*sum((tan(th).*sin(al))) ;
      
      n=length(fast) ;
      
      ap = 0;
      app = 0;
      
      for i=1:n-1
         for j=i+1:n
            ap = ap + (tan(th(i)).*tan(th(j)).*cos(al(i)-al(j))) ;
            app = app + (tan(th(i)).*tan(th(j)).*sin(al(i)-al(j))) ;
         end
      end
      
      ap = S.*(1-ap) ;
      app = S.*app ;
      ala = atan ( (app.^2.+Cs.^2.) ./ (app.*ap + Cs.*Cc) ) ;
      tha = atan ( (app) ./ (Cs.*cos(ala)-Cc.*sin(ala)) ) ;
   
      fast_eff = spol + (ala.*180./pi) ./ 2. ;
      tlag_eff = 2.*tha./w ;                                     
      
      fast_eff = unwind_pm_90(fast_eff) ;
   end
%**handle zero tlags      
   if (tlag_eff < 0)
      fast_eff = unwind_pm_90(fast_eff+90) ;
      tlag_eff = abs(tlag_eff) ;
   end   

return
%===============================================================================
