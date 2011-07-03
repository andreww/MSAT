% EFFECT_SPLITTING_N - N-layer effective splitting operator calculation.
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [fast_eff,tlag_eff]=MS_effective_splitting_N(f,spol,fast,tlag)
%  
%  Inputs:
%     f (scalar) : Dominant frequency of wave 
%     spol (scalar) : initial source polarisation
%     fast (scalar/vector) : fast direction(s) of layer(s) to be included
%     tlag (scalar/vector) : lag time(s) of layer(s) to be included
%
%  Fast and tlag can be scalars or vectors but must all be the same length.  
%
%  CAUTION! The method gives unstable results when
%  the source polarisation is near the effective fast direction. This would,
%  however probably be seen as a null result anyway.
%
%  Reference: Silver and Savage (GJI, v119 pp 949-963, 1994)  

% (C) James Wookey, 2005-2011
% (C) James Wookey and Andrew Walker, 2011

%===============================================================================
function [fast_eff,tlag_eff]=MS_effective_splitting(f,spol,fast,tlag)
%===============================================================================
%  check inputs are the same size   
   if ~isequal(length(fast),length(tlag))
      error('Input fast and tlag vectors must be of the same length')            
   end

%  remove any layers with 0 tlag time; these break the calculation
   ind = find(tlag~=0.0) ;
   fast2 = fast(ind) ;
   tlag2 = tlag(ind) ;
   fast = fast2 ;
   tlag = tlag2 ;
   
%  check for just one layer
   if length(fast)==1
      fast_eff = fast ;
      tlag_eff = tlag ;
      return
   end                  

%**first unwind the fast directions
   fast = MS_unwind_pm_90(fast) ;
   
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
