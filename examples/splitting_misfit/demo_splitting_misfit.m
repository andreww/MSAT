function demo_splitting_misfit()
% Function to demonstrate the misfit functions implemented in 
% MSAT to compare shear-wave splitting operator pairs. 


fast_ref = 0 ;
tlag_ref = 3 ;
%fast = 60 ;
%tlag = 3 ;
%
%MISFIT = ...
%   MS_splitting_misfit(fast,tlag, ....
%      fast_ref,tlag_ref,45,0.1,'mode','lam2','debug') 
%
%MISFIT = ...
%   MS_splitting_misfit(fast_ref,tlag_ref, ...
%      fast,tlag,45,0.1,'mode','lam2','debug') 

%return


figure

calc_misfit_surface(fast_ref,tlag_ref,-45,'lam2',0,1) ;

calc_misfit_surface(fast_ref,tlag_ref,-45,'lam2',1,2) ;

%calc_misfit_surface(fast_ref,tlag_ref,45,'lam2S',0,3) ;

%calc_misfit_surface(fast_ref,tlag_ref,45,'simple',0,4) ;


end

function []=calc_misfit_surface(fast_ref,tlag_ref,spol,modeStr,reverse,ip) ;

   fast = [-90:5:90] ;
   tlag = [0:0.1:4] ;

   [TLAG,FAST] = meshgrid(tlag,fast) ;
   MISFIT = FAST.*0 ;

   % generate misfit 
   for itlag=1:length(tlag)
      for ifast=1:length(fast)
         if reverse
            MISFIT(ifast,itlag) = ...
               MS_splitting_misfit(fast(ifast),tlag(itlag), ....
                  fast_ref,tlag_ref, spol,0.1,'mode',modeStr,'max_tlag',10) ;
         else
            MISFIT(ifast,itlag) = ...
               MS_splitting_misfit(fast_ref,tlag_ref, ...
                  fast(ifast),tlag(itlag),spol,0.1,'mode',modeStr,'max_tlag',10) ;
         end
      end
   end

   subplot(1,2,ip)
   contourf(TLAG,FAST,MISFIT,20) ; colorbar ;
   


end