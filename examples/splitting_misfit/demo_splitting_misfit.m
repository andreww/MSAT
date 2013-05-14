function demo_splitting_misfit()
% Function to demonstrate the misfit functions implemented in 
% MSAT to compare shear-wave splitting operator pairs. 


fast_ref = 0 ;
tlag_ref = 2 ;


figure('Position',[1 1 800 800]) ;

fprintf('Yellow diamond is the reference value.\n')


calc_misfit_surface(fast_ref,tlag_ref,-45,'lam2',0,1) ;
hold on
plot(tlag_ref,fast_ref,'kd','MarkerFaceColor','y','MarkerSize',12)
title('Mode = lam2')

calc_misfit_surface(fast_ref,tlag_ref,45,'lam2S',0,2) ;
hold on
plot(tlag_ref,fast_ref,'kd','MarkerFaceColor','y','MarkerSize',12)
title('Mode = lam2S')

calc_misfit_surface(fast_ref,tlag_ref,-45,'simple',0,3) ;
hold on
plot(tlag_ref,fast_ref,'kd','MarkerFaceColor','y','MarkerSize',12)
title('Mode = simple')

calc_misfit_surface(fast_ref,tlag_ref,45,'intensity',0,4) ;
hold on
plot(tlag_ref,fast_ref,'kd','MarkerFaceColor','y','MarkerSize',12)
title('Mode = intensity')


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

   subplot(2,2,ip)
   [C, H] = contourf(TLAG,FAST,MISFIT,20) ; colorbar ;
   set(H,'LineStyle','none') ;
   xlabel('TLAG (sec)')
   ylabel('FAST (deg)')


end