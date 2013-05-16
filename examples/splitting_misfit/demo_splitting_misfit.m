% DEMO_SPLITTING_MISFIT
% Function to demonstrate the misfit functions implemented in 
% MSAT to compare shear-wave splitting operator pairs. 
%
% This code plots a range of contour plots showing the
% misfit between a reference splitting operator and a range
% of other splitting operators, for different source
% polarisations and misfit modes in MS_SPLITTING_MISFIT.
%
% See also: MS_SPLITTING_MISFIT

% Copyright (c) 2013, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%==============================================================================
function demo_splitting_misfit()
%==============================================================================


fast_ref = 0 ;
tlag_ref = 2 ;

fprintf('Yellow diamond is the reference value.\n')


figure('Position',[1 1 800 800],'Name','Polarisation = 45') ;

   calc_misfit_surface(fast_ref,tlag_ref,45,'lam2',0,1) ;
   calc_misfit_surface(fast_ref,tlag_ref,45,'lam2S',0,2) ;
   calc_misfit_surface(fast_ref,tlag_ref,45,'simple',0,3) ;
   calc_misfit_surface(fast_ref,tlag_ref,45,'intensity',0,4) ;


figure('Position',[1 1 800 800],'Name','Polarisation = 30') ;

   calc_misfit_surface(fast_ref,tlag_ref,30,'lam2',0,1) ;
   calc_misfit_surface(fast_ref,tlag_ref,30,'lam2S',0,2) ;
   calc_misfit_surface(fast_ref,tlag_ref,30,'simple',0,3) ;
   calc_misfit_surface(fast_ref,tlag_ref,30,'intensity',0,4) ;

figure('Position',[1 1 800 800],'Name','Polarisation = 0') ;

   calc_misfit_surface(fast_ref,tlag_ref,0,'lam2',0,1) ;
   calc_misfit_surface(fast_ref,tlag_ref,0,'lam2S',0,2) ;
   calc_misfit_surface(fast_ref,tlag_ref,0,'simple',0,3) ;
   calc_misfit_surface(fast_ref,tlag_ref,0,'intensity',0,4) ;

end
%==============================================================================


%==============================================================================
function []=calc_misfit_surface(fast_ref,tlag_ref,spol,modeStr,reverse,ip) ;
%==============================================================================

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
   hold on
   plot(tlag_ref,fast_ref,'kd','MarkerFaceColor','y','MarkerSize',12)
   
   xlabel('TLAG (sec)')
   ylabel('FAST (deg)')
   title(['mode = ' modeStr]) ;

end
%==============================================================================
