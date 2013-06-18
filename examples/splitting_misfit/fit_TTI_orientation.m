% FIT_TTI_ORIENTATION - Example of inverting a (synthetic) set of shear 
%    wave splitting data for the best fitting tilted transversely-isotropic
%    (in this case, elliptical) medium. This 
%
%  In this example:
%
%  1. A synthetic TTI model is created, and shear wave splitting is predicted
%     for two raypaths across it.
%
%  2. The code then performs a grid search of orientations of a base tensor,
%     measuring the misfit between splitting associated with the same 
%     same raypaths and the original synthetics.  
%
%  3. The misfit surface is plotting, demonstrating the minimum misfit 
%     is at the correct orientation. 
%  
% See also: MS_SPLITTING_MISFIT, DEMO_SPLITTING_MISFIT

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

%==============================================================================
function fit_TTI_orientation()
%==============================================================================
   
   %% set the raypath parameters. ---------------------------------------------
   razi = [045 135 048] ;
   rinc = [000 000 032] ;
   dist = 50 ; % km
   
   %% set up the medium to be modelled. ---------------------------------------
   str = 60 ;
   dip = 25 ;

   vpi=2;
   vsi=1;
   rh=2000;
   
   [Craw]=MS_VTI(vpi,vsi,rh,0.05,0.05,0.05) ;
   
   %% generate the synthetic splitting data to be inverted. -------------------
   fprintf( ...
      'Creating dataset with orientation at strike = %2.0f, dip = %2.0f\n',...
         str,dip)
   [fast,tlag,ravs]=create_dataset(Craw,rh,str,dip,razi,rinc,dist) ;

   %% calculate a grid of misfits over a range of strikes and dips. -----------
   fprintf('Grid searching over strike and dip ... \n')
   [GSTR GDIP] = meshgrid([0:10:360],[0:5:90]) ; 
   GMIS=GSTR.*0 ;
   [n,m]=size(GSTR) ;
   for i=1:n
      for j=1:m
         str_dip = [GSTR(i,j) GDIP(i,j)] ;
         GMIS(i,j) = orientation_misfit(str_dip,Craw,rh,...
            0,0.1,fast,tlag,razi,rinc,dist)  ;
      end
   end

   %% find the minimum misfit in the grid. ------------------------------------
   [GMIS_minR,imin] = min(GMIS) ;
   [min_misfit,jmin] = min(GMIS_minR) ;
   imin = imin(jmin) ;
   
   %% plot the resulting misfit grid. -----------------------------------------
   figure
   contourf(GSTR,GDIP,GMIS,[0.00:0.01:0.1])
   colorbar        

   fprintf('Minimum (%10.3e) misfit occurs at: strike = %2.0f, dip = %2.0f\n',...
      min_misfit,GSTR(imin,jmin),GDIP(imin,jmin))
   
   hold on
   
   % plot the minimum
   plot(GSTR(imin,jmin),GDIP(imin,jmin),...
      'ks','MarkerFaceColor','w','MarkerSize',10)
      
   xlabel('Strike (deg)')
   ylabel('Dip (deg)')   
   title('Misfit of TTI orientation to splitting data')
   
   %% plot comparing obs with best fitting result
   % First get a best model
   fprintf('Plotting model and observations\n')
   Cbest = MS_rot3(Craw, 0, -GDIP(imin,jmin), GSTR(imin,jmin));
   MS_plot(Cbest,rh,'sdata',razi,rinc,fast',ravs','plotmap', {'avspol'}) ;

end   
%==============================================================================

%==============================================================================
function [misfit]=orientation_misfit(str_dip,Cbase,rh,...
   data_spol,data_dfreq,data_fast,data_tlag,razi,rinc,dist)
%==============================================================================
%
%  Calculate the misfit between the splitting predicted (for a set of raypaths
%  razi, rinc, dist) for a tensor (Cbase,rh) at an orientation (str_dip, 
%  a two-element vector) and a set of shear-wave splitting data (data_spol,
%  data_dfreq,data_fast,data_tlag). 
%

   % rotate the test matrix
   [Ctest]=MS_rot3(Cbase,0,-str_dip(2),str_dip(1)) ;
   
   % calculate associated splitting
   [pol,avs,vs1,vs2,~,~,~] = MS_phasevels(Ctest,rh,rinc,razi) ;
   fast = pol ;
   tlag = dist./vs2 - dist./vs1 ;
   
   
   % calculate the misfit (the average over all raypaths)
   misfit=0 ;
   for i=1:length(data_fast)
      misfit = misfit + MS_splitting_misfit(data_fast(i),data_tlag(i),...
                           fast(i),tlag(i),data_spol,data_dfreq,'mode','lam2') ;
   end
   misfit = misfit ./ length(data_fast) ;
   
end
%==============================================================================

%==============================================================================
function [fast,tlag,avs]=create_dataset(Craw,rh,str,dip,razi,rinc,dist)
%==============================================================================
%
%  Calculate a set of shear-wave splitting results (at raypaths razi,rinc,dist)
%  associated with an elastic tensor (Craw,rh) oriented at a specified
%  strike and dip (str,dip).
%

   [Cref]=MS_rot3(Craw,0,-dip,str) ;
   [pol,avs,vs1,vs2,~,~,~] = MS_phasevels(Cref,rh,rinc,razi) ;
   
   fast = pol ;
   tlag = dist./vs2 - dist./vs1 ;
   
   % display on the sphere.
   % MS_sphere(Cref,rh,'s','sdata',razi,rinc,pol,avs) ;
   
end
%==============================================================================