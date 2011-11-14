% SPLIT_MODEL.M - Example script modelling shear-wave splitting variation with
%                 backazimuth. 
%
% This script demonstrates simple shear-wave splitting modelling using the MSAT
% toolset. It predicts the SKS splitting variation with backazimuth associated 
% with two dipping, partially-aligned olivine layers. The splitting in each layer
% is calculated using the Christoffel equation (MS_phasvels), and combined using
% N-layer effective splitting equations (MS_effective_splitting_N; Silver and 
% Savage, GJI, 1994). Straight raypaths through the upper mantle are assumed. 
%
% See also: MS_EFFECTIVE_SPLITTING_N, MS_PHASEVELS

% (C) James Wookey and Andrew Walker, 2011
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

function split_model()

%  ** Setup model parameters
      S_slow = 4.814 ; % (SKS at 100 degrees from iasp91) ->
      aoi = 11.2262 ; % at 100 km depth, angle of incidence

%  ** Layer 1 (upper) parameters
      L1_depth = 0. ; L1_thick = 100. ; L1_dip = 0.0 ; % layer geometry
      L1_aaz = 30.0 ; % a-axis azimuth (rel. to down dip direction)
      L1_faln = 0.3 ; % fraction aligned
      
%  ** Layer 2 (lower) parameters
      L2_depth = 150. ; L2_thick = 60. ; L2_dip = 30.0 ; % layer geometry
      L2_aaz = 0.0 ;  % a-axis azimuth (rel. to down dip direction)
      L2_faln = 0.3 ; % fraction aligned 

      report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
                   L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) ; 

%  ** imaging parameters
      baz = [0:1:360] ; % 0 is down dip direction (perp. to strike)
      inc = -ones(size(baz)).*90 + aoi ;

%  ** calculate distances
      [dist1]=distance_in_dipping_layer(L1_dip,aoi,L1_thick,baz) ;
      [dist2]=distance_in_dipping_layer(L2_dip,aoi,L2_thick,baz) ;
      
%  ** load anisotropy, and generate an isotropic version of it. 
      [Cani,rh] = MS_elasticDB('olivine') ;
      [Ciso] = MS_decomp(MS_axes(Cani)) ;
      Cani = MS_rot3(Cani,90,0,0) ; % orientation for dry upper mantle

%  ** generate layer elasticities:
%        This is a Voigt-Reuss-Hill average of the appropriately rotated olivine
%        tensor and its isotropic equivalent.
      [L1_C,~] = MS_VRH([L1_faln 1-L1_faln],...
         MS_rot3(Cani,0,-L1_dip,L1_aaz,'order',[3 2 1]),rh, Ciso, rh) ;

      [L2_C,~] = MS_VRH([L2_faln 2-L2_faln],...
         MS_rot3(Cani,0,-L2_dip,L2_aaz,'order',[3 2 1]),rh, Ciso, rh) ;

%  ** interrogate elasticities to generate splitting parameters for each layer.
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L1_C, rh, inc, baz ) ;   
      fast1 = MS_unwind_pm_90((baz - pol')) ; % geog. reference frame
      tlag1 = dist1./vs2' - dist1./vs1' ;
            
      [ pol, ~, vs1, vs2, ~, ~, ~ ] = MS_phasevels( L2_C, rh, inc, baz ) ;   
      fast2 = MS_unwind_pm_90((baz - pol')) ; % geog. reference frame
      tlag2 = dist2./vs2' - dist2./vs1' ;

%  ** calculate the effective splitting between 2 layers    
      fast_eff = zeros(size(baz)) ;
      tlag_eff = fast_eff ;
      for i = 1:length(baz)
         [fast_eff(i),tlag_eff(i)] = ...
            MS_effective_splitting_N(0.125,baz(i), ...
            [fast2(i) fast1(i)],[tlag2(i) tlag1(i)]) ;
      end

%  ** make a figure of the results
      scrsz = get(0,'ScreenSize');
      figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)*0.9]) ;

%  ** lag times
      subplot(2,1,1)
      plot(baz,tlag1,'r--') ; hold on
      plot(baz,tlag2,'g--') ;
      plot(baz,tlag_eff,'k-','LineWidth',1.5)
      axis([0 360 0 4])
      xlabel('Polarisation (relative to downdip direction)')
      ylabel('Lag times (s)')
      legend('Upper layer','Lower layer','Total') ;

%  ** fast directions
      subplot(2,1,2)
      plot(baz,fast1,'r--') ; hold on
      plot(baz,fast2,'g--') ;
      plot(baz,fast_eff,'k-','LineWidth',1.5)
      axis([0 360 -90 90])
      xlabel('Polarisation (relative to downdip direction)')
      ylabel('Fast shear-wave orientation (degree)')      
      legend('Upper layer','Lower layer','Total') ;
      
return

function [dist]=distance_in_dipping_layer(dip,aoi,thick,baz)

%  ** calculate apparent dip
      alp = atand(tand(dip) .* sind(baz+90)) ;

%  ** calculate distances
      gam = 90 - alp - aoi ;
      bet = 90 + alp ;
      dist = thick .* sind(bet) ./ sind(gam) ;
return

function report_model(L1_depth, L1_thick, L1_dip, L1_aaz, L1_faln, ...
                      L2_depth, L2_thick, L2_dip, L2_aaz, L2_faln) 

%  ** Banner
      fprintf('\n** SKS Shear-wave splitting in two layer anisotropic models.\n\n')

%  ** Report layer details
      fprintf('Layer 1: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
         L1_depth, L1_thick, L1_dip)
      fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
         L1_aaz,L1_faln)
      fprintf('\n')

      fprintf('Layer 2: depth to top = %5.2fkm, thickness = %5.2f km, dip = %6.2f deg.\n', ...
         L2_depth, L2_thick, L2_dip)
      fprintf('         anisotropy: a-axis azimuth = %6.2f deg, v.f. aligned  = %4.1f\n', ...
         L2_aaz,L2_faln)

      fprintf('\n** Plotting backazimuthal variation of splitting parameters ...\n')


return