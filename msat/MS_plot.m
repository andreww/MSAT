% MS_PLOT - Plot phasevels/anisotropy on pole figures.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given an elasticity matrix and density, produce pole figures showing 
%     the P- and S-wave anisotrpy. 
%
%  % MS_plot(C, rh, ...)
%
% Usage: 
%     MS_plot(C, rh)                    
%         Produce three pole figures showing P-wave velocity, S-wave
%         anisotropy and fast S-wave polarisation direction
%
%     MS_plot(C, rh, ...)                    
%          Further arguments are EVALed (most useful to change various 
%          options, see source-code for details).
%
% See also: MS_SPHERE, MS_PHASEVELS

% Copyright (c) 2011, James Wookey and Andrew Walker
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


%=========================================================================
function MS_plot(C,rh,varargin)
%=========================================================================

%  ** Set defaults, these can be overriden in the function call
%  ** configure contouring options
      cvect = 10 ;     % number of contours
      VPcvect = NaN ;  % set these for different number of contours for VP
      AVScvect = NaN ; % and AVS plots

%  ** configure font options
      fntsz = 12;

%  ** configure colormap options
      cmap = jet(64) ;
      icmapflip = 1 ; % reverse the sense of the colourscale

%  ** process the optional arguments
      for i=1:length(varargin) 
         eval(varargin{i}) ; 
      end    

      if isnan(VPcvect)
         VPcvect = cvect ;
      end   
      if isnan(AVScvect)
         AVScvect = cvect ;
      end   
      
      % check the inputs: C
      assert(MS_checkC(C)==1, 'MS:PLOT:badC', 'MS_checkC error MS_plot') ;

%  ** buggy MATLAB contourf (version 7.1-7.3)
%
%     in these versions of MATLAB, the contourf routine doesn't seem able to
%     to cope with the spherical geometry surfaces. I therefore use 
%     contouf('v6',...) instead
%
      vstr = ver ;
      if vstr(1).Version >= 7.1 & vstr(1).Version <= 7.3
         buggyMATLAB = 1;
      else
         buggyMATLAB = 0;
      end   

%  ** set the viewing angle
      view_angle = [-90,90] ; % 1-North 2-West 3-Out of screen

      rad = pi./180 ;
      deg = 180./pi ;
      
      % Set up inc-az grids...
      [INC,AZ] = meshgrid([90:-6:0],[0:6:360]) ;
      
      % Invoke MS_phasevels to get wave velocs etc.
      [~,~,vs1,vs2,vp, S1P] = MS_phasevels(C,rh,...
        reshape(INC,61*16,1),reshape(AZ,61*16,1));
    
      % reverse so sph2cart() works properly
      AZ = -AZ;
      
      % Reshape results back to grids
      VS1 = reshape(vs1,61,16);
      VS2 = reshape(vs2,61,16);
      VP =  reshape(vp,61,16);
      VS1_x = reshape(S1P(:,1),61,16);
      VS1_y = reshape(S1P(:,2),61,16);
      VS1_z = reshape(S1P(:,3),61,16);
      
%  ** output average velocities
      VPiso=mean(mean(VP)) ;
      VSiso=mean([mean(mean(VS1)) mean(mean(VS2))]) ;

      fprintf('Isotropic average velocities: VP=%f, VS=%f\n',VPiso,VSiso) ;

%  ** Prepare window
      figure('Position',[1 1 1400 400]) ;

%  ** Setup a seismic colourmap (i.e. red->green->blue)
      if icmapflip
         cmap = flipud(cmap) ;
      end   
%  ** generate X/Y matrices for plotting
      [X,Y,Z] = sph2cart(AZ.*rad,INC.*rad,ones(size(AZ))) ;

%-------------------------------------------------------------------------------
%  ** VP anisotropy plot 
%-------------------------------------------------------------------------------
      subplot(1,3,1)

%  ** get VPmin and VPmax values and positions      
      [imin,jmin,VPmin,imax,jmax,VPmax] = minmax2d(VP) ;
      
      [VPmin_x,VPmin_y]=sph2cart(AZ(imin,jmin)*rad,INC(imin,jmin)*rad,1) ;   
      [VPmax_x,VPmax_y]=sph2cart(AZ(imax,jmax)*rad,INC(imax,jmax)*rad,1) ;   

      
      if (VPmin ~= VPmax)
         if buggyMATLAB
            [h1,h2]=contourf('v6',X,Y,VP,VPcvect) ;
            for j=1:length(h2)
               set(h2(j),'LineStyle','none')
            end
         else
            contourf(X,Y,VP,VPcvect,'LineStyle','none') ; 
         end
      else
         surf(X,Y,zeros(size(VP)),VP) ;
         shading flat ;
      end   
      view(view_angle)
      colormap(cmap) ;
      daspect([1 1 1]);
      colorbar('FontSize',fntsz) ;
      axis off
      title('V_P (km/s)','FontSize',fntsz+4,'FontWeight','bold') ;

      hold on

%  ** mark the max. min. values 
      h=plot(VPmax_x,VPmax_y,'ws') ;
      set(h,'MarkerFaceColor','black');
      h=plot(VPmin_x,VPmin_y,'wo') ;
      set(h,'MarkerFaceColor','black');

%  ** add some information to the plot            
      VPlabel1 = sprintf('Min. V_P =%6.2f, max. V_P =%6.2f',VPmin,VPmax) ;
      text(-1.15,0.8,VPlabel1,'FontSize',fntsz,'FontWeight','bold') ;
      VPmean = (VPmax+VPmin)./2.0 ;
      VPani = (VPmax-VPmin)/VPmean .* 100 ;
      VPlabel2 = sprintf('Anisotropy =%6.1f%%',VPani) ;
      text(-1.3,0.8,VPlabel2,'FontSize',fntsz,'FontWeight','bold') ;
      
%-------------------------------------------------------------------------------
%  ** dVS anisotropy plot 
%-------------------------------------------------------------------------------
      subplot(1,3,2)

      dVS =  (VS1-VS2) ;
      VSmean = (VS1+VS2)./2.0 ;
      AVS = 100.0*(dVS./VSmean) ;

%  ** get AVSmin and AVSmax values and positions      
      [imin,jmin,AVSmin,imax,jmax,AVSmax] = minmax2d(AVS) ;
      [AVSmin_x,AVSmin_y]=sph2cart(AZ(imin,jmin)*rad,INC(imin,jmin)*rad,1) ;   
      [AVSmax_x,AVSmax_y]=sph2cart(AZ(imax,jmax)*rad,INC(imax,jmax)*rad,1) ;   
      
      if (AVSmin ~= AVSmax)
         if buggyMATLAB
            contourf('v6',X,Y,AVS,AVScvect) ;
            [h1,h2]=contourf('v6',X,Y,AVS,AVScvect) ;
            for j=1:length(h2)
               set(h2(j),'LineStyle','none')
            end
         else
            contourf(X,Y,AVS,AVScvect,'LineStyle','none') ;
         end
      else
         surf(X,Y,zeros(size(AVS)),AVS) ;
         shading flat ;
      end   

      colormap(cmap) ;
      view(view_angle)
      daspect([1 1 1]);

      colorbar('FontSize',fntsz) ;
      caxis_save = caxis;
      axis tight; axis off
      title('dV_S (%)','FontSize',fntsz+4,'FontWeight','bold')
      
      hold on

%  ** mark the max. min. values 
      h=plot(AVSmax_x,AVSmax_y,'ws') ;
      set(h,'MarkerFaceColor','black');
      h=plot(AVSmin_x,AVSmin_y,'wo') ;
      set(h,'MarkerFaceColor','black');

%  ** add some information to the plot            
      AVSlabel1 = sprintf('V_S anisotropy, min =%6.2f, max =%6.2f',AVSmin,AVSmax) ;
      text(-1.15,1.0,AVSlabel1,'FontSize',fntsz,'FontWeight','bold') ;

%-------------------------------------------------------------------------------
%  ** Fast shear-wave polarisation plot
%-------------------------------------------------------------------------------
      subplot(1,3,3)
      
      if buggyMATLAB
         [h1,h2]=contourf('v6',X,Y,AVS,AVScvect) ;
         for j=1:length(h2)
            set(h2(j),'LineStyle','none')
         end
      else
         contourf(X,Y,AVS,AVScvect,'LineStyle','none') ;
      end
      
      colormap(cmap) ;
      colorbar('FontSize',fntsz) ;

      hold on ;

%  ** transform vectors
      [VS1_x,VS1_y,VS1_z] = vnormalise2(VS1_x,VS1_y,VS1_z) ;
      [XN,YN,ZN] = vnormalise2(X,Y,Z) ;

      A = zeros(3,61,16) ;
      B = zeros(3,61,16) ;
     
      A(1,:,:) = XN ;
      A(2,:,:) = YN ;
      A(3,:,:) = ZN ;
                  
      B(1,:,:) =  VS1_x ;
      B(2,:,:) =  VS1_y ;
      B(3,:,:) =  VS1_z ;
     
      C=cross(A,B) ;
      D=cross(A,C) ;
      
      VS1R_x(:,:) = D(1,:,:) ;
      VS1R_y(:,:) = D(2,:,:) ;
      VS1R_z(:,:) = D(3,:,:) ;
      
      [VS1R_x,VS1R_y,VS1R_z] = vnormalise2(VS1R_x,VS1R_y,VS1R_z) ;
 
%  ** define subset of polarisations to plot      
      cl = [1,3,5,7,9,11,15] ;
      drw = [60 10 5 3 2 2 2] ;

      ii=0 ;
      ip=0 ;
      np=136 ;
      ind1 = zeros(1,np) ;
      ind2 = zeros(1,np) ;
      
      for iinc=cl
         ii=ii+1 ;
         for iaz=1:drw(ii):61
            ip = ip + 1;
            ind1(ip) = iaz ;
            ind2(ip) = iinc ;
         end
      end   
      
%  ** rotate the particle motion vector so the normal to sphere is vertical
      for ip = 1:np
         [VS1R_xR(ind1(ip),ind2(ip)),VS1R_yR(ind1(ip),ind2(ip)),VS1R_zR(ind1(ip),ind2(ip))] = ...
                    rotate_pm_vector(...
         VS1R_x(ind1(ip),ind2(ip)),VS1R_y(ind1(ip),ind2(ip)),VS1R_z(ind1(ip),ind2(ip)),...
         AZ(ind1(ip),ind2(ip)),INC(ind1(ip),ind2(ip)));          
      end   

%  ** form the subsets
      X2 = zeros(1,np) ;
      Y2 = zeros(1,np) ;
      Z2 = zeros(1,np) ;
      U2 = zeros(1,np) ;
      V2 = zeros(1,np) ;
      W2 = zeros(1,np) ;

      for i=1:length(ind1) ;
         X2(i)=X(ind1(i),ind2(i)) ;
         Y2(i)=Y(ind1(i),ind2(i)) ;
         Z2(i)=Z(ind1(i),ind2(i)) ;
         U2(i)=VS1R_xR(ind1(i),ind2(i)) ;
         V2(i)=VS1R_yR(ind1(i),ind2(i)) ;
         W2(i)=VS1R_zR(ind1(i),ind2(i)) ;
      end


%  ** plot the vectors (changed to 2D to allow plotting with contourf surface)     
%      h=quiver3(X2,Y2,Z2,U2,V2,W2,0.18,'k.') ;
      h=quiver(X2,Y2,U2,V2,0.18,'w.') ;
      set(h,'LineWidth',3.0) ;

      h=quiver(X2,Y2,-U2,-V2,0.18,'w.') ;
      set(h,'LineWidth',3.0) ;

      h=quiver(X2,Y2,U2,V2,0.18,'k.') ;
      set(h,'LineWidth',2.0) ;

%      h=quiver3(X2,Y2,Z2,-U2,-V2,-W2,0.18,'k.') ;
      h=quiver(X2,Y2,-U2,-V2,0.18,'k.') ;
      set(h,'LineWidth',2.0) ;

%		plot(X2,Y2,'k.','MarkerSize',2.0) ;

      view(view_angle);
      daspect([1 1 1]);

      axis tight; axis off
      title('Fast-shear polarisation','FontSize',fntsz+4,'FontWeight','bold')
      
   return
%===============================================================================

%===============================================================================
   function [AN,BN,CN] = vnormalise2(A,B,C)
%===============================================================================
%
%     normalise a 3 vector
%
      VMAG = sqrt( A.^2 + B.^2 + C.^2 )  ; 

      AN = A./VMAG ;
      BN = B./VMAG ;
      CN = C./VMAG ;
      
   return
%===============================================================================

%===============================================================================
   function [xr,yr,zr] = rotate_pm_vector(x,y,z,az,in)
%===============================================================================
%           
%     Rotate the particle motion vector so propagation direction is vertical
%     (to plot in a 2D way)
%
      rad = pi./180;
      V=[x,y,z]' ;
      
%  ** build up rotation matrices      
      gamma = az .* rad ;
      
%  ** first rotation (about Z)
      R1 = [  cos(gamma), sin(gamma), 0 ;...
             -sin(gamma), cos(gamma), 0 ;...
                 0      ,    0      , 1 ] ;
                    
%  ** second rotation (about Y)
      beta = (90-in).*rad ;
      R2 = [  cos(beta), 0 , -sin(beta) ;...
                 0     , 1 ,     0      ;...
              sin(beta), 0 ,  cos(beta) ] ;

%  ** third rotation (about Z, reversing first rotation)
      gamma = -az .* rad ;
      R3 = [  cos(gamma), sin(gamma), 0 ;...
             -sin(gamma), cos(gamma), 0 ;...
                 0      ,    0      , 1 ] ;

%  ** apply rotation
      VR = R3 * R2 * (R1 * V) ;
      xr = VR(1) ; 
      yr = VR(2) ;
      zr = VR(3) ;

   return
%===============================================================================

%===============================================================================
   function [imin,jmin,Zmin,imax,jmax,Zmax] = minmax2d(Z)
%===============================================================================
%           
%     Find the maximum and minimum values in a (2d) matrix and return their
%     values and incidces
%
      [A,I] = min(Z) ; [Zmin,II] = min(A) ; imin = I(II) ; jmin = II ;
      [A,I] = max(Z) ; [Zmax,II] = max(A) ; imax = I(II) ; jmax = II ;

   return
%===============================================================================
