

function MS_polefigure(C,rh,varargin)

    [INC,AZ] = meshgrid([90:-6:0],[0:6:360]) ;
    [~,~,vs1,vs2,vp, S1P] = MS_phasevels(C,rh,...
        reshape(INC,61*16,1),reshape(AZ,61*16,1));
    %pol = reshape(pol,61,16);
    %avs = reshape(avs,61,16);
    vs1 = reshape(vs1,61,16);
    vs2 = reshape(vs2,61,16);
    vp =  reshape(vp,61,16);
    S1_X = reshape(S1P(:,1),61,16);
    S1_Y = reshape(S1P(:,2),61,16);
    S1_Z = reshape(S1P(:,3),61,16);
    
    if nargin>2
        plot_vpg(INC,AZ,vp,vs1,vs2,S1_X,S1_Y,S1_Z,varargin{:});
    else 
        plot_vpg(INC,AZ,vp,vs1,vs2,S1_X,S1_Y,S1_Z);
    end
return

% --------------------------
%  MATLAB function plot_vpg
% --------------------------
%
%  Usage: 
%     plot_vgp('file.vpg',...) - where file.vpg is a VPG file produced
%     by ematrix5 (or other program). Other arguments are EVALed (most
%     useful to change various options, see source-code for details). 
%
%  Purpose:    
%     Plot pole figures from EMATRIX5/6 output (VPG files)
%     At the moment phase velocities only are plotted
%
%  Bugs/known problems:       
%     None (workaround for contourf bug)
%
%  Written by James Wookey
%  Department of Earth Sciences, University of Bristol, UK
%  Version 1.2
%  Incept: March 2004
%  Last update: March 2007
%

%  Major changes log
%
%  v0.991 - fixed difference in plots
%  v0.992 - minor change to view options  
%  v0.993 - workaround for MATLAB v7.1 contourf bug
%  v0.994 - fixed a bug in the polarisation plotting
%  v0.995 - changed the subset of polarisations plotted to more even 
%           distribution on the sphere
%  v0.996 - Extended workaround to be standard, since 7.2 doesn't fix it.
%  v0.997 - Output average velocities also
%  v0.998 - Added a mechanism for configuring things; including plot orientation
%  v0.999 - Added a 'high resolution' option for publication quality graphics
%           and did some font tweaking
%  v1.000 - Restored contourf plots in buggyMATLAB using the contouf('v6'...)
%           syntax. Added configuration options for contouring.
%           Finally had to make code version 1. 
%  v1.1   - Finally got contourf plot behind the polarisation vectors working.
%           Other minor tweaks. This produces better vector image output for
%           EPS. Removed the highres and iport options, as these are a 
%           little redundant with the other changes. 
%  v1.2   - Reinstated proper contourf function since this got fixed in V7.4
%  
%-------------------------------------------------------------------------------
%  This software is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
   function plot_vpg(INC,AZ,VP,VS1,VS2,VS1_x,VS1_y,VS1_z,varargin)
%-------------------------------------------------------------------------------

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

%  ** check input            
%      if isempty(fname)
%         error('No filename specified') ;
%      end      

%  ** load the VPG file      
%      file = local_load_vpgfile(fname) ;

%  ** pull the bits out of the structure (for ease)
%      INC    = file.INC ;
%      AZ     = -file.AZ ; 
      AZ = -AZ;% reverse so sph2cart() works properly
%      VP     = file.VP ;
%      VS1    = file.VS1 ;
%      VS2    = file.VS2 ;
%      VP_x   = file.VP_x ;
%      VP_y   = file.VP_y ;
%      VP_z   = file.VP_z ;
 %     VS1_x  = file.VS1_x ;
 %     VS1_y  = file.VS1_y ;
 %     VS1_z  = file.VS1_z ;
 %     VS2_x  = file.VS2_x ;
 %     VS2_y  = file.VS2_y ;
 %     VS2_z  = file.VS2_z ;
 %     VPG_x  = file.VPG_x ;
 %     VPG_y  = file.VPG_y ;
 %     VPG_z  = file.VPG_z ;
 %     VS1G_x = file.VS1G_x ;
 %     VS1G_y = file.VS1G_y ;
 %     VS1G_z = file.VS1G_z ;
 %     VS2G_x = file.VS2G_x ;
 %     VS2G_y = file.VS2G_y ;
 %     VS2G_z = file.VS2G_z ;

%  ** clear the structure (to save memory)     
      clear vpgfile
      
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

%===============================================================================
   function [vpgfile] = local_load_vpgfile(fname)
%===============================================================================

%  ** check input            
      if isempty(fname)
         fprintf('ERROR: LOAD_VPGFILE: No filename specified\n') ;
         return
      end      
%  ** open and read file          
      fid=fopen(fname,'rt') ;
      vpgfile.TITLE=fscanf(fid,'%s',1) ;
      all_data=fscanf(fid,'%f',[21 inf])' ;
      [vpgfile.INC,vpgfile.AZ] = meshgrid([90:-6:0],[0:6:360]) ;

%  ** build AZI x INC matrices in structure          

%  ** phase velocity
      vpgfile.VP       = reshape(all_data(:,01),61,16) ;
      vpgfile.VS1      = reshape(all_data(:,02),61,16) ;
      vpgfile.VS2      = reshape(all_data(:,03),61,16) ;

%  ** particle motion (x-direction)
      vpgfile.VP_x     = reshape(all_data(:,04),61,16) ;
      vpgfile.VS1_x    = reshape(all_data(:,05),61,16) ;
      vpgfile.VS2_x    = reshape(all_data(:,06),61,16) ;

%  ** particle motion (y-direction)
      vpgfile.VP_y     = reshape(all_data(:,07),61,16) ;
      vpgfile.VS1_y    = reshape(all_data(:,08),61,16) ;
      vpgfile.VS2_y    = reshape(all_data(:,09),61,16) ;

%  ** particle motion (z-direction)
      vpgfile.VP_z     = reshape(all_data(:,10),61,16) ;
      vpgfile.VS1_z    = reshape(all_data(:,11),61,16) ;
      vpgfile.VS2_z    = reshape(all_data(:,12),61,16) ;

%  ** group velocities
      vpgfile.VPG_x    = reshape(all_data(:,13),61,16) ;
      vpgfile.VPG_y    = reshape(all_data(:,14),61,16) ;
      vpgfile.VPG_z    = reshape(all_data(:,15),61,16) ;
      vpgfile.VS1G_x   = reshape(all_data(:,16),61,16) ;
      vpgfile.VS1G_y   = reshape(all_data(:,17),61,16) ;
      vpgfile.VS1G_z   = reshape(all_data(:,18),61,16) ;
      vpgfile.VS2G_x   = reshape(all_data(:,19),61,16) ;
      vpgfile.VS2G_y   = reshape(all_data(:,20),61,16) ;
      vpgfile.VS2G_z   = reshape(all_data(:,21),61,16) ;


   return
%===============================================================================
