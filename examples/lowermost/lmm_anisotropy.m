% LMM_ANISOTROPY.M - Demonstrate the use of MSAT to explore complicated models. 
%
%  This code demonstrates the use of MSAT to calculate derived properties from
%  a set of elastic constants. The model used is a model of lowermost mantle
%  anisotropy from Walker et al (2011), specifically 'TX2008.V1.P010.dat'. See
%  http://www1.gly.bris.ac.uk/CoMITAC/software.htm for details. 
%
%  Reference:
%     A. M. Walker, A. M. Forte, J. Wookey, A. Nowacki and J.-M. Kendall (2011). 
%     Elastic anisotropy of D" predicted from global models of mantle flow. 
%     Geochem. Geophys. Geosyst., 12, Q10006, doi:10.1029/2011GC003732
%

% (C) James Wookey and Andrew Walker, 2011

function lmm_anisotropy()

%  ** Load one of the D" elastic models from Walker et al (G^3, 2011)
      [ n, rs, lats, lons, Cs ] = ...
         az_read_Si_file( 'TX2008.V1.P010.dat' ) ;

%  ** generate indices
      ila = 1+(lats - -90)/5 ;
      ilo = 1+(lons - -180)/5  ;
      
%  ** prepare grid
      [LON,LAT] = meshgrid(-180:5:180,-90:5:90) ;
      Cani = zeros(size(LON)) ;
      PercH  = zeros(size(LON)) ;

%  ** process the all elastic constants, calculate derived properties. 
      for i=1:length(lats)
         Chere = squeeze(Cs(:,:,i)) ;

%     ** amount of overall anisotropy
         Cani(ila(i),ilo(i))=MS_anisotropy(Chere) ;

%     ** proportion of hexagonality
         Chere = squeeze(Cs(:,:,i)) ;
         ChereR = MS_axes(Chere) ;
         [C_iso,C_hex,C_tet,C_ort,C_mon,C_tri] = MS_decomp(ChereR);
         P = MS_norms(ChereR,C_iso,C_hex,C_tet,C_ort,C_mon,C_tri) ;
         SumA = sum(P(2:end)) ;
         PercH(ila(i),ilo(i)) = P(2)./SumA ;
      end 

%  ** plot each of the derived properties on a globe
      plot_cprop(LON,LAT,Cani,[0 0.1],'Universal Elastic Anisotropy Index') ;
      plot_cprop(LON,LAT,PercH,[0 1.0],'Proportion of Hexagonality') ;

end

function plot_cprop(LON,LAT,Cprop,Cax,tstr,fn)
% Plot as a globe. 

%  ** fill out the halo points
      Cprop(:,1) = Cprop(:,end) ;
      Cprop(1,:) = mean(Cprop(2,:)) ;
      Cprop(end,:) = mean(Cprop(end-1,:)) ;

%  ** interpolate to smooth
      [ILON,ILAT] = meshgrid(-180:1:180,-90:1:90) ;

      CPI = interp2(LON,LAT,Cprop,ILON,ILAT,'cubic') ;
      
%  ** convert to spherical geometry
      RAD = zeros(size(CPI))+(6371-2891) ;
      [XX,YY,ZZ] = sph2cart(ILON*(pi/180),ILAT*(pi/180),RAD) ;

%  ** plot the data
      figure('Position',[1 1 1000 1000]) ;
      surf(XX,YY,ZZ,CPI)
      daspect([1 1 1]) ;
      shading interp ;
      lighting gouraud ;
      axis tight; axis off

      cmap = flipud(jet(64)) ;
      colormap(cmap) ;
      caxis(Cax) ;
      colorbar ;
      
%  ** add coastlines      
      load lmm_coasts ;
      [xc,yc,zc] = sph2cart(clon.*(pi/180),clat.*(pi/180), ...
                           zeros(size(clon))+(6371-2890)) ;
      hold on
      plot3(xc,yc,zc,'w-','LineWidth',1) ;
      view(80,20)   
      title(tstr,'FontSize',18)      
end

function [ n, rs, lats, lons, Cs ] = az_read_Si_file( filename )
% Function to read SI data from G^3 paper.

   all_data = dlmread(filename);
   n = size(all_data, 1);
   lats = all_data(:,1);
   lons = all_data(:,2);
   rs = all_data(:,3);
   Ci(n,1:6,1:6) = 0.0;
   Ci(:,1,:) = all_data(:,4:9);
   Ci(:,2,2:6) = all_data(:,10:14);
   Ci(:,3,3:6) = all_data(:,15:18);
   Ci(:,4,4:6) = all_data(:,19:21);
   Ci(:,5,5:6) = all_data(:,22:23);
   Ci(:,6,6) = all_data(:,24);

   Cs(1:6,1:6,n) = 0.0;
   for i=1:n
       thisC(1:6,1:6) = Ci(i,:,:);
       thisC(2,1) = thisC(1,2);
       thisC(3,1:2) = thisC(1:2,3);
       thisC(4,1:3) = thisC(1:3,4);
       thisC(5,1:4) = thisC(1:4,5);
       thisC(6,1:5) = thisC(1:5,6);
       Cs(1:6,1:6,i) = thisC(:,:);
   end

end