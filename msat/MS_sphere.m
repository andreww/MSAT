% CIJ_sphere : Plot a spherical plot of phasevels/anisotropy from a set of 
%              elastic constants.
%   
% Usage: CIJ_sphere(CC,rh,mode,...)
%     mode can be 'p'(-wave) or 's'(-wave).
%     extra arguments are EVALed; parameters that can be (re)set are:
%
%        cmap = jet(64) ; % colourmap
%        icmapflip = 1 ; % reverse the sense of the colourscale     
%        FSWTickLength=0.08 ; % Fast shear wave vector length
%
%        plotlabels = 1 ;
%        plotaxes = 1 ;
%        nofig = 0 ; % suppress the figure command, so plot can be sub-plotted
%        nocbar = 0 ; % suppress the colorbar
%        cax = NaN ; % colour axis, define to use - default is auto
%        dirs = [] ;  % directions of interest, this should be a vector of azi, 
%                       inc pairs, e.g., [az1 in1 az2 in2 az3 in3 az4 in4]
%        dlen = 1.5 ; % length of doi vectors 
%
%
% (C) James Wookey, December 2007.	

%-------------------------------------------------------------------------------
%  This software is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%-------------------------------------------------------------------------------

function MS_sphere(CC,rh,mode,varargin);

cmap = jet(64) ;
icmapflip = 1 ; % reverse the sense of the colourscale
FSWTickLength=0.08 ;
dirs = [] ;
dlen = 1.5 ;
plotlabels = 1 ;
plotaxes = 1 ;
nofig = 0 ;
nocbar = 0 ;
cax = NaN ; % define to use

% check the inputs: CC
if isnumeric(CC)
   if sum(size(CC)==6) ~= 2
      error('CC must be a 6x6 matrix')
   end
else   
   error('CC must be a 6x6 matrix')
end

[isgood] = MS_checkC(CC) ;

% check the inputs: rh
if isnumeric(rh)
   if sum(size(rh)==1) ~= 2
      error('rh must be a scalar')
   end
else
   error('rh must be a scalar')
end

if ischar(mode)
   if sum(strcmp({'p','s'},lower(mode)))~=1
      error('Mode must be ''S''(-wave) or ''P''(-wave) ')
   end
else
   error('Mode must be ''S''(-wave) or ''P''(-wave) ')
end     

% process the optional arguments, these are EVAL'ed in turn
for i=1:length(varargin) 
   eval(varargin{i}) ; 
end    

% check input pars
if mod(length(dirs),2)~=0
   error('Plot directions must be in (azi,inc) pairs') ;
end

%  ** Setup a seismic colourmap (i.e. red->green->blue)
if icmapflip
   cmap = flipud(cmap) ;
end   

% unpack the specified directions of interest
daz=dirs(1:2:length(dirs)).*pi/180 ;
din=dirs(2:2:length(dirs)).*pi/180 ;
ndir=length(daz) ;

if (nofig~=1)
   figure('position',[0 0 800 800])
end

%% load the triangulation
load SPHTR3.mat ;
[SF,avs,vs1,vs2,vp] = phasevels_local(CC,rh,inc,az) ;

if strcmp(lower(mode),'p')
   trisurf(faces,x,y,z,vp) ;
else
   trisurf(faces,x,y,z,avs) ;
end

shading interp
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]) ;
colormap(cmap)
daspect([1 1 1]) ;
hold on

if ~isnan(cax), caxis(cax), end

%% plot the shear-wave polarisation
if strcmp(lower(mode),'s') 
   load SPHTR2.mat ; % fewer sampled directions

%  find the closest points to each specified direction   
   [xd,yd,zd] = sph2cart(-daz,din,ones(size(daz))) ;

   for idir=1:ndir
%  ** calculate an array of distances to each point
      dist = sqrt((x-xd(idir)).^2+(y-yd(idir)).^2+(z-zd(idir)).^2) ;
      [val,ind] = min(dist) ;
      x(ind) = xd(idir) ;
      y(ind) = yd(idir) ;
      z(ind) = zd(idir) ;
      az(ind) = daz(idir)*180/pi ;
      inc(ind) = din(idir)*180/pi ;
   end
   
   [SF,avs,vs1,vs2,vp] = phasevels_local(CC,rh,inc,az) ;
   %% calculate PM vectors
   nv = length(x) ;
   for iv=1:nv
      XI=1.01.*[x(iv) y(iv) z(iv)] ;
      X1 = XI-FSWTickLength.*SF(iv,:);
      X2 = XI+FSWTickLength.*SF(iv,:);
      plot3(XI(1),XI(2),XI(3),'ko','MarkerSize',4,'MarkerFaceColor','k')
      plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)],'k-','LineWidth',2)
   end
end

%% draw the axes
if plotaxes
   plot3([-1.5 1.5],[0 0],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','--')
   plot3([0 0],[-1.5 1.5],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','--')
   plot3([0 0],[0 0],[-1.5 1.5],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','--')
   if plotlabels
      text(1.5,0,0,'X1','FontSize',12) ;
      text(0,1.5,0,'X2','FontSize',12) ;
      text(0,0,1.5,'X3','FontSize',12) ;
   end
end

for i=1:ndir
   [xd,yd,zd] = sph2cart(-daz(i),din(i),dlen) ;
   plot3([0 xd],[0 yd],[0 zd],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',3) ;
   if plotlabels
      text(xd*1.01,yd*1.01,zd*1.01,...
           sprintf('[%6.1f,%6.1f]', ...
           daz(i).*180./pi,din(i).*180./pi),'FontSize',12) ;
   end        
end   

%% final commands

axis off
daspect([1 1 1]) ;
if (nocbar~=1)
colorbar('NorthOutside','FontSize',14,'FontWeight','bold')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------ SUBROUTINES ----------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=======================================================================================  
function [VR] = V_rot3(V,alp,bet,gam)

%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;

VR = V * RR ;
 
return
%=======================================================================================  

%=======================================================================================  
	function [X] = cart2(irev,inc,azm)
%=======================================================================================  
%c convert from spherical to cartesian co-ordinates
%c north x=100  west y=010 up z=001
%c irev=+1 positive vector x
%c irev=-1 negative vector x
	caz=cosd(azm)  ;
   saz=sind(azm)  ;
   cinc=cosd(inc) ;
   sinc=sind(inc) ;
   X=[caz*cinc -saz*cinc sinc] ;
%c normalise to direction cosines
   r=sqrt(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)) ;
   
	X = X./r ;
	if(irev == -1), X = -X;, end
   return
%=======================================================================================  

%=======================================================================================  
	function [V,EIGVEC]=velo(X,rh,C)
%=======================================================================================  
%C PHASE-VELOCITY SURFACES IN AN ANISOTROPIC MEDIUM
%C revised April 1991
%C     X(3) - DIRECTION OF INTEREST
%C     RHO - DENSITY
%C     V - PHASE VELOCITIES (1,2,3= P,S,SS)
%C     EIGVEC(3,3) - eigenvectors stored by columns
%      implicit real*8(a-h,o-z)
%      dimension v(3),c(6,6),t(3,3),x(3),eval(3),eigvec(3,3),ei(3,3)
%      dimension ijkl(3,3)
%      data ((ijkl(i,j),j=1,3),i=1,3)/1,6,5,6,2,4,5,4,3/      
      
		ijkl = [1,6,5; ...
		        6,2,4; ...
		        5,4,3] ;
%c form symmetric matrix tik=cijkl*xj*xl
		for i=1:3
      	for k=1:3
	         T(i,k)=0.0 ;
      		for j=1:3
      			for l=1:3
      				m=ijkl(i,j) ;
      				n=ijkl(k,l) ;
      				T(i,k)=T(i,k)+C(m,n).*X(j).*X(l) ;
					end
				end
			end
		end
% determine the eigenvalues of symmetric tij
      [EIVEC EIVAL] = eig(T) ;

% calculate velocities and sort
		V_RAW = (sqrt([EIVAL(1,1) EIVAL(2,2) EIVAL(3,3)]./rh))*10. ;
		[V IND] = sort(V_RAW,2,'descend') ;
		EIGVEC = EIVEC ; % for dimensioning
		for i=1:3
			EIGVEC(:,i) = EIVEC(:,IND(i)) ;
		end

      return

%=======================================================================================  

% phasevels - calculate the phase velocity details for a set of elastic 
%                 constants - localised version
%
% (C) James Wookey, May 2007.	
%   
% Based on EMATRIX6 by D. Mainprice. Re-coded in MATLAB by JW
%
% Reference: Mainprice D. (1990). An efficient
%            FORTRAN program to calculate seismic anisotropy from
%            the lattice preferred orientation of minerals.
%            Computers & Gesosciences,vol16,pp385-393.
%
%
function [SF,avs,vs1,vs2,vp] = phasevels_local(C,rh,inc,azi)


		if (length(inc)~=length(azi))
			error('AZI and INC must be scalars or vectors of the same dimension');
		end	

%  ** convert GPa to MB file units (Mbars), density to g/cc
      C(:,:) = C(:,:) * 0.01 ;
      rh = rh ./ 1e3 ;
      
		avs = zeros(size(azi)) ;
		vp = zeros(size(azi)) ;
		vs1 = zeros(size(azi)) ;
		vs2 = zeros(size(azi)) ;
		pol = zeros(size(azi)) ;
		S1 = zeros(length(azi),3) ;
		S1P = zeros(length(azi),3) ;
      
%	** start looping
	for ipair = 1:length(inc)
		cazi = azi(ipair) ;
		cinc = inc(ipair) ;

%  ** create the cartesian vector
		XI = cart2(1,cinc,cazi) ;

%  ** compute phase velocities		
		[V,EIGVEC]=velo(XI,rh,C) ;
		
%  ** pull out the eigenvectors
		P  = EIGVEC(:,1) ;
      S1 = EIGVEC(:,2) ;

		if ~isreal(S1)
			S1
			fprintf('%f,%f\n',cinc,cazi)
			C_in
			error('bad') ;
		end
      S2 = EIGVEC(:,3) ;

%  ** calculate projection onto propagation plane      
      S1N = cross(XI,S1) ;
      S1P = cross(XI,S1N);

%  ** rotate into y-z plane to calculate angles
      [S1PR] = V_rot3(S1P,0,0,cazi) ;
		[S1PRR] = V_rot3(S1PR,0,cinc,0) ;

		ph = atan2(S1PRR(2),S1PRR(3)) .* 180/pi ;

%  ** transform angle to between -90 and 90
      if (ph < -90.), ph = ph + 180.;,end
      if (ph >  90.), ph = ph - 180.;,end

%	** calculate some useful values
		dVS =  (V(2)-V(3)) ;
      VSmean = (V(2)+V(3))/2.0 ;

      avs(ipair) = 100.0*(dVS/VSmean) ;
      vp(ipair) =  V(1) ;
      vs1(ipair) = V(2) ;
      vs2(ipair) = V(3) ;
		
		SF(ipair,:) = S1P(:) ;
	end % ipair = 1:length(inc_in)
		
return
%=======================================================================================  
