% MS_sphere : Plot a spherical plot of phasevels/anisotropy from a set of 
%              elastic constants.
%   
% Usage: MS_sphere(CC,rh,mode,...)
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
% (C) James Wookey, 2007-2011.	

%-------------------------------------------------------------------------------
%  This software is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%-------------------------------------------------------------------------------

function MS_sphere(CC,rh,mode,varargin)

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
   if sum(strcmpi({'p','s', 's1', 's2'},mode))~=1
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
[~,avs,vs1,vs2,vp] = MS_phasevels(CC,rh,inc,az) ;

if strcmpi(mode,'p')
   trisurf(faces,x,y,z,vp) ;
elseif strcmpi(mode,'s1')
   trisurf(faces,x,y,z,vs1) ;
elseif strcmpi(mode,'s2')
   trisurf(faces,x,y,z,vs2) ;
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
if sum(strcmpi({'s','s1','s2'},mode))==1 
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
   
   [~,~,~,~,~, SF, SS] = MS_phasevels(CC,rh,inc,az) ;
   %% calculate PM vectors
   nv = length(x) ;
   if sum(strcmpi(mode,{'s','s1'}))==1 
      for iv=1:nv
          XI=1.01.*[x(iv) y(iv) z(iv)] ;
          X1 = XI-FSWTickLength.*SF(iv,:);
          X2 = XI+FSWTickLength.*SF(iv,:);
          plot3(XI(1),XI(2),XI(3),'ko','MarkerSize',4,'MarkerFaceColor','k')
          plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)],'k-','LineWidth',2)
      end
   else
      for iv=1:nv
          XI=1.01.*[x(iv) y(iv) z(iv)] ;
          X1 = XI-FSWTickLength.*SS(iv,:);
          X2 = XI+FSWTickLength.*SS(iv,:);
          plot3(XI(1),XI(2),XI(3),'wo','MarkerSize',4,'MarkerFaceColor','w')
          plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)],'w-','LineWidth',2)
      end
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
cbax = colorbar('NorthOutside','FontSize',14,'FontWeight','bold');
if strcmpi(mode,'s')
    title(cbax,'S-wave anisotropy (%)');
elseif strcmpi(mode,'p')
    title(cbax,'P-wave velocity (m/s)');
elseif strcmpi(mode,'s1')
    title(cbax,'Fast S-wave velocity (m/s)');
elseif strcmpi(mode,'s2')
    title(cbax,'Slow S-wave velocity (m/s)');
end
end
return
