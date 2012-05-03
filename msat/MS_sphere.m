% MS_SPHERE - Plot a spherical plot of phasevels/anisotropy from a set of 
%              elastic constants.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Given an elasticity matrix and density, produce spherical figures showing 
%     the P- or S-wave anisotrpy. 
%
%  % MS_sphere(CC,rh,mode,...)
%
% Usage: 
%     MS_sphere(CC,rh,mode)                    
%         Mode can be 'p'(-wave), 's'(-wave), 's1' (fast s-wave) or 's2'
%         (slow s-wave).
%
%     MS_sphere(CC,rh,mode,...) 
%          Further arguments are optional, can be combined in any order and
%          change the behaviour of the function. Specific options are
%          listed below.
%
%     MS_sphere(..., 'reverse')
%          Reverse the sense of the colour map.
%
%     MS_sphere(..., 'FSWTickLength', length)
%          Set the length of the Fast shear wave direction markers, default
%          is 0.08.
%
%     MS_sphere(..., 'FSWMarkerSize', length)
%          Set the size of the Fast shear wave direction markers, default
%          is 4.
%
%     MS_sphere(..., 'cmap', CM)                    
%          Redefine the colormap. CM can either be a (nx3) matrix containing a 
%          colormap, or a string describing a function to generate such a 
%          matrix (such as the built-in MATLAB colormap functions). E.g.:
%          MS_sphere(C,rh,mode,'cmap','cool(64)') - uses the MATLAB 
%          function cool to generate a cyan-to-purple colourmap.
%          The default is 'jet', reversed so blue is fast/high, as is 
%          conventional for seismic velocity colorscales. 
%
%     MS_sphere(..., 'velmesh', levels)
%          Choose the number of points used to mesh the velocity contours
%          expressed in terms of the number of times the mesh spacing is 
%          halved. Default is 3.
%
%     MS_sphere(..., 'polmesh', levels)
%          Choose the number of points used to mesh the fast shear wave 
%          markers. Expressed in terms of the number of times the mesh spacing is 
%          halved. Default is 2.
%
%     MS_sphere(..., 'dirs', array)
%          Directions of interest, this should be a vector of azi, 
%          inc pairs, e.g., [az1 in1 az2 in2 az3 in3 az4 in4].
%
%     MS_sphere(..., 'dlen', lenght)
%          Set length of doi vectors. Default is 1.5
%
%     MS_sphere(..., 'nolabels')
%          Supress printing of the text labels.
%
%     MS_sphere(..., 'noaxes')
%          Supress printing of the direction axes.
%
%     MS_sphere(..., 'nocbar')
%          Supress printing of the colour bar.
%
%     MS_sphere(..., 'nofig')
%          Suppress the figure command, so plot can be sub-plotted.
%
% See also: MS_SPHERE, MS_PHASEVELS


% Copyright (c) 2007 - 2011, James Wookey
% Copyright (c) 2011, James Wookey and Andrew Walker
% Portions of this function are:
% Copyright (c) 1989 - 2002, Darren Weber, Jon Leech and Jim Buddenhagen.
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

% NB: redistributors must also comply with the GPL. See below.


function MS_sphere(CC,rh,mode,varargin)

cmap = jet(64) ;
icmapflip = 1 ; % reverse the sense of the colourscale
FSWTickLength=0.08 ;
FSWMarkerSize=4;
dirs = [] ;
dlen = 1.5 ;
plotlabels = 1 ;
plotaxes = 1 ;
nofig = 0 ;
nocbar = 0 ;
velmesh = 3;
polmesh = 2;
cax = NaN ; % define to use

% Changing defaults for case of slowness surfaces before the optional arguments are processed
if ischar(mode)
   if sum(strcmpi({'slowp','slows1','slows2'},mode))==1
	   cmap=copper ;
	   icmapflip = 1 ; % reverse the sense of the colourscale
	   velmesh = 5 ; % increase this to make it better
   end
end  

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'reverse'
               icmapflip = 1 ;
               iarg = iarg + 1 ;
            case 'fswticklength'
               FSWTickLength = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'fswmarkersize'
               FSWMarkerSize = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'dirs'
               dirs = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'dlen'
               dlen = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'velmesh'
               velmesh = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'polmesh'
               polmesh = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'nolabels'
               plotlabels = 0 ;
               iarg = iarg + 1 ;
            case 'noaxes'
               plotaxes = 0 ;
               iarg = iarg + 1 ;
            case 'nofig'
               nofig = 0 ;
               iarg = iarg + 1 ;
            case 'nocbar'
               nocbar = 0 ;
               iarg = iarg + 1 ;
            case 'cmap'
               cmarg = varargin{iarg+1} ;
               if isstr(cmarg)
                  eval(['cmap = ' cmarg ';']) ;
               else
                  cmap = cmarg ;
               end
               iarg = iarg + 2 ;
            otherwise 
               error(['Unknown option: ' varargin{iarg}]) ;   
         end   
      end 

% check the inputs: CC
assert(MS_checkC(CC)==1, 'MS:SPHERE:badC', 'MS_checkC error MS_sphere') ;

% check the inputs: rh
if isnumeric(rh)
   if sum(size(rh)==1) ~= 2
      error('rh must be a scalar')
   end
else
   error('rh must be a scalar')
end

if ischar(mode)
   if sum(strcmpi({'p','s', 's1', 's2', 'slowp', 'slows1', 'slows2'},mode))~=1
      error('MS:SPHERE:badmode', ...
          'Mode must be ''S''(-wave), ''S1'', ''S2'', or ''P''(-wave) ')
   end
else
   error('MS:SPHERE:badmode', ...
       'Mode must be ''S''(-wave) or ''P''(-wave) ')
end      

% check input pars
if mod(length(dirs),2)~=0
   error('MS:SPHERE:baddirs', 'Plot directions must be in (azi,inc) pairs') ;
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


% load the triangulation
[x, y, z, faces, az, inc] =  get_mesh(velmesh);

[~,avs,vs1,vs2,vp] = MS_phasevels(CC,rh,inc,az) ;

if strcmpi(mode,'p')
   trisurf(faces,x,y,z,vp) ;
elseif strcmpi(mode,'s1')
   trisurf(faces,x,y,z,vs1) ;
elseif strcmpi(mode,'s2')
   trisurf(faces,x,y,z,vs2) ;
elseif strcmpi(mode, 'slowp')
    slowness = 1./vp;
    nslow = slowness./max(slowness);
    xslow = x .* nslow;
    yslow = y .* nslow;
    zslow = z .* nslow;
    trisurf(faces, xslow, yslow, zslow, slowness);
elseif strcmpi(mode, 'slows1')
    slowness = 1./vs1;
    nslow = slowness./max(slowness);
    xslow = x .* nslow;
    yslow = y .* nslow;
    zslow = z .* nslow;
    trisurf(faces, xslow, yslow, zslow, slowness);
elseif strcmpi(mode, 'slows2')
    slowness = 1./vs2;
    nslow = slowness./max(slowness);
    xslow = x .* nslow;
    yslow = y .* nslow;
    zslow = z .* nslow;
    trisurf(faces, xslow, yslow, zslow, slowness);
else
   trisurf(faces,x,y,z,avs) ;
end

shading interp
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]) ;
colormap(cmap)
daspect([1 1 1]) ;
hold on

if ~isnan(cax), caxis(cax), end

% plot the shear-wave polarisation
if sum(strcmpi({'s','s1','s2','slows1','slows2'},mode))==1 
   [x, y, z, ~, az, inc] =  get_mesh(polmesh);

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
   
   [~,avs,vs1,vs2,vp, SF, SS] = MS_phasevels(CC,rh,inc,az) ;
   % calculate PM vectors
   nv = length(x) ;
   if sum(strcmpi(mode,{'s','s1','slows1'}))==1 
      for iv=1:nv
		  if strcmpi(mode, 'slows1')
			  slow = 1./vs1(iv);
			  nslow = slow./max(slowness);
		      xslow = x(iv) .* nslow;
		      yslow = y(iv) .* nslow;
		      zslow = z(iv) .* nslow;
			  XI=1.01.*[xslow yslow zslow] ;			  
		  else
			  XI=1.01.*[x(iv) y(iv) z(iv)] ;
		  end          
          X1 = XI-FSWTickLength.*SF(iv,:);
          X2 = XI+FSWTickLength.*SF(iv,:);
          plot3(XI(1),XI(2),XI(3),'ko','MarkerSize',FSWMarkerSize,'MarkerFaceColor','k')
          plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)],'k-','LineWidth',(FSWMarkerSize/2))
      end
   else
      for iv=1:nv
		  if strcmpi(mode, 'slows2')
			  slow = 1./vs2(iv);
			  nslow = slow./max(slowness);
		      xslow = x(iv) .* nslow;
		      yslow = y(iv) .* nslow;
		      zslow = z(iv) .* nslow;
			  XI=1.01.*[xslow yslow zslow] ;			  
		  else
			  XI=1.01.*[x(iv) y(iv) z(iv)] ;
		  end 
          X1 = XI-FSWTickLength.*SS(iv,:);
          X2 = XI+FSWTickLength.*SS(iv,:);
          plot3(XI(1),XI(2),XI(3),'wo','MarkerSize',FSWMarkerSize,'MarkerFaceColor','w')
          plot3([X1(1) X2(1)],[X1(2) X2(2)],[X1(3) X2(3)],'w-','LineWidth',(FSWMarkerSize/2))
      end
   end
end

% draw the axes
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

% final commands

axis off
daspect([1 1 1]) ;
if (nocbar~=1)
cbax = colorbar('NorthOutside','FontSize',14,'FontWeight','bold');
if strcmpi(mode,'s')
    title(cbax,'S-wave anisotropy (%)');
elseif strcmpi(mode,'p')
    title(cbax,'P-wave velocity (km/s)');
elseif strcmpi(mode,'s1')
    title(cbax,'Fast S-wave velocity (km/s)');
elseif strcmpi(mode,'s2')
    title(cbax,'Slow S-wave velocity (km/s)');
elseif strcmpi(mode,'slowp')
    title(cbax,'P-wave slowness (s/km)');
	light('Position',[-0.58674 -0.05336 0.80801],'Style','infinite');
elseif strcmpi(mode,'slows1')
    title(cbax,'Fast S-wave slowness (s/km)');
	light('Position',[-0.58674 -0.05336 0.80801],'Style','infinite');
elseif strcmpi(mode,'slows2')
    title(cbax,'Slow S-wave slowness (s/km)');
	light('Position',[-0.58674 -0.05336 0.80801],'Style','infinite');


end
end
return


% Tooling to build the triangulation... 

function [x, y, z, tri, az, inc ] = get_mesh(level)

    FV = sphere_tri('ico', level, 1, 0);
    x = FV.vertices(:,1);
    y = FV.vertices(:,2);
    z = FV.vertices(:,3);
    tri = FV.faces;
    
%   cartesian to spherical conversion.
    r = sqrt(x.^2 + y.^2 + z.^2);
    phi = atan2(y,x);
    theta = real(acos(z./r));

    az = -phi*(180.0/pi);
    inc = -(theta*(180.0/pi)-90.0);
    
return


function [FV] = sphere_tri(shape,maxlevel,r,winding)

% sphere_tri - generate a triangle mesh approximating a sphere
% 
% Usage: FV = sphere_tri(shape,Nrecurse,r,winding)
% 
%   shape is a string, either of the following:
%   'ico'   starts with icosahedron (most even, default)
%   'oct'   starts with octahedron
%   'tetra' starts with tetrahedron (least even)
%
%   Nrecurse is int >= 0, setting the recursions (default 0)
%
%   r is the radius of the sphere (default 1)
%
%   winding is 0 for clockwise, 1 for counterclockwise (default 0).  The
%   matlab patch command gives outward surface normals for clockwise
%   order of vertices in the faces (viewed from outside the surface).
%
%   FV has fields FV.vertices and FV.faces.  The vertices 
%   are listed in clockwise order in FV.faces, as viewed 
%   from the outside in a RHS coordinate system.
% 
% The function uses recursive subdivision.  The first
% approximation is an platonic solid, either an  icosahedron,
% octahedron or a tetrahedron.  Each level of refinement 
% subdivides each triangle face by a factor of 4 (see also 
% mesh_refine).  At each refinement, the vertices are 
% projected to the sphere surface (see sphere_project).
% 
% A recursion level of 3 or 4 is a good sphere surface, if
% gouraud shading is used for rendering.
% 
% The returned struct can be used in the patch command, eg:
% 
% % create and plot, vertices: [2562x3] and faces: [5120x3]
% FV = sphere_tri('ico',4,1);
% lighting phong; shading interp; figure;
% patch('vertices',FV.vertices,'faces',FV.faces,...
%       'facecolor',[1 0 0],'edgecolor',[.2 .2 .6]);
% axis off; camlight infinite; camproj('perspective');
% 
% See also: mesh_refine, sphere_project
%



% $Revision: 1.15 $ $Date: 2004/05/20 22:28:45 $

% Licence:  GNU GPL, no implied or express warranties
% Jon Leech (leech @ cs.unc.edu) 3/24/89
% icosahedral code added by Jim Buddenhagen (jb1556@daditz.sbc.com) 5/93
% 06/2002, adapted from c to matlab by Darren.Weber_at_radiology.ucsf.edu
% 05/2004, reorder of the faces for the 'ico' surface so they are indeed
% clockwise!  Now the surface normals are directed outward.  Also reset the
% default recursions to zero, so we can get out just the platonic solids.
%
%


if ~exist('shape','var') || isempty(shape),
    shape = 'ico';
end
% - AMW - fprintf('...creating sphere tesselation based on %s\n',shape);

% default maximum subdivision level
if ~exist('maxlevel','var') || isempty(maxlevel) || maxlevel < 0,
    maxlevel = 0;
end

% default radius
if ~exist('r','var') || isempty(r),
    r = 1;
end

if ~exist('winding','var') || isempty(winding),
    winding = 0;
end

% -----------------
% define the starting shapes

shape = lower(shape);

switch shape,
case 'tetra',
    
    % Vertices of a tetrahedron
    sqrt_3 = 0.5773502692;
    
    tetra.v = [  sqrt_3,  sqrt_3,  sqrt_3 ;   % +X, +Y, +Z  - PPP
                -sqrt_3, -sqrt_3,  sqrt_3 ;   % -X, -Y, +Z  - MMP
                -sqrt_3,  sqrt_3, -sqrt_3 ;   % -X, +Y, -Z  - MPM
                 sqrt_3, -sqrt_3, -sqrt_3 ];  % +X, -Y, -Z  - PMM
	
    % Structure describing a tetrahedron
    tetra.f = [ 1, 2, 3;
                1, 4, 2;
                3, 2, 4;
                4, 1, 3 ];
    
    FV.vertices = tetra.v;
    FV.faces    = tetra.f;
    
case 'oct',
    
    % Six equidistant points lying on the unit sphere
    oct.v = [  1,  0,  0 ;  %  X
              -1,  0,  0 ; 	% -X
               0,  1,  0 ;  %  Y
               0, -1,  0 ; 	% -Y
               0,  0,  1 ; 	%  Z
               0,  0, -1 ];	% -Z
	
    % Join vertices to create a unit octahedron
    oct.f = [ 1 5 3 ;    %  X  Z  Y  -  First the top half
              3 5 2 ;    %  Y  Z -X
              2 5 4 ;    % -X  Z -Y
              4 5 1 ;    % -Y  Z  X
              1 3 6 ;    %  X  Y -Z  -  Now the bottom half
              3 2 6 ;    %  Y  Z -Z
              2 4 6 ;    % -X  Z -Z
              4 1 6 ];   % -Y  Z -Z
    
    FV.vertices = oct.v;
    FV.faces    = oct.f;
    
case 'ico',
    
    % Twelve vertices of icosahedron on unit sphere
    tau = 0.8506508084; % t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
    one = 0.5257311121; % one=1/sqrt(1+t^2) , unit sphere
    
    ico.v( 1,:) = [  tau,  one,    0 ]; % ZA
    ico.v( 2,:) = [ -tau,  one,    0 ]; % ZB
    ico.v( 3,:) = [ -tau, -one,    0 ]; % ZC
    ico.v( 4,:) = [  tau, -one,    0 ]; % ZD
    ico.v( 5,:) = [  one,   0 ,  tau ]; % YA
    ico.v( 6,:) = [  one,   0 , -tau ]; % YB
    ico.v( 7,:) = [ -one,   0 , -tau ]; % YC
    ico.v( 8,:) = [ -one,   0 ,  tau ]; % YD
    ico.v( 9,:) = [   0 ,  tau,  one ]; % XA
    ico.v(10,:) = [   0 , -tau,  one ]; % XB
    ico.v(11,:) = [   0 , -tau, -one ]; % XC
    ico.v(12,:) = [   0 ,  tau, -one ]; % XD
    
    % Structure for unit icosahedron
    ico.f = [  5,  8,  9 ;
               5, 10,  8 ;
               6, 12,  7 ;
               6,  7, 11 ;
               1,  4,  5 ;
               1,  6,  4 ;
               3,  2,  8 ;
               3,  7,  2 ;
               9, 12,  1 ;
               9,  2, 12 ;
              10,  4, 11 ;
              10, 11,  3 ;
               9,  1,  5 ;
              12,  6,  1 ;
               5,  4, 10 ;
               6, 11,  4 ;
               8,  2,  9 ;
               7, 12,  2 ;
               8, 10,  3 ;
               7,  3, 11 ];
	
    FV.vertices = ico.v;
    FV.faces    = ico.f;
end


% -----------------
% refine the starting shapes with subdivisions
if maxlevel,
    
    % Subdivide each starting triangle (maxlevel) times
    for level = 1:maxlevel,
        
        % Subdivide each triangle and normalize the new points thus
        % generated to lie on the surface of a sphere radius r.
        FV = mesh_refine_tri4(FV);
        FV.vertices = sphere_project(FV.vertices,r);
        
        % An alternative might be to define a min distance
        % between vertices and recurse or use fminsearch
        
    end
end

if winding,
    % - AMW - fprintf('...returning counterclockwise vertex order (viewed from outside)\n');
    FV.faces = FV.faces(:,[1 3 2]);
% - AMW - else
    % - AMW - fprintf('...returning clockwise vertex order (viewed from outside)\n');
end

return

function [ FV ] = mesh_refine_tri4(FV)

% mesh_refine_tri4 - creates 4 triangle from each triangle of a mesh
%
% [ FV ] = mesh_refine_tri4( FV )
%
% FV.vertices   - mesh vertices (Nx3 matrix)
% FV.faces      - faces with indices into 3 rows
%                 of FV.vertices (Mx3 matrix)
% 
% For each face, 3 new vertices are created at the 
% triangle edge midpoints.  Each face is divided into 4
% faces and returned in FV.
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
% 
% It is assumed that the vertices are listed in clockwise order in
% FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
% system.
% 
% See also: mesh_refine, sphere_tri, sphere_project
% 


% ---this method is not implemented, but the idea here remains...
% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
% Alternatively, it could be done until some minimum mean 
% area of faces is achieved.  As is, it just refines once.


% $Revision: 1.12 $ $Date: 2004/05/10 21:01:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - AMW - tic;
% - AMW - fprintf('...refining mesh (tri4)...')

% NOTE
% The centroid is located one third of the way from each vertex to 
% the midpoint of the opposite side. Each median divides the triangle 
% into two equal areas; all the medians together divide it into six 
% equal parts, and the lines from the median point to the vertices 
% divide the whole into three equivalent triangles.

% Each input triangle with vertices labelled [A,B,C] as shown
% below will be turned into four new triangles:
%
% Make new midpoints
% a = (A+B)/2
% b = (B+C)/2
% c = (C+A)/2
%
%        B
%       /\
%      /  \
%    a/____\b       Construct new triangles
%    /\    /\       [A,a,c]
%   /  \  /  \      [a,B,b]
%  /____\/____\     [c,b,C]
% A	     c	   C    [a,b,c]
%

% Initialise a new vertices and faces matrix
Nvert = size(FV.vertices,1);
Nface = size(FV.faces,1);
V2 = zeros(Nface*3,3);
F2 = zeros(Nface*4,3);

for f = 1:Nface,
    
    % Get the triangle vertex indices
    NA = FV.faces(f,1);
    NB = FV.faces(f,2);
    NC = FV.faces(f,3);
    
    % Get the triangle vertex coordinates
    A = FV.vertices(NA,:);
    B = FV.vertices(NB,:);
    C = FV.vertices(NC,:);
    
    % Now find the midpoints between vertices
    a = (A + B) ./ 2;
    b = (B + C) ./ 2;
    c = (C + A) ./ 2;
    
    % Find the length of each median
    %A2blen = sqrt ( sum( (A - b).^2, 2 ) );
    %B2clen = sqrt ( sum( (B - c).^2, 2 ) );
    %C2alen = sqrt ( sum( (C - a).^2, 2 ) );
    
    % Store the midpoint vertices, while
    % checking if midpoint vertex already exists
    [FV, Na] = mesh_find_vertex(FV,a);
    [FV, Nb] = mesh_find_vertex(FV,b);
    [FV, Nc] = mesh_find_vertex(FV,c);
    
    % Create new faces with orig vertices plus midpoints
    F2(f*4-3,:) = [ NA, Na, Nc ];
    F2(f*4-2,:) = [ Na, NB, Nb ];
    F2(f*4-1,:) = [ Nc, Nb, NC ];
    F2(f*4-0,:) = [ Na, Nb, Nc ];
    
end

% Replace the faces matrix
FV.faces = F2;

% - AMW - t=toc; fprintf('done (%5.2f sec)\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, N] = mesh_find_vertex(FV,vertex)

    Vn = size(FV.vertices,1);
    Va = repmat(vertex,Vn,1);
    Vexist = find( FV.vertices(:,1) == Va(:,1) & ...
                   FV.vertices(:,2) == Va(:,2) & ...
                   FV.vertices(:,3) == Va(:,3) );
    if Vexist,
        if size(Vexist) == [1,1],
            N = Vexist;
        else,
            msg = sprintf('replicated vertices');
            error(msg);
        end
    else
        FV.vertices(end+1,:) = vertex;
        N = size(FV.vertices,1);
    end

return

function V = sphere_project(v,r,c)

% sphere_project - project point X,Y,Z to the surface of sphere radius r
% 
% V = sphere_project(v,r,c)
% 
% Cartesian inputs:
% v is the vertex matrix, Nx3 (XYZ)
% r is the sphere radius, 1x1 (default 1)
% c is the sphere centroid, 1x3 (default 0,0,0)
%
% XYZ are converted to spherical coordinates and their radius is
% adjusted according to r, from c toward XYZ (defined with theta,phi)
% 
% V is returned as Cartesian 3D coordinates
% 

% $Revision: 1.8 $ $Date: 2004/03/29 21:15:36 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('v','var'),
    msg = sprintf('SPHERE_PROJECT: No input vertices (X,Y,Z)\n');
    error(msg);
end

X = v(:,1);
Y = v(:,2);
Z = v(:,3);

if ~exist('c','var'),
    xo = 0;
    yo = 0;
    zo = 0;
else
    xo = c(1);
    yo = c(2);
    zo = c(3);
end

if ~exist('r','var'), r = 1; end

% alternate method is to use unit vector of V
% [ n = 'magnitude(V)'; unitV = V ./ n; ]
% to change the radius, multiply the unitV
% by the radius required.  This avoids the
% use of arctan functions, which have branches.


% Convert Cartesian X,Y,Z to spherical (radians)
theta = atan2( (Y-yo), (X-xo) );
phi   = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
% do not calc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);

%   Recalculate X,Y,Z for constant r, given theta & phi.
R = ones(size(phi)) * r;
x = R .* sin(phi) .* cos(theta);
y = R .* sin(phi) .* sin(theta);
z = R .* cos(phi);

V = [x y z];

return



