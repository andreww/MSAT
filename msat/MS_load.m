%-------------------------------------------------------------------------------
%                  MSAT - Matlab Seismic Anisotropy Toolkit 
%-------------------------------------------------------------------------------
% MS_load - Load a set of elastic constants from a file.
%-------------------------------------------------------------------------------
% 
% [ ... ] = MS_load( fname, ... )
%
% Usage: 
%     [C] = MS_load(fname)                    
%         Load the elastic constants in file 'fname'. See file format below. 
%
%     [C,r] = MS_load(fname,...)                    
%         Load the elastic constants and density in file 'fname'. See file 
%         format below.
%
%     [C,...] = MS_load(fname,...,'format',fmt)                    
%         Specify format for file. See below for descriptions. Available options 
%         are:
%            'simple' : Formatted text file (default)
%
%     [C,...] = MS_load(fname,...,'symmetry', mode) 
%         ** some parts not yet implemented **
%         Fill out elastic tensor based on symmetry, defined by mode. This can
%         take the following values:
%            'none' - nothing attempted, unspecified Cijs are zero (default)
%            'auto' - assume symmetry based on number of Cijs specified 
%            'iso' - isotropic (nCij=2) ; C33 and C66 must be specified.
%            'hex' - hexagonal (nCij=5) ; C11, C33, C44, C66 and C13 must be
%                       specified, x3 is symmetry axis
%            'vti' - synonym for hexagonal
%            'cubic' - cubic (nCij=3) ; C33, C66 and C12 must be specified
%
%     [C,...] = MS_load(fname,...,'eunits',unit)
%         Specify elasticity units; this is used to convert to GPa (msat's 
%         default operating unit)
%            'GPa' - elastic constants are already in GPa (default)
%            'mbar' - elastic constants are in mbar
%            'Pa' - elastic constants are in Pascals
%            'bar' - elastic constants are in bar
%
%     [C,...] = MS_load(fname,...,'dunits',unit)
%         Specify density units; this is used to convert to kg/m3 (msat's 
%         default operating unit)
%            'kgm3' - density is already in kg/m3 (default)
%            'gcc' - density is in g/cc
%
%     [C,...] = MS_load(fname,...,'Aij')
%         Elasticities are specified in density normalised form. This option
%         requires a density to be specified, the elasticities are returned with
%         the normalisation reversed (default format for msat). The density
%         normalisation is assumed to have been performed using the density 
%         specified in the file (so, obviously, one is required) without any
%         different unit scalings. It is *highly* advisable to check (for 
%         example using MS_phasevels) that sensible numbers result. The default
%         is to assume that the input is not density normalised. 
%
%     [C,...] = MS_load(fname,...,'force')
%         Disable post-load checking. You're on your own, buddy.
%
%  Supported file formats:
%  -----------------------
%
%  * 'simple' (default)
%        
%     In this mode MS_load expects an ascii file containing lines of the form:
%
%     I J C(I,J)
%
%     Elasticity tensor components are specified by their indices (i,j). Any 
%     number other than 1-6 in the i or j column denotes that a density is 
%     being given. 
%
%     So, for olivine, the file might look like:
%     
%     1 1 320.5
%     1 2 68.1
%     1 3 71.6
%     ...
%     6 6 78.7
%     0 0 3325
%     
%     Text following a '%' in any line is ignored (useful for header lines or 
%     comments)
%
%  * 'ematrix'
%
%    Elastic constants in the format read by EMATRIX (D. Mainprice). This 
%    consists of 2 header lines 
%
function [varargout] = MS_load(fname,varargin) ;

eunits = 'GPa' ;
dunits = 'kgm3' ;
aij = 0 ;
smode = 'none' ;
force = 0 ;
format = 'simple' ;

% check inputs here
% ...

%  ** process the optional arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
      case 'aij'
         aij = 1 ;
         iarg = iarg + 1 ;
      case 'force'
         force = 1 ;
         iarg = iarg + 1 ;
      case 'format'
         format = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'eunits'
         eunits = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'dunits'
         dunits = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      case 'symmetry'
         smode = varargin{iarg+1} ;
         iarg = iarg + 2 ;
      otherwise 
         error(['Unknown option: ' varargin{iarg}]) ;   
   end   
end 

switch lower(format)
case 'simple'
   [ ii, jj, ec, nec, rh ] = MS_load_simple( fname ) ;
otherwise
   error(['Specified file format "' format '" is not supported.']) ;
end   

% check if density is required but not provided.
if isnan(rh) & ( aij | nargout > 1 )
      error('No density is specified, but one is required.')
end

C = zeros(6,6) ;

for i=1:nec
   C(ii(i),jj(i)) = ec(i) ;
end

%  symmetry handling
switch lower(smode)
case 'none'
   % nothing required, just make symmetrical
   for i=1:6
      for j=(i+1):6
         C(j,i) = C(i,j) ;
      end
   end
otherwise
   % expand using MS_expand
   C = MS_expand(C,smode) ;
end

% unnormalise for density (if required)
if aij, C = C .* rh;, end

% perform unit conversions
switch lower(eunits)
case 'gpa'
case 'pa'
   C = C ./ 1e9 ;
case 'mbar'
   C = C .* 100 ;
case 'bar'
   C = C ./ 10e3 ;
otherwise
   error('Unsupported elasticity unit.')
end

switch lower(dunits)
case 'kgm3'
case 'gcc'
   rh = rh .* 1e3 ;
otherwise
   error('Unsupported density unit.')
end

% set the output arguments
switch nargout
case 1
   varargout(1) = {C} ;
case 2
   varargout(1) = {C} ;
   varargout(2) = {rh} ;
otherwise
   error('Requires 1 or 2 output arguments (Cij or [Cij,density]).')
end

% finally, check Cij matrix
if (~force)
   try
      MS_checkC(C) ;
   catch ME
      error(['Final matrix check failed with error message: ' ME.message]) ;
   end   
end   

return
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%  Loader functions
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%  'simple' format
%-------------------------------------------------------------------------------
function [ ii, jj, ec, nec, rh ] = MS_load_simple( fname )

% load the ascii file
try
   a = load(fname) ;
catch ME
   error(['File load failed with error: ' ME.message]) ;
end   

[nline ncol] = size(a) ;

if ncol~=3, error('File appears not to be in the specified format') ;, end

% remove any density lines
ind = find(a(:,1)>=1 & a(:,1)<=6 & a(:,2)>=1 & a(:,2)<=6) ;

% extract indices and constants
ii = a(ind,1) ;
jj = a(ind,2) ;
ec = a(ind,3) ;

% count 
[nec ndum] = size(ec) ;

% get density
if (nline-nec)==1
   % one density specified
   ind = find(a(:,1)<1 | a(:,1)>6 & a(:,2)<1 | a(:,2)>6) ;
   rh = a(ind,3) ;
elseif nline==nec
   % no density is specified
   rh = NaN ;
else   
   error('Multiple densities specified (only one permitted).')
end   
return
%-------------------------------------------------------------------------------
