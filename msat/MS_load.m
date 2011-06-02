% MS_load - Load a set of elastic constants from a file.
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
%     [C,...] = MS_load(fname,...,'eunit',unit)
%         Specify elasticity unit; this is used to convert to GPa (msat's 
%         default operating unit)
%            'GPa' - elastic constants are already in GPa (default)
%            'Mbar' - elastic constants are in Mbar
%            'Pa' - elastic constants are in Pascals
%            'bar' - elastic constants are in bar
%
%     [C,...] = MS_load(fname,...,'dunit',unit)
%         Specify density unit; this is used to convert to kg/m3 (msat's 
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
%    Elastic constants in the format read by EMATRIX (D. Mainprice). 
%

%===============================================================================
   function [varargout] = MS_load(fname,varargin) ;
%===============================================================================
      eunit = 'default' ;
      dunit = 'default' ;
      
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
            case 'eunit'
               eunit = varargin{iarg+1} ;
               iarg = iarg + 2 ;
               eunit_set = 1 ;
            case 'dunit'
               dunit = varargin{iarg+1} ;
               iarg = iarg + 2 ;
               dunit_set = 1 ;
            case 'symmetry'
               smode = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            otherwise 
               error(['Unknown option: ' varargin{iarg}]) ;   
         end   
      end 

%      
      switch lower(format)
      case 'simple'
         [ C, nec, rh ] = MS_load_simple( fname ) ;
         if nec<9 & strcmp(smode,'none')
            error('MS:LOAD:ExpansionRequired',...
               ['Fewer than 9 elastic constants were specified, and'  ...
                ' no symmetry based expansion' ' was requested']) ;
         end
      case 'ematrix'
%     ** forbid expansion from ematrix files. 
         if ~strcmp(smode,'none'), error('MS:LOAD:ExpandForbidden',...
            'Symmetry expansion is not supported from this file format.') ;,end
         [C,rh,eunit,dunit] = MS_load_ematrix(fname,eunit,dunit) ;
      otherwise
         error('MS:LOAD:BadFileFormat',...
            ['Specified file format "' format '" is not supported.']) ;
      end   

%  ** set default units if this hasn't already been done.
      switch lower(eunit)
      case 'default'
         eunit = 'GPa' ;
      otherwise
      % do nothing
      end
      
      switch lower(dunit)
      case 'default'
         dunit = 'kgm3' ;
      otherwise
      % do nothing
      end

%  ** check if density is required but not provided.
      if isnan(rh) & ( aij | nargout > 1 )
            error('MS:LOADdensityneeded',...
               'No density is specified, but one is required.')
      end

%  ** symmetry handling
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

%  ** unnormalise for density (if required)
      if aij, C = C .* rh;, end

%  ** perform unit conversions
      switch lower(eunit)
      case 'gpa'
      case 'pa'
         C = C ./ 1e9 ;
      case 'mbar'
         C = C .* 100 ;
      case 'bar'
         C = C ./ 10e3 ;
      otherwise
         error('MS:LOAD:BadEunit','Unsupported elasticity unit.')
      end
      
      switch lower(dunit)
      case 'kgm3'
      case 'gcc'
         rh = rh .* 1e3 ;
      otherwise
         error('MS:LOAD:BadDunit','Unsupported density unit.')
      end


%  ** finally, check Cij matrix
      if (~force)
         try
            MS_checkC(C) ;
         catch ME
            error('MS:LOAD:CheckFailed', ...
               ['Final matrix check failed with error message: ' ME.message]) ;
         end   
      end   

%  ** set the output arguments
      switch nargout
      case 0
         varargout(1) = {C} ;
      case 1
         varargout(1) = {C} ;
      case 2
         varargout(1) = {C} ;
         varargout(2) = {rh} ;
      otherwise
         error('MS:LOAD:BadOutputArgs', ...
            'Requires 1 or 2 output arguments (Cij or [Cij,density]).')
      end


   return
%===============================================================================

%===============================================================================
%-------------------------------------------------------------------------------
%  File loader functions
%-------------------------------------------------------------------------------
%===============================================================================

%===============================================================================
   function [ C , nec, rh ] = MS_load_simple( fname )
%===============================================================================
%  'simple' format

%  ** load the ascii file in a standard MATLAB way
      try
         a = load(fname) ;
      catch ME
         error(['File load failed with error: ' ME.message]) ;
      end   

      [nline ncol] = size(a) ;
      
      if ncol~=3, error('MS:LOAD_SIMPLE:BadFormat',...
         'File appears not to be in the specified format') ;, end

%  ** remove any density lines
      ind = find(a(:,1)>=1 & a(:,1)<=6 & a(:,2)>=1 & a(:,2)<=6) ;

%  ** extract indices and constants
      ii = a(ind,1) ;
      jj = a(ind,2) ;
      ec = a(ind,3) ;

%  ** count 
      [nec ndum] = size(ec) ;

%  ** get density
      if (nline-nec)==1
%     ** one density specified
         ind = find(a(:,1)<1 | a(:,1)>6 & a(:,2)<1 | a(:,2)>6) ;
         rh = a(ind,3) ;
      elseif nline==nec
%     ** no density is specified
         rh = NaN ;
      else   
         error('MS:LOAD_SIMPLE:MultipleDensities',...
         'Multiple densities specified (only one permitted).')
      end
      
      C = zeros(6,6) ;
      for i=1:nec
         C(ii(i),jj(i)) = ec(i) ;
      end
      
   return
%===============================================================================

%===============================================================================
%  'ematrix' format
%===============================================================================
   function [C, rh, eunit_out, dunit_out] = MS_load_ematrix(fname, eunit, dunit)
%===============================================================================
%  ** default units for ematrix files is Mbar & g/cc. Redefine defaults if user
%     has not explicitly defined units in call. 
      switch lower(eunit)
      case 'default'
         eunit_out = 'Mbar' ;
      otherwise
         eunit_out = eunit ;
      end

      switch lower(dunit)
      case 'default'
         dunit_out = 'gcc' ;
      otherwise
         dunit_out = dunit ;
      end

%  ** no density is specified in these files. 
      rh = NaN ;

%  ** read the file
      fid = fopen(fname,'rt') ;
%  ** header lines (these are ignored)
      TS = textscan(fid,'%*[^\n]',2) ;

%  ** rotations line, this is not used at the moment.
      TS = textscan(fid,'%f %f %f %f %f %f',1, 'CollectOutput', 1) ;
      ematrix_rotations = TS{1} ;

%  ** the good stuff
      TS = textscan(fid,'%f %f %f %f %f %f',6, 'CollectOutput', 1) ;
      C = TS{1} ;
      
      fclose(fid) ;
%===============================================================================
   return
%===============================================================================
