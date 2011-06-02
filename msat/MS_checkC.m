% MS_CHECKC - check consistency of a stiffness matrix against various criteria
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [ isgood ] = MS_checkC( C,... )
%
%	Inputs:
%     C = Stiffness tensor for checking
%
%  Outputs:
%     isgood = 1 for OK
%
%  Options:
%     Various checks can be disabled/enabled by specifying extra options:
%        'fast' - run a minimal set of checks (size and type only)
%        'warn' - warn only on more sophisticated errors than size/type
%        'nosymchk' - disable upper, lower symmetry checking
%        'nozerochk' - disable top left / trace zero checking 
%        'nopdefchk' - disable positive definiteness check 
% 
%     Various control parameters can be set by adding parameter,value pairs 
%     to the arguments
%        'thresh',x - set numerical threshold for zero to x (default = 1e-6)
%
%  See source code for further notes

% (C) James Wookey and Andrew Walker, 2011
function [isgood] = MS_checkC(C,varargin)

%  ** set defaults
      warn = 0 ;
      fast = 0 ;

      thresh = 1e-6 ;

      symchk = 1 ;
      zerochk = 1 ;
      pdefchk = 1 ;
      
%  ** process the optional arguments
      iarg = 1 ;
      while iarg<=(length(varargin))
         switch varargin{iarg}
            case 'fast'
               fast = 1 ;
               iarg = iarg + 1 ;
            case 'warn'
               warn = 1 ;
               iarg = iarg + 1 ;
            case 'nosymchk'
               symchk = 0 ;
               iarg = iarg + 1 ;
            case 'nozerochk'
               zerochk = 0 ;
               iarg = iarg + 1 ;
            case 'nopdef'
               magchk = 0 ;
               iarg = iarg + 1 ;
            case 'thresh'
               thresh = varargin{iarg+1}
               iarg = iarg + 2 ;
            otherwise 
               error('Unknown flag') ;   
         end   
      end 

      isgood = 1 ;

%  ** start with basic checks
      try % ismatrix is quite new, so suppress Undef func error.
          if ~ismatrix(C)
		    isgood = 0 ;
	       error('MS:CHECKCNotMatrix', ...
	         'Elasticity matrix error: Appears not to be a 2D matrix')
          end
      catch e
          if ~strcmp(e.identifier, 'MATLAB:UndefinedFunction')
              rethrow(e)
          end
      end
      [nr nc] = size(C) ;
      if (nr~=6 | nc~=6)
         isgood = 0;
		   error('MS:CHECKCnot6x6',...
		      'Elasticity matrix error: Appears not to be a 6x6 matrix')
		end

%  ** if fast checking is selected, we're done here.       
      if fast, return, end

%  ** otherwise, perform more elaborate tests.
      
%  ** check symmetry    
      if symchk
         D=C-C';
         if (length(find(abs(D)>thresh))>0)
            isgood = 0;
            if warn
               warning('MS:CHECKCnotsym_W',...
                  'Elasticity matrix warning: non-symmetric')
            else
               error('MS:CHECKCnotsym',...
                  'Elasticity matrix error: non-symmetric')
            end   
         end   
      end

%  ** check that enough values are set non-zero
      if zerochk
         isopart = [C(1,1) C(1,2) C(1,3) C(2,2) C(2,3) ...
                    C(3,3) C(4,4) C(5,5) C(6,6)] ;
         if (length(find(abs(isopart)<thresh))>0)
            isgood = 0;
            if warn
               warning('MS:CHECKCbadzeros_W',...
                  'Elasticity matrix warning: zeros in trace or top left')
            else
               error('MS:CHECKCbadzeros',...
                  'Elasticity matrix error: zeros in trace or top left')
            end   
         end
      end

%  ** check that matrix is positive definite
      if pdefchk
         try
            X=chol(C) ;
         catch ME   
            isgood = 0;
            if warn
               warning('MS:CHECKCnotposdef_W',...
                  'Elasticity matrix warning: not positive definite')
            else
               error('MS:CHECKCnotposdef',...
                  'Elasticity matrix error: not positive definite')
            end
         end               
      end
return
%=======================================================================================  
