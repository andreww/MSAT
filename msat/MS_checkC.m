%-------------------------------------------------------------------------------
%                  MSAT - Matlab Seismic Anisotropy Toolkit 
%-------------------------------------------------------------------------------
% MS_checkC - check consistency of a stiffness matrix, for use in MSAT codes
%-------------------------------------------------------------------------------
% [] = MS_checkC(C,...)
%
%	Inputs:
%     
%     C = Stiffness tensor for checking
%
%  Options:
%     Various checks can be disabled/enabled by specifying extra options:
%        'fast' - run a minimal set of checks (size and type only)
%        'warn' - warn only on more sophisticated errors than size/type
%        'nosymchk' - disable upper, lower symmetry checking
%        'nozerochk' - disable top left / trace zero checking 
% 
%     Various control parameters can be set by adding parameter,value pairs 
%     to the arguments
%        'thresh',x - set numerical threshold for zero to x (default = 1e-6)
%
%  See source code for further notes

% (C) James Wookey, 2011
function [isgood] = MS_checkC(C,varargin)


%  ** set defaults
      symchk = 1 ;
      zerochk = 1;
      thresh = 1e-6 ;
      warn = 0 ;
            
%  ** process the optional arguments
      iarg = 1 ;
      while iarg<=(length(varargin))
         switch varargin{iarg}
            case 'warn'
               warn = 1 ;
               iarg = iarg + 1 ;
            case 'nosymchk'
               symchk = 0 ;
               iarg = iarg + 1 ;
            case 'nozerochk'
               zerochk = 0 ;
               iarg = iarg + 1 ;
            case 'thresh'
               check_symmetry = varargin{iarg+1}
               iarg = iarg + 2 ;
            otherwise 
               error('Unknown flag') ;   
         end   
      end 

      isgood = 1 ;

%  ** start with basic checks
		if ~ismatrix(C)
		   isgood = 0 ;
		   error('Stiffness matrix error: Appears not to be a 2D matrix')
	   end
      [nr nc] = size(C) ;
      if (nr~=6 | nc~=6)
         isgood = 0;
		   error('Stiffness matrix error: Appears not to be a 2D matrix')
		end

%  ** check symmetry    
      if symchk
         D=C-C';
         if (length(find(abs(D)>thresh))>0)
            isgood = 0;
            if warn
               warning('Stiffness matrix warning: non-symmetric')
            else
               error('Stiffness matrix error: non-symmetric')
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
               warning('Stiffness matrix warning: zeros in trace or top left')
            else
               error('Stiffness matrix error: zeros in trace or top left')
            end   
         end
      end

%  ** check sensible magnitudes.	Make sure C(3,3). 
            
		
return
%=======================================================================================  
