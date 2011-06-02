% MS_INFO - Information about an elasticity matrix.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Prints out or returns a variety of information about an input elasticity
% matrix.
%
% [ ... ] = MS_info( C, ... )
%
% Usage: 
%     [] = MS_info(C)                    
%         Prints to the screen the default information about the matrix. 
%
%     [] = MS_info(C, 'mode', MODE)                    
%         Choose the set of outputs to return. Currently supported are:
%           'all' : all available information (default)
%
% Notes:
%
%

% Change log / copyright statements.
% 
%

%===============================================================================
   function [ ] = MS_info( C, varargin )
%===============================================================================

%  ** set defaults
      thresh = 0.01 ;
      imode = 'all' ;

      inf_vels = 0 ;
      inf_rot = 0 ;

%  ** check required_inputs here
      MS_checkC(C) ;

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'mode'  % parameter definition (value required)
               imode = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            otherwise 
               error('MS:INFO:UnknownOption',...
                  ['Unknown option: ' varargin{iarg}]) ;   
         end   
      end
      
      switch lower(imode)
         case 'all'
            inf_vels = 1 ;
            inf_rot = 1 ;
         otherwise
            error('MS:INFO:UnsupportedMode',...
               [ 'Unsupported mode: ' varargin{iarg}]) ;               
      end

%  ** rotational information
      if 
      fprintf(' * ROTATIONAL INFORMATION \n\n')
      try 
         [CR,RR] = MS_axes(C) ;
         ind = find(abs(RR) > thresh & abs(RR) < (1-thresh)) ;
         if ~isempty(ind)
            fprintf(['   Principle axes appear to be in a different', ...
                    ' orientation to this matrix.\n'])
          disp(' ');
          disp('    rotation matrix  =');
          disp(' ');
          disp(RR) ;       
         end           
      catch ME
         fprintf('   Principle axis check failed.\n\n')
      end
               
%  ** Velocity information
      

%===============================================================================
   return
%===============================================================================

%===============================================================================
%-------------------------------------------------------------------------------
%  SUBFUNCTIONS
%-------------------------------------------------------------------------------
%===============================================================================

%===============================================================================
   function [Y]=subfunction(X)
%===============================================================================
%  ** simple subfunction
      Y=X ;
%===============================================================================
   return
%===============================================================================


