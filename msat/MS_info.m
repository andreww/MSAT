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

      rh = NaN ;

%  ** check required_inputs here
      MS_checkC(C) ;

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         if isscalar(varargin{iarg}) % a density was specified
            rh = varargin{iarg} ;
            iarg = iarg + 1;
         else
            switch lower(varargin{iarg})
               case 'mode'  % parameter definition (value required)
                  imode = varargin{iarg+1} ;
                  iarg = iarg + 2 ;
               otherwise 
                  error('MS:INFO:UnknownOption',...
                     ['Unknown option: ' varargin{iarg}]) ;   
            end
         end   
      end
      
      switch lower(imode)
         case 'all'
            inf_vels = 1 ;
            
            if isnan(rh)
               warning('MS:INFO:NoVelocityInfo', ...
                  ['Velocity information is not available if a ', ...
                  'density is not specified']) ;
               inf_vels = 0 ;
            end   
            inf_rot = 1 ;
         otherwise
            error('MS:INFO:UnsupportedMode',...
               [ 'Unsupported mode: ' varargin{iarg}]) ;               
      end

%  ** rotational information
      if inf_rot
         fprintf('\n---ROTATIONAL INFORMATION \n\n')
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
            else
              fprintf(['   Elasticity matrix appears unrotated.\n'])
            end  
         catch ME
            fprintf('   Principle axis check failed.\n\n')
         end
      end
      
%  ** Velocity information
      if inf_vels
         fprintf('\n---VELOCITY INFORMATION \n\n')

         CR = MS_axes(C, 'nowarn') ;
         [Ciso] = MS_decomp(CR) ;
         
         fprintf('   Isotropic velocities: Vp = %6.3f ,  Vs = %6.3f\n\n', ...
            1e-3.*sqrt((Ciso(3,3)*1e9)./rh), 1e-3.*sqrt((Ciso(4,4)*1e9)./rh)) ; 
            
         
         azi = [  0  -90   0 ] ;
         inc = [  0    0   90] ;
         
         [pol,avs,vs1,vs2,vp] = MS_phasevels(C,rh,inc,azi) ;
         
         fprintf('   Direction      VP    VS1    VS2  VSAni \n') ;
         fprintf('          X1  %6.3f %6.3f %6.3f %6.3f \n', ...
            vp(1),vs1(1),vs2(1),avs(1)) ;
         fprintf('          X2  %6.3f %6.3f %6.3f %6.3f \n', ...
            vp(2),vs1(2),vs2(2),avs(2)) ;
         fprintf('          X3  %6.3f %6.3f %6.3f %6.3f \n', ...
            vp(3),vs1(3),vs2(3),avs(3)) ; 
         
         
         
      end
      

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


