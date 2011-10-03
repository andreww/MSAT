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
%     [] = MS_info(C, rho)                    
%         If the density is provided, velocity information is returned. 
%
%     [] = MS_info(C, rho, 'mode', MODE)
%         Choose the set of outputs to return. Currently supported are:
%           'all' : all available information (default)
%           'vel' : just report velocity information
%           'rot' : just report rotational information
%
%     [] = MS_info(C, 'mode', MODE)
%         As above, be default to 'rot'.
%
%     [report] = MS_info(C, ...)                    
%         As above, but return information in a string. Don't print to 
%         screen.
%

% Copyright (c) 2011, James Wookey and Andrew Walker
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

%===============================================================================
   function [ varargout ] = MS_info( C, varargin )
%===============================================================================

%  ** set defaults
      thresh = 0.01 ;
      imode = 'rot' ;
      inf_rot = 0;
      inf_vels = 0;
      rh = NaN ;
      OUT = '';

%  ** check required_inputs here
      MS_checkC(C) ;

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         if isscalar(varargin{iarg}) % a density was specified
            rh = varargin{iarg} ;
            iarg = iarg + 1;
            imode = 'all';
         else
            switch lower(varargin{iarg})
                case 'mode'  % parameter definition (value required)
                  imode = varargin{iarg+1} ;
                  iarg = iarg + 2 ;
                case 'vel'
                  imode = varargin{iarg+1} ;
                  iarg = iarg + 2 ;
                case 'rot'
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
          case 'vel'
            inf_vels = 1 ;
            if isnan(rh)
               error('MS:INFO:NoVelocityInfo', ...
                  ['Velocity information is not available if a ', ...
                  'density is not specified']) ;
            end   
          case 'rot'
            inf_rot = 1 ;
          otherwise
            error('MS:INFO:UnsupportedMode',...
               [ 'Unsupported mode: ' imode]) ;               
      end

%  ** rotational information
      if inf_rot
         OUT = sprintf([OUT, '\n---ROTATIONAL INFORMATION \n\n']);
         try 
            [CR,RR] = MS_axes(C) ;
            ind = find(abs(RR) > thresh & abs(RR) < (1-thresh)) ;
            if ~isempty(ind)
               OUT = sprintf([OUT, '   Principle axes appear to be in a different', ...
                       ' orientation to this matrix.\n']);
             OUT = sprintf([OUT, '\n    rotation matrix  =\n\n']);
             OUT = sprintf([OUT, '%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f\n'],... 
                  RR(1,1),RR(1,2),RR(1,3),RR(2,1),RR(2,2),RR(2,3),RR(3,1),RR(3,2),RR(3,3));       
            else
              OUT = sprintf([OUT, '   Elasticity matrix appears unrotated.\n']);
            end  
         catch ME
            OUT = sprintf([OUT, '   Principle axis check failed.\n\n']);
         end
      end
      
%  ** Velocity information
      if inf_vels
         OUT = sprintf([OUT, '\n---VELOCITY INFORMATION \n\n']);

         CR = MS_axes(C, 'nowarn') ;
         [Ciso] = MS_decomp(CR) ;
         
         OUT = sprintf([OUT, '   Isotropic velocities: Vp = %6.3f ,  Vs = %6.3f\n\n'], ...
            1e-3.*sqrt((Ciso(3,3)*1e9)./rh), 1e-3.*sqrt((Ciso(4,4)*1e9)./rh)) ; 
            
         
         azi = [  0  -90   0 ] ;
         inc = [  0    0   90] ;
         
         [pol,avs,vs1,vs2,vp] = MS_phasevels(C,rh,inc,azi) ;
         
         OUT = sprintf([OUT, '   Direction      VP    VS1    VS2  VSAni \n']) ;
         OUT = sprintf([OUT, '          X1  %6.3f %6.3f %6.3f %6.3f \n'], ...
            vp(1),vs1(1),vs2(1),avs(1)) ;
         OUT = sprintf([OUT, '          X2  %6.3f %6.3f %6.3f %6.3f \n'], ...
            vp(2),vs1(2),vs2(2),avs(2)) ;
         OUT = sprintf([OUT, '          X3  %6.3f %6.3f %6.3f %6.3f \n'], ...
            vp(3),vs1(3),vs2(3),avs(3)) ; 
         
         
         
      end
      switch nargout
          case 0
              fprintf(OUT)   
          case 1
              varargout(1) = {OUT} ;
          otherwise
              error('MS:BADOUTPUT', ...
                  'MS_info requires 0-1 output arguments.') ;
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


