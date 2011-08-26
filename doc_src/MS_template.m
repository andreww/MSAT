% MS_TEMPLATE - Template function for MSAT.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Short description of functionality. 
%
%  % [ required_outputs, ... ] = MS_template( required_inputs, ... )
%
% Usage: 
%     [required_outputs] = MS_template(required_inputs)                    
%         Basic behaviour 
%
%     [required_outputs] = MS_load(required_inputs,'option','value')                    
%         Description of input option
%     ...
%
%     [required_outputs, optional_output] = MS_load(required_inputs,...)
%         Description of optional output
%     ...
%
% Notes:
%     Notes on operation.
%
% See also: MS_NORMS, MS_INTERPOLATE

% Change log / copyright statements.
% 
%

%===============================================================================
   function [ varargout ] = MS_template( required_inputs, varargin )
%===============================================================================

%  ** set defaults
      isoption1 = 1 ;
      val2default = 1 ;

%  ** check required_inputs here
%  ...

%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case 'option1' % flag (i.e., no value required)
               aij = 1 ;
               iarg = iarg + 1 ;
            case 'value2'  % parameter definition (value required)
               val2default = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            otherwise 
               error('MS:TEMPLATE:UnknownOption',...
                  ['Unknown option: ' varargin{iarg}]) ;   
         end   
      end

%
%  ** MAIN CODE
%
      required_output1 = 0 ;
      optional_output1 = 0 ;
      optional_output2 = 0 ;

%  ** construct outputs
      switch nargout
      case 1
         varargout(1) = {required_output1} ;
      case 2
         varargout(1) = {required_output1} ;
         varargout(2) = {optional_output1} ;
      case 3
         varargout(1) = {required_output1} ;
         varargout(2) = {optional_output1} ;
         varargout(3) = {optional_output2} ;
      otherwise
         error('MS:TEMPLATE:BadOutputArgs','Requires 1-3 output arguments.')
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


