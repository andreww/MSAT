% MS_SAVE - Write a set of elastic constants from a file.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
% 
% MS_save( fname, C, rho, ... ) or MS_save( fid, C, rho, ... )
%
% Usage: 
%     MS_save(fname, C, rho)                    
%          Write the elastic constants, C, and density, rho, to file 
%          'fname'. C must be provided in GPa and rho must be provided in 
%          kg/m3. See file format below. 
%
%
%     MS_save(fid, C, rho)                    
%          Write the elastic constants, C, and density, rho, to already 
%          open file described by integer valued double 'fid'. Fid is
%          usually returned by a previous call to fopen. The file is not
%          closed on exit from this function and no checks are made that
%          the file is opened and can be used for output. C must be 
%          provided in GPa and rho must be provided in kg/m3. 
%          See file format below. 
%
%     MS_save(fname, C, rho, ..., 'format',fmt)                    
%         Specify format for file. See below for descriptions. Available 
%         options are:
%            'simple' : Formatted text file (default)
%            'ematrix' : Format read by EMATRIX (D. Mainprice).
%
%     MS_save(fname, C, rho, ..., 'eunit',unit)
%         Specify elasticity unit; this is used to convert from GPa (msat's 
%         default operating unit)
%            'GPa' - elastic constants are written in GPa (default)
%            'Mbar' - elastic constants are written in Mbar
%            'Pa' - elastic constants are written in Pascals
%            'bar' - elastic constants are written in bar
%         Note: when using a density normalised tensor, specify the unit 
%         *before* normalisation (i.e., a unit of pressure).
%
%     MS_save(fname, C, rho, ..., 'dunit',unit)
%         Specify density unit; this is used to convert from kg/m3 (msat's 
%         default operating unit)
%            'kgm3' - density is written in kg/m3 (default)
%            'gcc' - density is written g/cc
%
%     MS_save(fname, C, rho, ..., 'Aij')
%         Elasticities are written in density normalised form. This option
%         requires a density to be specified. The provided elastic
%         constants are first converted to the unit specified by the
%         'eunit' argument (if provided) and the density converted to the 
%         unit specified in the by the 'dunit' argument (if provided) then 
%         the elastic constants are divided by the density before being
%         written to the file.
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
% See also: MS_LOAD

% Copyright (c) 2012, James Wookey and Andrew Walker
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

function [varargout] = MS_save(fname, C, rho, varargin)

    eunit = 'default' ;
    dunit = 'default' ;

    aij = 0 ;
    smode = 'none' ;
    force = 0 ;
    format = 'simple' ;

    % Process the optional arguments
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
            otherwise
                error(['Unknown option: ' varargin{iarg}]) ;
        end
    end

    % Check input argument is sane
    if ~MS_checkC(C)
        error('MS:SAVE:badC', 'Argument C to MS_save is not valid')
    end

    % Set default units if this hasn't already been done.
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

    % Perform unit conversions
    switch lower(eunit)
        case 'gpa'
            % Do nothing
        case 'pa'
            C = C .* 1e9 ;
        case 'mbar'
            C = C ./ 100 ;
        case 'bar'
            C = C .* 10e3 ;
        otherwise
            error('MS:SAVE:BadEunit','Unsupported elasticity unit.')
    end

    switch lower(dunit)
        case 'kgm3'
        case 'gcc'
            rho = rho ./ 1e3 ;
        otherwise
            error('MS:SAVE:BadDunit','Unsupported density unit.')
    end

    % Normalise for density (if required)
    if aij
        C = C ./ rho;
    end

    % Save the file
    switch lower(format)
        case 'simple'
            MS_write_simple( fname, C, rho ) ;
        case 'ematrix'
            MS_write_ematrix( fname, C, rho ) ;
        otherwise
            error('MS:SAVE:BadFileFormat',...
                ['Specified file format "' format '" is not supported.']) ;
    end

end

function MS_write_simple( fname, C, rho)

    [fid, fid_state] = MS_open_write_fh(fname);
    
    for i = 1:6
        for j = i:6
            fprintf(fid, '%1i %1i %8e\n', i, j, C(i,j));
        end
    end
    fprintf(fid, '%1i %1i %f\n', 7, 7, rho);
    
    MS_close_write_fh(fid, fid_state);

end

function MS_write_ematrix( fname, C, rho)
      error('MS:SAVE:NotImplemented', 'Writing to ematrix format is not implemented')
end

function [fid, fid_state] = MS_open_write_fh(fname)
    % Open file handle and remember if it needs to be closed
    
    if ischar(fname) 
        % This is a string - open it and return fid
        fid = fopen(fname, 'wt'); % Should we use t here?
        fid_state = 1; % We will need to close this fif
    elseif isreal(fname) 
        % Not a string, but a number, so assume this is an fid
        % we should use, but not close later
        fid = fname;
        fid_state = 0; % We will not need to close this one
    else
        % This isn't going to work
        error('Cannot open file for write\n');
    end
end

function MS_close_write_fh(fid, fid_state)
    % If *we* opened the file handle, close it.
    if fid_state
        fclose(fid);
    end
end
      
      
      

