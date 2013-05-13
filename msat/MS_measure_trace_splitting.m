% MS_MEASURE_TRACE_SPLITTING - Measure SWS for a pair of orthogonal traces
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% [fast, tlag] = MS_measure_trace_splitting(time, T00, T90, d_fast, ...
%     d_tlag, max_tlag)
%
%  Inputs:
%     time (vector) : Series of time points where amplitudes are known (s)
%     T00 (vector)  : Series of amplitudes of transverse wave #1 before
%                     splitting applied. 
%     T90 (vector)  : Series of amplitudes of transverse wave #2 before
%                     splitting applied.  
%     d_fast (scalar) : Spacing of grid search for fast direction (degrees)
%     d_tlag (scalar) : Spacing of grid search lag time (s)
%     max_tlag (scalar) : Maximum lag time in grid search
%
%  Outputs:
%     fast (vector)  : Optimum fast direction (degrees).
%     tlag (vector)  : Optimum lag tim.
%
% Given a two traces (T00 and T90) on a common time axis (time) find
% the splitting operator (fast and tlag) that minimises the second
% eigenvalue of the covariance matrix. This is done in two steps. First
% using an approximate grid search (spacing and time limits set by 
% d_fast, d_tlag and max_tlag) then by using the Matlab built
% in simplex optimiser from the best grid point. 

% Copyright (c) 2013, James Wookey and Andrew Walker
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



function [fast, tlag] = MS_measure_trace_splitting(time, T00, T90, ...
        d_fast, d_tlag, max_tlag)
    % Given a two traces (T00 and T90) on a common time axis (time) find
    % the splitting operator (fast and tlag) that minimises the second
    % eigenvalue of the covariance matrix. This is done in two steps. First
    % using an approximate grid search (spacing and time limits set by 
    % d_fast, d_tlag and max_tlag) then by using the Matlab built
    % in simplex optimiser from the best grid point.

    % Do the grid search
    best_misfit = 100;
    for fast = 0.0:d_fast:180.0
        for tlag = 0.0:d_tlag:max_tlag
            split_op = [fast, tlag];
            misfit = splitting_function(split_op, time, T00, T90);
            if (misfit < best_misfit)
                best_split_op = split_op;
                best_misfit = misfit;
            end
        end
    end
    
    % Simplex minimisation
    [new_split_op, misfit] = fminsearch(...
                @(split_op) splitting_function(split_op, time, ...
                T00, T90), best_split_op, optimset('Display','notify'));
    if (misfit < best_misfit)
        best_split_op = new_split_op;
    end
            
    % Unpack results
    fast = best_split_op(1);
    tlag = best_split_op(2);        
end

function misfit = splitting_function(split_op, time, T00, T90)
    % Apply a spliting operator to two traces and calculate the misfit in a
    % way that can be used by measure_splitting_function
    
    % Unpack args
    fast = split_op(1);
    tlag = split_op(2);
    % Calculate the traces with the revered split
    [T00sp,T90sp] = MS_split_trace(time,T00,T90,fast,-tlag) ;
    % measure second eignvalue of the resulting wavelet
    % calculate the covariance matrix
    COVM = cov(T90sp,T00sp) ;
    % take the eigenvalues
    [~,D] = eig(COVM) ;
    % calculate normalised misfit.
    misfit = min([D(1,1) D(2,2)])./ max([D(1,1) D(2,2)]) ;     
end
