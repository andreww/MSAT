% MS_PLOT_TRACE - Plot one or more pairs of orthogonal traces
% 
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% MS_plot_trace(time, T00, T90, ...)
%
%  Inputs:
%     time (row vector)         : Series of time points where amplitudes
%                                 are known (s)  
%     T00 (row vector / matrix) : Series of amplitudes in first direction.
%     T90 (row vector / matrix) : Series of amplitudes in second 
%                                 direction.
%
%  Usage:
%      MS_plot_trace(time, T00, T90): If T00 and T90 are row vectors, plot %          amplitudes of a single trace as two curves on time axis. For 
%          multiple traces T00 and T90 for each trace can be 'stacked' to
%          form a matrix, e.g. MS_plot_trace(time,[T00a;T90b],[T00a;T90b]),
%          but note that these must be on a common time axis. Subplots are
%          then used to draw each trace.
%      MS_plot_trace(time, T00, T90, 'plots', plotstring): Alter the kind 
%          and/or number of plots for each trace. The variable plotstring 
%          is a series of single letter codes, each letter causes 
%          MS_plot_trace to generate a different kind of plot (again 
%          subplots may be used for multiple kinds of plot). Permitted 
%          letter codes are 'w': the default plot type described above; 
%          'd': plot of particle displacement, i.e. T00 vs T90.
%      MS_plot_trace(time, T00, T90, 'headings', {'first trace name'; ...
%          'second trace name'}): Give names to each trace in the plot.

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


function [] = MS_plot_trace(time, T00, T90, varargin) 
    % Plot some possibly useful diagnostics;

    assert(all(size(T00)==size(T90)), 'MS:PLOT_TRACE:invalidargs',...
        ['The sizes of the traces T00 and T90 must match in ' ...
         'MS_plot_trace']);
    assert(length(T00(1,:))==length(time), 'MS:PLOT_TRACE:invalidargs',...
        ['The length of T00(1,:) and the time axis must match in ' ...
         'MS_plot_trace']);
    
    % defaults
    plots = 'w';
    withheadings = 0;
    
    % process the optional arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'plots'
                plots = varargin{iarg+1};
                iarg = iarg + 2 ;
            case 'headings'
                headings = varargin{iarg+1};
                iarg = iarg + 2 ;
                withheadings = 1;
                assert(length(headings(:,1))==length(T00(:,1)),...
                    'MS:PLOT_TRACE:invalidargs', ...
                    ['The number of headings must match the number of '...
                     'traces in MS_plot_trace']);
            otherwise
                error('MS:PLOT_TRACE:UnknownOption',...
                    ['Unknown option in MS_plot_trace: ' varargin{iarg}]);
        end
    end
    
    % How many types of plot (width)
    n = size(plots);
    n = n(2);
    
    % How many pairs of traces (height)
    m = size(T00);
    m = m(1);
    
    figure
    pnum = 0;
    for i = 1:n
        for j = 1:m
            pnum = pnum + 1;
            subplot(n, m, pnum);
            switch lower(plots(1,i))
                case 'w'
                    plot_wavelet(time,T00(j,:),T90(j,:));
                case 'd'
                    plot_displace(T00(j,:), T90(j,:));
                otherwise
                    error('MS:PLOT_TRACE:UnknownType',...
                    ['Unknown plot type in MS_plot_trace: ' plots(1,i)]);
            end
            if (withheadings && (i == 1))
                title(headings(j,:));  
            end
        end
    end
  
end

function [] = plot_wavelet(time,T00,T90)
    plot(time,T00,'b-')
    hold on
    plot(time,T90,'r-')
end

function [] = plot_displace(T00, T90)
    plot(T00, T90, 'k-')
end
