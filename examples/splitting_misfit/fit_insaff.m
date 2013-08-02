% FIT_INSAFF - 

% Copyright (c) 2013, James Wookey and Andrew Walker
% Copyright (c) 2006, Joachim Vandekerckhove
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

function fit_insaff()

    % Fixed paramters - VP, VS rho (and h, which we need)

    vp = 3.600;
    vs = 1.850;
    rho = 920.0;
    [Ciso]=MS_iso(vp,vs,rho);
    Siso = inv(Ciso);
    nu = MS_poisson( Ciso, [1 0 0], [0 1 0]);
    h = (3.0*(1.0/Siso(1,1))*(2.0 - nu)) / ...
        (32.0*(1-nu^2.0));
        

    % read input
    dat = load('splitting_data_test.txt');
    azis = dat(:,1)';
    incs = 90.0-dat(:,2)';
    pols = dat(:,3)';
    avss = dat(:,4)';
    
    % setup model
    f = insaff_missfit_function(vp,vs,rho,h,azis, incs, pols, avss);
    
    options=anneal();
    options.Verbosity=2;
    
    [point mis] = anneal(@(p)f(p(1),p(2),p(3),p(4),p(5),p(6)),...
        [0.1 0.2 40.0 0.5 0.4 0.3],options) 
    zn = point(1);
    zt = point(2);
    alp = point(3);
    eps = point(4);
    del = point(5);
    gam = point(6);
    [C_model, ~] = MS_effective_medium('insaff', vp, vs, rho, ...
            abs(zn*zt*h), abs(zt*h), alp, eps, del, gam);
    MS_plot(C_model, rho, 'sdata', azis, incs, pols, avss);

end

function f = insaff_missfit_function(vp, vs, rho, h, azis, ...
    incs, pols, avss)

    obs_tlag = avss./vs; % tlags for unit distance
    dom_freq = 1.0/(max(obs_tlag)./4.0);
    f = @evaluate_misfit;
    function [misfit] = evaluate_misfit (zn, zt, alp, eps, del, gam)
        % Build the elastic model at this point
        [C_model, ~] = MS_effective_medium('insaff', vp, vs, rho, ...
            abs(zn*zt*h), abs(zt*h), alp, eps, del, gam);
        % Evaluate model splitting at data points
        [pol_model, avs, ~, ~, ~] = MS_phasevels(C_model, ...
            rho, incs, azis);
        % Turn vs into tlag...
        model_tlag = (avs./vs); %tlags for unit distance
        % Use vectorized splitting missfit to form missift.
        misfit = MS_splitting_misfit(pols, obs_tlag, pol_model, model_tlag,...
            0.0, dom_freq, 'mode','simple');
    end

end

function value = camel(x,y)

    value = (4-2.1*x.^2+x.^4/3).*x.^2+x.*y+4*(y.^2-1).*y.^2;

end



% Code below here is:
% Copyright (c) 2006, Joachim Vandekerckhove
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [minimum,fval] = anneal(loss, parent, options)
% ANNEAL  Minimizes a function with the method of simulated annealing
% (Kirkpatrick et al., 1983)
%
%  ANNEAL takes three input parameters, in this order:
%
%  LOSS is a function handle (anonymous function or inline) with a loss
%  function, which may be of any type, and needn't be continuous. It does,
%  however, need to return a single value.
%
%  PARENT is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the simulated annealing. If no
%  OPTIONS structure is provided, ANNEAL uses a default structure. OPTIONS
%  can contain any or all of the following fields (missing fields are
%  filled with default values):
%
%       Verbosity: Controls output to the screen.
%                  0 suppresses all output
%                  1 gives final report only [default]
%                  2 gives temperature changes and final report 
%       Generator: Generates a new solution from an old one.
%                  Any function handle that takes a solution as input and
%                  gives a valid solution (i.e. some point in the solution
%                  space) as output.
%                  The default function generates a row vector which
%                  slightly differs from the input vector in one element:
%                  @(x) (x+(randperm(length(x))==length(x))*randn/100)
%                  Other examples of possible solution generators:
%                  @(x) (rand(3,1)): Picks a random point in the unit cube
%                  @(x) (ceil([9 5].*rand(2,1))): Picks a point in a 9-by-5
%                                                 discrete grid
%                  Note that if you use the default generator, ANNEAL only
%                  works on row vectors. For loss functions that operate on
%                  column vectors, use this generator instead of the
%                  default:
%                  @(x) (x(:)'+(randperm(length(x))==length(x))*randn/100)'
%        InitTemp: The initial temperature, can be any positive number.
%                  Default is 1.
%        StopTemp: Temperature at which to stop, can be any positive number
%                  smaller than InitTemp. 
%                  Default is 1e-8.
%         StopVal: Value at which to stop immediately, can be any output of
%                  LOSS that is sufficiently low for you.
%                  Default is -Inf.
%       CoolSched: Generates a new temperature from the previous one.
%                  Any function handle that takes a scalar as input and
%                  returns a smaller but positive scalar as output. 
%                  Default is @(T) (.8*T)
%      MaxConsRej: Maximum number of consecutive rejections, can be any
%                  positive number.
%                  Default is 1000.
%        MaxTries: Maximum number of tries within one temperature, can be
%                  any positive number.
%                  Default is 300.
%      MaxSuccess: Maximum number of successes within one temperature, can
%                  be any positive number.
%                  Default is 20.
%
%
%  Usage:
%     [MINIMUM,FVAL] = ANNEAL(LOSS,NEWSOL,[OPTIONS]);
%          MINIMUM is the solution which generated the smallest encountered
%          value when input into LOSS.
%          FVAL is the value of the LOSS function evaluated at MINIMUM.
%     OPTIONS = ANNEAL();
%          OPTIONS is the default options structure.
%
%
%  Example:
%     The so-called "six-hump camelback" function has several local minima
%     in the range -3<=x<=3 and -2<=y<=2. It has two global minima, namely
%     f(-0.0898,0.7126) = f(0.0898,-0.7126) = -1.0316. We can define and
%     minimise it as follows:
%          camel = @(x,y)(4-2.1*x.^2+x.^4/3).*x.^2+x.*y+4*(y.^2-1).*y.^2;
%          loss = @(p)camel(p(1),p(2));
%          [x f] = ANNEAL(loss,[0 0])
%     We get output:
%               Initial temperature:     	1
%               Final temperature:       	3.21388e-007
%               Consecutive rejections:  	1027
%               Number of function calls:	6220
%               Total final loss:        	-1.03163
%               x =
%                  -0.0899    0.7127
%               f =
%                  -1.0316
%     Which reasonably approximates the analytical global minimum (note
%     that due to randomness, your results will likely not be exactly the
%     same).

%  Reference:
%    Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
%    Simulated Annealing. _Science, 220_, 671-680.

%   joachim.vandekerckhove@psy.kuleuven.be
%   $Revision: v5 $  $Date: 2006/04/26 12:54:04 $

def = struct(...
        'CoolSched',@(T) (.8*T),...
        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
        'InitTemp',1,...
        'MaxConsRej',1000,...
        'MaxSuccess',20,...
        'MaxTries',300,...
        'StopTemp',1e-8,...
        'StopVal',-Inf,...
        'Verbosity',1);

% Check input
if ~nargin %user wants default options, give it and stop
    minimum = def;
    return
elseif nargin<2, %user gave only objective function, throw error
    error('MATLAB:anneal:noParent','You need to input a first guess.');
elseif nargin<3, %user gave no options structure, use default
    options=def;
else %user gave all input, check if options structure is complete
    if ~isstruct(options)
        error('MATLAB:anneal:badOptions',...
            'Input argument ''options'' is not a structure')
    end
    fs = {'CoolSched','Generator','InitTemp','MaxConsRej',...
        'MaxSuccess','MaxTries','StopTemp','StopVal','Verbosity'};
    for nm=1:length(fs)
        if ~isfield(options,fs{nm}), options.(fs{nm}) = def.(fs{nm}); end
    end
end

% main settings
newsol = options.Generator;      % neighborhood space function
Tinit = options.InitTemp;        % initial temp
minT = options.StopTemp;         % stopping temp
cool = options.CoolSched;        % annealing schedule
minF = options.StopVal;
max_consec_rejections = options.MaxConsRej;
max_try = options.MaxTries;
max_success = options.MaxSuccess;
report = options.Verbosity;
k = 1;                           % boltzmann constant

% counters etc
itry = 0;
success = 0;
finished = 0;
consec = 0;
T = Tinit;
initenergy = loss(parent);
oldenergy = initenergy;
total = 0;
if report==2, fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy); end

while ~finished;
    itry = itry+1; % just an iteration counter
    current = parent; 
    
    % % Stop / decrement T criteria
    if itry >= max_try || success >= max_success;
        if T < minT || consec >= max_consec_rejections;
            finished = 1;
            total = total + itry;
            break;
        else
            T = cool(T);  % decrease T according to cooling schedule
            if report==2, % output
                fprintf(1,'  T = %7.5f, loss = %10.5f\n',T,oldenergy);
            end
            total = total + itry;
            itry = 1;
            success = 1;
        end
    end
    
    newparam = newsol(current);
    newenergy = loss(newparam);
    
    if (newenergy < minF),
        parent = newparam; 
        oldenergy = newenergy;
        break
    end
    
    if (oldenergy-newenergy > 1e-6)
        parent = newparam;
        oldenergy = newenergy;
        success = success+1;
        consec = 0;
    else
        if (rand < exp( (oldenergy-newenergy)/(k*T) ));
            parent = newparam;
            oldenergy = newenergy;
            success = success+1;
        else
            consec = consec+1;
        end
    end
end

minimum = parent;
fval = oldenergy;

if report;
    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', total);
    fprintf(1, '  Total final loss:        \t%g\n', fval);
end
end