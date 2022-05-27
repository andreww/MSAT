%% test_MS_effective_splitting_N_graphs

% try to check we get the solutions suggested by fig 2 (a) and (b)

expec = [-40 1.0];
[fe, tle] = MS_effective_splitting_N(1/20, 90, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res)
[fe, tle] = MS_effective_splitting_N(1/20, 0, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res)
[fe, tle] = MS_effective_splitting_N(1/5, 0, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res)
[fe, tle] = MS_effective_splitting_N(1/8, 0, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res)

expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

expec = [-20 0.8];
[fe, tle] = MS_effective_splitting_N(1/20, 10, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

%% test_MS_effective_splitting_N_graphs_2

% try to check we get the solutions suggested by fig 2 (c) and (d)

expec = [110-180 2.0];
[fe, tle] = MS_effective_splitting_N(1/20, 90, [140 90], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/20, 0, [140 90], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/5, 0, [140 90], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/8, 0, [140 90], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

expec = [-20 0.8];
[fe, tle] = MS_effective_splitting_N(1/20, 10, [90 140], [1.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

%% test_MS_effective_splitting_N_aggregation

% check that pre-aggregation of parameters works as expected.

% with an extra, cancelling pair of layers
expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 20 110 140], [1.0, 3.0, 3.0, 1.0]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

% last layer which stacks
expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140 140], [1.0, 0.5, 0.5]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

% last layer which partially cancels.
expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140 50], [1.0, 1.5, 0.5]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

% entirely cancelling stack
expec = [0 1.0];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [0 0 0 90 90 0 0], [1 1 1 2 2 1 1]);
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

%% test_MS_effective_splitting_N_errors

f = @()MS_effective_splitting_N(1/20, 10, [90], [1.0, 1.0]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80], [1.0]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80 70 60], [1 2 3 4 5]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80 70 60], [1 2 3]);
assertExceptionThrown(f, 'MS:ListsMustMatch')



% try to check we get the solutions suggested by fig 2 (a), (c) and (d)
% using the Gaussian Wavelet method. Seperate tests a bit because these 
% are slow.
%% test_MS_effective_splitting_N_graphs_1_GW
expec = [-40 1.0];
[fe, tle] = MS_effective_splitting_N(1/20, 90, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/20, 0, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/5, 0, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/8, 0, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1);

%% test_MS_effective_splitting_N_graphs_2_GW
% Errors (differences) become bigger as period gets longer
expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140], [1.0, 1.0], ...
     'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.15);

expec = [-20 0.8];
[fe, tle] = MS_effective_splitting_N(1/20, 10, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1);

%% test_MS_effective_splitting_N_graphs_3_GW

% try to check we get the solutions suggested by fig 2 (c) and (d)

expec = [110-180 2.0];
[fe, tle] = MS_effective_splitting_N(1/20, 90, [140 90], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
[fe, tle] = MS_effective_splitting_N(1/20, 0, [140 90], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)
% Errors (differences) become bigger as period gets longer
[fe, tle] = MS_effective_splitting_N(1/5, 0, [140 90], [1.0, 1.0], ...
     'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.25)
[fe, tle] = MS_effective_splitting_N(1/8, 0, [140 90], [1.0, 1.0], ...
     'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.25)

%% test_MS_effective_splitting_N_graphs_4_GW

% Errors (differences) become bigger as period gets longer
expec = [-25 0.8];
[fe, tle] = MS_effective_splitting_N(1/5, 10, [90 140], [1.0, 1.0], ...
     'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.15)

expec = [-20 0.8];
[fe, tle] = MS_effective_splitting_N(1/20, 10, [90 140], [1.0, 1.0], ...
    'mode', 'gaussianwavelet');
res = [fe tle];
assertElementsAlmostEqual(expec, res, 'relative', 0.1)

 