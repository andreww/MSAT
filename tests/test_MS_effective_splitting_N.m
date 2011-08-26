function test_suite = test_MS_effective_splitting_N
initTestSuite;
end

function test_MS_effective_splitting_N_graphs

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

end 

function test_MS_effective_splitting_N_graphs_2

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

end 

function test_MS_effective_splitting_N_errors

f = @()MS_effective_splitting_N(1/20, 10, [90], [1.0, 1.0]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80], [1.0]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80 70 60], [1 2 3 4 5]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

f = @()MS_effective_splitting_N(1/20, 10, [90 80 70 60], [1 2 3]);
assertExceptionThrown(f, 'MS:ListsMustMatch')

end