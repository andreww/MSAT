%% test_MS_build_isotropic_inputs

    [C r] = MS_elasticDB('ol');
    [K_vrh, G_vrh] = MS_polyaverage(C);
    [C, K, E, lam, mu, nu, M] = MS_build_isotropic('K', K_vrh,'mu',G_vrh);
    assertElementsAlmostEqual(MS_build_isotropic('lam', C(1,2),'mu',C(4,4)), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', C(4,4),'lam',C(1,2)), C);
    assertElementsAlmostEqual(MS_build_isotropic('M', C(1,1),'mu',C(4,4)), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', C(4,4),'M',C(1,1)), C);
    
    assertElementsAlmostEqual(MS_build_isotropic('e', E,'mu',mu), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', mu,'E',E), C);
    assertElementsAlmostEqual(MS_build_isotropic('K', K,'lam',lam), C);
    assertElementsAlmostEqual(MS_build_isotropic('lAm',lam, 'k', K), C);
    assertElementsAlmostEqual(MS_build_isotropic('k', K,'mu',mu), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', mu,'K',K), C);
    assertElementsAlmostEqual(MS_build_isotropic('nu', nu,'lam',lam), C);
    assertElementsAlmostEqual(MS_build_isotropic('laM', lam,'NU',nu), C);
    assertElementsAlmostEqual(MS_build_isotropic('nu', nu,'mu',mu), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', mu,'nu',nu), C);
    assertElementsAlmostEqual(MS_build_isotropic('e', E,'nu',nu), C);
    assertElementsAlmostEqual(MS_build_isotropic('Nu', nu,'E',E), C);
    assertElementsAlmostEqual(MS_build_isotropic('k', K,'nu',nu), C);
    assertElementsAlmostEqual(MS_build_isotropic('nu', nu,'K',K), C);
    assertElementsAlmostEqual(MS_build_isotropic('k', K,'E',E), C);
    assertElementsAlmostEqual(MS_build_isotropic('e', E,'K',K), C);
    assertElementsAlmostEqual(MS_build_isotropic('M', M,'mu',mu), C);
    assertElementsAlmostEqual(MS_build_isotropic('mu', mu,'M',M), C);



%% test_MS_build_isotropic_errors

    [C r] = MS_elasticDB('ol');
    [K_vrh, G_vrh] = MS_polyaverage(C);
    C = MS_build_isotropic('K', K_vrh,'mu',G_vrh);
    
    f = @()MS_build_isotropic('lam', C(1,2),'lam',C(1,2));
    assertExceptionThrown(f, 'MS:buildiso:wrongargs');
    
    f = @()MS_build_isotropic('bob', C(1,2),'lam',C(1,2));
    assertExceptionThrown(f, 'MS:buildiso:wrongargs');
    
    f = @()MS_build_isotropic('lam', C(1,2),'bob',C(1,2));
    assertExceptionThrown(f, 'MS:buildiso:wrongargs');
    
