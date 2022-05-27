%% test_MS_anisotropy_isotropic

  Ciso = [166.6667   66.6667   66.6667         0         0         0; ...
           66.6667  166.6667   66.6667         0         0         0; ...
           66.6667   66.6667  166.6667         0         0         0; ...
                 0         0         0   50.0000         0         0; ...
                 0         0         0         0   50.0000         0; ...
                 0         0         0         0         0   50.0000];

  [ uA, lma, zA, cbA ] = MS_anisotropy( Ciso );
    
  assertElementsAlmostEqual(uA, 0.0);
  assertElementsAlmostEqual(lma, 1.0);
  assertElementsAlmostEqual(zA, 1.0);
  assertElementsAlmostEqual(cbA, 0.0);


