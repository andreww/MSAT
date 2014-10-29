function test_suite = test_MS_rotRandom
initTestSuite;
end

function test_MS_rotRandom_isotropic

  % Check that a large number of randomly oriented crystals gives
  % an isotropic avarage. Can we think of better tests?
  
  C = MS_elasticDB('ol');

  % Calculate Universal anisotropy index of VRH mean of 10000 
  % random orientations using a fixed PRN stream.
  uA = MS_anisotropy(MS_VRH(ones(10000,1), ...
      MS_rotRandom(C, 'number', 10000, ...
                   'stream', RandStream('mrg32k3a', 'seed', 1)...
                   ), ones(10000,1)));
  assert(uA<1E-4, 'Not isotropic')
end
