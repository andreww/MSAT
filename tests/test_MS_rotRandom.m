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

function test_MS_rotRandom_euler

  % Check that the output Euler angles give CC and vica versa.
  C = MS_elasticDB('ol');
  [CC, eulers] = MS_rotRandom(C, 'number', 100);
  
  for i = 1:100
      CC_check = MS_rotEuler(C, eulers(1,i), eulers(2,i), eulers(3,i));
      C_check = MS_rotEuler(CC(:,:,i), eulers(1,i), eulers(2,i), ...
          eulers(3,i), 'sense', 'passive');
      
      assertElementsAlmostEqual(CC(:,:,i), CC_check)
      assertElementsAlmostEqual(C, C_check)
  end
  
end
