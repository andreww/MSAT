function test_suite = test_MS_rotM
initTestSuite;
end

function test_MS_rotM_simple

  I = [1 0 0 ; 0 1 0 ; 0 0 1];
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0), I); 
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [1 2 3]), I);
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [1 3 2]), I);
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [3 2 1]), I);
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [3 1 2]), I);
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [2 3 1]), I);
  assertElementsAlmostEqual(MS_rotM(0.0, 0.0, 0.0, 'order', [2 1 3]), I);
end

function test_MS_rotM_errors

    I = [1 0 0 ; 0 1 0 ; 0 0 1];
    f = @()MS_rotM(0.0, 0.0, 0.0, 'order', [1 2 4]);
    assertExceptionThrown(f, 'MS:ROTM:BadOrder')
    f = @()MS_rotM(0.0, 0.0, 0.0, 'order', [1 3 2 4]);
    assertExceptionThrown(f, 'MS:ROTM:BadOrder')
    f = @()MS_rotM(0.0, 0.0, 0.0, 'bob', [3 2 1]);
    assertExceptionThrown(f, 'MS:ROTM:UnknownOption')

end