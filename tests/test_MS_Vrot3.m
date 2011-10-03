function test_suite = test_MS_Vrot3
initTestSuite;
end

function test_MS_Vrot3_simple

  V = [1 0 0];
  assertElementsAlmostEqual(MS_Vrot3(V, 0.0, 0.0, 0.0), V); 
  assertElementsAlmostEqual(MS_Vrot3(V', 0.0, 0.0, 0.0), V');
  assertElementsAlmostEqual(MS_Vrot3(V, 0.0, 180.0, 0.0), [-1 0 0]); 
  assertElementsAlmostEqual(MS_Vrot3(V', 0.0, 180.0, 0.0), [-1 0 0]');
end

function test_MS_Vrot3_errors

    V = [1 0 0];
    f = @()MS_Vrot3(V, 0.0, 0.0, 0.0, 'order', [1 2 4]);
    assertExceptionThrown(f, 'MS:VROT3:BadOrder')
    f = @()MS_Vrot3(V, 0.0, 0.0, 0.0, 'order', [1 3 2 4]);
    assertExceptionThrown(f, 'MS:VROT3:BadOrder')
    f = @()MS_Vrot3(V, 0.0, 0.0, 0.0, 'bob', [1 2 3]);
    assertExceptionThrown(f, 'MS:VROT3:UnknownOption')

    f = @()MS_Vrot3([1 0 0 0], 0.0, 0.0, 0.0);
    assertExceptionThrown(f, 'MS:VROT3BadInputMatrix')

end