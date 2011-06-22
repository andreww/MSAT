function test_suite = test_MS_unwind_360
initTestSuite;
end

function test_MS_unwind_360_scalar
    assertElementsAlmostEqual(MS_unwind_360(30), 30);
    assertElementsAlmostEqual(MS_unwind_360(361), 1);
    assertElementsAlmostEqual(MS_unwind_360(-1), 359);
    assertElementsAlmostEqual(MS_unwind_360(0), 0);
    assertElementsAlmostEqual(MS_unwind_360(360), 0);
end

function test_MS_unwind_360_vector
    assertElementsAlmostEqual(MS_unwind_360([30, 40]), [30, 40]);
    assertElementsAlmostEqual(MS_unwind_360([0, 361]), [0, 1]);
    assertElementsAlmostEqual(MS_unwind_360([1, -1]), [1, 359]);
    assertElementsAlmostEqual(MS_unwind_360([0,0,0]), [0,0,0]);
    assertElementsAlmostEqual(MS_unwind_360([360]), [0]);
end

function test_MS_unwind_360_error
    f = @()MS_unwind_360(1E24);
    assertExceptionThrown(f, 'MS:UNWIND:toomanyits');
end