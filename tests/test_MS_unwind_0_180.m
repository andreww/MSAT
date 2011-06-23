function test_suite = test_MS_unwind_0_180
initTestSuite;
end

function test_MS_unwind_0_180_scalar
    assertElementsAlmostEqual(MS_unwind_0_180(30), 30);
    assertElementsAlmostEqual(MS_unwind_0_180(361), 1);
    assertElementsAlmostEqual(MS_unwind_0_180(-1), 179);
    assertElementsAlmostEqual(MS_unwind_0_180(0), 180);
    assertElementsAlmostEqual(MS_unwind_0_180(360), 180);
    assertElementsAlmostEqual(MS_unwind_0_180(180), 180);
end

function test_MS_unwind_0_180_vector
    assertElementsAlmostEqual(MS_unwind_0_180([30, 40]), [30, 40]);
    assertElementsAlmostEqual(MS_unwind_0_180([0, 361]), [180, 1]);
    assertElementsAlmostEqual(MS_unwind_0_180([1, -1]), [1, 179]);
    assertElementsAlmostEqual(MS_unwind_0_180([0,0,0]), [180,180,180]);
    assertElementsAlmostEqual(MS_unwind_0_180([360]), [180]);
end

function test_MS_unwind_0_180_error
    f = @()MS_unwind_0_180(1E24);
    assertExceptionThrown(f, 'MS:UNWIND:toomanyits');
end