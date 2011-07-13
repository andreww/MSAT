function test_suite = test_MS_unwind_pm_90
initTestSuite;
end

function test_MS_unwind_0_900_scalar
    assertElementsAlmostEqual(MS_unwind_pm_90(30), 30);
    assertElementsAlmostEqual(MS_unwind_pm_90(361), 1);
    assertElementsAlmostEqual(MS_unwind_pm_90(-1), -1);
    assertElementsAlmostEqual(MS_unwind_pm_90(0), 0);
    assertElementsAlmostEqual(MS_unwind_pm_90(360), 0);
    assertElementsAlmostEqual(MS_unwind_pm_90(180), 0);
    assertElementsAlmostEqual(MS_unwind_pm_90(-180), 0);
end

function test_MS_unwind_0_90_vector
    assertElementsAlmostEqual(MS_unwind_pm_90([30, 40]), [30, 40]);
    assertElementsAlmostEqual(MS_unwind_pm_90([0, 361]), [0, 1]);
    assertElementsAlmostEqual(MS_unwind_pm_90([1, -1]), [1, -1]);
    assertElementsAlmostEqual(MS_unwind_pm_90([0,0,0]), [0,0,0]);
    assertElementsAlmostEqual(MS_unwind_pm_90([360]), [0]);
end

function test_MS_unwind_0_90_error
    f = @()MS_unwind_pm_90('bob');
    assertExceptionThrown(f, 'MS:unwind_pm_90:BadInput');
end