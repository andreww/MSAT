function test_suite = test_MS_cij2cijkl
initTestSuite;
end

function test_MS_cij2cijkl_simple

    C = ones(6,6);
    CC = ones(3,3,3,3);
    assertElementsAlmostEqual(CC, MS_cij2cijkl(C));
    
end