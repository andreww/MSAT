%% test_MS_cijkl2cij_simple

    C = ones(6,6);
    CC = ones(3,3,3,3);
    assertElementsAlmostEqual(C, MS_cijkl2cij(CC));
    
