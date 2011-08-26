function test_suite = test_MS_phasevels
initTestSuite;
end

function test_MS_phasevels_stishovite_graph

    % Check we get the same results as those given in Mainprices review 
    % (Figure 3).
    
    
    [C, rh] = MS_elasticDB('stishovite');
    
    % [001], symmetry axis
    [pol, avs, vs1, vs2, vp] = MS_phasevels(C, rh, 90, 0);
    
    assert(vs1==vs2, 'VS1 ne VS2'); % By symmetry
    assert(avs==0, 'avs ne 0'); % By symmetry 
    assert((vp-13.5)^2<0.1^2, 'vp wrong'); % From graph
    assert((vs1-7.7)^2<0.1^2, 'vs wrong'); % From graph
    assert(isnan(pol), 'pol is wrong')
    
    % [100]
    [~, ~, vs1, vs2, vp] = MS_phasevels(C, rh, 0, 0);
    
    assert((vp-10.2)^2<0.1^2, 'vp wrong'); % From graph
    assert((vs1-8.4)^2<0.1^2, 'vs wrong'); % From graph
    assert((vs2-7.7)^2<0.1^2, 'vs wrong'); % From graph
    
    % [010]
    [~, ~, vs1, vs2, vp] = MS_phasevels(C, rh, 0, 90);
    
    assert((vp-10.2)^2<0.1^2, 'vp wrong'); % From graph
    assert((vs1-8.4)^2<0.1^2, 'vs wrong'); % From graph
    assert((vs2-7.7)^2<0.1^2, 'vs wrong'); % From graph

end 

function test_MS_phasevels_stishovite_list

    [C, rh] = MS_elasticDB('stishovite');
    
    % [001], symmetry axis
    [~, avs, vs1, vs2, vp] = MS_phasevels(C, rh, [90 90 90], [0 0 0]);
    
    assertElementsAlmostEqual(vs1, vs2); % By symmetry
    assertElementsAlmostEqual(avs, [0.0, 0.0, 0.0]'); % By symmetry 
    assertElementsAlmostEqual(vp, [13.5, 13.5, 13.5]', 'absolute', 0.5); % From graph
    assertElementsAlmostEqual(vs1, [7.7, 7.7, 7.7]', 'absolute', 0.5); % From graph
    

end

function test_MS_phasevels_stishovite_errors

    [C, rh] = MS_elasticDB('stishovite');
    
    f = @()MS_phasevels(C, rh, [90 90 90], [0 0]);
    assertExceptionThrown(f, 'MS:ListsMustMatch');
    
    f = @()MS_phasevels(C, rh, [90 90], [0 0 0]);
    assertExceptionThrown(f, 'MS:ListsMustMatch');
    
    % How shoould we test MS:PHASEVELS:vectornotreal
    
end