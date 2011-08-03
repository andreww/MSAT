function test_suite = test_MS_elasticDB
initTestSuite;
end

function test_MS_elasticDB_run

    infol = 'Single crystal olivine (Abramson et al, JGR, 1997; doi:10.1029/97JB00682)' ; 
    Cl = [320.5  68.1  71.6   0.0   0.0   0.0 ; ...
           68.1 196.5  76.8   0.0   0.0   0.0 ; ...
           71.6  76.8 233.5   0.0   0.0   0.0 ; ...
            0.0   0.0   0.0  64.0   0.0   0.0 ; ...
            0.0   0.0   0.0   0.0  77.0   0.0 ; ...
            0.0   0.0   0.0   0.0   0.0  78.7 ];
    rhl = 3355 ;

    [ C ] = MS_elasticDB('ol');
    assertElementsAlmostEqual(C, Cl)
    
    [ C, rh ] = MS_elasticDB('ol');
    assertElementsAlmostEqual(C, Cl)
    assertElementsAlmostEqual(rh, rhl)
    
    [ C, rh, info ] = MS_elasticDB('ol');
    assertElementsAlmostEqual(C, Cl)
    assertElementsAlmostEqual(rh, rhl)
    assert(all(info == infol))
    
end

function test_MS_elasticDB_errors

    f = @()MS_elasticDB('badname');
    assertExceptionThrown(f, 'MS:ELASTICDB:UNKNOWN')

    % This is open coded so that we can check the output
    % ... cannot work out how to do this by passing a function
    % handle into assertExceptionThrown. Maybe a patch to xunit
    % is in order
    val = 0;  
    try
        [out1, out2, out3, out4] = MS_elasticDB('ol');
    catch e
        assert(strcmp(e.identifier, ...
      'MS:BADOUTPUT'), 'Wrong error thrown')
       val = 1;
    end
    assert(val==1, 'No error thrown')

end