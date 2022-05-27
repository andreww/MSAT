%% test_MS_polyaverage_stishovite


   % Do we get the same numbers as WEIDNER et al. table 5?
   [C, rh] = MS_elasticDB('stishovite');
   [ K_vrh, G_vrh, K_v, G_v, K_r, G_r ] = MS_polyaverage(C);
   assertElementsAlmostEqual([K_vrh], [316], 'absolute', 0.5);
   assertElementsAlmostEqual([G_vrh], [220], 'absolute', 0.5);
   assertElementsAlmostEqual([K_r], [308], 'absolute', 0.5);
   assertElementsAlmostEqual([G_r], [208], 'absolute', 0.5);
   assertElementsAlmostEqual([K_v], [324], 'absolute', 0.5);
   assertElementsAlmostEqual([G_v], [232], 'absolute', 0.5);   


% No error cases to check!