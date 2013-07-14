
fprintf('Isotropic olivine for the matrix...')
[Col rol] = MS_elasticDB('ol');
[~, ~, Kol, Gol, ~, ~] = MS_polyaverage(Col);
Col_iso = build_isotropic(Kol, Gol)

fprintf('Isotropic diopside for the inclusion...')
[Cdi, rdi] = MS_elasticDB('di');
[~, ~, Kdi, Gdi, ~, ~] = MS_polyaverage(Cdi);
Cdi_iso = build_isotropic(Kdi, Gdi)

% triinc(Col_iso, Col_iso, [1 0 0; 0 1 0; 0 0 1], 1.0, 1.0, 1.0, 0.25)

fprintf('Checking we get the Eshelby tensor right...')
longaxis = [0.5 ; 0.75 ; 0.9; 0.95; 0.99999; 1.00001; 1.05; 1.1 ; 1.25; 1.5];
for i = 1:size(longaxis)
    Sanal = tandon_weng_num(vp_ol,vs1_ol,5000,longaxis(i),0.25,vp_di,vs1_di,5000);
    Snum = MS_cijkl2cij(eshelby(MS_cij2cijkl(Col_iso), longaxis(i), 1.0, 1.0));
    assertElementsAlmostEqual(Sanal, Snum, 'absolute', 0.005);
end
 

fprintf('Inclusion calculations...')



 [~, ~, vs1_ol, ~, vp_ol] = MS_phasevels(Col_iso, 5000, 0, 0);
 [~, ~, vs1_di, ~, vp_di] = MS_phasevels(Cdi_iso, 5000, 0, 0);
% %MS_tandon_and_weng(vp_ol,vs1_ol,5000,0.9999999,0.25,vp_di,vs1_di,5000)
% MS_TW_1 = MS_tandon_and_weng(vp_ol,vs1_ol,5000,1.0000001,0.25,vp_di,vs1_di,5000)
% MS_TW_2 = MS_tandon_and_weng(vp_ol,vs1_ol,5000,1.0000001,0.025,vp_di,vs1_di,5000)
% MS_TW_3 = MS_tandon_and_weng(vp_ol,vs1_ol,5000,1.0000001,0.0025,vp_di,vs1_di,5000)
% MS_TW_5 = MS_tandon_and_weng(vp_ol,vs1_ol,5000,1.0000001,0.000025,vp_di,vs1_di,5000)
% 
% TI_1 = triinc(Col_iso, Cdi_iso, [1 0 0; 0 1 0; 0 0 1], 1.0000001, 0.5, 1.0, 0.5)
MS_TW_3 = MS_effective_medium('t&w', Col_iso, rol,Cdi_iso, rdi, 1.1,0.025)
% TI_3 = triinc(Col_iso, Cdi_iso, [1 0 0; 0 1 0; 0 0 1], 1.1, 1.0, 1.0, 0.025)
[ Cout, Sout] = esh_inclusion(Col_iso, Cdi_iso, 1.0, 1.0, 1.1, 0.025)


