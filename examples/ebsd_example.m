% Script to make VRH avarage from EBSD data

[C, roh] = MS_elasticDB('ol')

% crystal symmetry
CS = {...
symmetry('m-3m'),...
symmetry('m-3m')};
% specimen symmetry
SS = symmetry('-1');
% specify file name
fname = fullfile(mtexDataPath,'aachen_ebsd','85_829grad_07_09_06.txt');
% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,SS,'interface','generic' ...
, 'ColumnNames', { 'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3' 'MAD' 'BC' 'BS' 'Bands' 'Error' 'ReliabilityIndex'}, 'Bunge', 'ignorePhase', 0)

All_Cs = MS_rot_mtex_rotation(C, get(ebsd,'orientations'));

voigt_ave = zeros(6,6) ;
reuss_ave = zeros(6,6) ;
for i = 1:49364
    voigt_ave = voigt_ave + cvec(:,:,i) .* (1/49364);
    reuss_ave = reuss_ave + inv(cvec(:,:,i)) .* (1/49364);
end

Cave = (voigt_ave + inv(reuss_ave)) ./ 2