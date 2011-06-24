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
, 'ColumnNames', { 'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3' 'MAD' 'BC' 'BS' 'Bands' 'Error' 'ReliabilityIndex'}, 'Bunge', 'ignorePhase', 0);

% Just phase 1...
[phi1s, thetas, phi2s] = Euler(get(ebsd(1),'orientations'), 'Bundge');
phi1s=phi1s*(180.0/pi);
thetas=thetas*(180.0/pi);
phi2s=phi2s*(180.0/pi);
num_xtals = length(phi1s);

%All_Cs = MS_rot_mtex_rotation(C, get(ebsd,'orientations'));
all_Cs = MS_rotEuler(C, phi1s, thetas, phi2s);

voigt_ave = zeros(6,6) ;
reuss_ave = zeros(6,6) ;
for i = 1:num_xtals
    voigt_ave = voigt_ave + all_Cs(:,:,i) .* (1/num_xtals);
    reuss_ave = reuss_ave + inv(all_Cs(:,:,i)) .* (1/num_xtals);
end

Cave = (voigt_ave + inv(reuss_ave)) ./ 2