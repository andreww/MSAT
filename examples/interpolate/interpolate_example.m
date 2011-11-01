
% Set up values we 'know' 
[C0, roh] = MS_elasticDB('olivine');
C30 = MS_rot3(C0, 0, 0, 30);
C60 = MS_rot3(C0, 0, 0, 60);
C90 = MS_rot3(C0, 0, 0, 90);
C120 = MS_rot3(C0, 0, 0, 120);

% This does not work:
C_interp_90 = MS_VRH([0.8, 0.2], C120, roh, C0, roh);

[Cx90 RR90] = MS_axes(C90)
[Cx120 RR120] = MS_axes(C120)
[Cx0 RR0] = MS_axes(C0)
[Cx30 RR30] = MS_axes(C30)
[Cx60 RR60] = MS_axes(C60)

