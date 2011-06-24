function [ CC ] = MS_rotEuler( C, phi1, theta, phi2 )
%MS_rot_mtex_rotation: rotate elastic matrix using mtex orentation object
%   The MTEX toolbox can represent cryatal orentations using a range of 
%   Euler angle conventions as instances of the orentation class. This 
%   function allows such instances to be used to rotate elatic constnats
%   tensors in MSAT.




    numEs = size(phi1, 1);
    
    % Convert all angles to radians
    phi1r = phi1*(pi/180.0);
    thetar = theta*(pi/180.0);
    phi2r = phi2*(pi/180.0);
    
    numdimCs = ndims(C);
    if ((numEs == 1) && (numdimCs == 2))
        % Two matrix scalar case
        CC = rot_Euler_scalar(C, phi1r, thetar, phi2r);
    elseif ((numEs > 1) && (numdimCs == 2))
        % Many rotations and one matrix case
        CC = zeros(6,6,numEs);
        for i = 1:numEs
            CC(:,:,i) = rot_Euler_scalar(C, phi1r(i), thetar(i), phi2r(i));
        end
    elseif ((numEs == 1) && (numdimCs == 3))
        % Many matrix and one rotation case
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C, 3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r, thetar, phi2r);
        end
    elseif ((numEs > 1) && (numdimCs == 3))
        % List of rotatiosn and matrices case
        assert((numEs == size(C, 3)), 'MS:ListsMustMatch', ...
            'The length of the lists of rotations and matrices must match');
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C,3)
            CC(:,:,i) = rot_Euler_scalar(C(:,:,i), phi1r(i), thetar(i), phi2r(i));
        end
    end
    
    
end


function [ CC ] = rot_Euler_scalar(C, phi1, theta, phi2)
    
    %FIXME - delegate this text to MS_rotR?
    assert(MS_checkC(C)==1, 'MS:Cinvalid', ...
        'Argument "C" must be a valid 6x6 elasticity matrix');
    
    % Extract Bundge convention Euler angles from rotation object
    % and pre-compute trig functions. phi1 etc are in radians.
    %[phi1, theta, phi2] = Euler(o, 'Bundge');
    cp1 = cos(phi1);
    sp1 = sin(phi1);
    cp2 = cos(phi2);
    sp2 = sin(phi2);
    cth = cos(theta);
    sth = sin(theta);

    % Form rotation matrix for Bunge convention of 
    % Euler angles, see eq 24 (pg81) of "Preferred
    % orientation in deformed meals and rocks... H-K Wenk (ed)"
    g = [cp1*cp2 - sp1*sp2*cth, sp1*cp2 + cp1*sp2*cth, sp2*sth ; ...
        -1*cp1*sp2 - sp1*cp2*cth, -1*sp1*sp2 + cp1*cp2*cth, cp2*sth; ...
        sp1*sth, -1*cp1*sth, cth];
    
    
    CC = MS_rotR(C, g);
   
end
