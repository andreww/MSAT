function [ CC ] = MS_rot_mtex_rotation( C, o )
%MS_rot_mtex_rotation: rotate elastic matrix using mtex orentation object
%   The MTEX toolbox can represent cryatal orentations using a range of 
%   Euler angle conventions as instances of the orentation class. This 
%   function allows such instances to be used to rotate elatic constnats
%   tensors in MSAT.


    try
       mtex_path;
    catch e
      if strcmp(e.identifier, 'MATLAB:UndefinedFunction')
          error('MS:NoMtex', ...
            'MTEX must be installed in order to use MS_rot_mtex_rotation');
      else
          error('MS:internaError', ...
            'Internal error in MS_rot_mtex_rotation');
      end
    end

    numOs = size(o, 1);
    numdimCs = ndims(C);
    if ((numOs == 1) && (numdimCs == 2))
        % Two matrix scalar case
        CC = rot_mtex_rotation_scalar(C, o);
    elseif ((numOs > 1) && (numdimCs == 2))
        % Many rotations and one matrix case
        CC = zeros(6,6,numOs);
        for i = 1:numOs
            CC(:,:,i) = rot_mtex_rotation_scalar(C, o(i));
        end
    elseif ((numOs == 1) && (numdimCs == 3))
        % Many matrix and one rotation case
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C, 3)
            CC(:,:,i) = rot_mtex_rotation_scalar(C(:,:,i), o);
        end
    elseif ((numOs > 1) && (numdimCs == 3))
        % List of rotatiosn and matrices case
        assert((numOs == size(C, 3)), 'MS:ListsMustMatch', ...
            'The length of the lists of rotations and matrices must match');
        CC = zeros(6,6,size(C, 3));
        for i = 1:size(C,3)
            CC(:,:,i) = rot_mtex_rotation_scalar(C(:,:,i), o(i));
        end
    end
    
    
end


function [ CC ] = rot_mtex_rotation_scalar(C, o)
    assert(isa(o, 'rotation'), 'MS:notarot', ...
        'Argument "o" must be a MTEX rotation object');
    
    %FIXME - delegate this text to MS_rotR?
    assert(MS_checkC(C)==1, 'MS:Cinvalid', ...
        'Argument "C" must be a valid 6x6 elasticity matrix');
    
    % Extract Bundge convention Euler angles from rotation object
    % and pre-compute trig functions. phi1 etc are in radians.
    [phi1, theta, phi2] = Euler(o, 'Bundge');
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

% Things to add...
% Make this work (JW's rot func)
% check g is a rotation matrix
% make it work for lists of g's an one matrix
% make it work for lists of C's and one g
% make it work for equal length lists?
% example for VRH av from VPSC output.

