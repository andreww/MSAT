% MS_PHASEVELS - Wave velocities in anisotropic media.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the phase velocity details for an elsticity matrix. 
%
%  [ pol, avs, vs1, vs2, vp, ...] = MS_phasevels( C, rh, inc, azi )
%
% Usage: 
%     [ pol, avs, vs1, vs2, vp ] = MS_phasevels( C, rh, inc, azi )                    
%         Calculate phase velocities from elasticity matrix C (in GPa) and
%         density rh (in kg/m^3) for a propogation direction defined by
%         an inclination and azimuth (both in degrees, see below). Output 
%         details are given below.
%
%     [ pol, avs, vs1, vs2, vp, SF, SS ] = MS_phasevels( C, rh, inc, azi )                    
%         Additionally output fast and slow S-wave polarisation in vector
%         form.
%
% Notes:
%     Azi is defined as the angle in degrees from the +ve 1-axis in x1-x2 
%     plane with +ve being clockwise when looking at origin from the
%     3-axis. Inc is defined as the angle in degrees from the x1-x2 plane
%     towards x3 with zero being in the x1-x2 plane. Inc and azi may be
%     scalars, or vectors of the same size. Outputs are:
%
%       'pol' = angle in plane normal to raypath of FSW                           
%              (deg, zero is x3 direction, +ve c'wise looking along             
%              raypath at origin)  
%       'avs' = shear-wave anisotropy
%       'vs1' = fast shear-wave velocity (m/s)
%       'vs2' = slow shear-wave velocity (m/s)
%       'vp'  = P-wave velocity (m/s)
%
%     and all are vectors of length equal to the input inc and azi vectors.
%     In the case of no S-wave splitting (vs1 and vs2 are equal to within
%     eps^1/2) pol is set to NaN. Optional outputs SF and SS are arrays of
%     size (length(inc),3), with each row corresponding to a polarisation
%     vector. This implementation is based on EMATRIX6 by D. Mainprice. 
%     Re-coded in MATLAB by James Wookey.
% 
% Reference: Mainprice D. (1990). An efficient
%            FORTRAN program to calculate seismic anisotropy from
%            the lattice preferred orientation of minerals.
%            Computers & Gesosciences, vol16, pp385-393.

% (C) James Wookey, 2007-2011
% (C) James Wookey and Andrew Walker, 2011

function [pol,avs,vs1,vs2,vp, S1P, S2P] = MS_phasevels(C,rh,inc,azi)

		if (length(inc)~=length(azi))
			error('AZI and INC must be scalars or vectors of the same dimension');
      end
      
%  ** convert inc, azi to column vectors if necessary      
      inc = reshape(inc,length(inc),1) ;
      azi = reshape(azi,length(azi),1) ;

      isotol = sqrt(eps); % Mbars
%  ** convert GPa to MB file units (Mbars), density to g/cc

      C(:,:) = C(:,:) * 0.01 ;
      rh = rh ./ 1e3 ;
      
		avs = zeros(size(azi)) ;
		vp = zeros(size(azi)) ;
		vs1 = zeros(size(azi)) ;
		vs2 = zeros(size(azi)) ;
		pol = zeros(size(azi)) ;
		S1 = zeros(length(azi),3) ;
		S1P = zeros(length(azi),3) ;
        S2P = zeros(length(azi),3) ;
        
%   ** Handle isotropic case quickly
     if isIsotropic(C, isotol)
         vp(:) = sqrt(( ((1.0/3.0)*(C(1,1)+2*C(1,2)))+ ...
                        ((4.0/3.0)*C(4,4)) )/rh)*10.0;
         vs1(:) = sqrt(C(4,4)/rh)*10.0; % Factor of 10 converts from
         vs2 = vs1;                     % Mbar to Pa.
         avs(:) = 0.0;
         pol(:) = NaN; % Both waves have same velocity... meaningless.
         S1P(:) = NaN;
         return
     end

%	** start looping
	for ipair = 1:length(inc)
		cazi = azi(ipair) ;
		cinc = inc(ipair) ;

%  ** create the cartesian vector
		XI = cart2(1,cinc,cazi) ;

%  ** compute phase velocities		
		[V,EIGVEC]=velo(XI,rh,C) ;
		
%  ** pull out the eigenvectors
		P  = EIGVEC(:,1) ;
      S1 = EIGVEC(:,2) ;

  		if ~isreal(S1)
            error_str = ['The S1 polarisation vector is not real!\n\n'...
                sprintf('inc = %f, azi = %f\n\n',cinc,cazi) ...
                sprintf('C = %f %f %f %f %f %f\n',C(1:6,1)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,2)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,3)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,4)) ...
                sprintf('    %f %f %f %f %f %f\n',C(1:6,5)) ...
                sprintf('    %f %f %f %f %f %f\n\n',C(1:6,6)) ...
                sprintf('S1 = %f %f %f\n',S1)];
  			error('MS_PHASEVELS:vectornotreal', error_str) ;
  		end
      S2 = EIGVEC(:,3) ;

%  ** calculate projection onto propagation plane      
      S1N = cross(XI,S1) ;
      S1P(ipair,:) = cross(XI,S1N);
      S2N = cross(XI,S2) ;
      S2P(ipair,:) = cross(XI,S2N);

%  ** rotate into y-z plane to calculate angles
      [S1PR] = V_rot3(S1P,0,0,cazi) ;
	  [S1PRR] = V_rot3(S1PR,0,cinc,0) ;

	   ph = atan2(S1PRR(2),S1PRR(3)) .* 180/pi ;

%  ** transform angle to between -90 and 90
      if (ph < -90.), ph = ph + 180.;,end
      if (ph >  90.), ph = ph - 180.;,end

%	** calculate some useful values
		dVS =  (V(2)-V(3)) ;
      VSmean = (V(2)+V(3))/2.0 ;

      avs(ipair) = 100.0*(dVS/VSmean) ;
      vp(ipair) =  V(1) ;
      vs1(ipair) = V(2) ;
      vs2(ipair) = V(3) ;
		
		pol(ipair) = ph ;
	end % ipair = 1:length(inc_in)

    % If any directions have zero avs (within machine accuracy)
    % set pol to NaN - array wise:
    isiso = real(avs > sqrt(eps)) ; % list of 1.0 and 0.0.
    pol = pol .* (isiso./isiso) ; % times by 1.0 or NaN. 

    S1P(:,1) = S1P(:,1) .* (isiso./isiso);
    S1P(:,2) = S1P(:,2) .* (isiso./isiso);
    S1P(:,3) = S1P(:,3) .* (isiso./isiso);
    
return
%=======================================================================================  

%=======================================================================================  
function [VR] = V_rot3(V,alp,bet,gam)

%  Make rotation matrix
a = alp * pi/180. ;
b = bet * pi/180. ;
g = gam * pi/180. ;

R1 = [ 1 0 0 ; 0 cos(a) sin(a) ; 0 -sin(a) cos(a) ] ;
R2 = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;
R3 = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;

RR =  R3 * R2 * R1;

VR = V * RR ;
 
return
%=======================================================================================  

%=======================================================================================  
	function [X] = cart2(irev,inc,azm)
%=======================================================================================  
%c convert from spherical to cartesian co-ordinates
%c north x=100  west y=010 up z=001
%c irev=+1 positive vector x
%c irev=-1 negative vector x
% NB: pre-converting azm and inc to radians and using
%     cos and sin directly (instead to cosd and sind) 
%     is ~10x faster making the function ~4x faster.
    azmr = azm.*(pi/180.0);
    incr = inc.*(pi/180.0);
    caz=cos(azmr)  ;
    saz=sin(azmr)  ;
    cinc=cos(incr) ;
    sinc=sin(incr) ;
    X=[caz*cinc -saz*cinc sinc] ;
%c normalise to direction cosines
    r=sqrt(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)) ;
   
	X = X./r ;
	if(irev == -1), X = -X;, end
   return
%=======================================================================================  

%=======================================================================================  
	function [V,EIGVEC]=velo(X,rh,C)
%=======================================================================================  
% PHASE-VELOCITY SURFACES IN AN ANISOTROPIC MEDIUM
% revised April 1991
%     X(3) - DIRECTION OF INTEREST
%     RHO - DENSITY
%     V - PHASE VELOCITIES (1,2,3= P,S,SS)
%     EIGVEC(3,3) - eigenvectors stored by columns
%
% Translated to MATLAB by James Wookey         
		ijkl = [1,6,5; ...
		        6,2,4; ...
		        5,4,3] ;
            
        % Form symmetric matrix Tik=Cijkl*Xj*Xl (summation convention)
        % note that this mostly-unrolled approach is ~5x faster than the
        % direct (four-looping) construction and allows us to use the 
        % symmetry of T to reduce the cost. I suspect there is a better
        % way, involving kron(X,X') and a single line for the summation, 
        % but I cannot see it.
        T = zeros(3,3);
        for j=1:3
            T(1,1)=T(1,1) + sum(C(ijkl(1,j),ijkl(1,1:3)).*X(j).*X(1:3)) ;
            T(1,2)=T(1,2) + sum(C(ijkl(1,j),ijkl(2,1:3)).*X(j).*X(1:3)) ;
            T(1,3)=T(1,3) + sum(C(ijkl(1,j),ijkl(3,1:3)).*X(j).*X(1:3)) ;
            T(2,2)=T(2,2) + sum(C(ijkl(2,j),ijkl(2,1:3)).*X(j).*X(1:3)) ;
            T(2,3)=T(2,3) + sum(C(ijkl(2,j),ijkl(3,1:3)).*X(j).*X(1:3)) ;
            T(3,3)=T(3,3) + sum(C(ijkl(3,j),ijkl(3,1:3)).*X(j).*X(1:3)) ;
        end
        % Impose the symmetry (build the lower-left corner).
        T(2,1) = T(1,2);
        T(3,1) = T(1,3);
        T(3,2) = T(2,3);
        
% determine the eigenvalues of symmetric tij
      [EIVEC EIVAL] = eig(T) ;

% calculate velocities and sort
		V_RAW = (sqrt([EIVAL(1,1) EIVAL(2,2) EIVAL(3,3)]./rh))*10. ;
		[V IND] = sort(V_RAW,2,'descend') ;
		EIGVEC = EIVEC ; % for dimensioning
		for i=1:3
			EIGVEC(:,i) = EIVEC(:,IND(i)) ;
		end

      return
%=======================================================================================  

function [ l ] = isIsotropic( C, tol )
    
    % Are we isotropic - assume matrix is symmetrical at this point.
    l = (abs(C(1,1)-C(2,2)) < tol) & (abs(C(1,1)-C(3,3)) < tol) & ...
        (abs(C(1,2)-C(1,3)) < tol) & (abs(C(1,2)-C(2,3)) < tol) & ...
        (abs(C(4,4)-C(5,5)) < tol) & (abs(C(4,4)-C(6,6)) < tol) & ...
        (abs(C(1,4)) < tol) & (abs(C(1,5)) < tol) & (abs(C(1,6)) < tol) & ...
        (abs(C(2,4)) < tol) & (abs(C(2,5)) < tol) & (abs(C(2,6)) < tol) & ...
        (abs(C(3,4)) < tol) & (abs(C(3,5)) < tol) & (abs(C(3,6)) < tol) & ...
        (abs(C(4,5)) < tol) & (abs(C(4,6)) < tol) & (abs(C(5,6)) < tol) & ...
        (((C(1,1)-C(1,2))/2.0)-C(4,4) < tol);
        
return
