function [PP  ,S1P ,S2P ,SNP ,SNS1,SNS2,VGP ,VGS1,VGS2] = MS_groupvels(C,rh,inc,azi)

      if (length(inc)~=length(azi))
			error('MS:ListsMustMatch', ...
                'AZI and INC must be scalars or vectors of the same dimension');
      end
      
%  ** convert inc, azi to column vectors if necessary      
      inc = reshape(inc,length(inc),1) ;
      azi = reshape(azi,length(azi),1) ;

 
      
      %  ** convert GPa to MB file units (Mbars), density to g/cc

      C(:,:) = C(:,:) * 0.01 ;
      rh = rh ./ 1e3 ;
      
      
      
      
    PP  = zeros(length(azi),3) ;
	S1P = zeros(length(azi),3) ;
    S2P = zeros(length(azi),3) ;
    SNP  = zeros(length(azi),3) ;
	SNS1 = zeros(length(azi),3) ;
    SNS2 = zeros(length(azi),3) ;
    VGP  = zeros(length(azi),3) ;
	VGS1 = zeros(length(azi),3) ;
    VGS2 = zeros(length(azi),3) ;
        

%	** start looping
	for ipair = 1:length(inc)
		cazi = azi(ipair) ;
		cinc = inc(ipair) ;

%  ** create the cartesian vector
		XI = cart2(cinc,cazi) ;

%  ** compute phase velocities		
		[V,EIGVEC]=velo(XI,rh,C) ;
		
%  ** pull out the eigenvectors
        PP(ipair,:)  = EIGVEC(:,1) ;
        S1P(ipair,:) = EIGVEC(:,2) ;
        S2P(ipair,:) = EIGVEC(:,3) ;
        
% ** slowness vectors
        SNP(ipair,:)  = XI/V(1)
        SNS1(ipair,:) = XI/V(2)
        SNS2(ipair,:) = XI/V(3)
        
% ** Group velocity vectors (need to convert density back first)
        VGP(ipair,:) = rayvel(C,SNP(ipair,:),rh*1e3)
        VGS1(ipair,:) = rayvel(C,SNS1(ipair,:),rh*1e3)
        VGS2(ipair,:) = rayvel(C,SNS2(ipair,:),rh*1e3)

        

	end % ipair = 1:length(inc_in)


    
return
%=======================================================================================  

function [c] = V_cross(a, b)
    % This is ~10 times quicker than the Matlab cross() function
    % for me (AMW). We assume the arguments are both 3-vectors and
    % avoid the checks and reshaping needed for the more general case.
    % (According to the prfiler, this moves cross from the most expensive
    % child function costing ~50% of the time to the third most expensive 
    % child costing ~10% of the time).

    c = [a(2).*b(3)-a(3).*b(2)
         a(3).*b(1)-a(1).*b(3)
         a(1).*b(2)-a(2).*b(1)];
    
return


%=======================================================================================
% Rather useing a general case rotation about three angles we use two 
% special case rotations around the c and b axes. This saves time as
% we don't need to convert 0 degrees to radians six times per call to 
% phasevels, or do two double matrix multiplications. This saves ~60%
% of the time spent in vector rotation (and, with the cross thing above)
% reduces the 10000 evaulations needed for MS_anisotropy lma from ~19
% secs to ~10 secs.

function [VR] = V_rot_gam(V,gam)

%  Make rotation matrix
g = gam * pi/180. ;

RR = [ cos(g) sin(g) 0 ; -sin(g) cos(g) 0 ; 0 0 1 ] ;
VR = V * RR ;
 
return

function [VR] = V_rot_bet(V,bet)

%  Make rotation matrix
b = bet * pi/180. ;

RR = [ cos(b) 0 -sin(b) ; 0 1 0 ; sin(b) 0 cos(b) ] ;

VR = V * RR ;
 
return

%=======================================================================================  

%=======================================================================================  
	function [X] = cart2(inc,azm)
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
% Revised, Andrew Walker 2012. 

% Form the 3x3 Christoffel tensor without converting the elasticity 
% from 6x6 to 3x3x3x3 form using the formula from page 1076 of 
% Winterstein 1990. This is > twice as fast as the quickest way I
% have found going via the full tensor form.
        gamma = [X(1) 0.0  0.0  0.0  X(3) X(2) ; ...
                 0.0  X(2) 0.0  X(3) 0.0  X(1) ; ...
                 0.0  0.0  X(3) X(2) X(1) 0.0 ];
        T = gamma * C * gamma';
         
% determine the eigenvalues of symmetric tij
        [EIVEC EIVAL] = eig(T) ;

% calculate velocities and sort
% note that we could get a significant speedup if
% we could avoid the sort - LAPACK usually does sort
% in order of incresing eigenvectors but I've found 
% cases where this does not happen (and we swap P and 
% S wave velocities). Note that MATLAB does not "guarantee
% that the eignevalues are not returned in sorted order"
% http://www.mathworks.com/matlabcentral/newsreader/view_original/394371 
% so we have to sort here...
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






function VG =rayvel(C,SN,rho)
% TO CALCULATE THE RAY-VELOCITY VECTOR CORRESPONDING TO A SLOWNESS VECTOR.
% Original fortran by David Mainprice as part of the EMATRIX code.
% Converted to Python by Alan Baird
%
% C: Stiffness tensor in Voigt Notation (6X6).
% SN: Slowness vector (3).
% rho: Density
%
% returns VG: Group velocity vector (3)

ijkl = [1,6,5; ...
        6,2,4; ...
        5,4,3] ;
        
        
gamma = [X(1) 0.0  0.0  0.0  X(3) X(2) ; ...
         0.0  X(2) 0.0  X(3) 0.0  X(1) ; ...
         0.0  0.0  X(3) X(2) X(1) 0.0 ];
         
F = gamma * C * gamma'-eye(3);
        
        
% Signed cofactors of F[i,k]
CF = zeros(3,3)

CF(1,1)=F(2,2)*F(3,3)-F(2,3)**2
CF(2,2)=F(1,1)*F(3,3)-F(1,3)**2
CF(3,3)=F(1,1)*F(2,2)-F(1,2)**2
CF(1,2)=F(2,3)*F(3,1)-F(2,1)*F(3,3)
CF(2,1)=CF(1,2)
CF(2,3)=F(3,1)*F(1,2)-F(3,2)*F(1,1)
CF(3,2)=CF(2,3)
CF(3,1)=F(1,2)*F(2,3)-F(1,3)*F(2,2)
CF(1,3)=CF(3,1)


% Derivatives of determinant elements
DF = zeros(3,3,3)
for i=1:3
    for j=1:3
        for k=1:3
            DF(i,j,k)=0.0
            for l=1:3
                DF(i,j,k) = DF(i,j,k) + (C(ijkl(i,j),ijkl(k,l))+ C(ijkl(k,j),ijkl(i,l)) ) * SN(l)                
            end
        end
    end
end

% Components of Gradient
DFD = zeros(3)
for k=1:3
    DFD(k) = 0.0
    for i=1:3
        for j=1:3
            DFD(k)=DFD(k)+DF(i,j,k)+CF(i,j)
        end        
    end
end

% Normalize to obtain group velocity
denom = 0.0
VG = zeros(3)
for i=1:3
    denom = denom+SN(i)*DFD(i)
end
for i=1:3
    VG(i) = DFD(i)/denom
end

end % function




