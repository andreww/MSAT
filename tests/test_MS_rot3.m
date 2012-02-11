function test_suite = test_MS_rot3
initTestSuite;
end

function test_MS_rot_3_triclinic

    [C,r] = MS_elasticDB('alb') ;
    C1 = MS_rot3(C,30,45,60) ;
    C2 = MS_rot3(C1,-30,-45,-60,'order',[3 2 1]) ;
    
    M = MS_rotM(30,45,60) ;
    
    C3=MS_rotR(C,M) ;
    
    assertElementsAlmostEqual(C, C2) ;
    assertElementsAlmostEqual(C1, C3) ;

end

function test_MS_rot_3_fuse

    % Check out matrix method works like the explicit sum 
    % ususally given for tensor rotation.
    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
          
    for i=1:1000
       
       % Generate random set of rotations (degs)
       a = rand(1)*360;
       b = rand(1)*360;
       g = rand(1)*360;
       
       % Calculate our rotated matrix
       C1 = MS_rot3(C_in,a,b,g);
       
       % Build the rotation matrix (NB - degrees)
       R = zeros(3,3,3) ;
       R(1,:,:) = [ 1 0 0 ; 0 cosd(a) sind(a) ; 0 -sind(a) cosd(a) ] ;
       R(2,:,:) = [ cosd(b) 0 -sind(b) ; 0 1 0 ; sind(b) 0 cosd(b) ] ;
       R(3,:,:) = [ cosd(g) sind(g) 0 ; -sind(g) cosd(g) 0 ; 0 0 1 ] ;
       orderV = [1 2 3] ;
       RR = squeeze(R(orderV(3),:,:)) * ...
             squeeze(R(orderV(2),:,:)) * ...
                squeeze(R(orderV(1),:,:));
       
       % 'Reference' rotated matrix     
       C2 = MS_rotRC4(C_in,RR) ;
       
       % Check the results match.
       assertElementsAlmostEqual(C1, C2)
            
    end
    
end

function test_MS_rot_3_scalar_vector

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs_in = zeros(6,6,3);
    Cs_in(:,:,1) = C_in;
    Cs_in(:,:,2) = C_in;
    Cs_in(:,:,3) = C_in;
    a = 90.0;
    b = 0.0;
    g = 0.0;
    as = [a ; a; a];
    bs = [b ; b ; b];
    gs = [g ; g ; g];
    
    % Scalar test
    assertElementsAlmostEqual(C_in, MS_rot3(C_in, a, b, g))
    % Vecor tests
    assertElementsAlmostEqual(Cs_in, MS_rot3(Cs_in, a, b, g))
    assertElementsAlmostEqual(Cs_in, MS_rot3(C_in, as, bs, gs))
    assertElementsAlmostEqual(Cs_in, MS_rot3(Cs_in, as, bs, gs))
end

function test_MS_rot_3_errors

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    badC_in = [  49.67  13.18   13.18   000.0   000.0   000.0
                 31.18  49.67   13.18   000.0   000.0   000.0
                 13.18  13.18   49.67   000.0   000.0   000.0
                 000.0  000.0   000.0   12.78   000.0   000.0
                 000.0  000.0   000.0   000.0   12.78   000.0
                 000.0  000.0   000.0   000.0   000.0   12.78];
    a = 90.0;
    b = 0.0;
    g = 0.0;
    
    f = @()MS_rot3(badC_in, a, b, g);
    assertExceptionThrown(f, 'MS:ROT3BadInputMatrix')
        
    f = @()MS_rot3(C_in, a, b, g, 'order',[1 1 1 1]);
    assertExceptionThrown(f, 'MS:ROT3:BadOrder')
    
    f = @()MS_rot3(C_in, a, b, g, 'order',[1 2 4]);
    assertExceptionThrown(f, 'MS:ROT3:BadOrder')
    
    f = @()MS_rot3(C_in, a, b, g, 'bob',[1 2 3]);
    assertExceptionThrown(f, 'MS:ROT3:UnknownOption')
    
   
end

function test_MS_rot_3_vecerrors

    C_in = [  49.67  13.18   13.18   000.0   000.0   000.0
              13.18  49.67   13.18   000.0   000.0   000.0
              13.18  13.18   49.67   000.0   000.0   000.0
              000.0  000.0   000.0   12.78   000.0   000.0
              000.0  000.0   000.0   000.0   12.78   000.0
              000.0  000.0   000.0   000.0   000.0   12.78];
    Cs_in = zeros(6,6,3);
    Cs_in(:,:,1) = C_in;
    Cs_in(:,:,2) = C_in;
    Cs_in(:,:,3) = C_in;
    Cs_in(:,:,4) = C_in;
    a = 90.0;
    b = 0.0;
    g = 0.0;
    as = [a ; a; a];
    bs = [b ; b ; b];
    gs = [g ; g ; g];
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, bs, gs);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, b, gs);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
    % 4 elastic matices and 3 rotation objects...
    f = @()MS_rot3(Cs_in, as, bs, g);
    assertExceptionThrown(f, 'MS:ListsMustMatch')
    
end


% Direct way to rotate (Cij) by a rotation matrix
% 
% [CR] = CIJ_rot3(C,R)
%
%  Inputs: the Cij elastic stiffness tensor, 3x3 rotation matrix
%
%  Output is the rotated Cij elastic stiffness tensor                  
%
function [CR] = MS_rotRC4(C,R)
RR =  R;
CCR = zeros(3,3,3,3);
[CC] = MS_cij2cijkl(C) ;
% rotate the elastic contants
for M=1:3
 for N=1:3
  for R=1:3
   for S=1:3
    CSUM = 0.0 ;
    for I=1:3
     for J=1:3
      for K=1:3
       for L=1:3
        AA = RR(M,I)*RR(N,J)*RR(R,K)*RR(S,L) ;
        CSUM = CSUM + AA * CC(I,J,K,L) ;
       end
      end
     end
    end
    CCR(M,N,R,S) = CSUM ;
   end
  end
 end
end
[CR] = MS_cijkl2cij(CCR) ;
end

