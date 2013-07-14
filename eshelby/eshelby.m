function [ esim, escr ] = eshelby (c4, ax1, ax2, ax3)

    % iopt should be 1
    % Assume ... about order of ax - do we need to check?
    
    esim = zeros(3,3,3,3);
    escr = zeros(3,3,3,3);
    
    % Init arrays that don't change between calls - could store these
    [alpha, aa1, aaww1, ww] = init_esh();
    
    
    abc = ax1*ax2*ax3;
    ratio1 = ax2/ax3;
    ratio2 = ax1/ax3;
    
    dte = [0.0 -0.7*ratio1+7 -ratio1+17 -ratio1+23 -ratio1+26 -ratio1+29.3 ...
           -ratio1+32 -ratio1+34.85 -ratio1+37 -ratio1+41.9 -ratio1+44.5];
       
    thecase=11;
    for i = 2:11
        if ((ratio2 >= dte(i-1)) && (ratio2 < dte(i)))
            thecase = i-1;
        end
    end
    
    npoints = 10 * 10; % Hard coded (Gauss-Lebot)
    
    pdil = 0.0;
    p = zeros(3,3);
    gamma2 = zeros(6,6);
    
    % Want C as 6,6 voigt here...
    c2 = MS_cijkl2cij(c4);
    
    for ny = 1:npoints
        
        aa1x = zeros(6);
        
        % A(i,j) = L(i,j,k,l)*a(j)*a(l)
        for i = 1:6
            aa1x(i) = aa1(thecase,i,ny); % Need to build aa1
        end
        
        a1 = esh_mult_voigt(c2,aa1x);
        
    
% c   If solving an elastic inclusion invert the system
% c   --> A(3,3) x X(3,3) = C(3,3)
% c   Inverts A(3,3) using explicit Voigt notation.
% c   Uses explicit form of C(3,3) to calculate solution in Voigt notation.
 
        a1inv = esh_inv3_voigt(a1);
        x1 = zeros(6);
        for i=1:6
            x1(i)=a1inv(i);
        end


        ro3=((alpha(thecase,1,ny)*ax1)^2+ ...
             (alpha(thecase,2,ny)*ax2)^2+ ...
             (alpha(thecase,3,ny)*ax3)^2)^1.5;
        abcoro3=abc/ro3;
    
          
% c   Compute the Eshelby integral Eq.B11 defining:
% c         Gamma(m,j,n,i)=T(m,n,i,j)=a(m)*a(j)*G(n,i)
% c   with the property:
% c         Gamma(m,j,n,i)=Gamma(j,m,n,i)=Gamma(m,j,i,n)

    
        for i = 1:6
            for j = 1:6
                gamma2(i,j) = gamma2(i,j)+aaww1(thecase,i,ny)*x1(j)*abcoro3;
            end
        end
    
    end
    
    % Back to 3*3*3*3 notation... for gamma
    
    gamma4 = MS_cij2cijkl_noMajSym(gamma2);
    
% c   Compute symmetric (distortion) Eshelby tensor from Eq.B9.
% c       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
% c   Compute anti-symmetric (rotation) Eshelby tensor from Eq.B9.
% c       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)

    for l = 1:3
        for k = 1:3
            for m = 1:3
                for n = 1:3
                    dumsim = 0.0;
                    dumscr = 0.0;
                    for j = 1:3
                        for i = 1:3
                            dumsim = dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*c4(i,j,k,l);
                            dumscr = dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*c4(i,j,k,l);
                        end
                    end
                    esim(n,m,k,l)=0.5*dumsim;
                    escr(n,m,k,l)=0.5*dumscr;
                end
            end
        end
    end
end

function [alpha, aa1, aaww1, ww] = init_esh( )

    punti = zeros(10,11);
    pesi = zeros(10,11);
    
    alpha = zeros(11,3,10*10);
    ww = zeros(11,10*10);
    aa1 = zeros(11,6,10*10);
    aaww1 = zeros(11,6,10*10);
    
% c****************************
% c  CASE 1
% c****************************
        punti(1,1)=4.71236594E-02;
        punti(2,1)=0.241774723;
        punti(3,1)=0.565131843;
        punti(4,1)=0.968887568;
        punti(5,1)=1.37937832;
        punti(6,1)=1.76221442;
        punti(7,1)=2.17270517;
        punti(8,1)=2.57646084;
        punti(9,1)=2.89981818;
        punti(10,1)=3.09446883;

        pesi(1,1)=0.120191820;
        pesi(2,1)=0.264987558;
        pesi(3,1)=0.373805553;
        pesi(4,1)=0.420841277;
        pesi(5,1)=0.390970200;
        pesi(6,1)=0.390970260;
        pesi(7,1)=0.420841366;
        pesi(8,1)=0.373805553;
        pesi(9,1)=0.264987499;
        pesi(10,1)=0.120192111;

% c****************************
% c  CASE 2
% c****************************
        punti(1,2)=1.57080423E-02;
        punti(2,2)=0.144995824;
        punti(3,2)=0.425559640;
        punti(4,2)=0.829968274;
        punti(5,2)=1.31460333;
        punti(6,2)=1.82698941;
        punti(7,2)=2.31162453;
        punti(8,2)=2.71603298;
        punti(9,2)=2.99659705;
        punti(10,2)=3.12588477;

        pesi(1,2)=5.41692823E-02;
        pesi(2,2)=0.207461149;
        pesi(3,2)=0.348739326;
        pesi(4,2)=0.452716887;
        pesi(5,2)=0.507709801;
        pesi(6,2)=0.507709682;
        pesi(7,2)=0.452716798;
        pesi(8,2)=0.348738998;
        pesi(9,2)=0.207461327;
        pesi(10,2)=5.41692935E-02;

% c****************************
% c  CASE 3
% c****************************
        punti(1,3)=3.76990959E-02;
        punti(2,3)=0.198626831;
        punti(3,3)=0.483041346;
        punti(4,3)=0.871647120;
        punti(5,3)=1.32964790;
        punti(6,3)=1.81194484;
        punti(7,3)=2.26994562;
        punti(8,3)=2.65855122;
        punti(9,3)=2.94296598;
        punti(10,3)=3.10389376;

        pesi(1,3)=9.68142375E-02;
        pesi(2,3)=0.224478707;
        pesi(3,3)=0.341134071;
        pesi(4,3)=0.430180043;
        pesi(5,3)=0.478189558;
        pesi(6,3)=0.478189170;
        pesi(7,3)=0.430180043;
        pesi(8,3)=0.341134191;
        pesi(9,3)=0.224478647;
        pesi(10,3)=9.68143344E-02;

% c****************************
% c  CASE 4
% c****************************
        punti(1,4)=3.45576368E-02;
        punti(2,4)=0.187556863;
        punti(3,4)=0.468425453;
        punti(4,4)=0.859980166;
        punti(5,4)=1.32527423;
        punti(6,4)=1.81631863;
        punti(7,4)=2.28161263;
        punti(8,4)=2.67316723;
        punti(9,4)=2.95403576;
        punti(10,4)=3.10703516;

        pesi(1,4)=8.95763785E-02;
        pesi(2,4)=0.217725381;
        pesi(3,4)=0.341026783;
        pesi(4,4)=0.435772508;
        pesi(5,4)=0.486694932;
        pesi(6,4)=0.486695170;
        pesi(7,4)=0.435772508;
        pesi(8,4)=0.341026902;
        pesi(9,4)=0.217725128;
        pesi(10,4)=8.95764604E-02;

% c****************************
% c  CASE 5
% c****************************
        punti(1,5)= 3.14158052E-02;
        punti(2,5)=0.177928671;
        punti(3,5)= 0.457155794;
        punti(4,5)= 0.851592362;
        punti(5,5)= 1.32222414;
        punti(6,5)= 1.81936860;
        punti(7,5)=2.29000044;
        punti(8,5)=2.68443704;
        punti(9,5)=2.96366405;
        punti(10,5)=3.11017680;

        pesi(1,5)=8.26927349E-02;
        pesi(2,5)=0.213228315;
        pesi(3,5)=0.342008322;
        pesi(4,5)=0.440196186;
        pesi(5,5)=0.492670894;
        pesi(6,5)=0.492670983;
        pesi(7,5)=0.440195888;
        pesi(8,5)=0.342008322;
        pesi(9,5)=0.213227972;
        pesi(10,5)=8.26930404E-02;

% c****************************
% c  CASE 6
% c****************************
        punti(1,6)= 2.98452154E-02;
        punti(2,6)=0.173592165;
        punti(3,6)=0.452448040;
        punti(4,6)=0.848216832;
        punti(5,6)=1.32101476;
        punti(6,6)=1.82057810;
        punti(7,6)= 2.29337597;
        punti(8,6)=2.68914461;
        punti(9,6)=2.96800065;
        punti(10,6)=3.11174774;

        pesi(1,6)=7.93928578E-02;
        pesi(2,6)=0.211627841;
        pesi(3,6)=0.342669785;
        pesi(4,6)=0.442057431;
        pesi(5,6)=0.495048553;
        pesi(6,6)=0.495048642;
        pesi(7,6)=0.442057490;
        pesi(8,6)=0.342670023;
        pesi(9,6)=0.211627468;
        pesi(10,6)=7.93929026E-02;    

% c****************************
% c  CASE 7
% c****************************
        punti(1,7)=2.67036632E-02;
        punti(2,7)=0.165752888;
        punti(3,7)=0.444431901;
        punti(4,7)=0.842614472;
        punti(5,7)=1.31902647;
        punti(6,7)= 1.82256627;
        punti(7,7)=2.29897833;
        punti(8,7)=2.69716072;
        punti(9,7)=2.97583985;
        punti(10,7)=3.11488938;

        pesi(1,7)=7.30879456E-02;
        pesi(2,7)=0.209402516;
        pesi(3,7)=0.344104946;
        pesi(4,7)=0.445234656;
        pesi(5,7)=0.498966068;
        pesi(6,7)= 0.498966306;
        pesi(7,7)=0.445234746;
        pesi(8,7)= 0.344104946;
        pesi(9,7)=0.209402665;
        pesi(10,7)=7.30878562E-02;

% c****************************
% c  CASE 8
% c****************************
        punti(1,8)=2.67036632E-02;
        punti(2,8)=0.165752888;
        punti(3,8)=0.444431901;
        punti(4,8)=0.842614472;
        punti(5,8)=1.31902647;
        punti(6,8)=1.82256627;
        punti(7,8)=2.29897833;
        punti(8,8)=2.69716072;
        punti(9,8)=2.97583985;
        punti(10,8)=3.11488938;

        pesi(1,8)=7.30879456E-02;
        pesi(2,8)=0.209402516;
        pesi(3,8)=0.344104946;
        pesi(4,8)=0.445234656;
        pesi(5,8)=0.498966068;
        pesi(6,8)=0.498966306;
        pesi(7,8)=0.445234746;
        pesi(8,8)=0.344104946;
        pesi(9,8)=0.209402665;
        pesi(10,8)= 7.30878562E-02;
% c****************************
% c  CASE 9
% c****************************
        punti(1,9)=2.43473575E-02;
        punti(2,9)=0.160516247;
        punti(3,9)=0.439386278;
        punti(4,9)=0.839168847;
        punti(5,9)=1.31781363;
        punti(6,9)=1.82377899;
        punti(7,9)=2.30242372;
        punti(8,9)=2.70220637;
        punti(9,9)=2.98107672;
        punti(10,9)=3.11724544;

        pesi(1,9)=6.86219111E-02;
        pesi(2,9)=0.208388865;
        pesi(3,9)=0.345189095;
        pesi(4,9)=0.447236270;
        pesi(5,9)=0.501360059;
        pesi(6,9)=0.501359940;
        pesi(7,9)=0.447236151;
        pesi(8,9)=0.345189214;
        pesi(9,9)=0.208388969;
        pesi(10,9)=6.86219335E-02;

% c****************************
% c  CASE 10
% c****************************
        punti(1,10)=2.19910536E-02;
        punti(2,10)=0.155757755;
        punti(3,10)=0.434985727;
        punti(4,10)=0.836206555;
        punti(5,10)=1.31677616;
        punti(6,10)= 1.82481658;
        punti(7,10)=2.30538607;
        punti(8,10)=2.70660710;
        punti(9,10)=2.98583508;
        punti(10,10)=3.11960149;

        pesi(1,10)=6.43825606E-02;
        pesi(2,10)=0.207786217;
        pesi(3,10)=0.346235514;
        pesi(4,10)=0.448981822;
        pesi(5,10)=0.503410578;
        pesi(6,10)= 0.503410578;
        pesi(7,10)=0.448981792;
        pesi(8,10)=0.346235693;
        pesi(9,10)=0.207785636;
        pesi(10,10)= 6.43827692E-02;
% c****************************
% c  CASE 11
% c****************************
        punti(1,11)=2.04204638E-02;
        punti(2,11)=0.152822554;
        punti(3,11)=0.432348520;
        punti(4,11)=0.834448099;
        punti(5,11)=1.31616223;
        punti(6,11)=1.82543063;
        punti(7,11)=2.30714464;
        punti(8,11)=2.70924401;
        punti(9,11)=2.98877001;
        punti(10,11)=3.12117243;

        pesi(1,11)=6.16818815E-02;
        pesi(2,11)=0.207559645;
        pesi(3,11)=0.346902698;
        pesi(4,11)=0.450027168;
        pesi(5,11)=0.504624724;
        pesi(6,11)= 0.504624426;
        pesi(7,11)=0.450027317;
        pesi(8,11)=0.346902847;
        pesi(9,11)=0.207559645;
        pesi(10,11)=6.16819337E-02;

        
        for tc = 1:11 % Loop over CASE
            xph(1:10) = punti(1:10,tc);
            xth(1:10) = punti(1:10,tc);
            wph(1:10) = pesi(1:10,tc);
            wth(1:10) = pesi(1:10,tc);
            
            for ith = 1:10 % Point in theta
                sinth=sin(xth(ith));
                costh=cos(xth(ith));
                simbtet=wth(ith)*sinth/(2.0*pi);
                
                for iph = 1:10 % point in phi
                    ny = iph+(ith-1)*10; % point in 2D int
                    ww(tc,ny)=simbtet*wph(iph);
                    alpha(tc,1,ny) = sinth*cos(xph(iph));
                    alpha(tc,2,ny) = sinth*sin(xph(iph));
                    alpha(tc,3,ny) = costh;
                    
                    aa2x = zeros(3,3);
                    aaww2x = zeros(3,3);
                    for i = 1:3
                        for j = 1:3
                            aa2x(i,j) = alpha(tc,i,ny)*alpha(tc,j,ny);
                            aaww2x(i,j) = aa2x(i,j)*ww(tc,ny);
                        end
                    end
                    % Turn aa2x and aaww2x second order tensors into
                    % 6-vectors and store.
                    % Quick symm check...
                    assert(aa2x(2,1)==aa2x(1,2))
                    assert(aa2x(3,1)==aa2x(1,3))
                    assert(aa2x(3,2)==aa2x(2,3))
                    aa1(tc,1,ny) = aa2x(1,1);
                    aa1(tc,2,ny) = aa2x(2,2);
                    aa1(tc,3,ny) = aa2x(3,3);
                    aa1(tc,4,ny) = aa2x(3,2);
                    aa1(tc,5,ny) = aa2x(3,1);
                    aa1(tc,6,ny) = aa2x(2,1);
                    aaww1(tc,1,ny) = aaww2x(1,1);
                    aaww1(tc,2,ny) = aaww2x(2,2);
                    aaww1(tc,3,ny) = aaww2x(3,3);
                    aaww1(tc,4,ny) = aaww2x(3,2);
                    aaww1(tc,5,ny) = aaww2x(3,1);
                    aaww1(tc,6,ny) = aaww2x(2,1);
                    
                    % aww array?
                    
                end 
            end
        end
end

function [ A ] = esh_mult_voigt(B, C)

    A = zeros(6);
    
      A(1)=B(1,1)*C(1)+B(6,6)*C(2)+B(5,5)*C(3) ... 
          +2*(B(5,6)*C(4)+B(1,5)*C(5)+B(1,6)*C(6));

      A(2)=B(6,6)*C(1)+B(2,2)*C(2)+B(4,4)*C(3) ...
          +2*(B(2,4)*C(4)+B(4,6)*C(5)+B(2,6)*C(6));

      A(3)=B(5,5)*C(1)+B(4,4)*C(2)+B(3,3)*C(3)...
          +2*(B(3,4)*C(4)+B(3,5)*C(5)+B(4,5)*C(6));

      A(4)=B(5,6)*C(1)+B(2,4)*C(2)+B(3,4)*C(3) ...
            +(B(2,3)+B(4,4))*C(4) ...
            +(B(3,6)+B(4,5))*C(5) ...
            +(B(4,6)+B(2,5))*C(6);

      A(5)=B(1,5)*C(1)+B(4,6)*C(2)+B(3,5)*C(3) ...
            +(B(3,6)+B(4,5))*C(4) ...
            +(B(1,3)+B(5,5))*C(5) ...
            +(B(1,4)+B(5,6))*C(6);

      A(6)=B(1,6)*C(1)+B(2,6)*C(2)+B(4,5)*C(3) ...
            +(B(4,6)+B(2,5))*C(4) ...
            +(B(1,4)+B(5,6))*C(5) ...
            +(B(1,2)+B(6,6))*C(6);

end

function [ AINV ] = esh_inv3_voigt( A )


      AINV = zeros(6);
      DET = A(1)*A(2)*A(3) + 2*A(4)*A(5)*A(6) - A(1)*A(4)*A(4) ...
          - A(2)*A(5)*A(5) - A(3)*A(6)*A(6);

      AINV(1) = ( A(2)*A(3) - A(4)*A(4))/DET;
      AINV(2) = ( A(1)*A(3) - A(5)*A(5))/DET;
      AINV(3) = ( A(1)*A(2) - A(6)*A(6))/DET;
      AINV(4) = (-A(1)*A(4) + A(5)*A(6))/DET;
      AINV(5) = ( A(4)*A(6) - A(2)*A(5))/DET;
      AINV(6) = (-A(3)*A(6) + A(4)*A(5))/DET;

end