%  [CC,rh]=MS_tandon_and_weng(vp,vs,rho,del,c,vpi,vsi,rhoi)
%
%  Calculate anisotropic elastic constants for an isotropic material containing
%  spheroidal inclusions with a symmetry axis aligned with the x1-direction.
%
%  From the theory of Tandon and Weng (1984), based on original FORTRAN code by Mike
%  Kendall.  
%
%  Input parameters:
%       vp,vs,rho : isotropic parameters of the matrix (km/s, kg/m3)
%             del : aspect ratio of spheroids:
%                      del < 1 = smarties
%                      del = 1 = spheres
%                      del > 1 = cigars
%               c : volume fractions of inclusions       
%  vpi, vsi, rhoi : isotropic parameters of the inclusions
%
%  Output parameters:
%     CC : Elastic constants (GPa)
%     rh : aggregate density (kg/m3)
%

 
function [CC,rh]=MS_tandon_and_weng(vp,vs,rho,del,c,vpi,vsi,rhoi) ;

%  weighted average density
   rh = (1.0-c)*rho + c*rhoi ;

   vp = vp * 1e3 ; % convert to m/s 
   vs = vs * 1e3 ; % convert to m/s

   amu = vs*vs*rho ;
   amui = vsi*vsi*rhoi ;
   alam = vp*vp*rho - 2.0*amu ;
   alami = vpi*vpi*rhoi - 2.0*amui ;
   bmi = alami + amui*2.0/3.0 ;
   bmps = alam + amu ;
%  Young's modulus for matrix
   E0 = amu*(3.0*alam + 2.0*amu)/(alam + amu) ;
%  Poisson's ratio of the matrix.
   anu = alam/(2.0*(alam + amu)) ;

%  Some time saving terms
   t1 = del*del - 1.0 ;
   t2 = 1.0 - anu ;
   t3 = 1.0 - 2.0*anu ;
   t4 = 3.0*del*del ;
   t5 = 1.0 - del*del ;
%
% D1, D2 and D3 from Tandon and Weng (1984) (just before equation (18)).
   D1 = 1.0 + 2.0*(amui - amu)/(alami - alam) ;
   D2 = (alam + 2.0*amu)/(alami - alam) ;
   D3 = alam/(alami-alam) ;
%
% g and g' terms (appendix of Tandon and Weng 1984). g is for spheroidal
% inclusions (del>1), whilst g' is for disc-like inclusions (del<1).
%
      if (del >= 1)
       acshdel = log(del + sqrt(t1)) ;
       g =(del*sqrt(t1) - acshdel)*del/(sqrt(t1)^3) ;
      else
%      g' below
       g =(acos(del) - del*sqrt(t5))*del/(sqrt(t5)^3) ;
      end
%
% Eshelby's Sijkl tensor (appendix of Tandon and Weng 1984).
%
       s11 = (t3 + (t4-1.0)/t1 - (t3 + t4/t1)*g)/(2.0*t2)                ;
       s22 = (t4/(t1*2.0) + (t3 - 9.0/(4.0*t1))*g)/(4.0*t2)              ;
       s33 = s22                                                         ;
       s23 = (del*del/(2.0*t1) - (t3 + 3.0/(4.0*t1))*g)/(4.0*t2)         ;
       s32 = s23                                                         ;
       s21 = (-2.0*del*del/t1 + (t4/t1 - t3)*g)/(4.0*t2)                 ;
       s31 = s21                                                         ;
       s12 = (-1.0*(t3 + 1.0/t1) + (t3 + 3.0/(2.0*t1))*g)/(2.0*t2)       ;
       s13 = s12                                                         ;
       s44 = (del*del/(2.0*t1) + (t3 - 3.0/(4.0*t1))*g)/(4.0*t2)         ;
       s66 = (t3 - (t1+2.0)/t1 - (t3 - 3.0*(t1+2.0)/t1)*g/2.0)/(4.0*t2)  ;
       s55 = s66                                                         ;
%
% Tandon and Weng's B terms (after equation 17).
       B1 = c*D1 + D2 + (1.0-c)*(D1*s11 + 2.0*s21)   ; 
       B2 = c + D3 + (1.0-c)*(D1*s12 + s22 + s23)    ; 
       B3 = c + D3 + (1.0-c)*(s11 + (1.0+D1)*s21)    ; 
       B4 = c*D1 + D2 + (1.0-c)*(s12 + D1*s22 + s23) ; 
       B5 = c + D3 + (1.0-c)*(s12 + s22 + D1*s23)    ; 
%
% Tandon and Weng's A terms (after equation 20).
       A1 = D1*(B4 + B5) - 2.0*B2      ;
       A2 = (1.0 + D1)*B2 - (B4 + B5)  ;
       A3 = B1 - D1*B3                 ;
       A4 = (1.0 + D1)*B1 - 2.0*B3     ;
       A5 = (1.0 - D1)/(B4 - B5)       ;
       A = 2.0*B2*B3 - B1*(B4+B5)      ;
%
% Tandon and Weng (1984) equations (25) (28) (31) (32) 
      E11 = E0 /(1.0+c*(A1+2.0*anu*A2)/A) ;
      E22 = E0 ...
      /(1.0+c*(-2.0*anu*A3 + (1.0-anu)*A4 + (1.0+anu)*A5*A)/(2.0*A)) ;
      amu12 = amu*(1.0 + c/(amu/(amui-amu) + 2.0*(1.0-c)*s66)) ;
      amu23 = amu*(1.0 + c/(amu/(amui-amu) + 2.0*(1.0-c)*s44)) ;
%
% Sayers equation (36)
      anu31 = anu - c*(anu*(A1+2.0*anu*A2)+(A3-anu*A4)) ...
                      /(A + c*(A1+2.0*anu*A2)) ;
%
% T&W equation (36)
%     aK12 term; bmps=plane strain bulk modulus
      anum = (1.0+anu)*(1.0-2.0*anu) ;
      denom = 1.0 - anu*(1.0+2.0*anu31) ...
         + c*(2.0*(anu31-anu)*A3 + (1.0-anu*(1.0+2.0*anu31))*A4)/A ;
      aK23 = bmps*anum/denom ;
      anu12tst = E11/E22 - (1.0/amu23 + 1.0/aK23)*E11/4.0 ;

%
% Cij - Sayers' (1992) equations (24)-(29).
% Conversion 
%
      CC(2,2) = amu23 + aK23            ;    
      CC(3,3) = CC(2,2)                 ;    
      CC(1,1) = E11 + 4.0*anu12tst*aK23 ;    
      CC(2,3) = -amu23 + aK23           ;    
      CC(1,2) = 2.0*anu31*aK23          ;    
      CC(1,3) = CC(1,2)                 ;    
      CC(5,5) = amu12                   ;    
      CC(6,6) = CC(5,5)                 ;    
      CC(4,4) = (CC(2,2)-CC(2,3))/2.0   ;    

% Fill out matrix by symmetry
% make symmetrical
for i=1:6
   for j=i:6
      CC(j,i) = CC(i,j) ;
   end
end
   

% convert to GPA
CC = CC./1e9 ; 
