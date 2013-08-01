% MS_effective_medium - generate elastic constants from various effective
%                       medium theories. 
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%  
%  Generate elastic constants for the effective medium parameters based on 
%  various theories, identified by the string 'theory'. Subsequent required 
%  arguments depend on the theory invoked. 
%
%  Currently available theories:
%
%Tandon and Weng, 84 ('tandon' or 't&w') Spheroids
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     (Isotropic host matrix with unidirectionally aligned isotropic spheroid 
%      inclusions) 
%
%
%     [Ceff,rh]=MS_effective_medium('t&w',vpm,vsm,rhm,vpi,vsi,rhi,del,f) or
%        Input parameters:
%           vpm,vsm,rhm : isotropic parameters of the matrix (km/s, kg/m3)
%                   del : aspect ratio of spheroids:
%                            del < 1 = smarties
%                            del = 1 = spheres
%                            del > 1 = cigars
%                    f : volume fractions of inclusions       
%        vpi, vsi, rhi : isotropic parameters of the inclusions
%
%     [Ceff,rh]=MS_effective_medium('t&w',Cm,rhm,Ci,rhi,del,f) or
%        Input parameters:
%             Cm,rh : elasticity and density of the matrix (GPa, kg/m3)
%                   del : aspect ratio of spheroids:
%                            del < 1 = smarties
%                            del = 1 = spheres
%                            del > 1 = cigars
%                    f : volume fractions of inclusions       
%         Ci, rhi : isotropic parameters of the inclusions
%
%       Output parameters:
%          Ceff : Elastic constants (GPa) (symmetry in X1 direction)
%          rh : aggregate density (kg/m3)
%
%
%Hudson, 1980, 1981 ('hudson' or 'crack') Cracks
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%     [Ceff,rh]=MS_effective_medium('hudson',vpm,vsm,rhm,vpc,vsc,rhc,ar,cd) or
%        Input parameters:
%           vpm,vsm,rhm : isotropic parameters of the matrix (km/s, kg/m3)
%                    ar : aspect ratio of cracks
%                    cd : crack density       
%        vpc, vsc, rhc : isotropic parameters of the fill material
%
%     [Ceff,rh]=MS_effective_medium('hudson',Cm,rhm,Ci,rhi,ar,cd) or
%        Input parameters:
%             Cm,rh : elasticity and density of the matrix (GPa, kg/m3)
%                    ar : aspect ratio of cracks
%                    cd : crack density     
%         Cc, rhc : isotropic parameters of the fill material
%
%       Output parameters:
%          Ceff : Elastic constants (GPa) (symmetry in X1 direction)
%          rh : aggregate density (kg/m3)
%
%Backus, 1962 ('backus') stack of thin horizontal layers
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%     [Ceff,rh]=MS_effective_medium('backus',thickness,vp,vs,rh) or
%        Input parameters:
%             thickness : layer thicknesses (vector)
%              vp,vs,rh : isotropic parameters of the layers (km/s, kg/m3) (vector)
%
%       Output parameters:
%          Ceff : Elastic constants (GPa) (symmetry in X1 direction)
%          rh : aggregate density (kg/m3)
%
%     [Ceff,rh]=MS_effective_medium('backus',thickness,C,rh) or
%        Input parameters:
%             thickness : layer thicknesses (vector)
%                     C : elasticity of the layers (GPa) (6x6xn tensor)
%                    rh : density of the layers (GPa, kg/m3) (vector)
%
%       Output parameters:
%          Ceff : Elastic constants (GPa) (symmetry in X3 direction)
%          rh : aggregate density (kg/m3)
%
%References
%~~~~~~~~~~
%
%  Tandon, GP and Weng, GJ. The Effect of Aspect Ratio of Inclusions on the 
%       Elastic Properties of Unidirectionally Aligned Composites. Polymer 
%       Composites, 5, pp 327-333, 1984.
%
%  Hudson, J. A., Overall properties of a cracked solid, Math. Proc. Camb. 
%       Phil. Soc., 88, pp 371-384. 1980
%
%  Hudson, J. A., Wave speeds and attenuation of elastic waves in material 
%       containing cracks, Geophys. J. R. Astr. Soc., 64, pp 133-150. 1981
%   
%  Crampin, S. Effective anisotropic elastic constants for wave propagation 
%       through cracked solids. Geophys. J. R. Astr. Soc., 76, pp 135-145, 1984
%
%  Backus, G. E., Long-Wave Elastic Anisotropy Produced by Horizontal Layering.  
%       J. Geophys. Res., pp 4427-4440
%
% See also: MS_elasticDB

% Copyright (c) 2011-2012, James Wookey and Andrew Walker
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [Ceff,rh]=MS_effective_medium(theory, varargin) ;

if ~ischar(theory)
    error('MS:EFFECTIVE_MEDIUM:BadTString', ...
        'A string specifying the theory to use is required.') ;
end

switch lower(theory)
    
    case {'tandon', 't&w'}
        if length(varargin)~=8 & length(varargin)~=6
            error('MS:EFFECTIVE_MEDIUM:TWWrongArgs', ...
                'Tandon and Weng (1984) requires 6 or 8 input parameters.') ;
        end
        if length(varargin)==6 % elasticity matrix form
            Cm = varargin{1} ; rhm = varargin{2} ;
            Ci = varargin{3} ; rhi = varargin{4} ;
            del = varargin{5} ; f = varargin{6} ;
            %     ** check the matrices
            MS_checkC(Cm) ;
            if MS_anisotropy( Cm) > 10*sqrt(eps) % not isotropic
                error('MS:EFFECTIVE_MEDIUM:TWBadC', ...
                    'Tandon and Weng (1984) requires isotropic inputs matrices.') ;
            end
            MS_checkC(Ci) ;
            if MS_anisotropy( Ci ) > 10*sqrt(eps) % not isotropic
                error('MS:EFFECTIVE_MEDIUM:TWBadC', ...
                    'Tandon and Weng (1984) requires isotropic input matrices.') ;
            end
            %     ** unload the
            vpm = sqrt(Cm(3,3)*1e3./rhm) ;  vsm = sqrt(Cm(6,6)*1e3./rhm) ;
            vpi = sqrt(Ci(3,3)*1e3./rhi) ;  vsi = sqrt(Ci(6,6)*1e3./rhi) ;
        else % velocity form
            vpm = varargin{1} ; vsm = varargin{2} ; rhm = varargin{3} ;
            vpi = varargin{4} ; vsi = varargin{5} ; rhi = varargin{6} ;
            del = varargin{7} ; f = varargin{8} ;
        end
        [Ceff,rh]=MS_tandon_and_weng(vpm,vsm,rhm,vpi,vsi,rhi,del,f) ;
        
        
    case {'hudson', 'crack'}
        if length(varargin)~=8 & length(varargin)~=6
            error('MS:EFFECTIVE_MEDIUM:HuWrongArgs', ...
                'Hudson (1981,1982) cracks require 6 or 8 input parameters.') ;
        end
        if length(varargin)==6 % elasticity matrix form
            Cm = varargin{1} ; rhm = varargin{2} ;
            Cc = varargin{3} ; rhc = varargin{4} ;
            aspr = varargin{5} ; cden = varargin{6} ;
            %     ** check the matrices
            MS_checkC(Cm) ;
            if MS_anisotropy( Cm) > 10*sqrt(eps) % not isotropic
                error('MS:EFFECTIVE_MEDIUM:HuBadC', ...
                    'Hudson (1981,1982) cracks require isotropic inputs matrices.') ;
            end
            MS_checkC(Cc) ;
            if MS_anisotropy( Cc ) > 10*sqrt(eps) % not isotropic
                error('MS:EFFECTIVE_MEDIUM:HuBadC', ...
                    'Hudson (1981,1982) cracks require isotropic input matrices.') ;
            end
            %     ** unload the velocities
            vpm = sqrt(Cm(3,3)*1e3./rhm) ;  vsm = sqrt(Cm(6,6)*1e3./rhm) ;
            vpc = sqrt(Cc(3,3)*1e3./rhc) ;  vsc = sqrt(Cc(6,6)*1e3./rhc) ;
        else % velocity form
            vpm = varargin{1} ; vsm = varargin{2} ; rhm = varargin{3} ;
            vpc = varargin{4} ; vsc = varargin{5} ; rhc = varargin{6} ;
            aspr = varargin{7} ; cden = varargin{8} ;
        end
        [Ceff,rh]=MS_hudson_cracks(vpm,vsm,rhm,vpc,vsc,rhc,aspr,cden) ;
        
        
    case {'backus'}
        if length(varargin)~=3 & length(varargin)~=4
            error('MS:EFFECTIVE_MEDIUM:BaWrongArgs', ...
                'Backus (1962) layering requires 3 or 4 input parameters.') ;
        end
        if length(varargin)==3 % elasticity matrix form
            h = varargin{1} ; C = varargin{2} ; rh = varargin{3} ;
            
            %     ** check the matrices
            [dum dum nC] = size(C) ;
            
            if length(rh)~=nC
                error('MS:EFFECTIVE_MEDIUM:BaVectorLengths', ...
                    'Backus (1962) layering requires equal length vectors.') ;
            end
            
            for iC = 1:nC
                Ctmp = C(:,:,iC) ;
                MS_checkC(Ctmp) ;
                if MS_anisotropy(Ctmp) > 10*sqrt(eps) % not isotropic
                    error('MS:EFFECTIVE_MEDIUM:BaBadC', ...
                        'Backus (1962) layering requires isotropic input matrices.') ;
                end
            end
            
            %     ** unload the velocities
            
            vp = zeros(1,nC) ;
            vs = zeros(1,nC) ;
            vp = sqrt(squeeze(C(3,3,:)).*1e3./rh)  ;
            vs = sqrt(squeeze(C(6,6,:)).*1e3./rh)  ;
        else % velocity form
            h = varargin{1} ; vp = varargin{2} ;
            vs = varargin{3} ; rh = varargin{4} ;
        end
        %  ** TODO: check vector lengths
        
        %  ** reshape so all vectors are the same orientation
        nl = length(h) ;
        h =  reshape(h,1,nl) ;
        vp = reshape(vp,1,nl) ;
        vs = reshape(vs,1,nl) ;
        rh = reshape(rh,1,nl) ;
        
        [Ceff,rh]=MS_backus_average(h,vp,vs,rh) ;
        
        
    otherwise
        error('MS:EFFECTIVE_MEDIUM:UnknownTheory', ...
            'Specified theory is not supported.') ;
        
end % of switch

end


     function [Ceff,rheff] = MS_backus_average(h,vp,vs,rh)
%-----------------------------------------------------------------------
%     Subroutine to perform Backus Averaging of a stack of horizontal
%     isotropic layers to form a homogenous VTI media.
%
%     Inputs:
%     
%     h(1..n)        : Individual thicknesses (not depths!) of the n layers
%     vp(1..n)       : Isotropic P-wave velocities of the n layers (km/s)
%     vs(1..n)       : Isotropic S-wave velocities of the n layers (km/s)
%     rh(1..n)       : Densities of the n layers (kg/m^3)
%
%     Outputs:
%
%     Ceff(6,6)      : Effective elasticity of the package. 
%     rheff          : Density of the package.
%-----------------------------------------------------------------------
   vp = vp.*1e3 ;
   vs = vs.*1e3 ;
   
%
%     Sum over all the layers, weighted by the layer thickness, to generate the
%     Backus parameters.
%
   a1 = sum(h.*(4.*rh.*vs.^2.*(vp.^2-vs.^2))./(vp.^2))./sum(h) ;
   a2 = sum(h.*(1./(rh.*vp.^2)))./sum(h) ;
   a3 = sum(h.*(vp.^2 - 2.*vs.^2)./(vp.^2))./sum(h) ;
   b1 = sum(h.*(2.*rh.*(vp.^2-2.*vs.^2).*vs.^2)./(vp.^2))./sum(h) ;
   b2 = sum(h.*(1./(rh.*vp.^2)))./sum(h) ;
   b3 = sum(h.*(vp.^2 - 2.*vs.^2)./(vp.^2))./sum(h) ;
   c1 = sum(h.*(1./(rh.*vp.^2)))./sum(h) ;
   f1 = sum(h.*(1./(rh.*vp.^2)))./sum(h) ;
   f2 = sum(h.*(vp.^2 - 2*vs.^2)./(vp.^2))./sum(h) ;
   l1 = sum(h.*(1./(rh.*vs.^2)))./sum(h) ;
   m1 = sum(h.*(rh.*vs.^2))./sum(h) ;
   rhave = sum(h.*rh)./sum(h) ;
%  
%  Form Backus coefficients
%   
   a = a1 + 1./a2 .* a3.^2 ;
   b = b1 + 1./b2 .* b3.^2 ;
   c = 1./c1 ;
   f = 1./f1 .* f2 ;
   l = 1./l1 ;
   m = m1 ;
%
%  Form the effective matrix.   
%   
   Ceff = [a b f 0 0 0 ; ...
           b a f 0 0 0 ; ...
           f f c 0 0 0 ; ...
           0 0 0 l 0 0 ; ...
           0 0 0 0 l 0 ; ...
           0 0 0 0 0 m ] ;
   rheff = rhave ;

%  convert to GPa
   Ceff = Ceff ./ 1e9 ;
   
end


function [CC,rh]=MS_tandon_and_weng(vp,vs,rho,vpi,vsi,rhoi,del,c)
% based on original FORTRAN code by Mike Kendall. 

%  weighted average density
   rh = (1.0-c)*rho + c*rhoi ;

   vp = vp * 1e3 ; % convert to m/s 
   vs = vs * 1e3 ; % convert to m/s
   vpi = vpi * 1e3 ; % convert to m/s 
   vsi = vsi * 1e3 ; % convert to m/s

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

% convert to GPa
CC = CC./1e9 ; 

end

function [cc,rhoeff]=MS_hudson_cracks(vp,vs,rho,vpc,vsc,rhoc,aspr,cden)
% based on original FORTRAN code by Mike Kendall. 

%  convert to m/s
   vp = vp * 1e3 ; 
   vs = vs * 1e3 ; 
   vpc = vpc * 1e3 ;
   vsc = vsc * 1e3 ;


%  calculate volume fraction       
   fv = cden*aspr*pi ;
   rhoeff = (1.0-fv)*rho + fv*rhoc ;
 
%  amu=shear modulus (mu), alam=lambda, bm=bulk modulus (kappa)
%  trailing "c" denotes crack property.
   amu = vs*vs*rho ;
   amuc = vsc*vsc*rhoc ;
   alam = vp*vp*rho - 2.0*amu ;
   alamc = vpc*vpc*rhoc - 2.0*amuc ;
   bmc = alamc + amuc*2.0/3.0 ;
 
%  Equation (4) of Crampin (1984, GJRAS)
   term1 = alam + 2.0*amu ;
   term2 = 3.0*alam + 4.0*amu ;
   term3 = alam + amu ;
   term4 = pi*aspr*amu ;
   AK = ((bmc+amuc*4.0/3.0)/term4)*(term1/term3) ;
   AM = (amuc*4.0/term4)*(term1/term2) ;
%
   U11 = (4.0/3.0)*(term1/term3)/(1.0+AK) ;
   U33 = (16.0/3.0)*(term1/term2)/(1.0+AM) ;
 
%  Equation (2) of Crampin (1984, GJRAS)
   term5= -cden/amu ;
   c111 = term1*term1*U11*term5 ;
   c122 = alam*alam*U11*term5 ;
   c133 = alam*alam*U11*term5 ;
   c144 = 0.0 ;
   c155 = amu*amu*U33*term5 ;
   c166 = amu*amu*U33*term5 ;
   c112 = alam*term1*U11*term5 ;
   c113 = alam*term1*U11*term5 ;
   c123 = alam*alam*U11*term5 ;
 
%  Equation (3) of Crampin (1984, GJRAS)
   term6= cden*cden/15.0 ;
   qterm= 15.0*alam*alam/amu/amu + 28.0*alam/amu + 28.0 ;
   Xterm= 2.0*amu*(3.0*alam + 8.0*amu)/term1 ;
   c211 = term1*qterm*U11*U11*term6 ;
   c222 = alam*alam*qterm*U11*U11*term6/term1 ;
   c233 = alam*alam*qterm*U11*U11*term6/term1 ;
   c244 = 0.0 ;
   c255 = Xterm*U33*U33*term6 ;
   c266 = Xterm*U33*U33*term6 ;
   c212 = alam*qterm*U11*U11*term6 ;
   c213 = alam*qterm*U11*U11*term6 ;
   c223 = alam*alam*qterm*U11*U11*term6/term1 ;

%  Host rock Cij
   c11 = vp*vp*rho    ;
   c22 = vp*vp*rho    ;
   c33 = vp*vp*rho    ;
   c44 = vs*vs*rho    ;
   c55 = vs*vs*rho    ;
   c66 = vs*vs*rho    ;
   c12 = c11-2.0*c44  ;
   c13 = c11-2.0*c44  ;
   c23 = c11-2.0*c44  ;

%  Build effective elastic tensor
   cc(1,1)=(c11+c111+c211) ;
   cc(2,2)=(c22+c122+c222) ;
   cc(3,3)=(c33+c133+c233) ;
   cc(4,4)=(c44+c144+c244) ;
   cc(5,5)=(c55+c155+c255) ;
   cc(6,6)=(c66+c166+c266) ;
   cc(1,2)=(c12+c112+c212) ;
   cc(1,3)=(c13+c113+c213) ;
   cc(2,3)=(c23+c123+c223) ;
 
   cc(2,1)=cc(1,2) ;
   cc(3,1)=cc(1,3) ;
   cc(3,2)=cc(2,3) ;
 
%  convert to GPa       
   cc = cc./1e9 ;

end
