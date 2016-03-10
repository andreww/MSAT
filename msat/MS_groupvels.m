% MS_PHASEVELS - Wave velocities in anisotropic media.
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Calculate the group velocity details for an elsticity matrix. 
%
%  [ VGP, VGS1, VGS2, ...] = MS_groupvels( C, rh, inc, azi )
%
% Usage: 
%     [ VGP, VGS1, VGS2, ...] = MS_groupvels( C, rh, inc, azi )                    
%         Calculate group velocity vectors from elasticity matrix C (in GPa) and
%         density rh (in kg/m^3) corresponding to a phase angle defined by
%         an inclination and azimuth (both in degrees). Output 
%         details are given below.
%
%     [ VGP, VGS1, VGS2, PE, S1E, S2E ] = MS_groupvels( C, rh, inc, azi )                 
%         Additionally output P, S1 and S2-wave polarisations in vector
%         form.
%
%     [ VGP, VGS1, VGS2, PE, S1E, S2E, SNP, SNS1, SNS2 ] = MS_groupvels( C, rh, inc, azi )                 
%         Additionally output P, S1 and S2 Slowness vectors
%

% Copyright (c) 2016, Alan Baird
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

function [ varargout ] = MS_groupvels(C,rh,inc,azi)
    
    [~,~,vs1,vs2,vp,~,~,PE,S1E,S2E,XIS] = MS_phasevels(C,rh,inc,azi);
    
    SNP  = zeros(length(azi),3) ;
	SNS1 = zeros(length(azi),3) ;
    SNS2 = zeros(length(azi),3) ;
    VGP  = zeros(length(azi),3) ;
	VGS1 = zeros(length(azi),3) ;
    VGS2 = zeros(length(azi),3) ;
    

        

%	** start looping
	for ipair = 1:length(inc)

% ** slowness vectors
        SNP(ipair,:)  = XIS(ipair,:)./vp(ipair);
        SNS1(ipair,:) = XIS(ipair,:)./vs1(ipair);
        SNS2(ipair,:) = XIS(ipair,:)./vs2(ipair);
        
% ** Group velocity vectors

        VGP(ipair,:) = rayvel(C,SNP(ipair,:),rh./1e3);
        VGS1(ipair,:) = rayvel(C,SNS1(ipair,:),rh./1e3);
        VGS2(ipair,:) = rayvel(C,SNS2(ipair,:),rh./1e3);

	end 
    
    switch nargout
    case 3
        varargout{1} = VGP ;
        varargout{2} = VGS1 ;
        varargout{3} = VGS2 ;


    case 6
        varargout{1} = VGP ;
        varargout{2} = VGS1 ;
        varargout{3} = VGS2 ;
        varargout{4} = PE ;
        varargout{5} = S1E ;
        varargout{6} = S2E ;

    
    case 9
        varargout{1} = VGP ;
        varargout{2} = VGS1 ;
        varargout{3} = VGS2 ;
        varargout{4} = PE ;
        varargout{5} = S1E ;
        varargout{6} = S2E ;
        varargout{7} = SNP ;
        varargout{8} = SNS1 ;
        varargout{9} = SNS2 ;
       
    otherwise   
       error('MS:GROUPVELS:BadOutputArgs','Requires 3, 6, or 9 output arguments.')
    end


    
return
%=======================================================================================  






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
        
        
gamma = [SN(1) 0.0  0.0  0.0  SN(3) SN(2) ; ...
         0.0  SN(2) 0.0  SN(3) 0.0  SN(1) ; ...
         0.0  0.0  SN(3) SN(2) SN(1) 0.0 ];
         
F = gamma * C * gamma'-eye(3).*rho;
        
        
% Signed cofactors of F[i,k]
CF = zeros(3,3);

CF(1,1)=F(2,2)*F(3,3)-F(2,3)^2;
CF(2,2)=F(1,1)*F(3,3)-F(1,3)^2;
CF(3,3)=F(1,1)*F(2,2)-F(1,2)^2;
CF(1,2)=F(2,3)*F(3,1)-F(2,1)*F(3,3);
CF(2,1)=CF(1,2);
CF(2,3)=F(3,1)*F(1,2)-F(3,2)*F(1,1);
CF(3,2)=CF(2,3);
CF(3,1)=F(1,2)*F(2,3)-F(1,3)*F(2,2);
CF(1,3)=CF(3,1);


% Derivatives of determinant elements
DF = zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            DF(i,j,k)=0.0;
            for l=1:3
                DF(i,j,k) = DF(i,j,k) + (C(ijkl(i,j),ijkl(k,l))+ C(ijkl(k,j),ijkl(i,l)) ) * SN(l);                
            end
        end
    end
end

% Components of Gradient
DFD = zeros(3,1);
for k=1:3
    DFD(k) = 0.0;
    for i=1:3
        for j=1:3
            DFD(k)=DFD(k)+DF(i,j,k)*CF(i,j);
        end        
    end
end

% Normalize to obtain group velocity
denom = 0.0;
VG = zeros(3,1);
for i=1:3
    denom = denom+SN(i)*DFD(i);
end
for i=1:3
    VG(i) = DFD(i)./denom;
end

return % function




