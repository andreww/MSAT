

function [ Cout ] = triinc( Cm, Ci, R, ax1, ax2, ax3, conc)

    
    % Build identity tensor in B-vector notation
    xid6 = zeros(6,6);
    % This trick assumes integer rounding a la Fortran...
    %for i = 1:6
    %    for j = 1:6
    %        xid6(i,j) = (i/j)*(j/i);
    %    end
    %end
    xid6(:,1) = 1.0;
    xid6(2,2) = 1.0;
    xid6(4,2) = 1.0;
    xid6(6,2) = 1.0;
    xid6(3,3) = 1.0;
    xid6(6,3) = 1.0;
    xid6(4,4) = 1.0;
    xid6(5,5) = 1.0;
    xid6(6,6) = 1.0;
    xid6 = xid6'

    
    % Rotate inclusion C matrix to inclusion ref frame
    Ci_r = MS_rotR(Ci, R);
    
    % Calculate Eshelby tensor 
    [esh, ~] = eshelby(MS_cij2cijkl(Ci_r), ax1, ax2, ax3);
    
    % Rotate Eshelby tensor back to sample axis system
    % Because Esh_ijkl != Esh_klij we cannot rotate this like the 
    % Cij matrix. Instead to it explicity.
    esh_sample_axis = zeros(3,3,3,3);
    RR = R';
    for i = 1:3
        for j = 1:3 
            for k = 1:3
                for l = 1:3
                    for ii = 1:3
                        for jj = 1:3
                            for kk = 1:3
                                for ll = 1:3
                                    esh_sample_axis(i,j,k,l) = ...
                                        esh_sample_axis(i,j,k,l) + ...
                                        RR(i,ii)*RR(j,jj)*RR(k,kk)*RR(l,ll)*...
                                        esh(ii,jj,kk,ll);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Ch basis - to 6*6 -- do I need this - YES!
    eshM = chg_basis_4(esh_sample_axis);
    
    % inverse Esh
    eshMi = inv(eshM);
    
    
    % Solve for C in a self consistent manner...
    
    % Ci is inclusion in sample axis and cij notation
    % we need B-notation too
    Ci4 = MS_cij2cijkl(Ci);
    CiM = chg_basis_4(Ci4);
    Cm4 = MS_cij2cijkl(Cm);
    CmM = chg_basis_4(Cm4);
    
    % Cub is Voigt average of Cm and Ci - use the cij notation
    % then put in B-notation
    Cub = zeros(6,6);
    for i = 1:6
        for j = 1:6
            Cub(i,j) = (1-conc)*Cm(i,j) + conc*Ci(i,j);
        end
    end
    CubM = chg_basis_4(MS_cij2cijkl(Cub));
    
    Cold = CubM;    
    mismatch = 1.0;
    
    while (mismatch > 0.00001)
             
        Ctilde = zeros(6,6);
        AC2 = zeros(6,6);
        for i = 1:6
            for j = 1:6
                for k = 1:6 %                   but this should be - ?
                    Ctilde(i,j) = Ctilde(i,j)+Cold(i,k)*(eshMi(k,j)-xid6(k,j));
                end
                AC2(i,j) = Cold(i,j)+Ctilde(i,j);
            end
        end
        MS_cijkl2cij(chg_basis_3(Ctilde))
    
        AC1 = zeros(1,6);
        for i = 1:6
            for j = 1:6
                AC1(i,j) = CiM(i,j)+Ctilde(i,j);
            end
        end
        AC1i = inv(AC1);
    
        AC = zeros(6,6);
        for i = 1:6
            for j = 1:6
                for k = 1:6
                    AC(i,j) = AC(i,j)+AC1i(i,k)*AC2(k,j);
                end
            end
        end
    
        Cnew = zeros(6,6);
        for i = 1:6
            for j = 1:6
                for k = 1:6
                    Cnew(i,j) = Cnew(i,j)+CiM(i,k)*AC(k,j);
                end
                % Voigt av with CmM
                Cnew(i,j) = Cnew(i,j)*conc + CmM(i,j)*(1-conc); 
            end
        end
    
        mismatch = get_mismatch(Cold, Cnew)    
        Cold = Cnew;
        
    end
  
    Cout = MS_cijkl2cij(chg_basis_3(Cnew));
    % Fix up symmetry - should I need to do this?
    %Cout = (0.5.*Cout + 0.5.*Cout');

end

function [ mismatch ] = get_mismatch(C1, C2)

    diff = zeros(6,6);
    av = zeros(6,6);
    diff_norm = 0.0;
    av_norm = 0.0;
    for i = 1:6
        for j = 1:6
            diff(i,j) = C1(i,j) - C2(i,j);
            av(i,j) = 0.5*(C1(i,j)+C2(i,j));
            diff_norm = diff_norm + diff(i,j)*diff(i,j);
            av_norm = av_norm + av(i,j)*av(i,j);
        end
    end
    diff_norm = sqrt(diff_norm);
    av_norm = sqrt(av_norm);
    
    mismatch = diff_norm / av_norm;
end

function [ M2 ] = chg_basis_4(T4) 

   B = get_B();
   
   M2 = zeros(6,6);
   for n = 1:6
       for m = 1:6
           for i = 1:3
               for j = 1:3
                   for k = 1:3
                       for l = 1:3
                           M2(n,m) = M2(n,m)+T4(i,j,k,l)*B(i,j,n)*B(k,l,m);
                       end
                   end
               end
           end
       end
   end
end

function [ T4 ] = chg_basis_3(M2)

    B = get_B();
    
    T4 = zeros(3,3,3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    for n = 1:6
                        for m = 1:6
                            T4(i,j,k,l) = T4(i,j,k,l)+M2(n,m)*B(i,j,n)*B(k,l,m);
                        end
                    end
                end
            end
        end
    end
end

function [ B ] = get_B()

    B = zeros(3,3,6);
    
    B(1,1,2) = -(sqrt(6)^-1);
    B(2,2,2) = -(sqrt(6)^-1);
    B(3,3,2) =  (sqrt(6)^-1)*2.0;
    
    B(1,1,1) = -(sqrt(2)^-1);
    B(2,2,1) =   sqrt(2)^-1;
    
    B(2,3,3) =   sqrt(2)^-1;
    B(3,2,3) =   sqrt(2)^-1;
    
    B(1,3,4) =   sqrt(2)^-1;
    B(3,1,4) =   sqrt(2)^-1;
   
    B(1,2,5) =   sqrt(2)^-1;
    B(2,1,5) =   sqrt(2)^-1;
  
    B(1,1,6) =   sqrt(3)^-1;
    B(2,2,6) =   sqrt(3)^-1;
    B(3,3,6) =   sqrt(3)^-1;
  
end
