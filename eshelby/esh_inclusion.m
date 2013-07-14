function [ Cout, Sout ] = esh_inclusion(Cmatrix, Cinclusion, a1, a2, a3, c)

    % Turn elasticity matricies into 3*3*3*3 tensors

    Cm = MS_cij2cijkl(Cmatrix);
    Sm = MS_cij2cijkl(inv(Cmatrix));
    Ci = MS_cij2cijkl(Cinclusion);
    Si = MS_cij2cijkl(inv(Cinclusion));
    C = zeros(3,3,3,3);
    % NB: should rotate Cm to elpsode axis system and then S back to 
    % sample reference frame - leave this for now (axis system of inclusion
    % parallel to that of sample.
    
    % Calculate Eshelby tensor, S.
    [ S, ~ ] = eshelby(Cm, a1, a2, a3);
    
    % Prep Sout for output
    Sout = MS_cijkl2cij(S);
    
    % Loop over stresses
    
    
    for i = 1:3
        for j = 1:3
    %i = 1
    %j = 1
            sigma = zeros(3,3);
            sigma(i,j) =  1.000;
            sigma(j,i) =  1.000;
            % Calculate reference strain from eq.2
            e0 = zeros(3,3);
            for k = 1:3
                for l = 1:3
                    %e0(k,l) = Sm(i,j,k,l) .* sigma(i,j);
                    if (Cm(i,j,k,l) ~= 0.0)
                        e0(k,l) = sigma(i,j)/Cm(i,j,k,l);
                    end
                end
            end
            
            % Init guess, estar = e0.
            %es = e0;
            es = zeros(3,3);
            es_v(1:6) = [es(1,1) es(2,2) es(3,3) ...
                         es(2,3) es(1,3) es(1,2)];
            f = @(es_v)trans_stress(es_v, Cm, Ci, e0, S, c);
            [es_v_o, squ] = fminsearch(f, es_v);
            
            es = [es_v_o(1) es_v_o(6) es_v_o(5);...
                  es_v_o(6) es_v_o(2) es_v_o(4);...
                  es_v_o(5) es_v_o(4) es_v_o(3)];
               
           
            
            % Calculate ebar (e.q. 11)
            
            ebar = e0 + c.*es ;
            
            for k = 1:3
                for l = 1:3
                    if (ebar(k,l) > 0.00001)
                        C(i,j,k,l) = sigma(i,j)/ebar(k,l);
                    else
                        C(i,j,k,l) = 0.0;
                    end
                end
            end
        end
    end
   
    Cout = MS_cijkl2cij(C) ;
   % Cout = (Cout + Cout')/2.0
   % Cout = inv(Cout)
    
end

function [ sum_sq_stress ] = trans_stress(es_v, Cm, Ci, e0, S, c)

    estar = [es_v(1) es_v(6) es_v(5);...
             es_v(6) es_v(2) es_v(4);...
             es_v(5) es_v(4) es_v(3)];

    stress = zeros(3,3);
    sum_sq_stress = 0.0;
    
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    estarstar = zeros(3,3);
                    for m = 1:3
                        for n = 1:3
                            estarstar(kk,ll) = estarstar(kk,ll) + ...
                                S(kk,ll,m,n)*estar(m,n);
                        end
                    end
                    stress(ii,jj) = stress(ii,jj) + (Ci(ii,jj,kk,ll)-Cm(ii,jj,kk,ll))*...
                       (e0(kk,ll)+(1-c)*estarstar(kk,ll)+c*estar(kk,ll)) + ...
                       Cm(ii,jj,kk,ll)*estar(kk,ll);
                end
            end
            
            sum_sq_stress = sum_sq_stress + (stress(ii,jj)^2);

        end
    end
    

end
