function test_glacial_splitting_C

   % This tests that we get the same elastic constants 
   % matrix as Geoff Lloyd for a given EBSD file and 
   % single crystal tensor. 
   
   % Cij for single crystal ice taken from 'Ice1h-Mb.txt'
   % sent by email, 17 June 2014. 
   % NB - set C(1:4,6) to 0.0 (from 0.0001).
   C_single = [  0.0139  0.0071  0.0058  0.0000  0.0000  0.0000; ...
                 0.0071  0.0139  0.0058  0.0000  0.0000  0.0000; ...
                 0.0058  0.0058  0.0150  0.0000  0.0000  0.0000; ...
                 0.0000  0.0000  0.0000  0.0030  0.0000  0.0000; ...
                 0.0000  0.0000  0.0000  0.0000  0.0030  0.0000; ...
                 0.0000  0.0000  0.0000  0.0000  0.0000  0.0034] ...
                 * 1000.0 ; % convert from Mbar to GPa.
   density = 0.9200 * 1000.0 ; % convert from g/cc to kg/m^3 
   
   % Results - from file 'V3329.txt' sent by email, 17 June 2014.
   C_vrh_OK =     [ .0134    .0067    .0066    .0000    .0000    .0000; ...
                    .0067    .0137    .0064    .0001   -.0001    .0000; ...
                    .0066    .0064    .0137   -.0001    .0001    .0000; ...
                    .0000    .0001   -.0001    .0032    .0000    .0000; ...
                    .0000   -.0001    .0001    .0000    .0037    .0000; ...
                    .0000    .0000    .0000    .0000    .0000    .0033] ...
                    * 1000.0 ; % convert from Mbar to GPa.
   density_OK = 0.9200 * 1000.0 ; % convert from g/cc to kg/m^3

   [C_poly, rho_poly] = Cs_from_EBSD_file(C_single, ...
                                          density, 'V3329g.txt');
                                      
   assertElementsAlmostEqual(density_OK, rho_poly);
   assertElementsAlmostEqual(C_vrh_OK, C_poly, 'absolute', 0.4);
                                      

end


function [C_poly, rho_poly] = Cs_from_EBSD_file(C_single, ...
                                                rho_single, filename)
    % Calculate poly-xtal elasticity of ice 
    % sample from EBSD data and single crystal 
    % elasticity
    
    % Read Eulers from data file
    [eulers, nxtls] = read_EBSD_txt(filename);
          
    % List of elasticities and densities for each crystal 
    Cs = zeros(6,6,nxtls);
    rhos = zeros(1,nxtls);
    for i = 1:nxtls
       Cs(:,:,i) = C_single(:,:);
       rhos(i) = rho_single;
    end

    % Rotate the Cs to give texture
    Cs = MS_rotEuler(Cs, eulers(1,:)', eulers(2,:)', eulers(3,:)', ...
        'sense', 'passive');
  
    % Form VRH mean - each grain has the same volume.
    [C_poly, rho_poly] = MS_VRH(ones(nxtls,1), Cs, rhos);
end


function [eulers, nxtl] = read_EBSD_txt(filename)
    % Given a file name, read the Euler angles 
    % given in the format used by the ice data
    % files. Output is a (3,nxtl) array of 
    % Euler angles (hopefully in Bunge convention
    % and degrees) and the number of crystals, nxtl.

    fid = fopen(filename); % Read - the default
    assert((fid~=-1), 'Could not open file %s', filename);
    fgetl(fid); % Header line - ignore
    data = fscanf(fid, '%f', [12 inf]);
    nxtl = length(data(1,:));
    eulers = zeros(3,nxtl);
    eulers(1,:) = data(3,:); % phi1
    eulers(2,:) = data(4,:); % phi1
    eulers(3,:) = data(5,:); % phi1
    fclose(fid);
end
