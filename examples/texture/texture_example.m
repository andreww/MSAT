% TEXTURE_EXAMPLE.M - Example script calculating elastic constants of a
%                     textured rock as measured by EBSD.
%
% This script demonstrates how to combine MTEX and MSAT to calculate the
% composite elastic constants of a rock consisiting of two phases where the
% crystals are partially alligned. Based on the MTEX Aachen EBSD example we
% read the EBSD data and assume that each data point represents a equal
% volume of the rock. MTEX is used to convert the EBSD measurment into a
% list of Euler angles. The MSAT function MS_rotEuler is then used to build
% a list of elastic constants tensors aligned with each measured crystal 
% orentation. Finally the tensors are averaged using MS_VRH and the seismic 
% properties of the sample reported. 
%
% For this script to work both MSAT and MTEX must be installed. It is worth
% noting that recent version of MTEX are able to calculate the composite
% elastic constants from an orientation distribution function (ODF) without
% having to sample a finite number of orientations. If the source data is
% from an EBSD run this implies the need to first create an ODF from the
% discrete orientation measurments. 
%
% See also: MS_VRH MS_ROTEULER LOADEBSD EULER 

% (C) James Wookey and Andrew Walker, 2011
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

function texture_example()

    fprintf('\nMSAT TEXTURE EXAMPLE SCRIPT\n\n');

    % Setup basic parameters - single crystal elastic constants and 
    % symmetry. Note that I don't know what the the sample is made of -
    % just that it's a two phase mixture with both phases belonging to the 
    % m-3m point group. I'm going to guess a NaCl / KCl mixture. 
    
    % Elastic constants and density
    [C1, rho1] = MS_elasticDB('NaCl');
    [C2, rho2] = MS_elasticDB('KCl');
    
    % crystal symmetry
    CS = {symmetry('m-3m'),...
          symmetry('m-3m')};
      
    % specimen symmetry
    SS = symmetry('-1');

    
    % Load in all the EBSD data using MTEX. NaCl data ends up in ebsd(1)
    % and KCl data ends up in ebsd(2). 
    tic; fprintf('Loading EBSD data ...');
    
    % specify file name. Use a local copy as the MTEX version seems to move
    % about between versions.
    fname = 'aachen_ebsd_85_829grad_07_09_06.txt';
    
    % create an MTEX EBSD object containing the data
    ebsd = loadEBSD(fname,CS,SS,'interface','generic' ...
        , 'ColumnNames', { 'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' ...
          'Euler3' 'MAD' 'BC' 'BS' 'Bands' 'Error' 'ReliabilityIndex'}, ...
          'Bunge', 'ignorePhase', 0);

    telap = toc; fprintf(' done (%4.2f secs)\n',telap);
    
    % For each EBSD measurment, create an elastic matrix with the measured 
    % orentation. Work on NaCl, then KCl before joining the elastic
    % matricies together
    
    tic; fprintf('Extracting Euler angles ...');
    % NaCl Euler angles
    [Nphi1s, Nthetas, Nphi2s] = Euler(get(ebsd(1),'orientations'), 'Bunge');
    % KCl Euler angles
    [Kphi1s, Kthetas, Kphi2s] = Euler(get(ebsd(2),'orientations'), 'Bunge');
    telap = toc; fprintf(' done (%4.2f secs)\n',telap);
  
    % rads to deg...
    
    Nphi1s=Nphi1s*(180.0/pi); Kphi1s=Kphi1s*(180.0/pi);
    Nthetas=Nthetas*(180.0/pi); Kthetas=Kthetas*(180.0/pi);
    Nphi2s=Nphi2s*(180.0/pi); Kphi2s=Kphi2s*(180.0/pi);
    num_nacl = length(Nphi1s); num_kcl = length(Kphi1s);
    num_xtals = length(Nphi1s) + length(Kphi1s);
    fprintf('%i5 NaCl and %i5 KCl measurments.\n', num_nacl, num_kcl);
    
    % Rotate all elastic constants.
    
    tic; fprintf('Creating list of rotated NaCl elasticity matrices ...');
    NaCl_Cs = MS_rotEuler(C1, Nphi1s, Nthetas, Nphi2s);
    NaCl_rhos = ones(num_xtals,1)*rho1;
    telap = toc; fprintf(' done (%4.2f secs)\n',telap);
    
    tic; fprintf('Creating list of rotated KCl elasticity matrices ...');
    KCl_Cs = MS_rotEuler(C2, Kphi1s, Kthetas, Kphi2s);
    KCl_rhos = ones(num_xtals,1)*rho2;
    telap = toc; fprintf(' done (%4.2f secs)\n',telap);

    
    % Build the input arguments for MS_VRH and calculate the
    % avarage.
    
    tic; fprintf('Calculating VRH average matrix ...');
    Cs = zeros(6,6,num_xtals);
    Cs(:,:,1:num_nacl) = NaCl_Cs(:,:,:);
    Cs(:,:,num_nacl+1:num_xtals) = KCl_Cs(:,:,:);
    vfs = ones(num_xtals,1); % Same volume fraction for each point.
    rhos = [NaCl_rhos KCl_rhos] ;
    
    [Cav, rhoav] = MS_VRH(vfs, Cs, rhos);
    telap = toc; fprintf(' done (%4.2f secs)\n\n',telap);

    
    % Finally, plot the S-wave anisotropy of the rock sample and report 
    % some data on elastic constants.
    plot(ebsd);
    MS_plot(Cav, rhoav);
    MS_info(Cav, rhoav);
    fprintf('\ndone.\n\n');
end