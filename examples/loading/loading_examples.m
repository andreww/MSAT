% LOADING_EXAMPLES.M - simple examples showing how to load data 
%                      into MSAT.
%
% This example shows how to read various file types into MSAT. 
% The function just prints text to the Matlab terminal. All 
% examples are for olivine data.
%
% See also: MS_LOAD, MS_ELASTICDB

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

function loading_examples()

      fprintf('\nReading files with MSAT...\n\n')

      % Load a file in default format assuming file units are GPa
      fprintf('Elasticity matrix from "Olivine.txt" - default format:')
      C = MS_load('Olivine.txt') 

      % Load a file in default format assuming file units are GPa and
      % density in kg.m^-3
      fprintf('Elasticity matrix and density from "Olivine.txt":')
      [C,rh] = MS_load('Olivine.txt')

      % Read density normalised elasticities, where the pre-normalised 
      % elasticities are in Pa and the density is in kg.m^-3 (i.e. units 
      % of the elasticity values in the file are m^2.s^-2). Note that 
      % unit conversion (from Pa to GPa) is performed after
      % denormalisation.
      fprintf('Elasticity matrix and density from "Olivine.Aij"')
      fprintf(' (file data is density normalised):')
      [C,rh] = MS_load('Olivine.Aij','Aij','eunit','Pa')
      
      % Load a file in 'ematrix' format    
      fprintf('Elasticity matrix from "Olivine.Ematrix",')
      fprintf(' file format of the ematrix program:')
      [C] = MS_load('Olivine.Ematrix','format','ematrix')

      % File with isotropic data - note symmetry handling. 
      fprintf('Elasticity matrix and density from "Olivine.iso",')
      fprintf(' this is an isotropic aggregate of olivine crystals:')
      [C,rh] = MS_load('Olivine.iso','symmetry','iso') 
      
      % Use the built in database. 
      fprintf('Elasticity matrix and density from the built in database:')
      [C,rh] = MS_elasticDB('Olivine') 
end
