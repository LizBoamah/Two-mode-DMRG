clear
clc
% Example script to run the two_site_dmrg function

% Parameters for the Hubbard model
N = 6;             % Number of sites
bd = 2;             % Local space dimension (spin up and spin down)
%D = 4;            % Bond dimension
U = 1;             % Onsite interaction strength
t = 1;             % Hopping parameter
max_sweeps = 20;   % Maximum number of sweeps
 tol = 1e-6;        % Convergence tolerance


% Run the two_site_dmrg function
[lowest_energy, energy_values, M, E_exact] = two_site_dmrg(N, bd, U,  t, max_sweeps, tol);

% % % Display the final energy
%  fprintf('Ground state energy: %.6f\n', energy);

