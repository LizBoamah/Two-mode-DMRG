function [Psi_exact, E_exact] = exact_diagonalization(Ho)
% Perform exact diagonalization on the Hamiltonian matrix
% Input:
%   H: Hamiltonian matrix
% Output:
%   E: Eigenvalues of the Hamiltonian matrix
%   psi: Eigenvectors of the Hamiltonian matrix

[Psi_exact,E_exact] = eig((Ho +Ho')/2);
E_exact = diag(E_exact);
% Find the lowest energy (minimum eigenvalue)
E_exact = min(E_exact);
end