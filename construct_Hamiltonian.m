function Ho = construct_Hamiltonian(t, U, N)
% The function construct the hamilotnian
% Input
% t = hopping term 
% U = strength of interaction
% N = Number of sites
% Output
% H = the Hamiltonian for either a spinless fermion or a tight bind.
% E = Eigenvalues / Energies
% Psi = wavefunction
    Ho = 0;
    
%     t = eye(2^N); % the hopping term can be given as a matrix(must be symmetric)
   % t = 1; % the hoppping term can calso be a constant
    for k = 1:N-1
        if U == 0 
            Ho = Ho -t*(creation_op(k,N) * annihilate_op(k+1,N))+ (-t*(creation_op(k+1,N) * annihilate_op(k,N)));
        else
            Ho =Ho -t*(creation_op(k,N) * annihilate_op(k+1,N))+ (-t*(creation_op(k+1,N) * annihilate_op(k,N)))+  ...
             + U*(num_op(k,N)*num_op(k+1,N));
        end
    end
end
    
    