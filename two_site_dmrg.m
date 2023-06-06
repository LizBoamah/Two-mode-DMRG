function [lowest_energy,energy_values, M,E_exact] = two_site_dmrg(N, bd, U, t, max_sweeps, tol)
% TWO_SITE_DMRG - Perform two-site DMRG for the Hubbard model
%
% N - number of sites
% bd - bond dimension
% D - maximum bond dimension
% U - onsite interaction strength
% t - hopping parameter
% max_sweeps - maximum number of sweeps
% tol - convergence tolerance

% energy - ground state energy
% M - MPS representation (cell array of tensors)
 E_left = [];
 E_right = [];
 energy_values = [];
%  energy_diffs = [];
% Initialize a basis states to form a initial wave function
Psi =mps_form_full_basis_new(bd,N,[0 1]',NaN);
%   Psi = mean(Psi, 2);
%  Psi = sum(Psi, 2);
  Psi= Psi(:, 1);
  Psi = Psi/norm(Psi);



% Set direction: right, left, mixed
dir = 'mixed';
 % Initialize the MPS representation
    M = mps_canonical(Psi,bd, N, dir,1);
    % Calculate the complex conjugate of each MPS tensor in the cell array
    M_ = cellfun(@conj, M, 'UniformOutput', false);


% Construct the MPO for the Hubbard Hamiltonian
 Ho = construct_Hamiltonian(t, U, N);
 H = MPOcompress(Ho,bd, N);
 [~,E_exact] = exact_diagonalization(Ho);


  

% Start the sweeping process
for sweep = 1:max_sweeps 
    % Sweep from left to right
    for s = 1:N-1       
       [left_env, right_env] = contract_environments(M, M_, H, s, N);

        % Combine the environments and the two-site tensor to build the effective Hamiltonian
        effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
    
        % Diagonalize the effective Hamiltonian
        [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
         energy_new = (min(diag(eig_val)));
         idx = find(diag(eig_val) == eig_val, 1);
          psi = eig_vec(:, idx);
          
           E_left = [E_left; energy_new];
 % Calculate the new energy after completing a left-to-right sweep
       % Reshape psi and perform SVD for s=1
         if s == 1
            psi_matrix = reshape(psi, [size(M{s}, 2), size(M{s+1}, 2)*size(M{s+1}, 3) ]);
        elseif s == N-1
            % Reshape psi and perform SVD for s=N
            psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
        else
            % Reshape psi and perform SVD for other values of s
            psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
        end
    
        [U, S, V] = svd(psi_matrix, 'econ');

if s == 1
            % Update M{s} and M{s+1}
            M{s} = reshape(U, [1, size(M{s}, 2), size(S, 1)]);
            M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 1), size(M{s+1}, 3)]);
        elseif s == N-1
            % Update M{s-1} and M{s}
            M{s} = reshape(U * S, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
            M{s+1} = reshape(V, [size(S, 1), size(M{s+1}, 2), 1]);
          else
             M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
             M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
end
        % Update the conjugate MPS tensor M_
        M_{s} = conj(M{s});
        M_{s+1} = conj(M{s+1});
   end

for s = N-1:-1:1
       [left_env, right_env] = contract_environments(M, M_, H, s, N);

        % Combine the environments and the two-site tensor to build the effective Hamiltonian
        effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
    
        % Diagonalize the effective Hamiltonian
        [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
        energy_new = (min(diag(eig_val)));
        idx = find(diag(eig_val) == eig_val, 1);
        psi = eig_vec(:, idx);
          E_right = [E_right; energy_new];
    if s == N-1
        % Reshape psi and perform SVD for s=N
        psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
    elseif s == 1
        % Reshape psi and perform SVD for s=1
        psi_matrix = reshape(psi, [size(M{s}, 2), size(M{s+1}, 1) * size(M{s+1}, 3)]);
    else
        % Reshape psi and perform SVD for other values of s
        psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
    end

    [U, S, V] = svd(psi_matrix, 'econ');

   if s == 1
        % Update M{s} and M{s+1}
        M{s} = reshape(U * S, [1, size(M{s}, 2), size(S, 1)]);
        M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
    elseif s == N-1
        % Update M{s-1} and M{s}
        M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
        M{s+1} = reshape(S * V', [size(S, 2), size(M{s+1}, 2), 1]);
    else
        % Update M{s} and M{s+1}
        M{s} = reshape(U * S, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
        M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
    end
  
    % Update the conjugate MPS tensor M_
    M_{s} = conj(M{s});
    M_{s+1} = conj(M{s+1});
end

% Store energy value after completing one full sweep
         E_full_sweep = (E_left(end) + E_right(end)/2);
        energy_values = [energy_values; E_full_sweep];
        lowest_energy = min(energy_values);
        disp(lowest_energy);

% % Check convergence after a full sweep
%     if sweep > 1  % Skip the first sweep as we don't have two energy values yet
%         energy_diff = abs(energy_values(end) - energy_values(end-1));
%         if energy_diff < tol
%             break;  % Stop the sweeping process if the solution has converged
%         end
%     end
% end
% 
% if converged
%     fprintf('DMRG converged after %d sweeps.\n', sweep);
% else
%     fprintf('DMRG did not converge within the maximum number of sweeps.\n');
% end
 
  
end






























































% for s = N-1:-1:1
% %      % Initialize the MPS representation
% %       M = mps_canonical(Psi,bd, N, dir,s);
% %     % Calculate the complex conjugate of each MPS tensor in the cell array
% %     M_ = cellfun(@conj, M, 'UniformOutput', false);
% 
%     % Perform the necessary updates and contractions for right-to-left sweep
%     if s == N-1
%         % Contract environment from the left
%         left_env = 1;
%         for j = 1:s-1
%             left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%         end
%         right_env = 1;
%     elseif s == 1
%         % Contract environment from the right
%         right_env = 1;
%         for j = N:-1:s+2
%             right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%         end
%         left_env = 1;
%     else
%         % Contract environment from the left
%         left_env = 1;
%         for j = 1:s-1
%             left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%         end
%         
%         % Contract environment from the right
%         right_env = 1;
%         for j = N:-1:s+2
%             right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%         end
%     end
% 
%     % Combine the environments and the two-site tensor to build the effective Hamiltonian
%     effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
% 
%     % Diagonalize the effective Hamiltonian
%     [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
%     energy_new = (min(diag(eig_val)));
%     idx = find(diag(eig_val) == eig_val, 1);
%     psi = eig_vec(:, idx);
%     E_right = [E_right; energy_new];
%    
% 
% %   fprintf('Sweep: %d, Sites: (%d, %d), Energy: %f\n', sweep, s, s+1, energy_new);
% %         
%     if s == N-1
%         % Reshape psi and perform SVD for s=N
%         psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     elseif s == 1
%         % Reshape psi and perform SVD for s=1
%         psi_matrix = reshape(psi, [size(M{s}, 2), size(M{s+1}, 1) * size(M{s+1}, 3)]);
%     else
%         % Reshape psi and perform SVD for other values of s
%         psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     end
% 
%     [U, S, V] = svd(psi_matrix, 'econ');
% 
%    if s == 1
%         % Update M{s} and M{s+1}
%         M{s} = reshape(U * S, [1, size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%     elseif s == N-1
%         % Update M{s-1} and M{s}
%         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(S * V', [size(S, 2), size(M{s+1}, 2), 1]);
%     else
%         % Update M{s} and M{s+1}
%         M{s} = reshape(U * S, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%         M{s+1} = reshape(V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%     end
% 
%     % Update the conjugate MPS tensor M_
%     M_{s} = conj(M{s});
%     M_{s+1} = conj(M{s+1});
% end
%   % Store energy value after completing one full sweep
%         E_full_sweep = (E_left(end) + E_right(end));
%         energy_values = [energy_values; E_full_sweep];
% 
% 
% end






































%         if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [1, size(M{s}, 2), size(S, 1)]);
%             M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%     %         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     %         M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%              M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%         end

%  % Perform the necessary updates and contractions for left-to-right sweep
%         if s == 1
%             % Contract environment from the right
%             right_env = 1;
%             for j = N:-1:s+2
%                 right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%             end
%             left_env = 1;
%         elseif s == N-1
%             % Contract environment from the left
%             left_env = 1;
%             for j = 1:s-1
%                 left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%             end
%             right_env = 1;
%         else
%             % Contract environment from the left
%             left_env = 1;
%             for j = 1:s-1
%                 left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%             end
%             
%             % Contract environment from the right
%             right_env = 1;
%             for j = N:-1:s+2
%                 right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%             end
%        end
% 
% 

% if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             end
% 
% 
%         % Update the conjugate MPS tensor M_
%         M_{s} = conj(M{s});
%         M_{s+1} = conj(M{s+1});




%         if s == 1
%             % Initialize the MPS representation
%             M = mps_canonical(Psi,bd, N, dir, s);
%             % Calculate the complex conjugate of each MPS tensor in the cell array
%             M_ = cellfun(@conj, M, 'UniformOutput', false);
%         
%         end
%        [left_env, right_env] = contract_environments(M, M_, H, s, N);
% 
%         % Combine the environments and the two-site tensor to build the effective Hamiltonian
%         effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env, s, N);
%     
%         % Diagonalize the effective Hamiltonian
%         [eig_vec, eig_val] = eig((effective_H + effective_H')/2);
%         energy_new = (min(diag(eig_val)));
%         idx = find(diag(eig_val) == eig_val, 1);
%         psi = eig_vec(:, idx);
%           
%          E_left = [E_left; energy_new];
%          
% 
%             % Calculate the new energy after completing a left-to-right sweep
%        
%       % Reshape psi and perform SVD for s=1
%          if s == 1
%             psi_matrix = reshape(psi, [size(M{s}, 2)* size(M{s},1), size(M{s+1}, 1) * size(M{s+1}, 2)]);
%         elseif s == N-1
%             % Reshape psi and perform SVD for s=N
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         else
%             % Reshape psi and perform SVD for other values of s
%             psi_matrix = reshape(psi, [size(M{s}, 2) * size(M{s}, 3), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         end
%     
%           [U, S, V] = svd(psi_matrix, 'econ');
% 
%        if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [ size(M{s}, 2), size(M{s}, 1), size(S, 1)]);
%             M{s+1} = reshape(S*V', [ size(M{s+1}, 3),size(M{s+1}, 1), size(M{s+1}, 2) ]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%              M{s} = reshape(U, [size(M{s}, 2), size(M{s}, 3), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 2)]);
%        end
% 
%         % Update the conjugate MPS tensor M_
%         M_{s} = conj(M{s});
%         M_{s+1} = conj(M{s+1});
% 
%        
% 
%       
%   
% 
% 
%      end
%     








%          if s == 1
%             psi_matrix = reshape(psi, [size(M{s}, 2)*size(M{s},3), size(M{s+1}, 1) * size(M{s+1}, 2)]);
%         elseif s == N-1
%             % Reshape psi and perform SVD for s=N
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         else
%             % Reshape psi and perform SVD for other values of s
%             psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%         end
%     



%        if s == 1
%             % Update M{s} and M{s+1}
%             M{s} = reshape(U, [size(S, 1),size(M{s}, 2), size(M{s}, 3)]);
%             M{s+1} = reshape(S*V', [size(M{s+1}, 1), size(M{s+1}, 2),size(S, 2)]);
%         elseif s == N-1
%             % Update M{s-1} and M{s}
%             M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), size(S, 2)]);
%             M{s} = reshape(V, [size(S, 1), size(M{s}, 2), 1]);
%         else
%             % Update M{s} and M{s+1}
%     %         M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     %         M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%              M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 1)]);
%              M{s+1} = reshape(S*V', [size(S, 2), size(M{s+1}, 2), size(M{s+1}, 3)]);
%        end




%         % Truncate singular values
%             truncation_idx = min(size(S, 1), D);
%             S = S(1:truncation_idx, 1:truncation_idx);
%             U = U(:, 1:truncation_idx);
%             V = V(:, 1:truncation_idx);
% 
%           if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%           end













% function[energy_new, M] = two_site_dmrg(N, bd, D, U, s, t, max_sweeps, tol)
% % TWO_SITE_DMRG - Perform two-site DMRG for the Hubbard model
% % 
% % N - number of sites
% % bd - bond dimension
% % D - maximum bond dimension
% % U - onsite interaction strength
% % t - hopping parameter
% % max_sweeps - maximum number of sweeps
% % tol - convergence tolerance
% 
% %energy - ground state energy
% %M - MPS representation (cell array of tensors)
% 
% %Construct the Hubbard Hamiltonian 
% % Ho = construct_Hamiltonian(t, U, N);
% 
% %Initialize a MPS 
% mps = init_random_mps(N, bd, D);
% 
% %Set direction: right, left, mixed
% dir = 'mixed';
% 
% %Construct the MPO for the Hubbard Hamiltonian
% % H = MPOcompress(Ho, bd, N);
%  H= hubbard_mpo_site(U, t, N, bd, D);
% %Initialize the MPS representation
% M = mps_canonicalM(mps, N,dir, s);
% 
% %Calculate the complex conjugate of each MPS tensor in the cell array
% 
%  M_ = cellfun(@conj, M,'UniformOutput', false);
% %Initialize variables for sweeping
%   energy_old = inf;
% 
% %Start the sweeping process
% for sweep = 1:max_sweeps
%     %Contract environment from the left
%     left_env = 1;
%     for j = 1:s-1
%         left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
%     end
% 
%     %Contract environment from the right
%     right_env = 1;
%     for j = N:-1:s+2
%         right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
%     end
% 
%     %Combine the environments and the two-site tensor to build the effective Hamiltonian
%     effective_H = build_effective_hamiltonian(left_env, H{s}, H{s+1}, right_env);
% 
%    % Diagonalize the effective Hamiltonian
%     [eig_vec, eig_val] = eig(effective_H);
%     energy_new = real(min(diag(eig_val)));
%     idx = find(diag(eig_val) == energy_new, 1);
%     psi = eig_vec(:, idx);
% 
%    % Reshape psi and perform SVD
%     psi_matrix = reshape(psi, [size(M{s}, 1) * size(M{s}, 2), size(M{s+1}, 2) * size(M{s+1}, 3)]);
%     [U, S, V] = svd(psi_matrix, 'econ');
% 
% %     %Truncate the bond dimension
% %     if size(S, 2) > D
% %         S = S(:, 1:D);
% %         U = U(:, 1:D);
% %         V = V(1:D, :);
% %     end
%     % Print the energy at the current pair of sites
%     fprintf('Sweep: %d, Sites: (%d, %d), Energy: %f\n', sweep, s, s+1, energy_new);
% 
% 
%     %Update MPS tensors
%     M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), size(S, 2)]);
%     M{s+1} = reshape(S*V, [size(S, 1), size(M{s+1}, 2), size(M{s+1}, 3)]);
% 
%  % Check for convergence
%     energy_diff = abs(energy_new - energy_old);
%     if energy_diff < tol
%         
%     end
%     energy_old = energy_new;
% end
% 
%    
% end
% 
% 



% % Check convergence
%     energy_diff = abs(energy_new - energy_old);
%     if energy_diff < tol
% %         converged = true;
%         break;
%     else
%         energy_old = energy_new;
%     end
% end
% energies = [energies, energy];
% 



%         % Truncate singular values
%             truncation_idx = min(size(S, 1), D);
%             S = S(1:truncation_idx, 1:truncation_idx);
%             U = U(:, 1:truncation_idx);
%             V = V(:, 1:truncation_idx);
% 
%           if s == 1
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [1, size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%             elseif s == N-1
%                 % Update M{s-1} and M{s}
%                 M{s-1} = reshape(U * S, [size(M{s-1}, 1), size(M{s-1}, 2), truncation_idx]);
%                 M{s} = reshape(V, [truncation_idx, size(M{s}, 2), 1]);
%             else
%                 % Update M{s} and M{s+1}
%                 M{s} = reshape(U, [size(M{s}, 1), size(M{s}, 2), truncation_idx]);
%                 M{s+1} = reshape(S*V', [truncation_idx, size(M{s+1}, 2), size(M{s+1}, 3)]);
%           end

































