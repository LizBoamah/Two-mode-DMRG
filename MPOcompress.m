% function [H] = MPOcompress( d, N)
% t=1;
% U=1;
% Ho = construct_Hamiltonian(t, U, N);
% 
%     ml = 1;
%    
% 
%     % Initialize MPO - tensor at site n
%     H = cell(1, N);
%     
% %     switch dir
% %         case 'left'
%             % Beginning to N-1
%             for l = 1:N-1
%                 W = reshape(Ho, [ml * (2 ^ d), d ^ (2 * (N - 1))]);
%                 [U, Lambda, V] = svd(W, 'econ');
% 
%                 % Calculate the bond dimension after truncation
%                 bd = size(U, 2);
% 
%                 % Update the MPO tensors
%                 H{l} = reshape(U, [ml,d, d, bd]);
% 
%                 % Update the Hamiltonian for the next iteration
%                 Ho = Lambda * V';
%                 ml = bd;
%             end
% 
%             % Last mode
%             H{N} = reshape(Ho, [ml, d ,d ,1]);

   % end

function H = MPOcompress(Ho,d, N)

    ml = 1;

    % Initialize MPO - tensor at site n
    H = cell(1, N);
    %--- beginning to N-1
    for l = 1:N - 1
        W = reshape(Ho, [ml * d^d, d^(2*(N - l))]); 
        [U, S, V] = svd(W, 'econ');
        new_ml = size(U, 2);

        U = reshape(U, [ml, d, d, new_ml]);
        V = S * V';

        H{l} = U;

        Ho = V;
        ml = new_ml;
    end

    %---last mode
    H{N} = reshape(Ho, [new_ml, d, d, 1]);
    

    
end
