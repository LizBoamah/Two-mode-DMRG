function effective_H = build_effective_hamiltonian(left_env, local_mpo1, local_mpo2, right_env, s, N)
% Build the effective Hamiltonian by combining the left and right environments and the two-site tensor

if s == 1
    % Contract the right environment with W2 and A2
    temp = tensorprod(right_env, local_mpo2, 2, 4);

    % Combine the left and right parts
    effective_H = tensorprod(local_mpo1, temp, 4, 1);
     effective_H = permute(effective_H, [1, 2, 3, 4, 5, 6, 7, 8]);
    
    siz1 = (size(effective_H, 1) * size(effective_H, 2) * size(effective_H, 3) * size(effective_H, 4));
    siz2 = (size(effective_H, 5) * size(effective_H, 6) * size(effective_H, 7) * size(effective_H, 8));


    effective_H = reshape(effective_H, siz1, siz2);

elseif s == N-1
    % Contract the left environment with W1 and A1
    temp = tensorprod(left_env, local_mpo1, 2, 1);

    % Combine the left and right parts
    effective_H = tensorprod(temp, local_mpo2, 2, 1);
    
    
    siz1 = (size(effective_H, 1) * size(effective_H, 2)* size(effective_H, 3));
    siz2 = (size(effective_H, 4) * size(effective_H, 5)* size(effective_H, 6));

    effective_H = reshape(effective_H, siz1, siz2);

else
    % Contract the left environment with W1 and A1
    temp1 = tensorprod(left_env, local_mpo1, 2, 1);

    % Contract the right environment with W2 and A2
    temp2 = tensorprod(right_env, local_mpo2, 2, 4);

    % Combine the left and right parts
    effective_H = tensorprod(temp1, temp2, 5, 3);
%     effective_H = permute(effective_H, [1, 3, 7, 5, 2, 4, 8, 6]);
     effective_H = permute(effective_H, [1, 2, 3, 4, 5, 6, 7,8]);


%     siz1 = (size(effective_H, 1) * size(effective_H, 3) * size(effective_H, 7) * size(effective_H, 5));
%     siz2 = (size(effective_H, 2) * size(effective_H, 4) * size(effective_H, 8) * size(effective_H, 6));

    siz1 = (size(effective_H, 1) * size(effective_H, 3) * size(effective_H, 4) * size(effective_H, 5));
    siz2 = (size(effective_H, 2) * size(effective_H, 5) * size(effective_H, 8) * size(effective_H, 7));

    effective_H = reshape(effective_H, siz1, siz2);
end
end










































% function effective_H = build_effective_hamiltonian(left_env,local_mpo1, local_mpo2,right_env)
% % Build the effective Hamiltonian by combining the left and right environments and the two-site tensor
% 
% % Contract the left environment with W1 and A1
% temp1 = tensorprod(left_env, local_mpo1, 2,1);
% 
% % Contract the right environment with W2 and A2
% temp2 = tensorprod(right_env, local_mpo2 ,2, 4);
% 
% % Combine the left and right parts
% effective_H = tensorprod(temp2, temp1, 5, 3);
% 
% effective_H = permute(effective_H, [1,3, 7, 5, 2, 4, 8,6]);
% 
% siz1 = (size(effective_H, 1)*size(effective_H, 3)*size(effective_H, 7)*size(effective_H, 5));
% siz2 = (size(effective_H, 2)*size(effective_H, 4)*size(effective_H, 8)* size(effective_H, 6));
% 
% effective_H = reshape(effective_H, siz1, siz2);
% end













































% function effective_H = build_effective_hamiltonian(left_env, local_mpo1, local_mpo2, right_env, s, N)
% % Build the effective Hamiltonian by combining the left and right environments and the two-site tensor
% 
% if s == 1
%     % Contract the right environment with W2 and A2
%     temp = tensorprod(right_env, local_mpo2, 2, 4);
% 
%     % Combine the left and right parts
%     effective_H = tensorprod(temp, local_mpo1, 5, 3);
%     
% elseif s == N-1
%     % Contract the left environment with W1 and A1
%     temp = tensorprod(left_env, local_mpo1, 2, 1);
% 
%     % Combine the left and right parts
%     effective_H = tensorprod(local_mpo2, temp, 1, 3);
%     
% else
%     % Contract the left environment with W1 and A1
%     temp1 = tensorprod(left_env, local_mpo1, 2, 1);
% 
%     % Contract the right environment with W2 and A2
%     temp2 = tensorprod(right_env, local_mpo2, 2, 4);
% 
%     % Combine the left and right parts
%     effective_H = tensorprod(temp2, temp1, 5, 3);
% end
% 
%  effective_H = permute(effective_H, [1,3,7,5,2,4,8,6]);
% % effective_H = permute(effective_H, [1,3,7,5,2,4,8,6]);
% 
% siz1 = (size(effective_H, 1)*size(effective_H, 3)*size(effective_H, 7)*size(effective_H, 6));
% siz2 = (size(effective_H, 2)*size(effective_H, 4)*size(effective_H, 8)*size(effective_H, 5));
% 
% effective_H = reshape(effective_H, siz1, siz2);
% end
