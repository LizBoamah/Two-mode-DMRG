function [left_env, right_env] = contract_environments(M, M_, H, s, N)
    if s == 1
        % Contract environment from the right
        right_env = 1;
        for j = N:-1:s+2
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
        left_env = 1;
    elseif s == N-1
        % Contract environment from the left
        left_env = 1;
        for j = 1:s-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end
        right_env = 1;
    else
        % Contract environment from the left
        left_env = 1;
        for j = 1:s-1
            left_env = contract_left_environment(left_env, M{j}, H{j}, M_{j});
        end
        
        % Contract environment from the right
        right_env = 1;
        for j =  N:-1:s+2
            right_env = contract_right_environment(right_env, M{j}, H{j}, M_{j});
        end
    end
end
