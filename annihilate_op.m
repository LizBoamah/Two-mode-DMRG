% An annihilation operator function
%-------------------------------------$
% Description
% This function returns the ith annihilation operator 

% Input
% N = the number of nodes
% i = ith annihilate operator

% Output
% ci = The ith annihilation operator 

% Define the basis operators
% [1 0 ] empty
% [ 0 1] occupied


function ci = annihilate_op(i,N)
    I = eye(2); %identity 
    Ph = diag([1 , -1 ]); % Phase operators
    c = [0 1 ; 0 0 ]; % annihilation operator 
    ci =1;
    for k = 1:N
        if k < i
            ci = kron(Ph ,ci );
        elseif k == i
            ci = kron(ci, c);
        else
        ci = kron(ci, I);
       end
    end
end




% function ci = annihilate_op(N)
%     I = eye(2); %identity 
%     Ph = diag([1 , -1 ]); % Phase operators
%     c = [0 1 ; 0 0 ]; % annihilation operator 
%     ci =cell(N, 1);
%     ci{1}= kron(c, I);
%     for k = 1:N-1
%         ci{k+1} = kron(Ph, ci{k});
%     end
% end



% Alternative to above
%     if k == 1 
%         ci{k} = kron(c, I);
%     else
%         ci{k} = kron(Ph, ci{k-1});
%     end
%      
 %identity 


    
