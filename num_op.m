% An number operator function for the down spin
%-------------------------------------$
% Description
% This function returns the ith number operator 

% Input
% N = the number of nodes
% i = ith number operator

% Output
% no = The ith number operator 

function  no = num_op(i , N)
    for k = 1:N
        no = creation_op(i, N) * annihilate_op(i, N);
    end    
end












% function  no = num_op(N)
%     cr = creation_op(N);
%     an = annihilate_op(N);
%     no = cell(1,N);
%     for k= 1:N
%         no{k} = cr{k} * an{k};  
%     end    
% end