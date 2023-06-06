% A creation operator function
%-------------------------------------$
% Description
% This function returns the ith creator operator 

% Input

% N = the number of nodes
% i = ith creation operator

% Output
% cpi = The ith creation operator 

% Define the basis operators
% [1 0 ] empty
% [0 1 ] occupied
function cpi = creation_op(i, N)
    cpi = annihilate_op(i,N)'; 
end






% function cpi = creation_op(N)
%     an = annihilate_op(N);
%     cpi = cell(1, N);
%     for k = 1:N
%         cpi{k}= an{k}';
%     end