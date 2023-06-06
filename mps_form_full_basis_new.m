%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% This generates a matrix of size [Dim,N] where Dim<=q^N subspace corresponding%
% to states with Qnumbers_tot                                                  %
%                                                                              %
% Inputs:                                                                      %
% q             @integer dimension of local states                             % 
%               (later generalize this for site dependent qi)                  %
% N             @integer number of sites                                       %
% Qnumbers_full @integer[q,1] quantum number for each local state              %
% Qnumbers_tot  @integer[1,SymNum] each element is total value of a symmetry   %
%               operator                                                       %
% Model         uses only Model.Operators and Model.SiteOperators              %
%                                                                              %
% Outputs:                                                                     %
% Basis         @integer [Dim, N]                                              %
%                                                                              %
%  mps_form_full_basis(q,N,Qnumbers_full,Qnumbers_tot)                         %
% [Basis]=mps_form_full_basis(q,N,Qnumbers_full,Qnumbers_tot, Model)           %
%                                                                              %
%                                                                              %
% Examples for spinless fermi q=2, 6 sites 3 fermions:                         %
%  mps_form_full_basis(2,6,[0 1]',3, Model)                                    %
% mps_form_full_basis(2,6,[0 1]',NaN, Model)                                  %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Basis]=mps_form_full_basis_new(q,N,Qnumbers_full,Qnumbers_tot)

%---determine stie dependent local dimensions
%Model
%Model.SiteOperators
%Model.Operators
%Qnumbers_full, q, N, error('aa')


q_i = zeros(1,N);
for ii=1:N
  q_i(ii) = q;
  Qnumbers_full_i{ii} = Qnumbers_full;
end
%%Qnumbers_full_i = cell(1, length(Model.SiteOperators));
%%for l=1:length(Model.SiteOperators)
%%  q_i(l) = Model.Operators.(Model.SiteOperators{l}).State.Dim;
%%  Qnumbers_full_i{l} = state_full(Model.Operators.(Model.SiteOperators{l}).State); 
%%end
%%%q_i
%%%Qnumbers_full_i 
%%%Model.SiteOperators
%%%Model
%%%error('aa')

%q,N,Qnumbers_full, Qnumbers_tot, error('ii')
%Basis = ones(q^N, N); %maximum possible size of the basis
Basis = ones(prod(q_i), N); %initialize Basis with states 1

d  = 1; %to include [1 1 1 ... 1] 
dd = 0;
%Qnumbers_tot
if find(isnan(Qnumbers_tot)) == length(Qnumbers_tot) dd=1; end%if %to keep [1 1 1 ..] for [NaN, NaN, ...]
tmp = Basis(1,:);
%while d<q^N
while d<prod(q_i)
  d=d+1;
  tmp(N) = tmp(N)+1;
  %---check if upper limit is reached for any register
  for l=N:-1:1
    %if tmp(l) > q
    if tmp(l) > q_i(l)
      tmp(l) = 1;
      tmp(l-1) = tmp(l-1) + 1;
    end%if
   end%l
   qnum_tmp = zeros(size(Qnumbers_tot));
   %tmp
   for l=1:N
     for ss=1:size(Qnumbers_full_i{l},2)
       %l,ss
       %Qnumbers_full_i{l}
       %tmp
       qnum_tmp(ss) = qnum_tmp(ss) + (Qnumbers_full_i{l}(tmp(l), ss)) ;
     end%ss
   end%ll
   %qnum_tmp
   %pause_t1
   %---keep only reduced hilbert space
   if ~any(any(qnum_tmp-Qnumbers_tot))
     dd=dd+1;
     Basis(dd,:) = tmp; %, pause_t1
   end
end%while
Basis = Basis(1:dd,:);

return

