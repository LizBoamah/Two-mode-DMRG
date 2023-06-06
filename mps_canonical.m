function [M] = mps_canonical(Psi,d, N, dir,s)

%Psi = Psi/(norm(Psi));
    ml = 1;
    m_mix =1;

    % Initialise MPS - tensor at site n
    M = cell(1,N);
%     for im = 1:N
%         M{im} = cell(1, d);
%     end

  switch dir

        case 'left'
            %--- beginning to N-1
            for l = 1:N - 1
                W = reshape(Psi, [ml * d, d ^ (N - l)]); 
                [U, S, V] = svd(W, 'econ');
                new_ml = size(U, 2);

                U = reshape(U, [ml, d, new_ml]);
                V = S * V';

%                 
                    M{l} = U(:, :, :);
%                 

                Psi = V;
                ml = new_ml;
            end
            %---last mode
            Psi = reshape(Psi, [ml, d, 1]);
%             for qq = 1:d
%                 M{N}{qq} = permute(reshape(Psi(:, qq), [ml, 1, d]), [1, 3, 2]);
                 M{N} = Psi(:,:,:);
%             end


         case 'right' % chk Ok
            %--- N to 2
            for l=N:-1:2
              W=reshape(Psi, d^(l-1),ml*d);
              [U,S,V]=svd(W, 'econ');
              V = ctranspose(V);
              new_ml = size(V, 1);
              V = reshape(V, [new_ml, d, ml]);
             
                M{l} = V(:,:,:);
              
              Psi = U*S;
              ml = new_ml;
            end
             %--- first mode
            
            Psi = reshape(Psi, [1, d, ml]);
           
              M{1} = Psi(:,:, :);
            
         case 'mixed' %chk ok
             %left canonical from 1 to s-1
                for l = 1:s-1
                    W = reshape(Psi, [ml * d, d ^ (N - l)]);
                    [U, S, V] = svd(W, 'econ');
                    new_ml = size(U, 2);

                    U = reshape(U, [ml, d, new_ml]);
                    V = S * V';

                    
                        M{l} = U(:, :, :);
                    
    
                    Psi = V;
                    ml = new_ml;
                end
 

                % right canonical from N to s+1
                for l=N:-1:s+2
                        W=reshape(Psi, d^(l-1),m_mix*d);
                      [U,S,V]=svd(W, 'econ');
                      V = ctranspose(V);
                      new_mlx = size(V, 1);
                      V = reshape(V, [new_mlx, d, m_mix]);
                     
                        M{l}= V(:,:,:);
                     
                      Psi = U*S;
                      m_mix = new_mlx; 

                end 
%                  %--- s-th mode
%                  %Psi = reshape(Psi, [m_mix, d, new_ml]);
%                  M{s} = reshape(Psi, [d, ml*d, new_mlx]);
%                    M{s} = Psi(:,:,:);
            % Sites s and s+1
            W = reshape(Psi, [ml * d, m_mix * d]);
            [U, S, V] = svd(W, 'econ');
            
            new_ml = size(U, 2);
            U = reshape(U, [ml, d, new_ml]);
            M{s} = U(:, :, :);
            
            new_mlx = size(V, 1);
            V = reshape(V, [new_ml, d, m_mix]);
            M{s+1} = V(:, :, :);
                

           
         

 end
