function [P,status,optval,scale] = cvxDC(E,G,normalization,mode) %#ok<INUSD,STOUT>

if normalization==1
    flag = [];
    counter = 0;
    for i = 1:length(E)
        if E(i)>1
            counter = counter+1;
            flag(counter) = i;
            G(i,:) = G(i,:)/E(i);
        end
    end
    %E = E(flag);
    E = ones(1,counter);
    G = G(flag,:);
end

[n_genes,n_cells] = size (G);
U = ones(n_cells+2,1);
U(n_cells+1,1) = -1;
U(n_cells+2,1) = 0;
Z = zeros(1,n_genes);
Y = [G';Z;-E];

cvx_begin quiet
eval(['cvx_solver ',mode]);
variable X(1,n_cells+2) nonnegative
minimize(X*Y*Y'*X');
subject to
X*U == 0;
X(1,n_cells+2) == 1;
cvx_end

P = X(1,1:n_cells);
scale = X(1,n_cells+1);
status = cvx_status;
optval = (P*G'-E)*(P*G'-E)';
