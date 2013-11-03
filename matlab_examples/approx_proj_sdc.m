clear all; close all
n = 100;
A = randn(n);A = A + A';

[V,E] = eig(A);
E = diag(E);

% test many steps
for k=1:100
    
    dA = 0.01*randn(n); dA = (dA + dA')/2;
    
    W = V'*dA*V;
    
    E_plus = E + diag(W);
    
    for i=1:n
        mEi = 1;
        for j=1:n
            if (i~=j)
                mEi = min(mEi,abs(E(i) - E(j)));
            end
        end
        if (mEi < 1e-3)
            for j=1:n
                if (i==j)
                    W(i,j) = 1;
                else
                    W(i,j) = 0;
                end
            end
        else
            for j=1:n
                if (i==j)
                    W(i,j) = 1;
                else
                    W(i,j) = W(i,j)/(E(i) - E(j));
                end
            end
        end
    end
    
    V_plus = V*W;
    V_plus = V_plus./(ones(n,1)*norms(V_plus));
    
    % [U,S] = eig(A+dA);
    % S = diag(S);
    
    del(k) = norm(dA,'fro')/norm(A + dA,'fro');
    
    err(k) = norm(A + dA - V_plus*diag(E_plus)*V_plus','fro')/norm(A + dA,'fro');
    
    A = A+dA;
    V = V_plus;
    E = E_plus;
    
end
plot(err); hold on; plot(del,'r')