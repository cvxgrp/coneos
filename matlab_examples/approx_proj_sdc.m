clear all; close all
randn('seed',0);rand('seed',0)
n = 500;
A = randn(n);A = A + A';

[V,E] = eig(A);
E = diag(E);

A_orig = A;
V_orig = V;
E_orig = E;
E_t = E;
% test many steps
for k=1:100
    
    dA = 0.001*randn(n); dA = (dA + dA')/2;
    
    W = V'*dA*V;
    
    E_plus = E + diag(W);
    
    for i=1:n
        %{
        mEi = 1;
        for j=1:n
            if (i~=j)
                mEi = min(mEi,abs(E(i) - E(j)));
            end
        end
        if (mEi < 1e-6)
            for j=1:n
                if (i==j)
                    W(i,j) = 1;
                else
                    W(i,j) = 0;
                end
            end
        else
        %}
            for j=1:n
                if (i==j)
                    W(i,j) = 1;
                else
                    W(i,j) = W(i,j)/(E(j) - E(i));
                end
            end
        %end
    end
    
    V_plus = V*W;
    V_plus = V_plus./(ones(n,1)*norms(V_plus));
    
    % [U,S] = eig(A+dA);
    % S = diag(S);
    
    del(k) = norm(dA,'fro')/norm(A,'fro');
    
    err(k) = norm(A + dA - V_plus*diag(E_plus)*V_plus','fro')/norm(A + dA,'fro');
    
    
    W = V_orig'*dA*V_orig;
    E_t = E_t + diag(W);

    err3(k) = norm(A + dA - V_orig*diag(E_t)*V_orig','fro')/norm(A,'fro');
    err4(k) = norm(A + dA - A_orig,'fro')/norm(A,'fro');

    
    A = A+dA;
    V = V_plus;
    E = E_plus;
    

    dA2 = A - A_orig;
    W = V_orig'*dA2*V_orig;
    
    E_plus = E_orig + diag(W);
    
    for i=1:n
        %{
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
        %}
            for j=1:n
                if (i==j)
                    W(i,j) = 1;
                else
                    W(i,j) = W(i,j)/(E(j) - E(i));
                end
            end
        %end
    end
        
    V_plus = V_orig*W;
    V_plus = V_plus./(ones(n,1)*norms(V_plus));

    err2(k) = norm(A - V_plus*diag(E_plus)*V_plus','fro')/norm(A,'fro');

    
end
plot(err); hold on; plot(err2,'r'); plot(err3,'g'); plot(err4,'k'); plot(del,'--')