clear all;
close all;
randn('seed',1234);rand('seed',1234)
N1 = 1000;
N2 = N1;

n = 2;m=N1+N2;

mu1 = 1*randn(1,n);
mu2 = 1*randn(1,n);
X1 = randn(N1,n)+ones(N1,1)*mu1;
X2 = randn(N2,n)+ones(N2,1)*mu2;

A = [X1;X2];
y = [-ones(N1,1);ones(N2,1)];

lam = 1;

cvx_begin
cvx_solver coneos_matlab
cvx_solver_settings('RELAX_X',0)
cvx_solver_settings('RHOX',1e-3)
cvx_solver_settings('ALPHA',1.8)
variables w_c(n) b_c s_c(m)
minimize(1/2*sum_square(w_c) + lam*norm(s_c,1))
y.*(A*w_c + b_c) >= 1 - s_c;
cvx_end


cvx_begin
variables w(n) b s(m)
minimize(1/2*sum_square(w) + lam*norm(s,1))
y.*(A*w + b) >= 1 - s;
cvx_end

if n==2
    figure();hold on;
    for i=1:N1
        plot(X1(i,1),X1(i,2),'bo');
    end
    for i=1:N2
        plot(X2(i,1),X2(i,2),'rx');
    end
    
    x1 = min(A(:,1)):0.1:max(A(:,1));

    x2 = (-b_c - x1.*w_c(1))./w_c(2);
    plot(x1,x2,'g');
    x2 = (-b - x1.*w(1))./w(2);
    plot(x1,x2,'k');
    
end

