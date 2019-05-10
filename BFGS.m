% BFGS Method
k = 0;
N = 500;
TOL1 = 1e-3;
TOL2 = 1e-3;
%beta0 = [94.9 0.009 90.1 113.0 20.0 73.8 140.0 20.0]';
beta0 = [96.0 0.0096 80.0 110.0 25.0 74.0 139.0 25.0]';
J = nabla_f(x0, beta);
B0 = J * J';
% B0 = Hessian(x0, y0, beta0);
tic
while k < N
    k = k + 1;
    d = - B0 \ I * nabla_g(x0, y0, beta0);
    %alpha = stepsize(x0, y0, beta0, d, TOL1);
    alpha = Armijo_stepsize(x0, y0, beta0, d, 0.2, 0.2, 1e-5, TOL2);
    beta = beta0 + alpha * d;
    if abs(g(x0, y0, beta0) - g(x0, y0, beta)) < TOL1
        break
    end
    s = beta - beta0;
    y = nabla_g(x0, y0, beta) - nabla_g(x0, y0, beta0);
    B0 = B0 + (y * y') / (y' * s) - (B0 * s * s' * B0) / (s' * B0 * s);
    beta0 = beta;
end
norm(nabla_g(x0, y0, beta))
g(x0, y0, beta)
toc