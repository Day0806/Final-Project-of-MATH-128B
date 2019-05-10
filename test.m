    m = 4;
    k = 0;
    y = zeros(m0, m+1);
    s = zeros(m0, m+1);
    a = zeros(m+1, 1);
    % Initial preparation
    %tic
    while k < m
        k = k + 1;
        d = - B0 \ I * nabla_g(x0, y0, beta0);
        alpha = stepsize(x0, y0, beta0, d, TOL2);
        %alpha = Armijo_stepsize(x0, y0, beta0, d, 1/10, 1/10, 1e-5, TOL);
        beta = beta0 + alpha * d;
        if abs(g(x0, y0, beta0) - g(x0, y0, beta)) < TOL1
            break
        end
        s(:, k+1) = beta - beta0;
        y(:, k+1) = nabla_g(x0, y0, beta) - nabla_g(x0, y0, beta0);
        B0 = B0 + (y(:, k+1) * y(:, k+1)') / (y(:, k+1)' * s(:, k+1)) - (B0 * s(:, k+1) * s(:, k+1)' * B0) / (s(:, k+1)' * B0 * s(:, k+1));
        beta0 = beta;
        output = g(x0, y0, beta0);
    end

    %limit implement
    z = B0 \ I * q;
    while k < N
        k = k + 1;
        d = -z;
        % Choice of step size
        alpha = stepsize(x0, y0, beta0, d, TOL2);
        %alpha = Armijo_stepsize(x0, y0, beta0, d, 1/10, 1/10, 1e-5, TOL);
        beta = beta0 + alpha * d;
        if abs(g(x0, y0, beta0) - g(x0, y0, beta)) < TOL1
            break
        end
        ss = beta - beta0;
        s = [s(:, 2:m+1) ss];
        yy = nabla_g(x0, y0, beta) - nabla_g(x0, y0, beta0);
        y = [y(:, 2:m+1) yy];

        % compute inv(B_{k+1}) * nabla_g(x_{k+1})
        q = nabla_g(x0, y0, beta0);
        for t = m + 1 : -1 : 1
            a(t) = (s(:, t)' * q) / (s(:, t)' * y(:, t));
            q = q - a(t) * y(:, t);
        end
        z = B0 \ I * q;
        for t = 1 : m + 1
            b = (y(:, t)' * z) / (s(:, t)' * y(:, t));
            z = z + (a(t) - b) * s(:, t);
        end
        beta0 = beta;
    end
    output = g(x0, y0, beta0);
    %toc