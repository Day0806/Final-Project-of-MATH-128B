function output = Armijo_stepsize( x0, y0, beta, d, s, gamma, sigma, TOL)
while g(x0, y0, beta) - g(x0, y0, beta + gamma * s * d) ...
        < -sigma * gamma * s * nabla_g(x0, y0, beta)' * d
    gamma = gamma * gamma;
    if gamma < TOL
        break
    end
end
output = gamma * s;
end