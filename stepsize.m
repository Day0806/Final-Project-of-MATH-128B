function output = stepsize( x0, y0, beta, d, TOL)
alpha = 1;
while g(x0, y0, beta + alpha * d) >= g(x0, y0, beta)
    alpha = alpha / 2;
    if alpha < TOL
        break
    end
end
output = alpha;
end

