function output = f( x, beta )
output = beta(1) * exp(-beta(2) * x) + beta(3) * exp(-(x - beta(4)).^2 / beta(5).^2) ...
    +beta(6) * exp(-(x - beta(7)).^2 / beta(8).^2);
end

