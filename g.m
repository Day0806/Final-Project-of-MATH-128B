function output = g( x, y, beta )
t = y - f(x, beta);
output = t' * t;
end

