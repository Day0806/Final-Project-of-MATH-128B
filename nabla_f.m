function output = nabla_f( x, beta )
vector2 = [                         exp(-beta(2).*x) ...
                              -beta(1).*x .*exp(-beta(2).*x) ...
                         exp(-(beta(4) - x).^2/beta(5).^2) ...
 -(beta(3).*exp(-(beta(4) - x).^2/beta(5).^2).*(2.*beta(4) - 2.*x))/beta(5).^2 ...
  (2.*beta(3).*exp(-(beta(4) - x).^2/beta(5).^2).*(beta(4) - x).^2)/beta(5).^3 ...
                         exp(-(beta(7) - x).^2/beta(8).^2) ...
 -(beta(6).*exp(-(beta(7) - x).^2/beta(8).^2).*(2.*beta(7) - 2.*x))/beta(8).^2 ...
  (2.*beta(6).*exp(-(beta(7) - x).^2/beta(8).^2).*(beta(7) - x).^2)/beta(8).^3];
output = vector2';
end

